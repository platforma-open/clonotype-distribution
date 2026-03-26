"""
In vivo compartment analysis: computes spatial restriction, temporal kinetics,
and cross-subject convergence metrics for clonal/cluster abundance data.

Input: CSV with columns [sampleId, elementId, abundance, <metadata columns>]
Output: Multiple CSV files with computed metrics.
"""

import argparse
import json
import math
import sys

import numpy as np
import polars as pl


def parse_args():
    parser = argparse.ArgumentParser(description="In vivo compartment analysis")
    parser.add_argument("input_file", help="Input CSV file")
    parser.add_argument("--calculation-mode", choices=["population", "intra-subject"],
                        default="population")
    parser.add_argument("--normalization", choices=["relative-frequency", "clr"],
                        default="relative-frequency")
    parser.add_argument("--compartment-columns", type=str, default="[]",
                        help="JSON array of compartment column names")
    parser.add_argument("--temporal-column", type=str, default="",
                        help="Temporal metadata column name")
    parser.add_argument("--timepoint-order", type=str, default="[]",
                        help="JSON array of timepoint values in order")
    parser.add_argument("--subject-column", type=str, required=True,
                        help="Subject metadata column name")
    parser.add_argument("--presence-threshold", type=float, default=0.0)
    parser.add_argument("--pseudo-count", type=float, default=1.0)
    parser.add_argument("--output-prefix", type=str, default="output")
    return parser.parse_args()


def compute_relative_frequency(df: pl.DataFrame) -> pl.DataFrame:
    """Compute per-sample relative frequency from abundance."""
    sample_totals = df.group_by("sampleId").agg(
        pl.col("abundance").sum().alias("sampleTotal")
    )
    df = df.join(sample_totals, on="sampleId")
    # Filter out zero-total samples
    df = df.filter(pl.col("sampleTotal") > 0)
    df = df.with_columns(
        (pl.col("abundance") / pl.col("sampleTotal")).alias("frequency")
    ).drop("sampleTotal")
    return df


def compute_clr(df: pl.DataFrame) -> pl.DataFrame:
    """Compute centered log-ratio transform per sample."""
    sample_totals = df.group_by("sampleId").agg(
        pl.col("abundance").sum().alias("sampleTotal")
    )
    df = df.join(sample_totals, on="sampleId")
    df = df.filter(pl.col("sampleTotal") > 0)
    df = df.with_columns(
        (pl.col("abundance") / pl.col("sampleTotal")).alias("frequency")
    ).drop("sampleTotal")

    # Multiplicative replacement for zeros
    min_nonzero = df.filter(pl.col("frequency") > 0)["frequency"].min()
    if min_nonzero is None or min_nonzero <= 0:
        min_nonzero = 1e-10
    delta = float(min_nonzero) * 0.65

    def clr_transform(group: pl.DataFrame) -> pl.DataFrame:
        freq = group["frequency"].to_numpy().astype(np.float64)
        # Replace zeros
        freq = np.where(freq == 0, delta, freq)
        freq = freq / freq.sum()
        # CLR
        log_freq = np.log(freq)
        geo_mean = np.mean(log_freq)
        clr_vals = log_freq - geo_mean
        return group.with_columns(pl.Series("frequency", clr_vals))

    df = df.group_by("sampleId", maintain_order=True).map_groups(clr_transform)
    return df


def shannon_entropy(p: np.ndarray) -> float:
    """Compute Shannon entropy H(p) = -sum(p_i * log2(p_i)) for p_i > 0."""
    p = p[p > 0]
    if len(p) == 0:
        return 0.0
    p = p / p.sum()
    return -float(np.sum(p * np.log2(p)))


def restriction_index(freq_by_compartment: np.ndarray) -> float:
    """RI = 1 - H(p) / log2(N) where N = number of compartments with nonzero presence."""
    nonzero = freq_by_compartment[freq_by_compartment > 0]
    n = len(nonzero)
    if n <= 1:
        return 1.0
    h = shannon_entropy(nonzero)
    return 1.0 - h / math.log2(n)


def compute_spatial_metrics(
    df: pl.DataFrame,
    compartment_col: str,
    subject_col: str,
    mode: str,
    presence_threshold: float,
) -> pl.DataFrame:
    """Compute RI, dominant compartment, and spatial breadth for one compartment variable."""
    categories = sorted(df[compartment_col].unique().to_list())
    n_categories = len(categories)

    if n_categories == 0:
        return pl.DataFrame()

    # Group by element and subject, get mean frequency per compartment
    per_subject_compartment = (
        df.group_by(["elementId", subject_col, compartment_col])
        .agg(pl.col("frequency").mean().alias("meanFreq"))
    )

    results = []

    if mode == "population":
        # Per-subject spatial metrics first
        per_subject_metrics = _compute_per_subject_spatial(
            per_subject_compartment, categories, subject_col, presence_threshold
        )

        # Consensus across subjects
        for element_id, group in per_subject_metrics.group_by("elementId"):
            eid = element_id[0] if isinstance(element_id, tuple) else element_id
            ris = group["ri"].drop_nulls().to_numpy()
            dominants = group["dominant"].to_list()
            subjects = group[subject_col].to_list()

            subject_prevalence = len(subjects)
            mean_ri = float(np.mean(ris)) if len(ris) > 0 else float("nan")
            std_ri = float(np.std(ris, ddof=1)) if len(ris) > 1 else float("nan")

            # Consensus dominant: mode, ties broken by highest mean frequency
            if dominants:
                dom_counts: dict[str, int] = {}
                dom_freq_sums: dict[str, float] = {}
                for i, d in enumerate(dominants):
                    if d is not None:
                        dom_counts[d] = dom_counts.get(d, 0) + 1
                        dom_freq_sums[d] = dom_freq_sums.get(d, 0.0)
                        row_freq = group.filter(pl.col(subject_col) == subjects[i])
                        if len(row_freq) > 0 and "dominant_freq" in row_freq.columns:
                            dom_freq_sums[d] += float(row_freq["dominant_freq"][0] or 0)

                max_count = max(dom_counts.values()) if dom_counts else 0
                tied = [k for k, v in dom_counts.items() if v == max_count]
                if len(tied) == 1:
                    consensus_dominant = tied[0]
                else:
                    consensus_dominant = max(tied, key=lambda k: dom_freq_sums.get(k, 0))
            else:
                consensus_dominant = None

            # Count dominant in each category
            count_dominant = {cat: 0 for cat in categories}
            for d in dominants:
                if d is not None and d in count_dominant:
                    count_dominant[d] += 1

            # Mean spatial breadth
            breadths = group["spatialBreadth"].to_numpy()
            mean_breadth = int(round(float(np.mean(breadths)))) if len(breadths) > 0 else 0

            row = {
                "elementId": eid,
                "ri": mean_ri,
                "dominant": consensus_dominant,
                "spatialBreadth": mean_breadth,
                "subjectPrevalence": subject_prevalence,
                "consensusDominant": consensus_dominant,
                "meanRi": mean_ri,
                "stdRi": std_ri,
            }
            for cat in categories:
                row[f"countDominantIn_{cat}"] = count_dominant[cat]
            results.append(row)

    else:
        # Intra-subject mode: per-subject metrics directly
        per_subject_metrics = _compute_per_subject_spatial(
            per_subject_compartment, categories, subject_col, presence_threshold
        )
        # Also compute aggregated summary
        for element_id, group in per_subject_metrics.group_by("elementId"):
            eid = element_id[0] if isinstance(element_id, tuple) else element_id
            ris = group["ri"].drop_nulls().to_numpy()
            dominants = group["dominant"].to_list()
            subjects = group[subject_col].to_list()

            subject_prevalence = len(subjects)
            mean_ri = float(np.mean(ris)) if len(ris) > 0 else float("nan")
            std_ri = float(np.std(ris, ddof=1)) if len(ris) > 1 else float("nan")

            # Consensus dominant
            if dominants:
                dom_counts: dict[str, int] = {}
                for d in dominants:
                    if d is not None:
                        dom_counts[d] = dom_counts.get(d, 0) + 1
                max_count = max(dom_counts.values()) if dom_counts else 0
                tied = [k for k, v in dom_counts.items() if v == max_count]
                consensus_dominant = tied[0] if tied else None
            else:
                consensus_dominant = None

            count_dominant = {cat: 0 for cat in categories}
            for d in dominants:
                if d is not None and d in count_dominant:
                    count_dominant[d] += 1

            breadths = group["spatialBreadth"].to_numpy()
            mean_breadth = int(round(float(np.mean(breadths)))) if len(breadths) > 0 else 0

            row = {
                "elementId": eid,
                "ri": mean_ri,
                "dominant": consensus_dominant,
                "spatialBreadth": mean_breadth,
                "subjectPrevalence": subject_prevalence,
                "consensusDominant": consensus_dominant,
                "meanRi": mean_ri,
                "stdRi": std_ri,
            }
            for cat in categories:
                row[f"countDominantIn_{cat}"] = count_dominant[cat]
            results.append(row)

    if not results:
        return pl.DataFrame()
    return pl.DataFrame(results)


def _compute_per_subject_spatial(
    per_subject_compartment: pl.DataFrame,
    categories: list[str],
    subject_col: str,
    presence_threshold: float,
) -> pl.DataFrame:
    """Compute per-subject RI, dominant compartment, spatial breadth."""
    results = []

    for (element_id, subject), group in per_subject_compartment.group_by(
        ["elementId", subject_col]
    ):
        freq_map = dict(zip(group[group.columns[2]].to_list(), group["meanFreq"].to_list()))
        freq_arr = np.array([freq_map.get(cat, 0.0) for cat in categories])

        ri = restriction_index(freq_arr)
        dominant_idx = int(np.argmax(freq_arr))
        dominant = categories[dominant_idx] if freq_arr[dominant_idx] > 0 else None
        dominant_freq = float(freq_arr[dominant_idx])
        breadth = int(np.sum(freq_arr > presence_threshold))

        results.append({
            "elementId": element_id,
            subject_col: subject,
            "ri": ri,
            "dominant": dominant,
            "dominant_freq": dominant_freq,
            "spatialBreadth": breadth,
        })

    if not results:
        return pl.DataFrame()
    return pl.DataFrame(results)


def compute_temporal_metrics(
    df: pl.DataFrame,
    temporal_col: str,
    timepoint_order: list[str],
    subject_col: str,
    mode: str,
    pseudo_count: float,
) -> pl.DataFrame:
    """Compute peak timepoint, temporal shift index, and log2 kinetic delta."""
    if len(timepoint_order) < 2:
        return pl.DataFrame()

    rank_map = {tp: i for i, tp in enumerate(timepoint_order)}
    t_count = len(timepoint_order)

    # Filter to known timepoints
    df = df.filter(pl.col(temporal_col).is_in(timepoint_order))

    if mode == "population":
        # Average frequency across subjects at each timepoint
        per_tp = (
            df.group_by(["elementId", temporal_col])
            .agg(pl.col("frequency").mean().alias("meanFreq"))
        )
    else:
        # Per subject per timepoint
        per_tp = (
            df.group_by(["elementId", subject_col, temporal_col])
            .agg(pl.col("frequency").mean().alias("meanFreq"))
        )

    results = []

    if mode == "population":
        for element_id, group in per_tp.group_by("elementId"):
            eid = element_id[0] if isinstance(element_id, tuple) else element_id
            tp_freq = dict(zip(group[temporal_col].to_list(), group["meanFreq"].to_list()))
            row = _compute_temporal_for_element(eid, tp_freq, timepoint_order, rank_map, t_count, pseudo_count)
            results.append(row)
    else:
        # Average per-subject temporal metrics
        per_subject_temporal = []
        for (element_id, subject), group in per_tp.group_by(["elementId", subject_col]):
            eid = element_id if not isinstance(element_id, tuple) else element_id
            tp_freq = dict(zip(group[temporal_col].to_list(), group["meanFreq"].to_list()))
            row = _compute_temporal_for_element(eid, tp_freq, timepoint_order, rank_map, t_count, pseudo_count)
            per_subject_temporal.append(row)

        if not per_subject_temporal:
            return pl.DataFrame()

        ps_df = pl.DataFrame(per_subject_temporal)
        for element_id, group in ps_df.group_by("elementId"):
            eid = element_id[0] if isinstance(element_id, tuple) else element_id
            row = {
                "elementId": eid,
                "peakTimepoint": group["peakTimepoint"].mode().to_list()[0] if len(group) > 0 else None,
                "temporalShiftIndex": float(group["temporalShiftIndex"].mean()) if len(group) > 0 else float("nan"),
                "log2KineticDelta": float(group["log2KineticDelta"].mean()) if len(group) > 0 else float("nan"),
            }
            results.append(row)

    if not results:
        return pl.DataFrame()
    return pl.DataFrame(results)


def _compute_temporal_for_element(
    element_id: str,
    tp_freq: dict[str, float],
    timepoint_order: list[str],
    rank_map: dict[str, int],
    t_count: int,
    pseudo_count: float,
) -> dict:
    """Compute temporal metrics for a single element."""
    freqs = np.array([tp_freq.get(tp, 0.0) for tp in timepoint_order])

    # Peak timepoint
    peak_idx = int(np.argmax(freqs))
    peak_tp = timepoint_order[peak_idx]

    # Temporal shift index
    total_freq = freqs.sum()
    if total_freq > 0 and t_count > 1:
        ranks = np.arange(t_count, dtype=np.float64)
        tsi = float((np.sum(ranks * freqs) / total_freq) / (t_count - 1))
    else:
        tsi = float("nan")

    # Log2 kinetic delta
    # Find first and last timepoints with nonzero frequency
    nonzero_indices = np.where(freqs > 0)[0]
    if len(nonzero_indices) >= 2:
        first_freq = float(freqs[nonzero_indices[0]])
        last_freq = float(freqs[nonzero_indices[-1]])
        log2kd = math.log2((last_freq + pseudo_count) / (first_freq + pseudo_count))
    elif len(nonzero_indices) == 1:
        log2kd = 0.0
    else:
        log2kd = 0.0

    return {
        "elementId": element_id,
        "peakTimepoint": peak_tp,
        "temporalShiftIndex": tsi,
        "log2KineticDelta": log2kd,
    }


def compute_subject_prevalence(df: pl.DataFrame, subject_col: str) -> pl.DataFrame:
    """Count distinct subjects per element."""
    return (
        df.filter(pl.col("abundance") > 0)
        .group_by("elementId")
        .agg(pl.col(subject_col).n_unique().alias("subjectPrevalence"))
    )


def build_heatmap_data(
    df: pl.DataFrame,
    compartment_col: str,
    subject_col: str,
    mode: str,
) -> pl.DataFrame:
    """Build heatmap data: element x compartment -> normalized frequency."""
    if mode == "population":
        # Average across subjects
        heatmap = (
            df.group_by(["elementId", compartment_col])
            .agg(pl.col("frequency").mean().alias("normalizedFrequency"))
        )
    else:
        heatmap = (
            df.group_by(["elementId", compartment_col])
            .agg(pl.col("frequency").mean().alias("normalizedFrequency"))
        )
    return heatmap.rename({compartment_col: "compartmentCategory"})


def build_temporal_line_data(
    df: pl.DataFrame,
    temporal_col: str,
    timepoint_order: list[str],
    mode: str,
) -> pl.DataFrame:
    """Build temporal line plot data: element x timepoint -> frequency."""
    df = df.filter(pl.col(temporal_col).is_in(timepoint_order))

    if mode == "population":
        line_data = (
            df.group_by(["elementId", temporal_col])
            .agg(pl.col("frequency").mean().alias("frequency"))
        )
    else:
        line_data = (
            df.group_by(["elementId", temporal_col])
            .agg(pl.col("frequency").mean().alias("frequency"))
        )
    return line_data.rename({temporal_col: "timepoint"})


def build_prevalence_histogram(prevalence_df: pl.DataFrame) -> pl.DataFrame:
    """Build prevalence histogram: prevalenceCount -> cloneCount."""
    return (
        prevalence_df.group_by("subjectPrevalence")
        .agg(pl.len().alias("cloneCount"))
        .rename({"subjectPrevalence": "prevalenceCount"})
        .sort("prevalenceCount")
    )


def main():
    args = parse_args()

    # Read input
    df = pl.read_csv(args.input_file)

    # Ensure abundance is numeric (may arrive as string from XSV export, with empty strings for missing values)
    if df["abundance"].dtype == pl.Utf8 or df["abundance"].dtype == pl.String:
        df = df.with_columns(
            pl.when(pl.col("abundance") == "")
            .then(pl.lit(0.0))
            .otherwise(pl.col("abundance").cast(pl.Float64))
            .alias("abundance")
        )
    elif not df["abundance"].dtype.is_float():
        df = df.with_columns(pl.col("abundance").cast(pl.Float64))

    compartment_columns = json.loads(args.compartment_columns)
    temporal_column = args.temporal_column if args.temporal_column else None
    timepoint_order = json.loads(args.timepoint_order)
    subject_col = args.subject_column
    mode = args.calculation_mode
    prefix = args.output_prefix

    # Normalize
    if args.normalization == "clr":
        df = compute_clr(df)
    else:
        df = compute_relative_frequency(df)

    # Subject prevalence (always computed)
    prevalence = compute_subject_prevalence(df, subject_col)
    prevalence.write_csv(f"{prefix}_prevalence.csv")

    # Prevalence histogram
    histogram = build_prevalence_histogram(prevalence)
    histogram.write_csv(f"{prefix}_prevalence_histogram.csv")

    # Spatial metrics per compartment variable
    for comp_col in compartment_columns:
        spatial = compute_spatial_metrics(df, comp_col, subject_col, mode, args.presence_threshold)
        if len(spatial) > 0:
            spatial.write_csv(f"{prefix}_spatial_{comp_col}.csv")

        # Heatmap data
        heatmap = build_heatmap_data(df, comp_col, subject_col, mode)
        if len(heatmap) > 0:
            heatmap.write_csv(f"{prefix}_heatmap_{comp_col}.csv")

    # Temporal metrics
    if temporal_column and len(timepoint_order) >= 2:
        temporal = compute_temporal_metrics(
            df, temporal_column, timepoint_order, subject_col, mode, args.pseudo_count
        )
        if len(temporal) > 0:
            temporal.write_csv(f"{prefix}_temporal.csv")

        # Temporal line data
        line_data = build_temporal_line_data(df, temporal_column, timepoint_order, mode)
        if len(line_data) > 0:
            line_data.write_csv(f"{prefix}_temporal_line.csv")

    # Write combined main table
    # Start with prevalence
    main_table = prevalence

    # Join spatial metrics
    for comp_col in compartment_columns:
        spatial_file = f"{prefix}_spatial_{comp_col}.csv"
        try:
            spatial = pl.read_csv(spatial_file)
            # Rename columns to avoid conflicts
            rename_map = {}
            for col in spatial.columns:
                if col != "elementId":
                    rename_map[col] = f"{col}_{comp_col}"
            spatial = spatial.rename(rename_map)
            main_table = main_table.join(spatial, on="elementId", how="left")
        except Exception:
            pass

    # Join temporal metrics
    if temporal_column and len(timepoint_order) >= 2:
        temporal_file = f"{prefix}_temporal.csv"
        try:
            temporal = pl.read_csv(temporal_file)
            main_table = main_table.join(temporal, on="elementId", how="left")
        except Exception:
            pass

    main_table.write_csv(f"{prefix}_main.csv")

    print(f"Analysis complete. Mode: {mode}, Normalization: {args.normalization}")
    print(f"Compartment variables: {compartment_columns}")
    if temporal_column:
        print(f"Temporal variable: {temporal_column}, Order: {timepoint_order}")
    print(f"Results written with prefix: {prefix}")


if __name__ == "__main__":
    main()
