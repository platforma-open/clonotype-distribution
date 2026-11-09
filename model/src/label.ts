import type { BlockData } from "./types";

/** Settings part of the block label: the calculation mode, and CLR when chosen. */
export function labelFromParams(data: BlockData): string {
  return [calculationModeLabel(data), normalizationLabel(data)].filter(Boolean).join(", ");
}

/**
 * Subtitle shown when the user has not set a custom block label: the calculation mode, the chosen
 * columns' result-pool labels, then CLR when chosen.
 */
export function defaultSubtitle(data: BlockData): string {
  return [calculationModeLabel(data), data.subtitleLabel, normalizationLabel(data)]
    .filter(Boolean)
    .join(", ");
}

// Internals

function calculationModeLabel(data: BlockData): string {
  return data.calculationMode === "intra-subject" ? "Intra-Subject" : "Population";
}

/**
 * CLR is the non-default normalization; it is named so blocks that differ only in normalization
 * stay distinguishable downstream.
 */
function normalizationLabel(data: BlockData): string | undefined {
  return data.normalization === "clr" ? "CLR" : undefined;
}
