# Spatiotemporal Analysis

Platforma block that quantifies how clonotypes or clusters distribute across tissues, timepoints, and subjects.

## Metrics

**Spatial** (requires compartment variable):
- Restriction Index (RI) — Shannon-entropy-based measure of compartment specificity (0 = uniform, 1 = exclusive)
- Dominant Compartment — category with highest abundance
- Spatial Breadth — number of compartments with presence above threshold

**Temporal** (requires timepoint variable with 2+ ordered values):
- Peak Timepoint — timepoint with highest mean frequency
- Temporal Shift Index (TSI) — frequency-weighted position on the time axis (0 = early, 1 = late)
- Log2 Kinetic Delta — log2 fold-change between first and last observed timepoints

**Cross-subject** (population-level mode):
- Subject Prevalence — number of subjects where a clone is detected
- Consensus Dominant — most common dominant compartment across subjects
- Mean RI / RI StdDev — cross-subject restriction index statistics

## Inputs

- **Abundance column** — clonotype/cluster count or UMI column from MiXCR clonotyping or similar
- **Subject variable** — metadata column identifying biological subjects (required)
- **Compartment variable** — metadata column for tissue/organ/site (optional)
- **Temporal variable** — metadata column for timepoints (optional)

## Modes

- **Population-Level** — averages metrics across subjects, reports consensus statistics
- **Intra-Subject** — computes per-subject metrics, then summarizes

## Development

```bash
pnpm install
pnpm build        # production build
pnpm build:dev    # dev build (local packages)
pnpm test         # integration tests (requires running platforma backend)
```
