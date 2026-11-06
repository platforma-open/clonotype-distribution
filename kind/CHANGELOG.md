# @platforma-open/milaboratories.spatiotemporal-analysis.kind

## 2.0.0

### Major Changes

- 4fa5937: Columns are stored by their result-pool column ids: the abundance, the subset it is narrowed to, and the subject, grouping and temporal columns. The timepoints of the timepoint order are stored as column value references of the temporal column, so the block that produced the metadata knows which of its values this block uses. The model finds columns and builds the options of every column picker through the column collection API instead of the deprecated result pool API, and the workflow takes every column through a column bundle. The settings panel reads the timepoints from a frame of the temporal column alone, and shows why when they cannot be read.
  Opening a project created before this change keeps its abundance and subset, and clears its subject, grouping and temporal columns and the timepoint order, since they were stored by ids anchored to the abundance column and as plain values, which name no producing block; they have to be picked again before the block runs. The initialization params of the kind carry the same ids and references, so templates made for an earlier version of the kind cannot be applied.
  Timepoints that are no longer among the values of the temporal column stay in the order and are marked as not available instead of being removed silently, with a tooltip naming the block that was to provide them, as a link to it, and what went wrong. Picking another temporal column still resets the order to the values of that column.

## 1.1.0

### Minor Changes

- 0e45d30: Update the SDK and accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) with the abundance. The abundance list is built with buildDatasetOptions (requireEnrichments: false), with the same labels as before. The subset restricts the clonotypes before the analysis, so frequencies, restriction indices and prevalence describe the subset only. Exported columns of a subset run carry pl7.app/inputSubset and a trace rooted in the subset. Unfiltered runs are unchanged.
