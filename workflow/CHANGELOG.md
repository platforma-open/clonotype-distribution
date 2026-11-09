# @platforma-open/milaboratories.spatiotemporal-analysis.workflow

## 1.0.1

### Patch Changes

- b6fe3c6: Compute the block label without the UI

  The trace label is now built by the workflow from the native labels of the chosen grouping,
  temporal and subject columns plus the calculation mode and normalization, so opening the block
  no longer changes its arguments and a block that was never opened still gets a meaningful label.
  The subtitle shows the result-pool labels of the chosen columns. Existing blocks show as stale
  once after the update.

## 1.0.0

### Major Changes

- 4fa5937: Columns are stored by their result-pool column ids: the abundance, the subset it is narrowed to, and the subject, grouping and temporal columns. The timepoints of the timepoint order are stored as column value references of the temporal column, so the block that produced the metadata knows which of its values this block uses. The model finds columns and builds the options of every column picker through the column collection API instead of the deprecated result pool API, and the workflow takes every column through a column bundle. The settings panel reads the timepoints from a frame of the temporal column alone, and shows why when they cannot be read.
  Opening a project created before this change keeps its abundance and subset, and clears its subject, grouping and temporal columns and the timepoint order, since they were stored by ids anchored to the abundance column and as plain values, which name no producing block; they have to be picked again before the block runs. The initialization params of the kind carry the same ids and references, so templates made for an earlier version of the kind cannot be applied.
  Timepoints that are no longer among the values of the temporal column stay in the order and are marked as not available instead of being removed silently, with a tooltip naming the block that was to provide them, as a link to it, and what went wrong. Picking another temporal column still resets the order to the values of that column.

## 0.4.0

### Minor Changes

- 0e45d30: Update the SDK and accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) with the abundance. The abundance list is built with buildDatasetOptions (requireEnrichments: false), with the same labels as before. The subset restricts the clonotypes before the analysis, so frequencies, restriction indices and prevalence describe the subset only. Exported columns of a subset run carry pl7.app/inputSubset and a trace rooted in the subset. Unfiltered runs are unchanged.

## 0.3.3

### Patch Changes

- e5bf7da: Release every package through npm trusted publishing (OIDC) so each package's trusted publisher is used before it expires.
- Updated dependencies [e5bf7da]
  - @platforma-open/milaboratories.spatiotemporal-analysis.software@0.3.2

## 0.3.2

### Patch Changes

- Updated dependencies [57d43ca]
  - @platforma-open/milaboratories.spatiotemporal-analysis.software@0.3.1

## 0.3.1

### Patch Changes

- 21f4cc8: Stamp the output trace with the per-instance block label so two Clonotype Distribution blocks on the same dataset stay distinguishable downstream.

  Previously every instance emitted the same static trace label (`"Clonotype Distribution - <mode>"`), so downstream pickers (Lead Selection, etc.) showed duplicate options with no suffix — the label derivation had no differing trace entry to work from. The workflow now resolves the label as `customBlockLabel || defaultBlockLabel || "Clonotype Distribution"` and injects it into every exported column, so the derivation can append a distinguishing suffix. The model already carried the `customBlockLabel`/`defaultBlockLabel` fields and the subtitle binding; the workflow now consumes them, and the UI default label folds in the CLR normalization so blocks differing only in normalization stay distinct too.

## 0.3.0

### Minor Changes

- 55a62fb: Modernize the block on BlockModelV3, vectorize the Python pipeline with Polars, and publish PColumn specs before computation completes.

  - **Early spec export.** Every `xsv.importFile` call now passes `splitDataAndSpec: true`; a shared `buildTracedPf` helper wraps each result through `pframes.pFrameBuilder` with trace injection. Downstream blocks (e.g. Lead Selection) can resolve anchored queries and populate their UI the moment specs land, instead of waiting for the full workflow.
  - **Remove `wf.prepare`.** The outer `wf.body` now runs immediately and delegates to a new ephemeral template `process.tpl.tengo` that awaits only `PColumnBundle` state (specs, not data). Output specs publish to the result pool shortly after block start instead of after the full Python analysis — `wf.prepare`'s `Final`-state await previously gated the whole body on upstream data availability, cancelling the benefit of `splitDataAndSpec`.
  - **BlockModelV3 migration.** Replace `BlockModel` V1 with `BlockModelV3` + `DataModelBuilder`, merge `args` and `uiState` into a single `BlockData` type, and keep existing projects working via `upgradeLegacy`. UI switches to `defineAppV3` and reads `app.model.data.*`.
  - **Polars vectorization.** Rewrite `compartment_analysis.py` hot paths — grouping, temporal kinetics, CLR normalization — as Polars expressions, eliminating per-row Python loops.
  - **CID stability for reruns.** Resolve spec gaps 2, 3, 4, 5: the heatmap, temporal-line, and prevalence-histogram PFrames now carry `blockId` in their axis and column domains, matching the main PFrame and preventing content-ID collisions across block instances.
  - **Python test suite.** Add 101 behavioral tests covering grouping restriction, temporal kinetics, cross-subject convergence, CLR normalization, and empty/single-element edge cases — 77 % coverage on `compartment_analysis.py`. Shared fixtures, strengthened assertions.
  - **UI polish.** Reorder sidebar to surface subject prevalence first, hide graph tabs whose driving variable is unset, include the subject name in chart subtitles, and guard conditional output resolves.
  - **Discrete-filter derivation.** Prefer the upstream `pl7.app/discreteValues` annotation on the grouping column (new samples-and-data attaches it) — this fast path requires no data scan and preserves early spec export. For legacy upstreams that lack the annotation, fall back to scanning the grouping column via `col.data.getDataAsJson()`; this blocks the inner template on data availability and delays spec publication for those projects, but keeps the dropdown working. Peak Timepoint uses `args.timepointOrder` directly. Lead Selection renders dropdowns listing the actual categories (e.g. `spleen`, `lung`, `mediastinal lymph node`) instead of free-text fields.

### Patch Changes

- Updated dependencies [55a62fb]
  - @platforma-open/milaboratories.spatiotemporal-analysis.software@0.3.0

## 0.2.0

### Minor Changes

- 72b5956: Minor bug fix in python code

### Patch Changes

- Updated dependencies [72b5956]
  - @platforma-open/milaboratories.spatiotemporal-analysis.software@0.2.0
