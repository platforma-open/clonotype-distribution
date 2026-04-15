---
"@platforma-open/milaboratories.spatiotemporal-analysis.software": minor
"@platforma-open/milaboratories.spatiotemporal-analysis.workflow": minor
"@platforma-open/milaboratories.spatiotemporal-analysis.model": minor
"@platforma-open/milaboratories.spatiotemporal-analysis.test": minor
"@platforma-open/milaboratories.spatiotemporal-analysis.ui": minor
"@platforma-open/milaboratories.spatiotemporal-analysis": minor
---

Modernize the block on BlockModelV3, vectorize the Python pipeline with Polars, and publish PColumn specs before computation completes.

- **Early spec export.** Every `xsv.importFile` call now passes `splitDataAndSpec: true`; a shared `buildTracedPf` helper wraps each result through `pframes.pFrameBuilder` with trace injection. Downstream blocks (e.g. Lead Selection) can resolve anchored queries and populate their UI the moment specs land, instead of waiting for the full workflow.
- **BlockModelV3 migration.** Replace `BlockModel` V1 with `BlockModelV3` + `DataModelBuilder`, merge `args` and `uiState` into a single `BlockData` type, and keep existing projects working via `upgradeLegacy`. UI switches to `defineAppV3` and reads `app.model.data.*`.
- **Polars vectorization.** Rewrite `compartment_analysis.py` hot paths — grouping, temporal kinetics, CLR normalization — as Polars expressions, eliminating per-row Python loops.
- **CID stability for reruns.** Resolve spec gaps 2, 3, 4, 5: the heatmap, temporal-line, and prevalence-histogram PFrames now carry `blockId` in their axis and column domains, matching the main PFrame and preventing content-ID collisions across block instances.
- **Python test suite.** Add 101 behavioral tests covering grouping restriction, temporal kinetics, cross-subject convergence, CLR normalization, and empty/single-element edge cases — 77 % coverage on `compartment_analysis.py`. Shared fixtures, strengthened assertions.
- **UI polish.** Reorder sidebar to surface subject prevalence first, hide graph tabs whose driving variable is unset, include the subject name in chart subtitles, and guard conditional output resolves.
- **Discrete-filter propagation.** Read `pl7.app/discreteValues` from the upstream Grouping and Temporal Variable metadata columns and forward it to the `Dominant Group`, `Consensus Dominant`, and `Peak Timepoint` output specs. `Peak Timepoint` falls back to `args.timepointOrder` when the temporal column lacks the annotation. Lead Selection now renders these columns as dropdowns listing the actual categories (e.g. `spleen`, `lung`, `mediastinal lymph node`) instead of free-text fields.
