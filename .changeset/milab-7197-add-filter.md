---
"@platforma-open/milaboratories.spatiotemporal-analysis": minor
"@platforma-open/milaboratories.spatiotemporal-analysis.kind": minor
"@platforma-open/milaboratories.spatiotemporal-analysis.model": minor
"@platforma-open/milaboratories.spatiotemporal-analysis.ui": minor
"@platforma-open/milaboratories.spatiotemporal-analysis.workflow": minor
---

Update the SDK and accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) with the abundance. The abundance list is built with buildDatasetOptions (requireEnrichments: false), with the same labels as before. The subset restricts the clonotypes before the analysis, so frequencies, restriction indices and prevalence describe the subset only. Exported columns of a subset run carry pl7.app/inputSubset and a trace rooted in the subset. Unfiltered runs are unchanged.
