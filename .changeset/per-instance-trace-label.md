---
'@platforma-open/milaboratories.spatiotemporal-analysis.workflow': patch
'@platforma-open/milaboratories.spatiotemporal-analysis.ui': patch
---

Stamp the output trace with the per-instance block label so two Clonotype Distribution blocks on the same dataset stay distinguishable downstream.

Previously every instance emitted the same static trace label (`"Clonotype Distribution - <mode>"`), so downstream pickers (Lead Selection, etc.) showed duplicate options with no suffix — the label derivation had no differing trace entry to work from. The workflow now resolves the label as `customBlockLabel || defaultBlockLabel || "Clonotype Distribution"` and injects it into every exported column, so the derivation can append a distinguishing suffix. The model already carried the `customBlockLabel`/`defaultBlockLabel` fields and the subtitle binding; the workflow now consumes them, and the UI default label folds in the CLR normalization so blocks differing only in normalization stay distinct too.
