---
"@platforma-open/milaboratories.spatiotemporal-analysis.model": patch
"@platforma-open/milaboratories.spatiotemporal-analysis.ui": patch
"@platforma-open/milaboratories.spatiotemporal-analysis.workflow": patch
---

Compute the block label without the UI

The trace label is now built by the workflow from the native labels of the chosen grouping,
temporal and subject columns plus the calculation mode and normalization, so opening the block
no longer changes its arguments and a block that was never opened still gets a meaningful label.
The subtitle shows the result-pool labels of the chosen columns, refreshed when the user picks a
column. Existing blocks show as stale once after the update.
