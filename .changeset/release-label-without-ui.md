---
"@platforma-open/milaboratories.spatiotemporal-analysis": patch
---

Release the block label computed without the UI, so opening the block no longer changes its arguments

In 2.0.0 the settings page wrote `defaultBlockLabel` -- the calculation mode, the result-pool labels
of the chosen columns and CLR -- into the block data, and that field was part of the args. A block
created from a template starts with an empty label, so the first time its window was opened after
the run, the args changed and the block and every block below it went stale. The fix (the settings
label computed in the model, the column labels added by the workflow, the subtitle text kept out of
the args) was merged after 2.0.0 but changed only the model, UI and workflow packages, so no block
version carrying it was published. This releases it. A block calculated with 2.0.0 shows as stale
once after the update, because its args no longer carry `defaultBlockLabel`.
