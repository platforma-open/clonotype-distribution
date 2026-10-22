---
"@platforma-open/milaboratories.spatiotemporal-analysis": major
"@platforma-open/milaboratories.spatiotemporal-analysis.kind": major
"@platforma-open/milaboratories.spatiotemporal-analysis.model": major
"@platforma-open/milaboratories.spatiotemporal-analysis.ui": major
"@platforma-open/milaboratories.spatiotemporal-analysis.workflow": major
---

Columns are stored by their result-pool column ids: the abundance and the subject, grouping and temporal columns. The timepoints of the timepoint order are stored as column value references of the temporal column, so the block that produced the metadata knows which of its values this block uses. The model finds columns and builds the options of every column picker through the column collection API instead of the deprecated result pool API, and the workflow takes every column through a column bundle. The settings panel reads the timepoints from a frame of the temporal column alone, and shows why when they cannot be read.
Opening a project created before this change keeps its abundance, and clears its subject, grouping and temporal columns and the timepoint order, since they were stored by ids anchored to the abundance column and as plain values, which name no producing block; they have to be picked again before the block runs. The initialization params of the kind carry the same ids and references, so templates made for an earlier version of the kind cannot be applied.
Timepoints that are no longer among the values of the temporal column stay in the order and are marked as not available instead of being removed silently, with a tooltip naming the block that was to provide them, as a link to it, and what went wrong. Picking another temporal column still resets the order to the values of that column.
