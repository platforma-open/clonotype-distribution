---
"@platforma-open/milaboratories.spatiotemporal-analysis.kind": patch
"@platforma-open/milaboratories.spatiotemporal-analysis": patch
---

Publish the kind with its corrected doc comment as a new version, so the block release goes through

The registry keeps a published kind version immutable and identifies it by a hash of the kind's
`src/`. A comment in `kind/src/index.ts` was corrected after kind 2.0.0 had been published (it
described `defaultBlockLabel`, which is no longer a field), so publishing block 2.0.1 failed with
"Immutable kind version republished with different content" and 2.0.1 was never released. The kind
goes to 2.0.1; its params contract is the same as 2.0.0. The block release also carries the fix from
2.0.1: opening the block no longer changes its arguments.
