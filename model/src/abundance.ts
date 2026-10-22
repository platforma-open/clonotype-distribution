import type { ColumnUniversalId, PlRef, RelaxedColumnSelector } from "@platforma-sdk/model";
import { createGlobalPObjectId, deriveColumnOptions } from "@platforma-sdk/model";
import type { ColumnPickerOption } from "./columns";
import { exact, resultPool } from "./columns";

/**
 * Primary raw abundances per sample of upstream blocks, as abundance picker options labelled
 * with the column's own label (e.g. "Number of Reads"). Every option id is a result-pool leaf,
 * so it names the block that produced the column.
 */
export function getAbundanceOptions(): ColumnPickerOption[] {
  const columns = resultPool().filter({ include: PRIMARY_ABUNDANCE });
  return deriveColumnOptions(columns, { includeNativeLabel: true }).map(({ id, label }) => ({
    value: id,
    label,
  }));
}

/** The result-pool id of the upstream output a ref points to: it names the same block and column. */
export function toColumnId(ref: PlRef): ColumnUniversalId {
  return createGlobalPObjectId(ref.blockId, ref.name);
}

// Internals

/** Primary raw abundances per sample: the first of exactly two axes is the sample. */
const PRIMARY_ABUNDANCE: RelaxedColumnSelector = {
  axes: [{ name: [exact("pl7.app/sampleId")] }, {}],
  partialAxesMatch: false,
  annotations: {
    "pl7.app/isAbundance": [exact("true")],
    "pl7.app/abundance/normalized": [exact("false")],
    "pl7.app/abundance/isPrimary": [exact("true")],
  },
};
