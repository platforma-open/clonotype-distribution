import type { ColumnUniversalId, ColumnValueOption, RenderCtxBase } from "@platforma-sdk/model";
import {
  deriveColumnOptions,
  deriveColumnValueOptions,
  readColumnValues,
} from "@platforma-sdk/model";
import type { ColumnPickerOption } from "./columns";
import { axisSelector, exact, findColumn, resultPool } from "./columns";

/** Options of the timepoint order: the values of the temporal column. */
export type TimepointOptions = {
  /**
   * The temporal column the options were built for: they describe the picks only while it stays
   * picked.
   */
  computedFor: ColumnUniversalId;
  options: ColumnValueOption[];
};

/**
 * Sample-metadata columns of upstream blocks whose only axis is the abundance column's sample
 * axis (axis 0), as options of the subject, grouping and temporal pickers. Every option id is a
 * result-pool leaf, so it names the block that produced the column. `undefined` only while the
 * result pool is still resolving; an abundance column that is gone from a settled pool offers no
 * options.
 */
export function getMetadataOptions(
  ctx: RenderCtxBase,
  abundanceRef: ColumnUniversalId,
): ColumnPickerOption[] | undefined {
  const sampleAxis = findColumn(ctx, abundanceRef)?.getSpec().axesSpec[0];
  if (sampleAxis === undefined) return resultPool(ctx).isFinal() ? [] : undefined;

  const columns = resultPool(ctx).filter({
    include: {
      name: [exact(METADATA_COLUMN_NAME)],
      axes: [axisSelector(sampleAxis)],
      partialAxesMatch: false,
    },
  });
  return deriveColumnOptions(columns).map(({ id, label }) => ({ value: id, label }));
}

/**
 * Values of the temporal column as options of the timepoint order, each a column value id naming
 * the temporal column, and through it the block that produced the value. `undefined` until a
 * temporal column is picked. A column whose values cannot be read now — still resolving, not
 * computed, gone, or not stored inline — offers no values, so its picks show as not available.
 */
export function getTimepointOptions(
  ctx: RenderCtxBase,
  temporalColumnRef: ColumnUniversalId | undefined,
): TimepointOptions | undefined {
  if (temporalColumnRef === undefined) return undefined;
  const column = findColumn(ctx, temporalColumnRef);
  const entries = column === undefined ? undefined : readColumnValues(column);
  return { computedFor: temporalColumnRef, options: deriveColumnValueOptions(entries ?? []) };
}

// Internals

const METADATA_COLUMN_NAME = "pl7.app/metadata";
