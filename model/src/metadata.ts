import type { ColumnUniversalId, RenderCtxBase } from "@platforma-sdk/model";
import { deriveColumnOptions } from "@platforma-sdk/model";
import type { ColumnPickerOption } from "./columns";
import { axisSelector, exact, findColumn, resultPool } from "./columns";

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

// Internals

const METADATA_COLUMN_NAME = "pl7.app/metadata";
