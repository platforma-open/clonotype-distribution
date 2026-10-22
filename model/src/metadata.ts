import type {
  ColumnUniversalId,
  PFrameHandle,
  PObjectId,
  RenderCtxBase,
} from "@platforma-sdk/model";
import { ColumnRecipe, deriveColumnOptions, extractPObjectId } from "@platforma-sdk/model";
import type { ColumnPickerOption } from "./columns";
import { axisSelector, exact, findColumn, resultPool } from "./columns";

/**
 * Where the UI reads the values of a picked metadata column from: a frame of the column alone,
 * while the column is there. `column` is the pick the frame was made for, so the UI can tell one
 * that still describes an earlier pick; `frameColumn` is the column as the frame holds it.
 */
export type MetadataColumnFrame =
  | { readonly type: "resolving"; readonly column: ColumnUniversalId }
  /** The column is gone: a stored id may name a column its upstream block no longer produces. */
  | { readonly type: "absent"; readonly column: ColumnUniversalId }
  | {
      readonly type: "present";
      readonly column: ColumnUniversalId;
      readonly frame: PFrameHandle;
      readonly frameColumn: PObjectId;
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

/** A frame of the metadata `column` alone, while the column and its frame are there. */
export function getMetadataColumnFrame(
  ctx: RenderCtxBase,
  column: ColumnUniversalId,
): MetadataColumnFrame {
  const status = ColumnRecipe.getStatus(column, { ctx: ctx.ctx });
  if (status !== "present") return { type: status, column };
  const frame = ctx.createPFrame([column]);
  if (frame === undefined) return { type: "resolving", column };
  // The frame holds a column by the bare id of its leaf.
  return { type: "present", column, frame, frameColumn: extractPObjectId(column) };
}

// Internals

const METADATA_COLUMN_NAME = "pl7.app/metadata";
