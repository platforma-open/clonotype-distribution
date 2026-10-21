import type {
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  RenderCtxBase,
  TreeNodeAccessor,
} from "@platforma-sdk/model";
import { BlockModelV3, createPFrameForGraphs, createPlDataTableV3 } from "@platforma-sdk/model";
import { kind } from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";
import { getAbundanceOptions } from "./abundance";
import { findColumn, getOutputColumns, toIdsAndSpecs, toPColumns } from "./columns";
import { blockDataModel } from "./dataModel";
import { getMetadataOptions, getTimepointOptions } from "./metadata";
import type { BlockArgs } from "./types";

export { blockDataModel } from "./dataModel";
export type { ColumnPickerOption } from "./columns";
export type { TimepointOptions } from "./metadata";
export type { BlockArgs, BlockData, LegacyUiState } from "./types";
export type * from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";

export const platforma = BlockModelV3.create({ dataModel: blockDataModel, kind })

  // The run gate. Throwing surfaces the reason in the UI, which a disabled Run
  // button would not.
  .args<BlockArgs>((data): BlockArgs => {
    // Strip UI state; everything else maps 1:1 to BlockArgs.
    const {
      tableState: _tableState,
      heatmapState: _heatmapState,
      temporalLineState: _temporalLineState,
      prevalenceHistogramState: _prevalenceHistogramState,
      ...args
    } = data;
    const {
      abundanceRef,
      groupingColumnRef,
      temporalColumnRef,
      timepointOrder,
      subjectColumnRef,
      calculationMode,
    } = args;
    if (abundanceRef === undefined) throw new Error("Abundance ref required");
    const hasGrouping = groupingColumnRef !== undefined;
    const hasTemporal = temporalColumnRef !== undefined && timepointOrder.length >= 2;
    if (!hasGrouping && !hasTemporal)
      throw new Error("At least grouping or temporal variable required");
    if (calculationMode === "intra-subject" && subjectColumnRef === undefined)
      throw new Error("Subject required in intra-subject mode");
    return args;
  })

  // Inverse of the kind's init-params contract: every field a user sets by
  // hand. `defaultBlockLabel` is derived by a watchEffect in
  // ui/src/pages/MainPage.vue from the chosen columns' option labels, so it is
  // projected into args (the workflow reads it for the trace) but never
  // templated.
  .templateParams((data) => ({
    abundanceRef: data.abundanceRef,
    calculationMode: data.calculationMode,
    subjectColumnRef: data.subjectColumnRef,
    groupingColumnRef: data.groupingColumnRef,
    temporalColumnRef: data.temporalColumnRef,
    timepointOrder: data.timepointOrder,
    normalization: data.normalization,
    presenceThreshold: data.presenceThreshold,
    minAbundanceThreshold: data.minAbundanceThreshold,
    minSubjectCount: data.minSubjectCount,
    topN: data.topN,
    customBlockLabel: data.customBlockLabel,
  }))

  // Abundance column options
  .output("abundanceOptions", (ctx) => getAbundanceOptions(ctx))

  // Metadata column options
  .output("metadataOptions", (ctx) => {
    const abundanceRef = ctx.data.abundanceRef;
    if (abundanceRef === undefined) return undefined;
    return getMetadataOptions(ctx, abundanceRef);
  })

  // Dataset spec for detecting cluster vs clonotype
  .output("datasetSpec", (ctx) => {
    const abundanceRef = ctx.data.abundanceRef;
    if (abundanceRef === undefined) return undefined;
    return findColumn(ctx, abundanceRef)?.getSpec();
  })

  // Options of the timepoint order: the values of the temporal column
  .output("timepointOptions", (ctx) => getTimepointOptions(ctx, ctx.data.temporalColumnRef))

  // Main output table
  .outputWithStatus("mainTable", (ctx) => {
    const columns = getOutputColumns(ctx, ctx.outputs?.resolve("mainPf"));
    if (columns === undefined) return undefined;
    return createPlDataTableV3(ctx, {
      // The block's own result columns, all stored as they are, full-joined on the element.
      primaryColumns: columns.getColumns(),
      columns: [],
      tableState: ctx.data.tableState,
    });
  })

  // Heatmap PFrame + raw columns for graph defaults (requires grouping variable)
  .outputWithStatus("heatmapPf", (ctx): PFrameHandle | undefined => {
    if (ctx.data.groupingColumnRef === undefined) return undefined;
    return createGraphFrame(ctx, resolveOptional(ctx, "heatmapPf"));
  })
  .output("heatmapPCols", (ctx) => {
    if (ctx.data.groupingColumnRef === undefined) return undefined;
    return getGraphColumns(ctx, resolveOptional(ctx, "heatmapPf"));
  })

  // Temporal line PFrame + raw columns for graph defaults (requires temporal variable)
  .outputWithStatus("temporalLinePf", (ctx): PFrameHandle | undefined => {
    if (ctx.data.temporalColumnRef === undefined) return undefined;
    return createGraphFrame(ctx, resolveOptional(ctx, "temporalLinePf"));
  })
  .output("temporalLinePCols", (ctx) => {
    if (ctx.data.temporalColumnRef === undefined) return undefined;
    return getGraphColumns(ctx, resolveOptional(ctx, "temporalLinePf"));
  })

  // Prevalence histogram PFrame + raw columns for graph defaults (requires subject variable)
  .outputWithStatus("prevalenceHistogramPf", (ctx): PFrameHandle | undefined => {
    if (ctx.data.subjectColumnRef === undefined) return undefined;
    return createGraphFrame(ctx, resolveOptional(ctx, "prevalenceHistogramPf"));
  })
  .output("prevalenceHistogramPCols", (ctx) => {
    if (ctx.data.subjectColumnRef === undefined) return undefined;
    return getGraphColumns(ctx, resolveOptional(ctx, "prevalenceHistogramPf"));
  })

  // R3: Per-subject detail PFrame (intra-subject mode with subject variable)
  .outputWithStatus("perSubjectPf", (ctx): PFrameHandle | undefined => {
    if (ctx.data.calculationMode !== "intra-subject" || ctx.data.subjectColumnRef === undefined)
      return undefined;
    return createGraphFrame(ctx, resolveOptional(ctx, "perSubjectPf"));
  })
  .output("perSubjectPCols", (ctx) => {
    if (ctx.data.calculationMode !== "intra-subject" || ctx.data.subjectColumnRef === undefined)
      return undefined;
    return getGraphColumns(ctx, resolveOptional(ctx, "perSubjectPf"));
  })

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .title(() => "Clonotype Distribution")

  .subtitle((ctx) => ctx.data.customBlockLabel || ctx.data.defaultBlockLabel)

  .sections((ctx) => {
    const sections: { type: "link"; href: `/${string}`; label: string }[] = [
      { type: "link", href: "/", label: "Main" },
    ];
    if (ctx.data.subjectColumnRef !== undefined) {
      sections.push({ type: "link", href: "/prevalence", label: "Subject Prevalence" });
    }
    if (ctx.data.groupingColumnRef !== undefined) {
      sections.push({ type: "link", href: "/heatmap", label: "Distribution Heatmap" });
    }
    if (ctx.data.temporalColumnRef !== undefined) {
      sections.push({ type: "link", href: "/temporal", label: "Temporal Trajectory" });
    }
    return sections;
  })

  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;

// Internals

/**
 * One of this block's graph outputs. The workflow produces each only for the settings it ran
 * with, so the field may be absent for good.
 */
function resolveOptional(ctx: RenderCtxBase, field: string): TreeNodeAccessor | undefined {
  return ctx.outputs?.resolve({ field, allowPermanentAbsence: true });
}

/** A graph frame of the columns of one of this block's outputs, once they have all arrived. */
function createGraphFrame(
  ctx: RenderCtxBase,
  accessor: TreeNodeAccessor | undefined,
): PFrameHandle | undefined {
  const columns = getOutputColumns(ctx, accessor)?.getColumns();
  if (columns === undefined) return undefined;
  return createPFrameForGraphs(ctx, toPColumns(columns));
}

/** Ids and specs of the columns of one of this block's outputs, for graph defaults. */
function getGraphColumns(
  ctx: RenderCtxBase,
  accessor: TreeNodeAccessor | undefined,
): PColumnIdAndSpec[] | undefined {
  const columns = getOutputColumns(ctx, accessor)?.getColumns();
  if (columns === undefined) return undefined;
  return toIdsAndSpecs(columns);
}
