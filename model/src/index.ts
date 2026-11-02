import type {
  DatasetOption,
  InferOutputsType,
  PFrameHandle,
  PObjectSpec,
} from "@platforma-sdk/model";
import {
  BlockModelV3,
  buildDatasetOptions,
  createGlobalPObjectId,
  createPFrameForGraphs,
  createPlDataTableV2,
  isPColumnSpec,
  plRefsEqual,
} from "@platforma-sdk/model";
import { kind } from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";
import { blockDataModel } from "./dataModel";
import type { BlockArgs } from "./types";

export { blockDataModel } from "./dataModel";
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
      filterRef,
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
    return {
      ...args,
      ...(filterRef !== undefined && {
        inputFilter: createGlobalPObjectId(filterRef.blockId, filterRef.name),
      }),
    };
  })

  // Inverse of the kind's init-params contract: every field a user sets by
  // hand. `defaultBlockLabel` is derived by a watchEffect in
  // ui/src/pages/MainPage.vue from the chosen columns' option labels, so it is
  // projected into args (the workflow reads it for the trace) but never
  // templated.
  .templateParams((data) => ({
    abundanceRef: data.abundanceRef,
    filterRef: data.filterRef,
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

  // Primary raw abundances, each with its subset columns (e.g. Repertoire Labeling labels) as filters.
  .output("abundanceOptions", (ctx): DatasetOption[] => {
    const isPrimaryAbundance = (spec: PObjectSpec) =>
      isPColumnSpec(spec) &&
      spec.axesSpec.length === 2 &&
      spec.axesSpec[0].name === "pl7.app/sampleId" &&
      spec.annotations?.["pl7.app/isAbundance"] === "true" &&
      spec.annotations?.["pl7.app/abundance/normalized"] === "false" &&
      spec.annotations?.["pl7.app/abundance/isPrimary"] === "true";

    // Labels as before, with the column's own label (e.g. "Number of Reads"): buildDatasetOptions
    // derives them without it and takes no label options for the datasets.
    const labeled = ctx.resultPool.getOptions(isPrimaryAbundance, {
      label: { includeNativeLabel: true },
    });

    const options =
      buildDatasetOptions(ctx, {
        primary: isPrimaryAbundance,
        // Only subsets keyed by the element axis alone: a sample subset would drop whole samples.
        filter: (spec) =>
          isPColumnSpec(spec) &&
          spec.axesSpec.length === 1 &&
          spec.axesSpec[0].name !== "pl7.app/sampleId",
        requireEnrichments: false,
      }) ?? [];

    return options.map((option) => {
      const label = labeled.find((o) => plRefsEqual(o.ref, option.primary.ref, true))?.label;
      return label === undefined ? option : { ...option, primary: { ...option.primary, label } };
    });
  })

  // Metadata column options
  .output("metadataOptions", (ctx) => {
    const anchor = ctx.data.abundanceRef;
    if (anchor === undefined) return undefined;
    return ctx.resultPool.getCanonicalOptions({ main: anchor }, [
      {
        axes: [{ anchor: "main", idx: 0 }],
        name: "pl7.app/metadata",
      },
    ]);
  })

  // Dataset spec for detecting cluster vs clonotype
  .output("datasetSpec", (ctx) => {
    if (ctx.data.abundanceRef) return ctx.resultPool.getPColumnSpecByRef(ctx.data.abundanceRef);
    return undefined;
  })

  // PFrame containing temporal column data (for fetching unique values in UI)
  .output("temporalColumnPframe", (ctx) => {
    const { temporalColumnRef, abundanceRef } = ctx.data;
    if (!temporalColumnRef || !abundanceRef) return undefined;

    const cols = ctx.resultPool.getAnchoredPColumns(
      { main: abundanceRef },
      JSON.parse(temporalColumnRef) as never,
    );
    if (!cols || cols.length === 0) return undefined;
    return ctx.createPFrame(cols);
  })

  // Column ID for the temporal column (needed by getSingleColumnData in UI)
  .output("temporalColumnId", (ctx) => {
    const { temporalColumnRef, abundanceRef } = ctx.data;
    if (!temporalColumnRef || !abundanceRef) return undefined;

    const cols = ctx.resultPool.getAnchoredPColumns(
      { main: abundanceRef },
      JSON.parse(temporalColumnRef) as never,
    );
    return cols?.[0]?.id;
  })

  // Main output table
  .outputWithStatus("mainTable", (ctx) => {
    const pCols = ctx.outputs?.resolve("mainPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPlDataTableV2(ctx, pCols, ctx.data.tableState);
  })

  // Heatmap PFrame + raw columns for graph defaults (requires grouping variable)
  .outputWithStatus("heatmapPf", (ctx): PFrameHandle | undefined => {
    if (ctx.data.groupingColumnRef === undefined) return undefined;
    try {
      const pCols = ctx.outputs?.resolve("heatmapPf")?.getPColumns();
      if (pCols === undefined) return undefined;
      return createPFrameForGraphs(ctx, pCols);
    } catch {
      return undefined;
    }
  })
  .output("heatmapPCols", (ctx) => {
    if (ctx.data.groupingColumnRef === undefined) return undefined;
    try {
      return ctx.outputs?.resolve("heatmapPf")?.getPColumns();
    } catch {
      return undefined;
    }
  })

  // Temporal line PFrame + raw columns for graph defaults (requires temporal variable)
  .outputWithStatus("temporalLinePf", (ctx): PFrameHandle | undefined => {
    if (ctx.data.temporalColumnRef === undefined) return undefined;
    try {
      const pCols = ctx.outputs?.resolve("temporalLinePf")?.getPColumns();
      if (pCols === undefined) return undefined;
      return createPFrameForGraphs(ctx, pCols);
    } catch {
      return undefined;
    }
  })
  .output("temporalLinePCols", (ctx) => {
    if (ctx.data.temporalColumnRef === undefined) return undefined;
    try {
      return ctx.outputs?.resolve("temporalLinePf")?.getPColumns();
    } catch {
      return undefined;
    }
  })

  // Prevalence histogram PFrame + raw columns for graph defaults (requires subject variable)
  .outputWithStatus("prevalenceHistogramPf", (ctx): PFrameHandle | undefined => {
    if (ctx.data.subjectColumnRef === undefined) return undefined;
    try {
      const pCols = ctx.outputs?.resolve("prevalenceHistogramPf")?.getPColumns();
      if (pCols === undefined) return undefined;
      return createPFrameForGraphs(ctx, pCols);
    } catch {
      return undefined;
    }
  })
  .output("prevalenceHistogramPCols", (ctx) => {
    if (ctx.data.subjectColumnRef === undefined) return undefined;
    try {
      return ctx.outputs?.resolve("prevalenceHistogramPf")?.getPColumns();
    } catch {
      return undefined;
    }
  })

  // R3: Per-subject detail PFrame (intra-subject mode with subject variable)
  .outputWithStatus("perSubjectPf", (ctx): PFrameHandle | undefined => {
    if (ctx.data.calculationMode !== "intra-subject" || ctx.data.subjectColumnRef === undefined)
      return undefined;
    try {
      const pCols = ctx.outputs?.resolve("perSubjectPf")?.getPColumns();
      if (pCols === undefined) return undefined;
      return createPFrameForGraphs(ctx, pCols);
    } catch {
      return undefined;
    }
  })
  .output("perSubjectPCols", (ctx) => {
    if (ctx.data.calculationMode !== "intra-subject" || ctx.data.subjectColumnRef === undefined)
      return undefined;
    try {
      return ctx.outputs?.resolve("perSubjectPf")?.getPColumns();
    } catch {
      return undefined;
    }
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
