import {
  DataModelBuilder,
  createPlDataTableStateV2,
  isGlobalPObjectId,
  isColumnValueId,
} from "@platforma-sdk/model";
import { kind } from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";
import { toColumnId } from "./abundance";
import type { BlockArgsV1, BlockDataV1, BlockDataV2 } from "./stored_data";
import type { BlockData, LegacyUiState } from "./types";

export const blockDataModel = new DataModelBuilder({ kind })
  .from<BlockDataV1>("v1")
  .upgradeLegacy<BlockArgsV1, LegacyUiState>(({ args, uiState }) => ({
    ...args,
    ...uiState,
  }))
  // Columns become result-pool ids and the timepoints column value ids
  .migrate<BlockDataV2>("v2", toColumnIds)
  // The old label mixed the chosen columns' result-pool labels with settings; the settings part
  // is recomputed from data, and the column part returns once the block is opened.
  .migrate<BlockData>("v3", ({ defaultBlockLabel: _, ...rest }) => ({
    ...rest,
    subtitleLabel: "",
  }))
  .init(({ params }) => ({
    customBlockLabel: params?.customBlockLabel ?? "",
    // The columns' result-pool labels are out of reach here; the UI fills them in once it opens.
    subtitleLabel: "",
    abundanceRef: params?.abundanceRef,
    filterRef: params?.filterRef,
    calculationMode: params?.calculationMode ?? ("population" as const),
    groupingColumnRef: params?.groupingColumnRef,
    temporalColumnRef: params?.temporalColumnRef,
    timepointOrder: params?.timepointOrder ?? [],
    subjectColumnRef: params?.subjectColumnRef,
    normalization: params?.normalization ?? ("relative-frequency" as const),
    presenceThreshold: params?.presenceThreshold ?? 0,
    minAbundanceThreshold: params?.minAbundanceThreshold ?? 0,
    minSubjectCount: params?.minSubjectCount ?? 1,
    topN: params?.topN ?? 20,

    tableState: createPlDataTableStateV2(),
    heatmapState: {
      title: "Distribution heatmap",
      template: "heatmap",
      currentTab: null,
    },
    temporalLineState: {
      title: "Temporal frequency trajectory",
      template: "curve_dots",
      currentTab: null,
      layersSettings: {
        curve: {
          smoothing: false,
        },
      },
    },
    prevalenceHistogramState: {
      title: "Subject prevalence distribution",
      template: "bar",
      currentTab: null,
      layersSettings: {
        bar: { fillColor: "#5b9bd5" },
      },
    },
  }));

// Internals

/**
 * Moves stored data to column ids. A ref to an upstream output becomes the result-pool id of
 * that output, which names the same block and column. The subject, grouping and temporal
 * columns were picked by ids anchored against the abundance column and the timepoints stored as
 * plain values; neither names the block that produced the metadata, so they cannot become
 * column ids and column value ids and are dropped, to be picked again.
 */
function toColumnIds(data: BlockDataV1): BlockDataV2 {
  const column = (stored: string | undefined) => (isGlobalPObjectId(stored) ? stored : undefined);
  return {
    ...data,
    abundanceRef: data.abundanceRef === undefined ? undefined : toColumnId(data.abundanceRef),
    filterRef: data.filterRef === undefined ? undefined : toColumnId(data.filterRef),
    subjectColumnRef: column(data.subjectColumnRef),
    groupingColumnRef: column(data.groupingColumnRef),
    temporalColumnRef: column(data.temporalColumnRef),
    timepointOrder: data.timepointOrder.filter(isColumnValueId),
  };
}
