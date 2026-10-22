import {
  DataModelBuilder,
  createPlDataTableStateV2,
  isColumnUniversalId,
  isColumnValueId,
} from "@platforma-sdk/model";
import { kind } from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";
import { toColumnId } from "./abundance";
import type { BlockArgsV1, BlockDataV1 } from "./stored_data";
import type { BlockData, LegacyUiState } from "./types";

export const blockDataModel = new DataModelBuilder({ kind })
  .from<BlockDataV1>("v1")
  .upgradeLegacy<BlockArgsV1, LegacyUiState>(({ args, uiState }) => ({
    ...args,
    ...uiState,
  }))
  // Columns become result-pool ids and the timepoints column value ids
  .migrate<BlockData>("v2", toColumnIds)
  .init(({ params }) => ({
    // Left empty on purpose. The default label is a join of the chosen columns'
    // human option labels, and those come from the result pool, which `init`
    // cannot reach -- so a block created from a template traces the static
    // "Clonotype Distribution" fallback until the settings panel derives the
    // real label. Anything computable here (the mode, the CLR marker) would be
    // a prefix of that label, not the label, so it would still be replaced.
    defaultBlockLabel: "",
    customBlockLabel: params?.customBlockLabel ?? "",
    abundanceRef: params?.abundanceRef,
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
function toColumnIds(data: BlockDataV1): BlockData {
  const column = (stored: string | undefined) => (isColumnUniversalId(stored) ? stored : undefined);
  return {
    ...data,
    abundanceRef: data.abundanceRef === undefined ? undefined : toColumnId(data.abundanceRef),
    subjectColumnRef: column(data.subjectColumnRef),
    groupingColumnRef: column(data.groupingColumnRef),
    temporalColumnRef: column(data.temporalColumnRef),
    timepointOrder: data.timepointOrder.filter(isColumnValueId),
  };
}
