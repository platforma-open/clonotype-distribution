import { DataModelBuilder, createPlDataTableStateV2 } from "@platforma-sdk/model";
import { kind } from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";
import type { BlockArgs, BlockData, LegacyUiState } from "./types";

export const blockDataModel = new DataModelBuilder({ kind })
  .from<BlockData>("v1")
  .upgradeLegacy<BlockArgs, LegacyUiState>(({ args, uiState }) => ({
    ...args,
    ...uiState,
  }))
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
