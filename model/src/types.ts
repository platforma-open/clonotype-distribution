import type { GraphMakerState } from "@milaboratories/graph-maker";
import type { PlDataTableStateV2, PlRef, SUniversalPColumnId } from "@platforma-sdk/model";
// The calculation-mode and normalization vocabularies live in the kind: its
// init-params contract names them and a kind cannot import from the model.
import type {
  CalculationMode,
  Normalization,
} from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";

/**
 * Workflow-facing args — every user-committed analysis decision. Both label
 * fields are here because the workflow reads them for the provenance trace.
 */
export type BlockArgs = {
  defaultBlockLabel: string;
  customBlockLabel: string;
  abundanceRef?: PlRef;
  calculationMode: CalculationMode;
  groupingColumnRef?: SUniversalPColumnId;
  temporalColumnRef?: SUniversalPColumnId;
  timepointOrder: string[];
  subjectColumnRef?: SUniversalPColumnId;
  normalization: Normalization;
  presenceThreshold: number;
  minAbundanceThreshold: number;
  minSubjectCount: number;
  topN: number;
};

/** Unified V3 data — the args above plus the four view states the workflow never sees. */
export type BlockData = BlockArgs & {
  tableState: PlDataTableStateV2;
  heatmapState: GraphMakerState;
  temporalLineState: GraphMakerState;
  prevalenceHistogramState: GraphMakerState;
};

/** Legacy on-disk UI channel, read only by `.upgradeLegacy`. */
export type LegacyUiState = {
  tableState: PlDataTableStateV2;
  heatmapState: GraphMakerState;
  temporalLineState: GraphMakerState;
  prevalenceHistogramState: GraphMakerState;
};
