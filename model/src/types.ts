import type { GraphMakerState } from "@milaboratories/graph-maker";
import type { ColumnUniversalId, ColumnValueId, PlDataTableStateV2 } from "@platforma-sdk/model";
// The calculation-mode and normalization vocabularies live in the kind: its
// init-params contract names them and a kind cannot import from the model.
import type {
  CalculationMode,
  Normalization,
} from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";

/**
 * Workflow-facing args — every user-committed analysis decision. Both label
 * fields are here because the workflow builds the provenance trace label from
 * them and the native labels of the chosen columns.
 * The abundance and metadata columns are result-pool column ids; the timepoints
 * are column value ids of the temporal column.
 */
export type BlockArgs = {
  customBlockLabel: string;
  /** Settings part of the trace label: the calculation mode, and CLR when chosen. */
  paramsLabel: string;
  abundanceRef?: ColumnUniversalId;
  /** Canonical id of the subset column; set only when a filter is picked. */
  inputFilter?: string;
  calculationMode: CalculationMode;
  groupingColumnRef?: ColumnUniversalId;
  temporalColumnRef?: ColumnUniversalId;
  timepointOrder: ColumnValueId[];
  subjectColumnRef?: ColumnUniversalId;
  normalization: Normalization;
  presenceThreshold: number;
  minAbundanceThreshold: number;
  minSubjectCount: number;
  topN: number;
};

/** Unified V3 data — the args above plus the four view states the workflow never sees. */
export type BlockData = Omit<BlockArgs, "inputFilter" | "paramsLabel"> & {
  /**
   * Labels of the chosen subject, grouping and temporal columns as the result pool shows them.
   * Kept in sync by the UI; feeds the subtitle only, never the workflow.
   */
  subtitleLabel?: string;
  /** Subset column picked with the abundance in PlDatasetSelector (e.g. a Repertoire Labeling label). */
  filterRef?: ColumnUniversalId;
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
