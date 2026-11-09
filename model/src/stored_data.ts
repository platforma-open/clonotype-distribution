import type { PlRef } from "@platforma-sdk/model";
import type { BlockData, LegacyUiState } from "./types";

/**
 * Args as stored at `v1` and before it: the abundance column as a ref to an upstream output,
 * the subject, grouping and temporal columns as ids anchored against the abundance column, and
 * the timepoints as plain values of the temporal column.
 */
export type BlockArgsV1 = {
  defaultBlockLabel: string;
  customBlockLabel: string;
  abundanceRef?: PlRef;
  calculationMode: "population" | "intra-subject";
  groupingColumnRef?: string;
  temporalColumnRef?: string;
  timepointOrder: string[];
  subjectColumnRef?: string;
  normalization: "relative-frequency" | "clr";
  presenceThreshold: number;
  minAbundanceThreshold: number;
  minSubjectCount: number;
  topN: number;
};

/**
 * Data as stored at `v1`: the args above, the subset the abundance is narrowed to as a ref to an
 * upstream output, and the four view states.
 */
export type BlockDataV1 = BlockArgsV1 & LegacyUiState & { filterRef?: PlRef };

/**
 * Data as stored at `v2`: the block label the UI derived from the chosen columns' result-pool
 * labels and the settings, in place of the subtitle text.
 */
export type BlockDataV2 = Omit<BlockData, "subtitleLabel"> & { defaultBlockLabel: string };
