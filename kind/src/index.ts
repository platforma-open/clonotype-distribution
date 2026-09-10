import { assertParamsObject, defineBlockKind } from "@platforma-sdk/block-kind";
import type { PlRef, SUniversalPColumnId } from "@platforma-sdk/model";
import {
  isAnchoredPColumnId,
  isColumnUniversalId,
  isPlRef,
  parseJsonSafely,
} from "@platforma-sdk/model";
import { name, version } from "../package.json" with { type: "json" };

/** Whether distribution metrics are pooled across the cohort or computed within each subject. */
export type CalculationMode = "population" | "intra-subject";

/** How per-sample abundances are scaled before the metrics are computed. */
export type Normalization = "relative-frequency" | "clr";

/**
 * This block's init-params contract — everything a user sets by hand: the
 * abundance column to profile, the calculation mode, the three metadata
 * columns the metrics are computed over, the chronological order of the
 * timepoints, the normalization, the four numeric thresholds, and the subtitle
 * they type.
 *
 * `defaultBlockLabel` is absent: a `watchEffect` in `ui/src/pages/MainPage.vue`
 * derives it from the chosen columns' option labels, which only exist once the
 * result pool has resolved them.
 *
 * Every field is optional. A half-configured block is ordinary state the UI
 * reaches -- the abundance column alone leaves every metric column unset, and
 * the block stays that way until the user picks one -- and the projection hands
 * that state back untouched, so a required field would break the export/apply
 * round trip.
 */
export type BlockParams = {
  abundanceRef?: PlRef;
  calculationMode?: CalculationMode;
  subjectColumnRef?: SUniversalPColumnId;
  groupingColumnRef?: SUniversalPColumnId;
  temporalColumnRef?: SUniversalPColumnId;
  timepointOrder?: string[];
  normalization?: Normalization;
  presenceThreshold?: number;
  minAbundanceThreshold?: number;
  minSubjectCount?: number;
  topN?: number;
  customBlockLabel?: string;
};

// Identity (`name`/`version`) comes from this package's own `package.json`, so
// the on-wire `{name}@{version}` reference can never drift from what npm
// publishes; the bundler inlines the JSON import.
export const kind = defineBlockKind<BlockParams>({
  name,
  version,
  parseInitializationParams,
});

// Internals

const CALCULATION_MODES: readonly string[] = ["population", "intra-subject"];
const NORMALIZATIONS: readonly string[] = ["relative-frequency", "clr"];

/** The same contract at runtime, for params arriving from a template file rather than typed code. */
function parseInitializationParams(value: unknown): BlockParams {
  assertParamsObject(value);

  const {
    abundanceRef,
    calculationMode,
    subjectColumnRef,
    groupingColumnRef,
    temporalColumnRef,
    timepointOrder,
    normalization,
    presenceThreshold,
    minAbundanceThreshold,
    minSubjectCount,
    topN,
    customBlockLabel,
  } = value;

  if (abundanceRef !== undefined && !isPlRef(abundanceRef)) {
    throw new Error(
      "'abundanceRef' must be a reference to an upstream column, written as { block, name }.",
    );
  }
  if (calculationMode !== undefined && !CALCULATION_MODES.includes(calculationMode as string)) {
    throw new Error(`'calculationMode' must be one of: ${CALCULATION_MODES.join(", ")}.`);
  }
  if (subjectColumnRef !== undefined && !isColumnId(subjectColumnRef)) {
    throw new Error("'subjectColumnRef' must be a metadata column identifier.");
  }
  if (groupingColumnRef !== undefined && !isColumnId(groupingColumnRef)) {
    throw new Error("'groupingColumnRef' must be a metadata column identifier.");
  }
  if (temporalColumnRef !== undefined && !isColumnId(temporalColumnRef)) {
    throw new Error("'temporalColumnRef' must be a metadata column identifier.");
  }
  // The values of the temporal column in chronological order. They are the
  // column's own cell values, so anything the user's metadata holds is legal --
  // only the array-of-strings envelope is checked.
  if (timepointOrder !== undefined) {
    if (!Array.isArray(timepointOrder) || !timepointOrder.every((v) => typeof v === "string")) {
      throw new Error("'timepointOrder' must be an array of timepoint values.");
    }
  }
  if (normalization !== undefined && !NORMALIZATIONS.includes(normalization as string)) {
    throw new Error(`'normalization' must be one of: ${NORMALIZATIONS.join(", ")}.`);
  }
  // A clone's frequency within a group, so a fraction of 1.
  if (presenceThreshold !== undefined && !isFraction(presenceThreshold)) {
    throw new Error("'presenceThreshold' must be a number between 0 and 1.");
  }
  // An abundance floor compared against raw per-sample abundance, which need
  // not be a whole number.
  if (minAbundanceThreshold !== undefined && !isNonNegativeNumber(minAbundanceThreshold)) {
    throw new Error("'minAbundanceThreshold' must be a number greater than or equal to 0.");
  }
  // Both are counts: subjects a clone must appear in, and clones to plot.
  if (minSubjectCount !== undefined && !isPositiveInteger(minSubjectCount)) {
    throw new Error("'minSubjectCount' must be an integer greater than or equal to 1.");
  }
  if (topN !== undefined && !isPositiveInteger(topN)) {
    throw new Error("'topN' must be an integer greater than or equal to 1.");
  }
  if (customBlockLabel !== undefined && typeof customBlockLabel !== "string") {
    throw new Error("'customBlockLabel' must be a string.");
  }

  return {
    abundanceRef,
    calculationMode: calculationMode as CalculationMode | undefined,
    subjectColumnRef: subjectColumnRef as SUniversalPColumnId | undefined,
    groupingColumnRef: groupingColumnRef as SUniversalPColumnId | undefined,
    temporalColumnRef: temporalColumnRef as SUniversalPColumnId | undefined,
    timepointOrder: timepointOrder as string[] | undefined,
    normalization: normalization as Normalization | undefined,
    presenceThreshold,
    minAbundanceThreshold,
    minSubjectCount,
    topN,
    customBlockLabel,
  };
}

/**
 * A column identifier as this block stores it: a canonically serialized JSON key.
 * `isColumnUniversalId` covers the key forms the SDK's id encoding uses, but every
 * column id here comes from `resultPool.getCanonicalOptions`, which mints an
 * *anchored* key -- a shape none of those recognizes even though the SDK types it
 * `SUniversalPColumnId`. Both forms are accepted, or the kind would refuse the ids
 * the block itself writes into a template.
 */
function isColumnId(value: unknown): value is SUniversalPColumnId {
  if (typeof value !== "string") return false;
  return isColumnUniversalId(value) || isAnchoredPColumnId(parseJsonSafely(value));
}

function isFraction(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value) && value >= 0 && value <= 1;
}

function isNonNegativeNumber(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value) && value >= 0;
}

function isPositiveInteger(value: unknown): value is number {
  return typeof value === "number" && Number.isInteger(value) && value >= 1;
}
