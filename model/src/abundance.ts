import type {
  ColumnUniversalId,
  DatasetOption,
  DatasetSelection,
  PColumnSpec,
  PlRef,
  PObjectSpec,
  RelaxedColumnSelector,
  RenderCtxBase,
} from "@platforma-sdk/model";
import {
  buildDatasetOptions,
  createDatasetSelection,
  createGlobalPObjectId,
  createPrimaryRef,
  deriveColumnOptions,
  isPColumnSpec,
  isPlRef,
  parseColumnIdSafely,
} from "@platforma-sdk/model";
import { exact, resultPool } from "./columns";

/** The abundance and the subset it is narrowed to, as the block stores them. */
export type AbundancePick = {
  abundanceRef?: ColumnUniversalId;
  filterRef?: ColumnUniversalId;
};

/**
 * The abundance picker's options: primary raw abundances per sample of upstream blocks,
 * labelled with the column's own label (e.g. "Number of Reads"), each with its subset columns
 * (e.g. Repertoire Labeling labels) as filters.
 */
export function getAbundanceOptions(ctx: RenderCtxBase): DatasetOption[] {
  // Only the filters are taken from here: the abundances below are matched to them by column id.
  const withFilters =
    buildDatasetOptions(ctx, {
      primary: isPrimaryAbundance,
      // Only subsets keyed by the element axis alone: a sample subset would drop whole samples.
      filter: (spec) =>
        isPColumnSpec(spec) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name !== "pl7.app/sampleId",
      requireEnrichments: false,
    }) ?? [];

  const columns = resultPool().filter({ include: PRIMARY_ABUNDANCE });
  return deriveColumnOptions(columns, { includeNativeLabel: true }).flatMap(
    ({ id, label }): DatasetOption[] => {
      const ref = toPlRef(id);
      if (ref === undefined) return [];
      const primary = { ref, label };
      const filters = withFilters.find((o) => toColumnId(o.primary.ref) === id)?.filters;
      return [filters === undefined ? { primary } : { primary, filters }];
    },
  );
}

/** The abundance picker's selection of the stored abundance and subset. */
export function toAbundanceSelection(pick: AbundancePick): DatasetSelection | undefined {
  const column = pick.abundanceRef === undefined ? undefined : toPlRef(pick.abundanceRef);
  if (column === undefined) return undefined;
  const filter = pick.filterRef === undefined ? undefined : toPlRef(pick.filterRef);
  return createDatasetSelection(createPrimaryRef(column, filter));
}

/** The abundance and subset an abundance picker selection stores. */
export function fromAbundanceSelection(selection: DatasetSelection | undefined): AbundancePick {
  const primary = selection?.primary;
  return {
    abundanceRef: primary === undefined ? undefined : toColumnId(primary.column),
    filterRef: primary?.filter === undefined ? undefined : toColumnId(primary.filter),
  };
}

/** The result-pool id of the upstream output a ref points to: it names the same block and column. */
export function toColumnId(ref: PlRef): ColumnUniversalId {
  return createGlobalPObjectId(ref.blockId, ref.name);
}

// Internals

/** Primary raw abundances per sample: the first of exactly two axes is the sample. */
const PRIMARY_ABUNDANCE: RelaxedColumnSelector = {
  axes: [{ name: [exact("pl7.app/sampleId")] }, {}],
  partialAxesMatch: false,
  annotations: {
    "pl7.app/isAbundance": [exact("true")],
    "pl7.app/abundance/normalized": [exact("false")],
    "pl7.app/abundance/isPrimary": [exact("true")],
  },
};

/** The same columns as `PRIMARY_ABUNDANCE`, told by their spec. */
function isPrimaryAbundance(spec: PObjectSpec): spec is PColumnSpec {
  if (!isPColumnSpec(spec)) return false;
  if (spec.axesSpec.length !== 2 || spec.axesSpec[0].name !== "pl7.app/sampleId") return false;
  const a = spec.annotations ?? {};
  return (
    a["pl7.app/isAbundance"] === "true" &&
    a["pl7.app/abundance/normalized"] === "false" &&
    a["pl7.app/abundance/isPrimary"] === "true"
  );
}

/** The upstream output a result-pool id names, as the ref the abundance picker works with. */
function toPlRef(id: ColumnUniversalId): PlRef | undefined {
  const key = parseColumnIdSafely(id);
  return isPlRef(key) ? key : undefined;
}
