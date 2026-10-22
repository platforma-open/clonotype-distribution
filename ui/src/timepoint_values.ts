import type { MetadataColumnFrame } from "@platforma-open/milaboratories.spatiotemporal-analysis.model";
import type { ColumnUniversalId, ColumnValueOption } from "@platforma-sdk/model";
import { deriveColumnValueOptions, PFrameImpl, pTableValue } from "@platforma-sdk/model";
import type { ShallowRef } from "vue";
import { shallowRef, watch } from "vue";

/** The values of the picked temporal column, as far as they are read. */
export type TimepointValuesState =
  | { readonly type: "loading" }
  /** The column is gone: a stored id may name a column its upstream block no longer produces. */
  | { readonly type: "absent" }
  | { readonly type: "read"; readonly options: ColumnValueOption[] }
  | { readonly type: "failed"; readonly message: string };

/**
 * The values of the picked temporal `column` as options of the timepoints, read through the
 * `source` the model gives for it. `undefined` while no column is picked. Loading while the source
 * is resolving or still describes an earlier pick. Called from `setup`.
 */
export function useTimepointValues(params: {
  source: () => MetadataColumnFrame | undefined;
  column: () => ColumnUniversalId | undefined;
}): Readonly<ShallowRef<TimepointValuesState | undefined>> {
  const state = shallowRef<TimepointValuesState | undefined>();
  let request = 0;

  watch(
    // An output recomputed after its producer reruns keeps its object and changes the frame
    // inside it; copying it makes the read follow the frame.
    [params.column, () => copyOf(params.source())],
    ([column, source]) => {
      const current = ++request;
      if (column === undefined) {
        state.value = undefined;
        return;
      }
      if (source === undefined || source.column !== column || source.type === "resolving") {
        state.value = LOADING;
        return;
      }
      if (source.type === "absent") {
        state.value = ABSENT;
        return;
      }
      state.value = LOADING;
      readTimepointOptions(source).then(
        (options) => {
          if (current === request) state.value = { type: "read", options };
        },
        (error: unknown) => {
          if (current === request) state.value = { type: "failed", message: messageOf(error) };
        },
      );
    },
    { immediate: true },
  );

  return state;
}

/**
 * The distinct values of the column `source` holds, read in the UI through its frame, as options
 * of the timepoints, ordered naturally ("day 2" before "day 10"). Missing and empty values are not
 * offered. Rejects when the column has more distinct values than the timepoints can list.
 */
export async function readTimepointOptions(
  source: Extract<MetadataColumnFrame, { type: "present" }>,
): Promise<ColumnValueOption[]> {
  const response = await new PFrameImpl(source.frame).getUniqueValues({
    columnId: source.frameColumn,
    filters: [],
    limit: VALUES_LIMIT,
  });
  if (response.overflow) {
    throw new Error(`The column has more than ${VALUES_LIMIT} distinct values`);
  }

  const values: string[] = [];
  for (let row = 0; row < response.values.data.length; row++) {
    const value = pTableValue(response.values, row, { na: null });
    if (value !== null) values.push(String(value));
  }
  return deriveColumnValueOptions(source.column, values);
}

// Internals

/** More distinct values than the timepoints can usefully list. */
const VALUES_LIMIT = 100_000;

const LOADING: TimepointValuesState = { type: "loading" };

const ABSENT: TimepointValuesState = { type: "absent" };

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}

function copyOf<T extends object>(source: T | undefined): T | undefined {
  return source === undefined ? undefined : { ...source };
}
