import type { ColumnValueOption } from "@platforma-sdk/model";
import { readColumnValueOptions } from "@platforma-sdk/model";
import type { PlColumnValueSettings } from "@platforma-sdk/ui-vue";
import type { ShallowRef, WatchSource } from "vue";
import { shallowRef, watch } from "vue";

/** The values of the picked temporal column, as far as they are read. */
export type TimepointOptions =
  /** No temporal column is picked. */
  | { readonly type: "none" }
  /** The values are not read yet. */
  | { readonly type: "loading" }
  /** The temporal column's values, in natural order ("Day 2" before "Day 10"). */
  | { readonly type: "read"; readonly options: ColumnValueOption[] }
  /** The temporal column is gone: its upstream block no longer produces it. */
  | { readonly type: "absent" }
  /** The values cannot be read. */
  | { readonly type: "failed"; readonly message: string };

/**
 * Reads the values of the temporal column through the frame its `settings` hold. A read made for
 * earlier settings is dropped. Called from `setup`.
 */
export function useTimepointOptions(
  settings: WatchSource<PlColumnValueSettings>,
): Readonly<ShallowRef<TimepointOptions>> {
  const state = shallowRef<TimepointOptions>(NONE);
  let request = 0;

  watch(
    settings,
    (current) => {
      const read = ++request;
      switch (current.type) {
        case "none":
          state.value = NONE;
          return;
        case "pending":
          state.value = LOADING;
          return;
        case "absent":
          state.value = ABSENT;
          return;
        case "present":
          state.value = LOADING;
          readColumnValueOptions(current.model).then(
            (options) => {
              if (read === request) state.value = { type: "read", options };
            },
            (error: unknown) => {
              if (read === request) state.value = { type: "failed", message: messageOf(error) };
            },
          );
      }
    },
    { immediate: true },
  );

  return state;
}

// Internals

const NONE: TimepointOptions = { type: "none" };
const LOADING: TimepointOptions = { type: "loading" };
const ABSENT: TimepointOptions = { type: "absent" };

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
