import {
  BlockStorageFacadeCallbacks,
  createColumnValueId,
  createGlobalPObjectId,
  parseColumnValueId,
  stringifyJson,
} from "@platforma-sdk/model";
import { beforeAll, describe, expect, test } from "vitest";

const callbackRegistry: Record<string, (...args: string[]) => unknown> = {};

beforeAll(async () => {
  Object.assign(globalThis, { cfgRenderCtx: { callbackRegistry } });
  await import("./index");
});

// Entry ids of the upstream blocks as a template file holds them, and the ids the project being
// built gave those blocks.
const blockIds = { "template-mixcr": "project-mixcr", "template-samples": "project-samples" };

const column = (blockId: string, name: string) => createGlobalPObjectId(blockId, name);
const timepoints = (temporal: string) =>
  ["Day 0", "Day 7", "Day 14", "Day 28"].map((value) =>
    createColumnValueId({ column: temporal as never, value }),
  );

/** Params of a calculated block as a template made by export holds them. */
const params = {
  abundanceRef: column("template-mixcr", "pf.abundance0"),
  calculationMode: "population",
  customBlockLabel: "",
  groupingColumnRef: column("template-samples", "pf.metadata.tissue"),
  minAbundanceThreshold: 0,
  minSubjectCount: 1,
  normalization: "relative-frequency",
  presenceThreshold: 0,
  subjectColumnRef: column("template-samples", "pf.metadata.donor"),
  temporalColumnRef: column("template-samples", "pf.metadata.timepoint"),
  timepointOrder: timepoints(column("template-samples", "pf.metadata.timepoint")),
  topN: 20,
};

describe("a block created from template params", () => {
  test("points its columns and timepoints at the blocks of the project", () => {
    const { data } = parseStorage(initialStorage(params));
    const temporal = column("project-samples", "pf.metadata.timepoint");
    expect(data).toMatchObject({
      abundanceRef: column("project-mixcr", "pf.abundance0"),
      groupingColumnRef: column("project-samples", "pf.metadata.tissue"),
      subjectColumnRef: column("project-samples", "pf.metadata.donor"),
      temporalColumnRef: temporal,
      timepointOrder: timepoints(temporal),
    });
    // The timepoints name the temporal column itself, so none of them shows as not available.
    for (const id of data.timepointOrder as string[]) {
      expect(parseColumnValueId(id as never).column).toBe(temporal);
    }
  });

  test("exports the params it was created from", () => {
    expect(templateParams(initialStorage(params))).toEqual(relocated(params));
  });

  test("keeps its args when its window fills in the subtitle", () => {
    const storage = parseStorage(initialStorage(params));
    const before = deriveArgs(storage);
    expect(before).toMatchObject({ issues: [] });

    // What the window writes once it opens: the result-pool labels of the chosen columns.
    const opened = {
      ...storage,
      data: { ...storage.data, subtitleLabel: "Tissue, Timepoint, Donor" },
    };
    expect(deriveArgs(opened)).toEqual(before);
  });

  test("has args that no field outside its params reaches", () => {
    // A field the params do not carry is one only the window can fill in -- a label it reads
    // from the result pool, a view state. Were it in the args, opening the window of a
    // calculated block would make the block and every block below it stale.
    const storage = parseStorage(initialStorage(params));
    const exported = templateParams(initialStorage(params)) as Record<string, unknown>;
    const filledByWindow = Object.keys(storage.data).filter((key) => !(key in exported));
    expect(filledByWindow).not.toHaveLength(0);

    const before = deriveArgs(storage);
    for (const key of filledByWindow) {
      const opened = { ...storage, data: { ...storage.data, [key]: "written by the window" } };
      expect(deriveArgs(opened), key).toEqual(before);
    }
  });
});

// Internals

type Storage = Record<string, unknown> & { data: Record<string, unknown> };

/** The params a template written for the project would hold: every reference relocated. */
function relocated(value: unknown): unknown {
  let text = JSON.stringify(value);
  for (const [from, to] of Object.entries(blockIds)) text = text.replaceAll(from, to);
  return JSON.parse(text);
}

function initialStorage(entryParams: unknown): string {
  const result = callback(BlockStorageFacadeCallbacks.StorageInitialFromParams)(
    stringifyJson(entryParams),
    stringifyJson(blockIds),
  ) as { error?: string; storageJson?: string };
  if (result.error !== undefined) throw new Error(result.error);
  return result.storageJson!;
}

/** Storage with its block data taken out as `data`, so a test can change it. */
function parseStorage(storageJson: string): Storage {
  const { __data, ...rest } = JSON.parse(storageJson) as Record<string, unknown>;
  return { ...rest, data: __data as Record<string, unknown> };
}

function deriveArgs({ data, ...rest }: Storage): unknown {
  return callback(BlockStorageFacadeCallbacks.ArgsDerive)(stringifyJson({ ...rest, __data: data }));
}

function templateParams(storageJson: string): unknown {
  const result = callback(BlockStorageFacadeCallbacks.InitializationParamsDerive)(storageJson) as {
    error?: string;
    value?: unknown;
  };
  if (result.error !== undefined) throw new Error(result.error);
  return result.value;
}

function callback(name: string): (...args: string[]) => unknown {
  const registered = callbackRegistry[name];
  if (registered === undefined) throw new Error(`The ${name} callback is not registered`);
  return registered;
}
