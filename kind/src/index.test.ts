import { createColumnValueId, createGlobalPObjectId } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import { kind } from "./index";

const parse = (params: unknown) => kind.parseInitializationParams(params);

/**
 * The anchored id `resultPool.getCanonicalOptions` minted for the metadata
 * columns before they were stored by result-pool ids. It names no producing
 * block, so the contract no longer accepts it.
 */
const ANCHORED_ID =
  '{"axes":[{"anchor":"main","idx":0}],"domain":{"pl7.app/metadata/id":"tissue"},"name":"pl7.app/metadata"}';

/** The serialized key forms `isColumnUniversalId` recognizes. */
const GLOBAL_ID = '{"__isRef":true,"blockId":"b1","name":"pl7.app/metadata"}';
const LOCAL_ID = '{"name":"pl7.app/metadata","resolvePath":["a","b"]}';
const FILTERED_ID = `{"__isFiltered":true,"axisFilters":[[0,"lung"]],"source":${JSON.stringify(GLOBAL_ID)}}`;
const OVERRIDDEN_ID = `{"__isOverridden":true,"source":${JSON.stringify(GLOBAL_ID)},"specOverrides":{"annotations":{"pl7.app/label":"Tissue"}}}`;

const ABUNDANCE_ID = '{"__isRef":true,"blockId":"b1","name":"pf/abundance"}';

/** Timepoints picked from the temporal column, as column value ids. */
const timepoint = (value: string) =>
  createColumnValueId({ column: createGlobalPObjectId("b1", "pl7.app/metadata"), value });

const COLUMN_FIELDS = ["subjectColumnRef", "groupingColumnRef", "temporalColumnRef"] as const;

describe.each(COLUMN_FIELDS)("%s", (field) => {
  it.each([
    ["a global key id", GLOBAL_ID],
    ["a local key id", LOCAL_ID],
    ["a filtered key id", FILTERED_ID],
    ["an overridden key id", OVERRIDDEN_ID],
  ])("accepts %s", (_label, id) => {
    expect(parse({ [field]: id })).toEqual({ [field]: id });
  });

  it.each([
    ["an anchored id, as getCanonicalOptions minted it", ANCHORED_ID],
    ["a string that is not JSON", "pl7.app/metadata"],
    ["malformed JSON", '{"name":"pl7.app/metadata"'],
    ["JSON that is not a column key", '{"foo":1}'],
    ["JSON that is not an object", '"pl7.app/metadata"'],
    ["a number", 42],
    ["null", null],
    ["the key form rather than its serialization", { name: "pl7.app/metadata", resolvePath: [] }],
  ])("rejects %s", (_label, id) => {
    expect(() => parse({ [field]: id })).toThrow(
      `'${field}' must be a metadata column, written as a column id.`,
    );
  });
});

describe("abundanceRef", () => {
  it("accepts a column id", () => {
    expect(parse({ abundanceRef: ABUNDANCE_ID })).toEqual({ abundanceRef: ABUNDANCE_ID });
  });

  it.each([
    ["a PlRef", { __isRef: true, blockId: "b1", name: "pf/abundance" }],
    ["an anchored id", ANCHORED_ID],
    ["a number", 42],
  ])("rejects %s", (_label, abundanceRef) => {
    expect(() => parse({ abundanceRef })).toThrow(
      "'abundanceRef' must be an upstream column, written as a column id.",
    );
  });
});

describe("the enumerated fields", () => {
  it.each(["population", "intra-subject"])("accepts calculationMode %s", (calculationMode) => {
    expect(parse({ calculationMode })).toEqual({ calculationMode });
  });

  it.each(["relative-frequency", "clr"])("accepts normalization %s", (normalization) => {
    expect(parse({ normalization })).toEqual({ normalization });
  });

  it("rejects a calculationMode outside the vocabulary", () => {
    expect(() => parse({ calculationMode: "per-sample" })).toThrow("'calculationMode' must be one");
  });

  it("rejects a normalization outside the vocabulary", () => {
    expect(() => parse({ normalization: "log" })).toThrow("'normalization' must be one");
  });
});

describe("timepointOrder", () => {
  it.each([
    ["an empty array, as a block with no temporal column has", []],
    ["values of the temporal column", [timepoint("d0"), timepoint("d7"), timepoint("d28")]],
    ["numeric values of the temporal column, as text", [timepoint("0"), timepoint("7")]],
  ])("accepts %s", (_label, timepointOrder) => {
    expect(parse({ timepointOrder })).toEqual({ timepointOrder });
  });

  it.each([
    ["a bare column value id", timepoint("d0")],
    ["plain values", ["d0", "d7"]],
    ["an array of numbers", [0, 7]],
    ["an array with a hole", [timepoint("d0"), null]],
  ])("rejects %s", (_label, timepointOrder) => {
    expect(() => parse({ timepointOrder })).toThrow("'timepointOrder' must be an array");
  });
});

describe("the numeric thresholds", () => {
  it.each([0, 0.0001, 0.5, 1])("accepts presenceThreshold %s", (presenceThreshold) => {
    expect(parse({ presenceThreshold })).toEqual({ presenceThreshold });
  });

  it.each([-0.1, 1.1, Number.NaN, Number.POSITIVE_INFINITY, "0.5"])(
    "rejects presenceThreshold %s",
    (presenceThreshold) => {
      expect(() => parse({ presenceThreshold })).toThrow("'presenceThreshold' must be a number");
    },
  );

  it.each([0, 2.5, 1000])(
    "accepts minAbundanceThreshold %s -- raw abundance need not be whole",
    (minAbundanceThreshold) => {
      expect(parse({ minAbundanceThreshold })).toEqual({ minAbundanceThreshold });
    },
  );

  it.each([-1, Number.NaN, "0"])("rejects minAbundanceThreshold %s", (minAbundanceThreshold) => {
    expect(() => parse({ minAbundanceThreshold })).toThrow("'minAbundanceThreshold' must be");
  });

  it.each(["minSubjectCount", "topN"] as const)("accepts %s as a count", (field) => {
    expect(parse({ [field]: 1 })).toEqual({ [field]: 1 });
    expect(parse({ [field]: 20 })).toEqual({ [field]: 20 });
  });

  it.each(["minSubjectCount", "topN"] as const)("rejects a non-count %s", (field) => {
    for (const bad of [0, -1, 1.5, Number.NaN, "5"]) {
      expect(() => parse({ [field]: bad })).toThrow(`'${field}' must be an integer`);
    }
  });
});

describe("the params envelope", () => {
  it("accepts an empty object -- every field is optional", () => {
    expect(parse({})).toEqual({});
  });

  it("accepts a fully configured block", () => {
    const params = {
      abundanceRef: ABUNDANCE_ID,
      calculationMode: "intra-subject",
      subjectColumnRef: GLOBAL_ID,
      groupingColumnRef: GLOBAL_ID,
      temporalColumnRef: GLOBAL_ID,
      timepointOrder: [timepoint("d0"), timepoint("d7")],
      normalization: "clr",
      presenceThreshold: 0.01,
      minAbundanceThreshold: 5,
      minSubjectCount: 3,
      topN: 50,
      customBlockLabel: "Lung vs blood",
    };
    expect(parse(params)).toEqual(params);
  });

  it("drops keys the contract does not name", () => {
    expect(parse({ topN: 20, notAParam: "x" })).toEqual({ topN: 20 });
  });

  it("rejects params that are not an object", () => {
    expect(() => parse(null)).toThrow();
    expect(() => parse([ANCHORED_ID])).toThrow();
    expect(() => parse(5)).toThrow();
  });

  it("rejects a customBlockLabel that is not a string", () => {
    expect(() => parse({ customBlockLabel: 42 })).toThrow("'customBlockLabel' must be a string.");
  });
});
