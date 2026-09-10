import { describe, expect, it } from "vitest";
import { kind } from "./index";

const parse = (params: unknown) => kind.parseInitializationParams(params);

/**
 * The shape `resultPool.getCanonicalOptions` mints for `metadataOptions`, which is
 * what the settings panel writes into all three metadata column fields. Built from
 * the anchored selector the model asks for -- `axes: [{ anchor: "main", idx: 0 }]`,
 * `name: "pl7.app/metadata"` -- with the domain the pool fills in.
 *
 * `isColumnUniversalId` does not recognize this form even though the SDK types it
 * `SUniversalPColumnId`. A kind that refused it would refuse the ids the block itself
 * writes, so applying a template exported from this block would fail. Anything that
 * stops this test passing has broken the export/apply round trip.
 */
const ANCHORED_ID =
  '{"axes":[{"anchor":"main","idx":0}],"domain":{"pl7.app/metadata/id":"tissue"},"name":"pl7.app/metadata"}';

/** The other four serialized key forms `isColumnUniversalId` recognizes. */
const GLOBAL_ID = '{"__isRef":true,"blockId":"b1","name":"pl7.app/metadata"}';
const LOCAL_ID = '{"name":"pl7.app/metadata","resolvePath":["a","b"]}';
const FILTERED_ID = `{"__isFiltered":true,"axisFilters":[[0,"lung"]],"source":${JSON.stringify(GLOBAL_ID)}}`;
const OVERRIDDEN_ID = `{"__isOverridden":true,"source":${JSON.stringify(GLOBAL_ID)},"specOverrides":{"annotations":{"pl7.app/label":"Tissue"}}}`;

const COLUMN_FIELDS = ["subjectColumnRef", "groupingColumnRef", "temporalColumnRef"] as const;

describe.each(COLUMN_FIELDS)("%s", (field) => {
  it.each([
    ["an anchored id, as getCanonicalOptions mints it", ANCHORED_ID],
    ["a global key id", GLOBAL_ID],
    ["a local key id", LOCAL_ID],
    ["a filtered key id", FILTERED_ID],
    ["an overridden key id", OVERRIDDEN_ID],
  ])("accepts %s", (_label, id) => {
    expect(parse({ [field]: id })).toEqual({ [field]: id });
  });

  it.each([
    ["a string that is not JSON", "pl7.app/metadata"],
    ["malformed JSON", '{"name":"pl7.app/metadata"'],
    ["JSON that is not a column key", '{"foo":1}'],
    ["JSON that is not an object", '"pl7.app/metadata"'],
    ["a number", 42],
    ["null", null],
    ["the key form rather than its serialization", { name: "pl7.app/metadata", resolvePath: [] }],
  ])("rejects %s", (_label, id) => {
    expect(() => parse({ [field]: id })).toThrow(
      `'${field}' must be a metadata column identifier.`,
    );
  });
});

describe("abundanceRef", () => {
  it("accepts a PlRef", () => {
    const abundanceRef = { __isRef: true as const, blockId: "b1", name: "pf/abundance" };
    expect(parse({ abundanceRef })).toEqual({ abundanceRef });
  });

  it.each([
    ["a column id", ANCHORED_ID],
    ["an object missing the marker", { blockId: "b1", name: "pf/abundance" }],
    ["a number", 42],
  ])("rejects %s", (_label, abundanceRef) => {
    expect(() => parse({ abundanceRef })).toThrow("'abundanceRef' must be a reference");
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
    ["the column's own values", ["d0", "d7", "d28"]],
  ])("accepts %s", (_label, timepointOrder) => {
    expect(parse({ timepointOrder })).toEqual({ timepointOrder });
  });

  it.each([
    ["a bare string", "d0"],
    ["an array of numbers", [0, 7]],
    ["an array with a hole", ["d0", null]],
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
      abundanceRef: { __isRef: true as const, blockId: "b1", name: "pf/abundance" },
      calculationMode: "intra-subject",
      subjectColumnRef: ANCHORED_ID,
      groupingColumnRef: ANCHORED_ID,
      temporalColumnRef: ANCHORED_ID,
      timepointOrder: ["d0", "d7"],
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
