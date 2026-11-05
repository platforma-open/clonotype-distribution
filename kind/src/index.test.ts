import {
  createColumnValueId,
  createGlobalPObjectId,
  createLocalPObjectId,
} from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import { kind } from "./index";

const parse = (params: unknown) => kind.parseInitializationParams(params);

/** A metadata column of an upstream block, as the settings panel stores it: its result-pool id. */
const GLOBAL_ID = createGlobalPObjectId("b1", "pl7.app/metadata");

/**
 * The anchored form the settings panel stored before the columns became result-pool ids. It
 * names a column relative to the abundance column, which the dataset picker cannot show.
 */
const ANCHORED_ID =
  '{"axes":[{"anchor":"main","idx":0}],"domain":{"pl7.app/metadata/id":"tissue"},"name":"pl7.app/metadata"}';

/** Other column id forms, which name no single output of an upstream block. */
const LOCAL_ID = createLocalPObjectId(["a", "b"], "pl7.app/metadata");
const FILTERED_ID = `{"__isFiltered":true,"axisFilters":[[0,"lung"]],"source":${JSON.stringify(GLOBAL_ID)}}`;
const OVERRIDDEN_ID = `{"__isOverridden":true,"source":${JSON.stringify(GLOBAL_ID)},"specOverrides":{"annotations":{"pl7.app/label":"Tissue"}}}`;

const COLUMN_FIELDS = [
  "abundanceRef",
  "filterRef",
  "subjectColumnRef",
  "groupingColumnRef",
  "temporalColumnRef",
] as const;

describe.each(COLUMN_FIELDS)("%s", (field) => {
  it("accepts the result-pool id of an upstream output", () => {
    expect(parse({ [field]: GLOBAL_ID })).toEqual({ [field]: GLOBAL_ID });
  });

  it.each([
    ["an anchored id", ANCHORED_ID],
    ["a local id", LOCAL_ID],
    ["a filtered id", FILTERED_ID],
    ["an overridden id", OVERRIDDEN_ID],
    ["a PlRef object rather than its id", { __isRef: true, blockId: "b1", name: "pf/column" }],
    ["a string that is not JSON", "pl7.app/metadata"],
    ["malformed JSON", '{"name":"pl7.app/metadata"'],
    ["a number", 42],
    ["null", null],
  ])("rejects %s", (_label, id) => {
    expect(() => parse({ [field]: id })).toThrow(`'${field}' must be`);
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
  const timepoint = (value: string) => createColumnValueId({ column: GLOBAL_ID, value });

  it.each([
    ["an empty array, as a block with no temporal column has", []],
    ["value ids of the temporal column", [timepoint("d0"), timepoint("d7"), timepoint("d28")]],
  ])("accepts %s", (_label, timepointOrder) => {
    expect(parse({ timepointOrder })).toEqual({ timepointOrder });
  });

  it.each([
    ["a bare string", "d0"],
    ["plain values, as stored before value ids", ["d0", "d7"]],
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
      abundanceRef: createGlobalPObjectId("b1", "pf/abundance"),
      filterRef: createGlobalPObjectId("b2", "labels.igg"),
      calculationMode: "intra-subject",
      subjectColumnRef: createGlobalPObjectId("b0", "metadata.donor"),
      groupingColumnRef: createGlobalPObjectId("b0", "metadata.tissue"),
      temporalColumnRef: GLOBAL_ID,
      timepointOrder: [
        createColumnValueId({ column: GLOBAL_ID, value: "d0" }),
        createColumnValueId({ column: GLOBAL_ID, value: "d7" }),
      ],
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
    expect(() => parse([GLOBAL_ID])).toThrow();
    expect(() => parse(5)).toThrow();
  });

  it("rejects a customBlockLabel that is not a string", () => {
    expect(() => parse({ customBlockLabel: 42 })).toThrow("'customBlockLabel' must be a string.");
  });
});
