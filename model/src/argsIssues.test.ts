import {
  BlockStorageFacadeCallbacks,
  createBlockStorage,
  createColumnValueId,
  createGlobalPObjectId,
  stringifyJson,
} from "@platforma-sdk/model";
import { beforeAll, describe, expect, test } from "vitest";

const callbackRegistry: Record<string, (storageJson: string) => unknown> = {};

beforeAll(async () => {
  Object.assign(globalThis, { cfgRenderCtx: { callbackRegistry } });
  await import("./index");
});

const abundance = createGlobalPObjectId("mixcr", "abundance");
const grouping = createGlobalPObjectId("samples", "metadata.tissue");
const temporal = createGlobalPObjectId("samples", "metadata.timepoint");
const subject = createGlobalPObjectId("samples", "metadata.donor");
const timepoint = (value: string) => createColumnValueId({ column: temporal, value });

const configured = {
  defaultBlockLabel: "Tissue",
  customBlockLabel: "",
  abundanceRef: abundance,
  calculationMode: "population",
  groupingColumnRef: grouping,
  timepointOrder: [],
  normalization: "relative-frequency",
  presenceThreshold: 0,
  minAbundanceThreshold: 0,
  minSubjectCount: 1,
  topN: 20,
  tableState: { version: 2 },
  heatmapState: { title: "Distribution heatmap", template: "heatmap" },
  temporalLineState: { title: "Temporal frequency trajectory", template: "curve_dots" },
  prevalenceHistogramState: { title: "Subject prevalence distribution", template: "bar" },
};

describe("args issues", () => {
  test("a block with a grouping column derives runnable args without its view states", () => {
    const result = deriveArgs(configured);
    expect(result).toMatchObject({
      value: { abundanceRef: abundance, groupingColumnRef: grouping, topN: 20 },
      issues: [],
    });
    expect(result).not.toHaveProperty("value.tableState");
  });

  test("a temporal column with two timepoints is enough to compare along", () => {
    const data = {
      ...configured,
      groupingColumnRef: undefined,
      temporalColumnRef: temporal,
      timepointOrder: [timepoint("day 0"), timepoint("day 7")],
    };
    expect(deriveArgs(data)).toMatchObject({ issues: [] });
  });

  test("intra-subject mode runs once the subject column is chosen", () => {
    const data = { ...configured, calculationMode: "intra-subject", subjectColumnRef: subject };
    expect(deriveArgs(data)).toMatchObject({ issues: [] });
  });

  test.for([
    {
      why: "no abundance column",
      data: { ...configured, abundanceRef: undefined },
      message: "Abundance ref required",
    },
    {
      why: "neither a grouping nor a temporal column",
      data: { ...configured, groupingColumnRef: undefined },
      message: "At least grouping or temporal variable required",
    },
    {
      why: "a temporal column with a single timepoint and no grouping column",
      data: {
        ...configured,
        groupingColumnRef: undefined,
        temporalColumnRef: temporal,
        timepointOrder: [timepoint("day 0")],
      },
      message: "At least grouping or temporal variable required",
    },
    {
      why: "intra-subject mode without a subject column",
      data: { ...configured, calculationMode: "intra-subject" },
      message: "Subject required in intra-subject mode",
    },
  ])("$why blocks Run with its message", ({ data, message }) => {
    expect(deriveArgs(data)).toMatchObject({
      error: message,
      issues: [{ severity: "error", message }],
    });
  });

  test("every unmet requirement is reported", () => {
    const data = {
      ...configured,
      abundanceRef: undefined,
      groupingColumnRef: undefined,
      calculationMode: "intra-subject",
    };
    expect(deriveArgs(data)).toMatchObject({
      issues: [
        { severity: "error", message: "Abundance ref required" },
        { severity: "error", message: "At least grouping or temporal variable required" },
        { severity: "error", message: "Subject required in intra-subject mode" },
      ],
    });
  });
});

// Internals

function deriveArgs(data: unknown): unknown {
  const derive = callbackRegistry[BlockStorageFacadeCallbacks.ArgsDerive];
  if (derive === undefined) throw new Error("The args derive callback is not registered");
  return derive(stringifyJson(createBlockStorage(data)));
}
