import { createColumnValueId, createGlobalPObjectId } from "@platforma-sdk/model";
import { describe, expect, test } from "vitest";
import { kind } from "@platforma-open/milaboratories.spatiotemporal-analysis.kind";
import { blockDataModel } from "./dataModel";

const abundance = { __isRef: true, blockId: "mixcr", name: "abundance" };

/** A metadata column as the v1 pickers stored it: an id anchored against the abundance column. */
const anchored = (metadataName: string) =>
  JSON.stringify({
    axes: [{ anchor: "main", idx: 0 }],
    domain: { "pl7.app/metadata/name": metadataName },
    name: "pl7.app/metadata",
  });

const viewStates = {
  tableState: { version: 2 },
  heatmapState: { title: "Distribution heatmap", template: "heatmap" },
  temporalLineState: { title: "Temporal frequency trajectory", template: "curve_dots" },
  prevalenceHistogramState: { title: "Subject prevalence distribution", template: "bar" },
};

const storedAtV1 = {
  defaultBlockLabel: "Donor, Timepoint",
  customBlockLabel: "",
  abundanceRef: abundance,
  calculationMode: "intra-subject",
  subjectColumnRef: anchored("Donor"),
  groupingColumnRef: anchored("Tissue"),
  temporalColumnRef: anchored("Timepoint"),
  timepointOrder: ["day 0", "day 7", "day 30"],
  normalization: "clr",
  presenceThreshold: 0.1,
  minAbundanceThreshold: 2,
  minSubjectCount: 3,
  topN: 15,
  ...viewStates,
};

describe("v2 migration", () => {
  const { data } = blockDataModel.migrate({ version: "v1", data: storedAtV1 });

  test("the abundance ref becomes the result-pool id of the same upstream output", () => {
    expect(data.abundanceRef).toBe(createGlobalPObjectId("mixcr", "abundance"));
  });

  test("a subset picked with the abundance becomes its result-pool id", () => {
    const filterRef = { __isRef: true, blockId: "labeling", name: "lead" };
    const migrated = blockDataModel.migrate({ version: "v1", data: { ...storedAtV1, filterRef } });
    expect(migrated.data.filterRef).toBe(createGlobalPObjectId("labeling", "lead"));
  });

  test("anchored subject, grouping and temporal columns and plain timepoints are reset", () => {
    expect(data.subjectColumnRef).toBeUndefined();
    expect(data.groupingColumnRef).toBeUndefined();
    expect(data.temporalColumnRef).toBeUndefined();
    expect(data.timepointOrder).toEqual([]);
  });

  test("the analysis settings, custom label and view states are kept", () => {
    expect(data).toMatchObject({
      customBlockLabel: "",
      calculationMode: "intra-subject",
      normalization: "clr",
      presenceThreshold: 0.1,
      minAbundanceThreshold: 2,
      minSubjectCount: 3,
      topN: 15,
      ...viewStates,
    });
  });

  test("columns that already are result-pool ids and their value ids are kept", () => {
    const temporal = createGlobalPObjectId("samples", "metadata.timepoint");
    const subject = createGlobalPObjectId("samples", "metadata.donor");
    const day0 = createColumnValueId({ column: temporal, value: "day 0" });
    const migrated = blockDataModel.migrate({
      version: "v1",
      data: {
        ...storedAtV1,
        subjectColumnRef: subject,
        temporalColumnRef: temporal,
        timepointOrder: [day0, "day 7"],
      },
    }).data;
    expect(migrated.subjectColumnRef).toBe(subject);
    expect(migrated.groupingColumnRef).toBeUndefined();
    expect(migrated.temporalColumnRef).toBe(temporal);
    expect(migrated.timepointOrder).toEqual([day0]);
  });

  test("a block without an abundance column keeps it unset", () => {
    const { abundanceRef: _, ...withoutAbundance } = storedAtV1;
    const migrated = blockDataModel.migrate({ version: "v1", data: withoutAbundance }).data;
    expect(migrated.abundanceRef).toBeUndefined();
  });

  test("the migrated picks are what the kind accepts", () => {
    const filterRef = { __isRef: true, blockId: "labeling", name: "lead" };
    const migrated = blockDataModel.migrate({ version: "v1", data: { ...storedAtV1, filterRef } });
    const { abundanceRef, subjectColumnRef, groupingColumnRef, temporalColumnRef, timepointOrder } =
      migrated.data;
    expect(() =>
      kind.parseInitializationParams({
        abundanceRef,
        filterRef: migrated.data.filterRef,
        subjectColumnRef,
        groupingColumnRef,
        temporalColumnRef,
        timepointOrder,
      }),
    ).not.toThrow();
  });

  test("a migrated block is left as is when it is loaded again", () => {
    expect(blockDataModel.migrate({ version: "v3", data }).data).toEqual(data);
  });
});

describe("v3 migration", () => {
  const storedAtV2 = {
    ...storedAtV1,
    defaultBlockLabel: "Intra-Subject, Tissue, Day, Donor, CLR",
    customBlockLabel: "My distribution",
    abundanceRef: createGlobalPObjectId("mixcr", "abundance"),
    subjectColumnRef: undefined,
    groupingColumnRef: undefined,
    temporalColumnRef: undefined,
    timepointOrder: [],
  };
  const { data } = blockDataModel.migrate({ version: "v2", data: storedAtV2 });

  test("the derived label is dropped and the subtitle text starts empty", () => {
    expect(data).not.toHaveProperty("defaultBlockLabel");
    expect(data.subtitleLabel).toBe("");
  });

  test("the custom label and the settings are kept", () => {
    expect(data).toMatchObject({
      customBlockLabel: "My distribution",
      calculationMode: "intra-subject",
      normalization: "clr",
    });
  });
});
