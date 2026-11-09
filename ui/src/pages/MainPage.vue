<script setup lang="ts">
import {
  defaultSubtitle,
  fromAbundanceSelection,
  toAbundanceSelection,
} from "@platforma-open/milaboratories.spatiotemporal-analysis.model";
import type { ColumnUniversalId, DatasetSelection } from "@platforma-sdk/model";
import {
  PlAccordion,
  PlAccordionSection,
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlBtnGroup,
  PlDatasetSelector,
  PlDropdown,
  PlElementList,
  PlIcon24,
  PlNumberField,
  PlSlideModal,
  PlTooltip,
  PlValueLabel,
  usePlColumnValueSettings,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, ref, watch, watchEffect } from "vue";
import { useApp } from "../app";
import { useTimepointOptions } from "../timepoints";

const app = useApp();

// The selector picks the abundance, or the abundance narrowed by one of its subset columns.
const abundanceSelection = computed<DatasetSelection | undefined>({
  get: () => toAbundanceSelection(app.model.data),
  set: (selection) => {
    const { abundanceRef, filterRef } = fromAbundanceSelection(selection);
    app.model.data.abundanceRef = abundanceRef;
    app.model.data.filterRef = filterRef;
  },
});

const settingsAreShown = ref(app.model.outputs.datasetSpec === undefined);
const showSettings = () => {
  settingsAreShown.value = true;
};

const calculationModeOptions = [
  { label: "Population-Level", value: "population" },
  { label: "Intra-Subject", value: "intra-subject" },
];

const normalizationOptions = [
  { label: "Relative Frequency", value: "relative-frequency" },
  { label: "CLR Transform", value: "clr" },
];

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.mainTable,
});

type MetadataColumnField = "subjectColumnRef" | "groupingColumnRef" | "temporalColumnRef";

function setMetadataColumn(field: MetadataColumnField, ref: ColumnUniversalId | undefined) {
  app.model.data[field] = ref;
  refreshSubtitleLabel();
}

// The subtitle names the chosen grouping, temporal and subject columns as the result pool labels
// them. Called only from the pickers' setters, so reopening the block never rewrites it.
function refreshSubtitleLabel() {
  const metadataOptions = app.model.outputs.metadataOptions ?? [];
  const { groupingColumnRef, temporalColumnRef, subjectColumnRef } = app.model.data;
  app.model.data.subtitleLabel = [groupingColumnRef, temporalColumnRef, subjectColumnRef]
    .filter((ref) => ref !== undefined)
    .map((ref) => metadataOptions.find((o) => o.value === ref)?.label)
    .filter(Boolean)
    .join(", ");
}

// Subject is required only in intra-subject mode
const subjectRequired = computed(() => app.model.data.calculationMode === "intra-subject");

/** The values of the temporal column picked now, as far as they are read. */
const timepointValuesState = useTimepointOptions(
  usePlColumnValueSettings({
    source: () => app.model.outputs.timepointOptions,
    valueColumn: () => app.model.data.temporalColumnRef,
  }),
);

/**
 * Values of the temporal column as options of the timepoints. `undefined` while they are read,
 * so the timepoints are shown unchecked until then. A column whose values are gone or cannot be
 * read offers none; a picked timepoint missing from them stays picked and is shown as not
 * available.
 */
const timepointOptions = computed(() => {
  const state = timepointValuesState.value;
  if (state.type === "none" || state.type === "loading") return undefined;
  return state.type === "read" ? state.options : [];
});

/** Why the values of the temporal column could not be read. */
const timepointValuesError = computed(() => {
  const state = timepointValuesState.value;
  return state.type === "failed" ? state.message : undefined;
});

/** Every value of the temporal column, once read. */
const timepointValues = computed(() => {
  const state = timepointValuesState.value;
  return state.type === "read" ? state.options.map((o) => o.value) : [];
});

/**
 * Whether the timepoint order differs in membership from the values of the temporal column: some
 * value is left out, or some picked timepoint is no longer among them.
 */
const canResetTimepointOrder = computed(() => {
  const values = new Set(timepointValues.value);
  const picked = new Set(app.model.data.timepointOrder);
  if (values.size === 0) return false;
  return [...values].some((v) => !picked.has(v)) || [...picked].some((v) => !values.has(v));
});

const resetTimepointOrder = () => {
  if (timepointValues.value.length) {
    app.model.data.timepointOrder = [...timepointValues.value];
  }
};

// The timepoints are values of the temporal column, so picking another column (or none)
// invalidates them outright. Only a change of the column clears them: timepoints persisted with
// the block survive a remount of the page.
watch(
  () => app.model.data.temporalColumnRef,
  () => {
    app.model.data.timepointOrder = [];
  },
);

// An empty timepoint order takes every value of the temporal column, in natural order, once
// they are read.
watchEffect(() => {
  const vals = timepointValues.value;
  if (vals.length > 0 && app.model.data.timepointOrder.length === 0) {
    app.model.data.timepointOrder = [...vals];
  }
});

const isTimepointOrderOpen = ref(true);
const isAdvancedOpen = ref(false);
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="defaultSubtitle(app.model.data)"
    title="Clonotype Distribution"
  >
    <template #append>
      <PlBtnGhost @click.stop="showSettings">
        Settings
        <template #append>
          <PlIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlAgDataTableV2
      v-model="app.model.data.tableState"
      :settings="tableSettings"
      not-ready-text="Data is not computed"
      show-export-button
      show-search-field
    />
  </PlBlockPage>

  <PlSlideModal v-model="settingsAreShown">
    <template #title>Settings</template>

    <PlDatasetSelector
      v-model="abundanceSelection"
      :options="app.model.outputs.abundanceOptions"
      label="Select abundance"
      clearable
      required
    />

    <PlBtnGroup
      v-model="app.model.data.calculationMode"
      :options="calculationModeOptions"
      label="Calculation mode"
    />

    <PlDropdown
      :model-value="app.model.data.subjectColumnRef"
      :options="app.model.outputs.metadataOptions"
      label="Subject variable"
      :required="subjectRequired"
      clearable
      @update:model-value="(ref) => setMetadataColumn('subjectColumnRef', ref)"
    >
      <template #tooltip>
        Required for Intra-Subject mode. Optional in Population-Level mode — when omitted, all
        samples are pooled and cross-subject metrics are skipped.
      </template>
    </PlDropdown>

    <PlDropdown
      :model-value="app.model.data.groupingColumnRef"
      :options="app.model.outputs.metadataOptions"
      label="Grouping variable"
      clearable
      @update:model-value="(ref) => setMetadataColumn('groupingColumnRef', ref)"
    >
      <template #tooltip>
        Categorical metadata column (e.g. tissue, organ) used to compute distribution metrics:
        restriction index, dominant group, and breadth.
      </template>
    </PlDropdown>

    <PlDropdown
      :model-value="app.model.data.temporalColumnRef"
      :options="app.model.outputs.metadataOptions"
      label="Temporal variable"
      clearable
      @update:model-value="(ref) => setMetadataColumn('temporalColumnRef', ref)"
    >
      <template #tooltip>
        Ordered metadata column (e.g. timepoint, day) used to compute temporal expansion metrics:
        Log2 Peak Delta, Temporal Shift Index, and Log2 Kinetic Delta.
      </template>
    </PlDropdown>

    <PlAccordion multiple>
      <PlAccordionSection
        v-if="app.model.data.temporalColumnRef"
        v-model="isTimepointOrderOpen"
        label="Timepoint order"
      >
        <div style="display: flex; margin-bottom: -15px">
          Define timepoint order
          <PlTooltip class="info">
            <template #label>Define timepoint order</template>
            <template #tooltip>
              <div>
                Drag to reorder timepoints chronologically. The order determines how temporal
                metrics (TSI, Log2 Peak Delta, Log2 Kinetic Delta) are computed. First = earliest,
                last = latest.
              </div>
            </template>
          </PlTooltip>
        </div>
        <PlAlert v-if="timepointValuesError !== undefined" type="error">
          {{ timepointValuesError }}
        </PlAlert>
        <PlElementList v-model:items="app.model.data.timepointOrder">
          <template #item-title="{ item }">
            <PlValueLabel :value="item" :options="timepointOptions" />
          </template>
        </PlElementList>
        <PlBtnGhost v-if="canResetTimepointOrder" @click="resetTimepointOrder">
          Reset to default
          <template #append>
            <PlIcon24 name="reverse" />
          </template>
        </PlBtnGhost>
      </PlAccordionSection>

      <PlAccordionSection v-model="isAdvancedOpen" label="Advanced settings">
        <PlBtnGroup
          v-model="app.model.data.normalization"
          :options="normalizationOptions"
          label="Normalization"
        />

        <PlNumberField
          v-model="app.model.data.minAbundanceThreshold"
          label="Min abundance threshold"
          :min-value="0"
          :step="1"
        >
          <template #tooltip>
            Filter clones with abundance below this value in ALL samples before computation. Default
            0 includes everything.
          </template>
        </PlNumberField>

        <PlNumberField
          v-model="app.model.data.minSubjectCount"
          label="Min subject count"
          :min-value="1"
          :step="1"
        >
          <template #tooltip>
            Averaged cross-subject metrics (Mean RI, Mean Log2PD, etc.) are set to NaN when a clone
            is present in fewer subjects than this threshold. Default: 1 (no filtering).
          </template>
        </PlNumberField>

        <PlNumberField
          v-model="app.model.data.topN"
          label="Top N clones (temporal plot)"
          :min-value="1"
          :step="1"
        >
          <template #tooltip>
            Number of top clones to show in the temporal trajectory plot, ranked by absolute Log2
            Peak Delta. Default: 20.
          </template>
        </PlNumberField>

        <PlNumberField
          v-model="app.model.data.presenceThreshold"
          label="Presence threshold"
          :min-value="0"
          :max-value="1"
          :step="0.0001"
        >
          <template #tooltip>
            Minimum frequency for a clone to be considered present in a group. Default 0 means any
            detection counts.
          </template>
        </PlNumberField>
      </PlAccordionSection>
    </PlAccordion>
  </PlSlideModal>
</template>
