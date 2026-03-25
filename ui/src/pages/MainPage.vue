<script setup lang="ts">
import {
  PlBlockPage,
  PlAgDataTableV2,
  PlDropdownRef,
  PlBtnGhost,
  PlBtnGroup,
  PlNumberField,
  PlSlideModal,
  PlMaskIcon24,
  usePlDataTableSettingsV2,
} from '@platforma-sdk/ui-vue';
import { computed, ref } from 'vue';
import { useApp } from '../app';

const app = useApp();
const settingsOpen = ref(true);

const calculationModeOptions = [
  { label: 'Population-Level', value: 'population' },
  { label: 'Intra-Subject', value: 'intra-subject' },
];

const normalizationOptions = [
  { label: 'Relative Frequency', value: 'relative-frequency' },
  { label: 'CLR Transform', value: 'clr' },
];

const abundanceOptions = computed(() => app.model.outputs.abundanceOptions ?? []);
const metadataOptions = computed(() => app.model.outputs.metadataOptions ?? []);

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.mainTable,
});
</script>

<template>
  <PlBlockPage>
    <template #title>In Vivo Compartment Analysis</template>

    <template #append>
      <PlBtnGhost @click.stop="() => (settingsOpen = true)">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlAgDataTableV2
      v-model="app.model.ui.tableState"
      :settings="tableSettings"
    />

    <PlSlideModal v-model="settingsOpen" title="Settings">
      <div style="display: flex; flex-direction: column; gap: 16px; padding: 16px;">
        <PlDropdownRef
          v-model="app.model.args.abundanceRef"
          :options="abundanceOptions"
          label="Abundance column"
        />

        <PlBtnGroup
          v-model="app.model.args.calculationMode"
          :options="calculationModeOptions"
          label="Calculation mode"
        />

        <PlDropdownRef
          v-model="app.model.args.subjectVariable"
          :options="metadataOptions"
          label="Subject variable"
          clearable
        />

        <div>
          <h4 style="margin: 0 0 8px 0;">Compartment variables</h4>
          <div
            v-for="(compVar, idx) in app.model.args.compartmentVariables"
            :key="idx"
            style="display: flex; gap: 8px; align-items: center; margin-bottom: 8px;"
          >
            <PlDropdownRef
              v-model="app.model.args.compartmentVariables[idx].columnRef"
              :options="metadataOptions"
              :label="'Compartment ' + (idx + 1)"
              clearable
              style="flex: 1;"
            />
            <PlBtnGhost @click="app.model.args.compartmentVariables.splice(idx, 1)">
              Remove
            </PlBtnGhost>
          </div>
          <PlBtnGhost @click="app.model.args.compartmentVariables.push({ columnRef: undefined as any, label: '' })">
            Add compartment variable
          </PlBtnGhost>
        </div>

        <div>
          <h4 style="margin: 0 0 8px 0;">Temporal variable</h4>
          <PlDropdownRef
            v-model="app.model.args.temporalVariable.columnRef"
            :options="metadataOptions"
            label="Timepoint column"
            clearable
          />
        </div>

        <PlBtnGroup
          v-model="app.model.args.normalization"
          :options="normalizationOptions"
          label="Normalization"
        />

        <PlNumberField
          v-model="app.model.args.presenceThreshold"
          label="Presence threshold"
          :min="0"
          :max="1"
          :step="0.0001"
        />

        <PlNumberField
          v-model="app.model.args.pseudoCount"
          label="Pseudo-count (for Log2 Kinetic Delta)"
          :min="0"
          :max="100"
          :step="1"
        />
      </div>
    </PlSlideModal>
  </PlBlockPage>
</template>
