<script setup lang="ts">
import { PlBlockPage } from '@platforma-sdk/ui-vue';
import { GraphMaker } from '@milaboratories/graph-maker';
import { useApp } from '../app';

const app = useApp();
</script>

<template>
  <PlBlockPage>
    <template #title>Subject Prevalence Distribution</template>
    <GraphMaker
      v-model="app.model.ui.prevalenceHistogramState"
      chartType="histogram"
      :data-state-key="app.model.outputs.prevalenceHistogramPf"
      :p-frame="app.model.outputs.prevalenceHistogramPf"
    />
  </PlBlockPage>
</template>
