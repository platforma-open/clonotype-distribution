<script setup lang="ts">
import { PlBlockPage } from '@platforma-sdk/ui-vue';
import { GraphMaker } from '@milaboratories/graph-maker';
import { useApp } from '../app';

const app = useApp();
</script>

<template>
  <PlBlockPage>
    <template #title>Compartment Heatmap</template>
    <GraphMaker
      v-model="app.model.ui.heatmapState"
      chartType="heatmap"
      :data-state-key="app.model.outputs.heatmapPf"
      :p-frame="app.model.outputs.heatmapPf"
    />
  </PlBlockPage>
</template>
