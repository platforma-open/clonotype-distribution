<script setup lang="ts">
import { PlBlockPage } from '@platforma-sdk/ui-vue';
import { GraphMaker } from '@milaboratories/graph-maker';
import { useApp } from '../app';

const app = useApp();
</script>

<template>
  <PlBlockPage>
    <template #title>Temporal Frequency Trajectory</template>
    <GraphMaker
      v-model="app.model.ui.temporalLineState"
      chartType="scatterplot"
      :data-state-key="app.model.outputs.temporalLinePf"
      :p-frame="app.model.outputs.temporalLinePf"
    />
  </PlBlockPage>
</template>
