import { platforma } from "@platforma-open/milaboratories.spatiotemporal-analysis.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import { computed, watch } from "vue";
import MainPage from "./pages/MainPage.vue";
import HeatmapPage from "./pages/HeatmapPage.vue";
import TemporalLinePage from "./pages/TemporalLinePage.vue";
import PrevalenceHistogramPage from "./pages/PrevalenceHistogramPage.vue";

export const sdkPlugin = defineAppV3(platforma, (app) => {
  app.model.data.customBlockLabel ??= "";

  syncSubtitleLabel(app.model);

  return {
    progress: () => {
      return app.model.outputs.isRunning;
    },
    routes: {
      "/": () => MainPage,
      "/heatmap": () => HeatmapPage,
      "/temporal": () => TemporalLinePage,
      "/prevalence": () => PrevalenceHistogramPage,
    },
  };
});

export const useApp = sdkPlugin.useApp;

type AppModel = ReturnType<typeof useApp>["model"];

// The subtitle names the chosen grouping, temporal and subject columns as the result pool labels
// them: undefined while the options are loading; a chosen column the loaded options no longer
// list is left out.
function syncSubtitleLabel(model: AppModel) {
  const subtitleLabel = computed<string | undefined>(() => {
    const { groupingColumnRef, temporalColumnRef, subjectColumnRef } = model.data;
    const refs = [groupingColumnRef, temporalColumnRef, subjectColumnRef].filter(
      (ref) => ref !== undefined,
    );
    if (refs.length === 0) return "";
    const metadataOptions = model.outputs.metadataOptions;
    if (metadataOptions === undefined) return undefined;
    return refs
      .map((ref) => metadataOptions.find((o) => o.value === ref)?.label)
      .filter(Boolean)
      .join(", ");
  });

  watch(
    subtitleLabel,
    (label) => {
      if (label !== undefined && label !== model.data.subtitleLabel) {
        model.data.subtitleLabel = label;
      }
    },
    { immediate: true },
  );
}
