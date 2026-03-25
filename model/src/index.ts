import type { GraphMakerState } from '@milaboratories/graph-maker';
import type {
  InferOutputsType,
  PFrameHandle,
  PlDataTableStateV2,
  PlRef,
  SUniversalPColumnId,
} from '@platforma-sdk/model';
import {
  BlockModel,
  createPFrameForGraphs,
  createPlDataTableStateV2,
  createPlDataTableV2,
} from '@platforma-sdk/model';

export type CompartmentVariableConfig = {
  columnRef: SUniversalPColumnId;
  label: string;
};

export type TemporalVariableConfig = {
  columnRef?: SUniversalPColumnId;
  timepointOrder: string[];
};

export type BlockArgs = {
  defaultBlockLabel: string;
  customBlockLabel: string;

  // Input
  abundanceRef?: PlRef;

  // Mode
  calculationMode: 'population' | 'intra-subject';

  // Variable assignments
  compartmentVariables: CompartmentVariableConfig[];
  temporalVariable: TemporalVariableConfig;
  subjectVariable?: SUniversalPColumnId;

  // Normalization
  normalization: 'relative-frequency' | 'clr';

  // Thresholds
  presenceThreshold: number;
  pseudoCount: number;
};

export type UiState = {
  tableState: PlDataTableStateV2;
  heatmapState: GraphMakerState;
  temporalLineState: GraphMakerState;
  prevalenceHistogramState: GraphMakerState;
  settingsOpen: boolean;
};

export function getDefaultBlockArgs(): BlockArgs {
  return {
    defaultBlockLabel: '',
    customBlockLabel: '',
    calculationMode: 'population',
    compartmentVariables: [],
    temporalVariable: { timepointOrder: [] },
    normalization: 'relative-frequency',
    presenceThreshold: 0,
    pseudoCount: 1,
  };
}

export const model = BlockModel.create()

  .withArgs<BlockArgs>(getDefaultBlockArgs())

  .withUiState<UiState>({
    tableState: createPlDataTableStateV2(),
    heatmapState: {
      title: 'Compartment heatmap',
      template: 'heatmap',
      currentTab: null,
    },
    temporalLineState: {
      title: 'Temporal frequency trajectory',
      template: 'curve_dots',
      currentTab: null,
      layersSettings: {
        curve: {
          smoothing: false,
        },
      },
    },
    prevalenceHistogramState: {
      title: 'Subject prevalence distribution',
      template: 'bins',
      currentTab: null,
      layersSettings: {
        bins: { fillColor: '#99e099' },
      },
      axesSettings: {
        axisY: {
          axisLabelsAngle: 90,
          scale: 'log',
        },
      },
    },
    settingsOpen: true,
  })

  .argsValid((ctx) => {
    const { abundanceRef, compartmentVariables, temporalVariable, subjectVariable } = ctx.args;
    if (abundanceRef === undefined) return false;
    // At least one variable must be configured
    const hasCompartment = compartmentVariables.length > 0
      && compartmentVariables.every((v) => v.columnRef !== undefined);
    const hasTemporal = temporalVariable.columnRef !== undefined
      && temporalVariable.timepointOrder.length >= 2;
    if (!hasCompartment && !hasTemporal) return false;
    // Subject variable required for both modes
    if (subjectVariable === undefined) return false;
    return true;
  })

  // Abundance column options (same pattern as enrichment block)
  .output('abundanceOptions', (ctx) =>
    ctx.resultPool.getOptions([{
      axes: [
        { name: 'pl7.app/sampleId' },
        {},
      ],
      annotations: {
        'pl7.app/isAbundance': 'true',
        'pl7.app/abundance/normalized': 'false',
        'pl7.app/abundance/isPrimary': 'true',
      },
    }], { includeNativeLabel: true }),
  )

  // Metadata column options (anchored to abundance sampleId axis)
  .output('metadataOptions', (ctx) => {
    const anchor = ctx.args.abundanceRef;
    if (anchor === undefined) return undefined;
    return ctx.resultPool.getCanonicalOptions({ main: anchor },
      [{
        axes: [
          { anchor: 'main', idx: 0 },
        ],
        name: 'pl7.app/metadata',
      }],
    );
  })

  // Main output table
  .outputWithStatus('mainTable', (ctx) => {
    const pCols = ctx.outputs?.resolve('mainPf')?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPlDataTableV2(ctx, pCols, ctx.uiState.tableState);
  })

  // Heatmap PFrame
  .outputWithStatus('heatmapPf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('heatmapPf')?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPFrameForGraphs(ctx, pCols);
  })

  // Temporal line PFrame
  .outputWithStatus('temporalLinePf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('temporalLinePf')?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPFrameForGraphs(ctx, pCols);
  })

  // Prevalence histogram PFrame
  .outputWithStatus('prevalenceHistogramPf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('prevalenceHistogramPf')?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPFrameForGraphs(ctx, pCols);
  })

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .title(() => 'In Vivo Compartment Analysis')

  .subtitle((ctx) => ctx.args.customBlockLabel || ctx.args.defaultBlockLabel)

  .sections((_ctx) => [
    { type: 'link', href: '/', label: 'Main' },
    { type: 'link', href: '/heatmap', label: 'Compartment Heatmap' },
    { type: 'link', href: '/temporal', label: 'Temporal Trajectory' },
    { type: 'link', href: '/prevalence', label: 'Subject Prevalence' },
  ])

  .done(2);

export type BlockOutputs = InferOutputsType<typeof model>;
