import FlowPicker from './FlowPicker.jsx';
import {
  Database,
  BadgeCheck,
  LifeBuoy,
  PlayCircle,
  PauseCircle,
  Rocket,
  Cpu,
  Eye,
  EyeOff,
  Upload,
  FileSpreadsheet,
} from 'lucide-react';

const TopBar = ({
  flows,
  activeFlowId,
  onSelectFlow,
  onShowData,
  narrativeActions = [],
  activeNarrativeId,
  onNarrative,
  preflight,
  onStartTour,
  showTourCta,
  demoMode,
  onToggleDemoMode,
  demoSteps = [],
  onDemoStep,
  onShowHelp,
  flowBanner,
  successMetrics = [],
  macToggle = false,
  onToggleMac = () => {},
  showTourToggle,
  tourRunning,
  onToggleTour,
  children,
  onUploadScenario = () => {},
  onLoadSampleScenario = () => {},
}) => {
  const handleScenarioFile = (event) => {
    const file = event.target.files?.[0];
    if (file) {
      onUploadScenario(file);
      event.target.value = '';
    }
  };

  return (
    <header className="relative flex flex-col gap-3 rounded-xl border border-orange-200 bg-white p-4 shadow-sm">
      <div className="flex flex-wrap items-start justify-between gap-4">
        <div className="max-w-[60%]">
          <h1 className="text-lg font-semibold text-orange-700">
            Reagent Robustness Simulator
          </h1>
          <p className="text-sm text-orange-600/70">
            Compare candidate stability across temperature and pH windows.
          </p>
          {flowBanner && (
            <div className="mt-2 rounded-lg border border-orange-200 bg-orange-50 px-3 py-2 text-sm text-orange-700">
              {flowBanner}
            </div>
          )}
        </div>
        <div className="flex flex-wrap items-center gap-2">
          <button
            id="preflight-indicator"
            type="button"
            className="flex items-center gap-2 rounded-full border border-orange-200 bg-orange-50 px-3 py-1 text-xs font-medium text-orange-700"
            title={preflight?.summary ?? 'Chem preflight status'}
          >
            <BadgeCheck size={14} />
            {preflight?.summary
              ? `Azure ML Preflight: ${preflight.summary}`
              : 'Azure ML Preflight: Mocked ChemOps'}
          </button>
          <div className="flex items-center gap-2">
            <button
              type="button"
              onClick={onShowHelp}
              className="flex items-center gap-2 rounded-full border border-orange-200 bg-white px-3 py-1 text-xs font-medium text-orange-700 hover:bg-orange-50"
            >
              <LifeBuoy size={14} />
              How to use
            </button>
            <span className="rounded-full border border-orange-200 bg-orange-50 px-3 py-1 text-xs text-orange-600/70">
              Enable `USE_REAL_PROPka` / `USE_REAL_DOCKING` to run the full chem stack.
            </span>
          </div>
          <button
            id="show-data"
            type="button"
            onClick={onShowData}
            className="flex items-center gap-2 rounded-lg border border-orange-200 bg-white px-3 py-2 text-xs font-medium text-orange-700 hover:bg-orange-50"
          >
            <Database size={14} />
            Show Data
          </button>
          <label
            htmlFor="scenario-upload"
            className="flex cursor-pointer items-center gap-2 rounded-lg border border-orange-200 bg-white px-3 py-2 text-xs font-medium text-orange-700 hover:bg-orange-50"
          >
            <Upload size={14} />
            Upload Scenario
            <input
              id="scenario-upload"
              type="file"
              accept=".csv"
              onChange={handleScenarioFile}
              className="hidden"
            />
          </label>
          <button
            type="button"
            onClick={onLoadSampleScenario}
            className="flex items-center gap-2 rounded-lg border border-orange-200 bg-white px-3 py-2 text-xs font-medium text-orange-700 hover:bg-orange-50"
          >
            <FileSpreadsheet size={14} />
            Sample Scenario
          </button>
          {showTourCta && (
            <button
              type="button"
              onClick={onStartTour}
              className="flex items-center gap-2 rounded-full border border-orange-200 bg-white px-3 py-1 text-xs font-medium text-orange-700 hover:bg-orange-50"
            >
              <Rocket size={14} />
              Start tour
            </button>
          )}
        </div>
      </div>

      <div className="flex flex-wrap items-center justify-between gap-2">
        <div className="flex flex-wrap items-center gap-2">
          <div className="flex items-center gap-2 rounded-full border border-orange-200 bg-orange-50 px-3 py-2 text-xs font-medium text-orange-700">
            <span>Demo mode</span>
            <button
              type="button"
              onClick={onToggleDemoMode}
              className={`flex items-center gap-1 rounded-full border px-2 py-[2px] text-xs ${
                demoMode
                  ? 'border-orange-300 bg-white text-orange-700'
                  : 'border-orange-200 bg-orange-50 text-orange-500'
              }`}
            >
              {demoMode ? <PauseCircle size={12} /> : <PlayCircle size={12} />}
              {demoMode ? 'On' : 'Off'}
            </button>
            {demoMode && (
              <div className="flex items-center gap-1 text-xs">
                {demoSteps.map((step) => (
                  <button
                    type="button"
                    key={step.id}
                    onClick={() => onDemoStep(step.id)}
                    className="rounded-md border border-orange-200 bg-white px-2 py-[1px] text-xs text-orange-700 hover:bg-orange-50"
                  >
                    {step.label}
                  </button>
                ))}
              </div>
            )}
          </div>
          <FlowPicker id="flow-tabs" flows={flows} activeFlowId={activeFlowId} onSelect={onSelectFlow} />
        </div>
        <div id="narrative-buttons" className="flex flex-wrap items-center gap-2">
          {narrativeActions.map((action) => (
            <button
              key={action.id}
              type="button"
              onClick={() => onNarrative(action.id)}
              className={`rounded-lg border px-3 py-2 text-xs font-medium transition ${
                activeNarrativeId === action.id
                  ? 'border-orange-300 bg-orange-50 text-orange-700'
                  : 'border-transparent bg-white text-orange-500 hover:border-orange-200 hover:text-orange-700'
              }`}
              title={action.description}
            >
              {action.label}
            </button>
          ))}
        </div>
        <div className="flex flex-wrap items-center gap-2">
          {successMetrics.map((metric) => (
            <div
              key={metric.label}
              className="rounded-lg border border-orange-200 bg-white px-3 py-[6px] text-xs text-orange-700"
              title={metric.tooltip}
            >
              <div className="text-[11px] font-medium text-orange-600/70">{metric.label}</div>
              <div className="text-sm font-semibold text-orange-700">{metric.value}</div>
            </div>
          ))}
          <button
            type="button"
            onClick={onToggleMac}
            className={`flex items-center gap-1 rounded-full border px-3 py-[6px] text-xs font-medium transition ${
              macToggle ? 'border-orange-300 bg-white text-orange-700' : 'border-orange-200 bg-orange-50 text-orange-500 hover:bg-orange-100'
            }`}
            title="Toggle full chem stack to estimate Azure MAC usage"
          >
            <Cpu size={14} />
            MAC
          </button>
          {showTourToggle && (
            <button
              type="button"
              onClick={onToggleTour}
              className={`flex items-center gap-1 rounded-full border px-3 py-[6px] text-xs font-medium transition ${
                tourRunning
                  ? 'border-orange-300 bg-white text-orange-700'
                  : 'border-orange-200 bg-orange-50 text-orange-500 hover:bg-orange-100'
              }`}
              title={tourRunning ? 'Stop guided tour' : 'Start guided tour'}
            >
              {tourRunning ? <EyeOff size={14} /> : <Eye size={14} />}
              Tour
            </button>
          )}
        </div>
        <div className="text-xs text-orange-600/70">
          {macToggle
            ? 'Full chem stack assumes docking and PROPKA integrations are enabled.'
            : 'Demo mode uses mocked chem operations for instant interactions.'}
        </div>
      </div>
      {children && <div className="mt-3 w-full">{children}</div>}
    </header>
  );
};

export default TopBar;
