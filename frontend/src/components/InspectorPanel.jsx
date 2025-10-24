import {
  Thermometer,
  Droplets,
  Activity,
  ClipboardList,
  Lightbulb,
  ArrowRightLeft,
  FileCode,
} from 'lucide-react';

const formatPercent = (value) => `${Math.round(value * 100)}%`;
const formatDelta = (value) => `${value >= 0 ? '+' : ''}${value.toFixed(2)}`;
const safeFixed = (value, digits = 2) =>
  typeof value === 'number' && Number.isFinite(value) ? value.toFixed(digits) : '—';

const InspectorPanel = ({
  candidate,
  stats,
  previousStats,
  rationale,
  narrativeNote,
  scenarioActive,
  simulationDetails,
  breakdown,
  curveData = [],
  cellMatrix = [],
}) => {
  if (!candidate) {
    return (
      <div className="rounded-xl border border-orange-200 bg-white p-4 text-orange-600/80">
        Select a candidate to inspect stability metrics.
      </div>
    );
  }

  return (
    <div className="flex h-full flex-col gap-3">
      <div className="rounded-xl border border-orange-200 bg-white p-4 shadow-sm">
        <div className="flex items-center justify-between text-sm font-medium text-orange-600/80">
          <span className="text-orange-500/80">{candidate.id}</span>
          {scenarioActive && (
            <span className="rounded-full border border-orange-200 bg-orange-50 px-2 py-[2px] text-xs text-orange-700">
              Scenario window
            </span>
          )}
        </div>
        <h2 className="mt-1 text-xl font-semibold text-orange-700">{candidate.name}</h2>
        <p className="mt-1 text-sm text-orange-600/70">{candidate.notes}</p>
      </div>

      <div className="flex-1 min-h-0 space-y-3 overflow-y-auto pr-1">
        <div className="grid grid-cols-2 gap-3 text-sm text-orange-600/80">
          <MetricCard
            icon={Thermometer}
            label="Max stable temp"
            value={`${stats.maxStableTemp}°C`}
            detail="highest temp >= 0.82"
          />
          <MetricCard
            icon={Droplets}
            label="Stable pH window"
            value={stats.phWindow}
            detail="pH range with >=0.82"
          />
          <MetricCard
            icon={Activity}
            label="Area above 0.82"
            value={formatPercent(stats.areaAbove)}
            detail="fraction of grid"
          />
          <MetricCard
            icon={ClipboardList}
            label="Scenario avg"
            value={stats.windowAvg.toFixed(2)}
            detail="current window mean"
          />
        </div>

        {breakdown?.totals?.count ? (
          <div className="rounded-xl border border-slate-200 bg-white p-4 text-sm text-slate-600">
            <div className="flex items-center justify-between text-xs font-medium text-slate-500">
              <span>Score breakdown</span>
              <span>{breakdown.totals.count} scenario cells</span>
            </div>
            <div className="mt-2 grid grid-cols-2 gap-3 font-mono">
              {['base', 'tempBonus', 'tempPenalty', 'phBonus', 'phPenalty', 'scenarioBias', 'modifierBoost', 'jitter', 'clampLoss'].map(
                (field) => {
                  const labelMap = {
                    base: 'Base grid',
                    tempBonus: 'Temp +',
                    tempPenalty: 'Temp -',
                    phBonus: 'pH +',
                    phPenalty: 'pH -',
                    scenarioBias: 'Scenario focus',
                    modifierBoost: 'Modifier',
                    jitter: 'Jitter',
                    clampLoss: 'Clamp',
                  };
                  const avg = (breakdown.totals[field] ?? 0) / breakdown.totals.count;
                  return (
                    <div key={field} className="space-y-[2px]">
                      <div className="flex items-center justify-between text-xs text-slate-500">
                        <span>{labelMap[field] ?? field}</span>
                        <span
                          className={
                            avg > 0
                              ? 'text-emerald-600'
                              : avg < 0
                              ? 'text-rose-600'
                              : 'text-slate-500'
                          }
                        >
                          {formatDelta(avg)}
                        </span>
                      </div>
                    </div>
                  );
                }
              )}
            </div>
            <div className="mt-3 rounded-lg border border-slate-200 bg-slate-50 px-3 py-2 font-mono text-sm text-slate-700">
              Scenario window score: {safeFixed(breakdown.totals.final / breakdown.totals.count)}
            </div>
          </div>
        ) : null}

        {curveData.length > 0 && (
          <div className="rounded-xl border border-slate-200 bg-white p-4 text-sm text-slate-600">
            <div className="flex items-center justify-between text-xs font-medium text-slate-500">
              <span>Stability vs. temperature</span>
            </div>
            <Sparkline data={curveData} />
            <div className="mt-2 flex items-center justify-between text-xs text-slate-500">
              <span>Raw mean</span>
              <span>Adjusted mean</span>
              <span>Scenario mean</span>
            </div>
          </div>
        )}

        {cellMatrix.length > 0 && (
          <div className="rounded-xl border border-slate-200 bg-white p-4 text-sm text-slate-600">
            <div className="flex items-center gap-2 text-xs font-medium text-slate-500">
              <ArrowRightLeft size={14} />
              Temp × pH breakdown
            </div>
            <div className="mt-3 space-y-3">
              {cellMatrix.map((row) => (
                <div key={row.temp} className="space-y-2">
                  <div className="text-xs font-medium text-slate-500">
                    {row.temp}°C
                  </div>
                  <div className="grid grid-cols-2 gap-2">
                    {row.cells.map((cell) => {
                      const delta = cell.adjusted - cell.raw;
                      const deltaClass =
                        delta > 0 ? 'text-emerald-600' : delta < 0 ? 'text-rose-600' : 'text-slate-500';
                      return (
                        <div
                          key={`${cell.tempIndex}:${cell.phIndex}`}
                          className={`rounded-lg border px-3 py-2 ${
                            cell.inScenario
                              ? 'border-blue-200 bg-blue-50'
                              : 'border-slate-200 bg-slate-50'
                          }`}
                        >
                          <div className="flex items-center justify-between text-xs text-slate-500">
                            <span>pH {cell.ph}</span>
                            {cell.inScenario && <span className="text-blue-600">Window</span>}
                          </div>
                          <div className="mt-1 text-base font-semibold text-slate-800">
                            {safeFixed(cell.adjusted)}
                          </div>
                          <div className={`text-xs font-mono ${deltaClass}`}>
                            {safeFixed(cell.raw)} → {safeFixed(cell.adjusted)} ({formatDelta(delta)})
                          </div>
                        </div>
                      );
                    })}
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {previousStats && (
          <div className="rounded-xl border border-slate-200 bg-white p-4 text-sm text-slate-600">
            <div className="flex items-center gap-2 text-xs font-medium text-slate-500">
              <ArrowRightLeft size={14} />
              What changed
            </div>
            <div className="mt-2 space-y-1 font-mono">
              <ChangeRow label="Scenario Avg" prev={previousStats.windowAvg} next={stats.windowAvg} />
              <ChangeRow label="Area Above" prev={previousStats.areaAbove} next={stats.areaAbove} formatter={(v) => `${Math.round(v * 100)}%`} />
            </div>
          </div>
        )}

        {rationale && (
          <div className="rounded-xl border border-slate-200 bg-white p-4 text-sm text-slate-600">
            <div className="flex items-center gap-2 font-semibold text-slate-700">
              <Lightbulb size={14} />
              Why in top-5
            </div>
            <p className="mt-2 text-slate-600">{rationale.why}</p>
            <p className="mt-1 text-slate-500">Watch: {rationale.watch}</p>
          </div>
        )}

        {simulationDetails && (
          <div className="rounded-xl border border-slate-200 bg-white p-4 text-sm text-slate-600">
            <div className="flex items-center gap-2 text-xs font-medium text-slate-500">
              <FileCode size={14} />
              Runtime chem simulate
            </div>
            <div className="mt-2 space-y-1 font-mono">
              <div>Dock score: {simulationDetails.dock_score}</div>
              <div>RMSD: {simulationDetails.rmsd} Å</div>
              <div>Hydrogen bonds: {simulationDetails.hbonds}</div>
            </div>
          </div>
        )}

        {narrativeNote && (
          <div className="rounded-xl border border-blue-200 bg-blue-50 p-3 text-sm text-blue-800">
            {narrativeNote}
          </div>
        )}
      </div>
    </div>
  );
};

const MetricCard = ({ icon: Icon, label, value, detail }) => (
  <div className="rounded-lg border border-slate-200 bg-white p-3">
    <div className="flex items-center gap-2 text-xs font-medium text-slate-500">
      <Icon size={14} className="text-slate-500" />
      {label}
    </div>
    <div className="mt-2 text-lg font-semibold text-slate-800">{value}</div>
    <div className="text-xs text-slate-400">{detail}</div>
  </div>
);

const ChangeRow = ({ label, prev, next, formatter = (v) => v.toFixed(2) }) => (
  <div className="flex items-center justify-between">
    <span>{label}</span>
    <span>
      {formatter(prev)} → <strong className="text-slate-700">{formatter(next)}</strong>
    </span>
  </div>
);

export default InspectorPanel;

const Sparkline = ({ data }) => {
  if (!data || data.length === 0) return null;
  const width = 260;
  const height = 80;
  const maxValue = Math.max(
    1,
    ...data.map((point) => Math.max(point.raw ?? 0, point.adjusted ?? 0, point.scenario ?? 0))
  );
  const minValue = Math.min(
    0,
    ...data.map((point) => Math.min(point.raw ?? 0, point.adjusted ?? 0, point.scenario ?? 0))
  );
  const range = Math.max(0.0001, maxValue - minValue);

  const toPoint = (value, index) => {
    const x = data.length > 1 ? (index / (data.length - 1)) * width : width / 2;
    const normalized = (value - minValue) / range;
    const y = height - normalized * height;
    return `${x},${y}`;
  };

  const adjustedPoints = data.map((point, index) => toPoint(point.adjusted, index)).join(' ');
  const rawPoints = data.map((point, index) => toPoint(point.raw, index)).join(' ');
  const scenarioPoints = data
    .map((point, index) =>
      point.scenario !== null && point.scenario !== undefined
        ? toPoint(point.scenario, index)
        : null
    )
    .filter(Boolean)
    .join(' ');

  return (
    <svg viewBox={`0 0 ${width} ${height}`} className="mt-3 w-full">
      <polyline
        fill="none"
        stroke="#94a3b8"
        strokeWidth="2"
        strokeDasharray="4 3"
        points={rawPoints}
      />
      <polyline fill="none" stroke="#2563eb" strokeWidth="2.5" points={adjustedPoints} />
      {scenarioPoints && (
        <polyline
          fill="none"
          stroke="#10b981"
          strokeWidth="2.5"
          strokeDasharray="5 4"
          points={scenarioPoints}
        />
      )}
    </svg>
  );
};
