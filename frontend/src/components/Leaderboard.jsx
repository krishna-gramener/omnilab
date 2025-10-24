import { useEffect, useMemo, useState } from 'react';
import { ShieldCheck, AlertTriangle, Info } from 'lucide-react';

const flagColors = {
  shipping_friendly: 'bg-orange-50 border-orange-200 text-orange-700',
  wide_window: 'bg-amber-50 border-amber-200 text-amber-700',
  buffer_c_ok: 'bg-orange-100 border-orange-200 text-orange-700',
  heat_tolerant: 'bg-orange-50 border-orange-200 text-orange-700',
  cold_chain_ready: 'bg-amber-50 border-amber-200 text-amber-700',
  mg2_plus: 'bg-orange-100 border-orange-200 text-orange-700',
  cold_boost: 'bg-amber-50 border-amber-200 text-amber-700',
};

const flagTooltips = {
  shipping_friendly: 'Performs during transit; low risk in shipping studies.',
  wide_window: 'Maintains activity across the broadest pH/temperature window.',
  buffer_c_ok: 'Compatible with Buffer C formulations.',
  heat_tolerant: 'Retains stability at elevated temperatures.',
  cold_chain_ready: 'Designed for cold-chain integrity checks.',
  mg2_plus: 'Demo effect +0.03 for Mg-compatible probes; actual magnitude varies by assay.',
  cold_boost: 'Enhanced low-temp performance via formulation tweak.',
  legacy_reference: 'Reference candidate for comparison.',
};

const formatPercent = (value) => `${Math.round(value * 100)}%`;
const formatDelta = (value) => `${value >= 0 ? '+' : ''}${value.toFixed(2)}`;

const Leaderboard = ({
  items = [],
  activeId,
  onSelect,
  scenarioActive,
  artifactStatus = [],
  breakdowns = new Map(),
}) => {
  const [showBreakdown, setShowBreakdown] = useState(false);
  const activeBreakdown = breakdowns?.get(activeId);
  useEffect(() => {
    if (!activeBreakdown) {
      setShowBreakdown(false);
    }
  }, [activeId, activeBreakdown]);

  const breakdownSummary = useMemo(() => {
    if (!activeBreakdown?.totals) return null;
    const { totals } = activeBreakdown;
    const count = totals.count || 0;
    if (count === 0) {
      return { count: 0, averages: {} };
    }
    const avg = (field) => (totals[field] ?? 0) / count;
    const averages = {
      base: avg('base'),
      tempBonus: avg('tempBonus'),
      tempPenalty: avg('tempPenalty'),
      phBonus: avg('phBonus'),
      phPenalty: avg('phPenalty'),
      scenarioBias: avg('scenarioBias'),
      modifierBoost: avg('modifierBoost'),
      jitter: avg('jitter'),
      clampLoss: avg('clampLoss'),
      final: avg('final'),
    };
    const nets = {
      tempNet: averages.tempBonus - averages.tempPenalty,
      phNet: averages.phBonus - averages.phPenalty,
      clampNet: -averages.clampLoss,
    };
    return { count, averages, nets };
  }, [activeBreakdown]);

  return (
    <div className="flex h-full flex-col gap-3">
      <div className="flex items-center justify-between">
        <h3 className="text-sm font-semibold text-orange-700">
          {scenarioActive ? 'Scenario Ranking' : 'Top-5 Robustness'}
        </h3>
        {scenarioActive && (
          <span className="rounded-full border border-orange-200 bg-orange-50 px-3 py-1 text-xs text-orange-700">
            Scenario window applied
          </span>
        )}
        {activeBreakdown && (
          <button
            type="button"
            onClick={() => setShowBreakdown((prev) => !prev)}
            className="flex items-center gap-1 rounded-full border border-orange-200 bg-white px-3 py-1 text-xs text-orange-700 hover:bg-orange-50"
          >
            <Info size={12} className="text-orange-500" />
            {showBreakdown ? 'Hide breakdown' : 'Score breakdown'}
          </button>
        )}
      </div>

      {artifactStatus.length > 0 && (
        <div className="flex flex-wrap items-center gap-2 text-xs text-orange-600/80">
          {artifactStatus.map((badge) => (
            <span
              key={badge.id}
              className="flex items-center gap-1 rounded-full border border-orange-200 bg-orange-50 px-2 py-[2px] text-xs text-orange-700"
              title={badge.tooltip}
            >
              <span>{badge.label}</span>
              <span>{badge.symbol}</span>
            </span>
          ))}
        </div>
      )}

      {showBreakdown && breakdownSummary && (
        <div className="rounded-lg border border-orange-200 bg-orange-50 p-3 text-sm text-orange-700">
          <div className="flex items-center justify-between text-xs font-medium text-orange-600/80">
            <span>Scenario breakdown</span>
            <span>{breakdownSummary.count} cells</span>
          </div>
          <div className="mt-2 grid grid-cols-2 gap-3 font-mono text-xs">
            <BreakdownLine
              label="Base grid"
              value={breakdownSummary.averages.base}
              variant="neutral"
            />
            <BreakdownLine
              label="Temp net"
              value={breakdownSummary.nets.tempNet}
              detail={`+${breakdownSummary.averages.tempBonus.toFixed(2)} / -${breakdownSummary.averages.tempPenalty.toFixed(2)}`}
            />
            <BreakdownLine
              label="pH net"
              value={breakdownSummary.nets.phNet}
              detail={`+${breakdownSummary.averages.phBonus.toFixed(2)} / -${breakdownSummary.averages.phPenalty.toFixed(2)}`}
            />
            <BreakdownLine
              label="Scenario focus"
              value={breakdownSummary.averages.scenarioBias}
            />
            <BreakdownLine
              label="Modifier"
              value={breakdownSummary.averages.modifierBoost}
            />
            <BreakdownLine
              label="Clamp"
              value={breakdownSummary.nets.clampNet + breakdownSummary.averages.jitter}
              detail={`jitter ${formatDelta(breakdownSummary.averages.jitter)} / clamp ${formatDelta(
                breakdownSummary.nets.clampNet
              )}`}
            />
          </div>
          <div className="mt-3 flex items-center justify-between rounded-lg border border-orange-200 bg-white px-3 py-2 font-mono text-sm text-orange-700">
            <span>Scenario window score</span>
            <span>{breakdownSummary.averages.final.toFixed(2)}</span>
          </div>
        </div>
      )}

      <div className="mt-2 flex-1 min-h-0 flex flex-col gap-2 overflow-y-auto pr-1">
        {items.map((item, index) => {
          const active = item.id === activeId;
          const progress = scenarioActive ? item.scenarioScore : item.robustness;
          return (
            <button
              key={item.id}
              onClick={() => onSelect(item.id)}
              className={`rounded-xl border px-3 py-3 text-left transition ${
                active
                  ? 'border-blue-400 bg-blue-50 shadow-sm'
                  : 'border-slate-200 bg-white hover:border-blue-200'
              }`}
            >
              <div className="flex items-center justify-between text-sm font-medium text-orange-700">
                <span>
                  #{scenarioActive ? index + 1 : item.rank} • {item.name}
                </span>
                <div className="flex items-center gap-2">
                  {item.delta !== undefined && item.delta !== 0 && (
                    <span
                      className={`rounded-full px-2 py-[1px] text-xs font-medium ${
                        item.delta < 0 ? 'bg-rose-100 text-rose-700' : 'bg-emerald-100 text-emerald-700'
                      }`}
                    >
                      {item.delta > 0 ? '▲ ' : '▼ '}
                      {Math.abs(item.delta)}
                    </span>
                  )}
                  <span>{formatPercent(progress)}</span>
                </div>
              </div>
              <div className="mt-2 h-2 w-full overflow-hidden rounded-full bg-orange-100">
                <div
                  className="h-full rounded-full bg-orange-500"
                  style={{ width: `${Math.min(100, progress * 100)}%` }}
                />
              </div>
              <div className="mt-3 space-y-1 text-sm leading-snug text-orange-600/80">
                <div className="flex items-center gap-2 text-orange-600/80">
                  <ShieldCheck size={14} className="text-orange-500" />
                  <span className="truncate" title={item.why}>
                    {item.why}
                  </span>
                </div>
                <div className="flex items-center gap-2 text-orange-500">
                  <AlertTriangle size={14} className="text-orange-400" />
                  <span className="truncate" title={item.watch}>
                    {item.watch}
                  </span>
                </div>
              </div>
              {item.flags?.length > 0 && (
                <div className="mt-2 flex flex-wrap gap-2">
                  {item.flags.map((flag) => (
                    <span
                      key={flag}
                      className={`rounded-full border px-2 py-[2px] text-xs ${
                        flagColors[flag] || 'border-slate-200 bg-slate-100 text-slate-600'
                      }`}
                      title={flagTooltips[flag] || 'Chem heuristic flag for this candidate.'}
                    >
                      {flag.replace(/_/g, ' ')}
                    </span>
                  ))}
                </div>
              )}
            </button>
          );
        })}
      </div>
    </div>
  );
};

export default Leaderboard;

const BreakdownLine = ({ label, value, detail }) => {
  const variant = value > 0 ? 'pos' : value < 0 ? 'neg' : 'neutral';
  const colorClass =
    variant === 'pos'
      ? 'text-orange-700'
      : variant === 'neg'
      ? 'text-rose-600'
      : 'text-orange-600/80';
  return (
    <div className="space-y-[2px]">
      <div className="flex items-center justify-between">
        <span className="text-xs font-medium text-orange-600/80">{label}</span>
        <span className={`text-sm font-semibold ${colorClass}`}>{formatDelta(value)}</span>
      </div>
      {detail && <div className="text-[11px] text-orange-500/70">{detail}</div>}
    </div>
  );
};
