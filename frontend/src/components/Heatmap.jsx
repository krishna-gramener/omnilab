import { useEffect, useMemo, useRef, useState } from 'react';
import {
  ResponsiveContainer,
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  ZAxis,
  Tooltip,
} from 'recharts';
import { Play, Pause, Info } from 'lucide-react';

const COLOR_STOPS = [
  { stop: 0, color: [127, 29, 29] },
  { stop: 0.2, color: [220, 38, 38] },
  { stop: 0.4, color: [249, 115, 22] },
  { stop: 0.6, color: [251, 191, 36] },
  { stop: 0.72, color: [217, 119, 6] },
  { stop: 0.82, color: [132, 204, 22] },
  { stop: 0.9, color: [34, 197, 94] },
  { stop: 0.97, color: [21, 128, 61] },
  { stop: 1, color: [6, 95, 70] },
];

const clamp = (value, min = 0, max = 1) => Math.min(max, Math.max(min, value));

const interpolateColor = (value) => {
  const clampedValue = clamp(value);
  if (clampedValue <= COLOR_STOPS[0].stop) {
    return COLOR_STOPS[0].color;
  }
  if (clampedValue >= COLOR_STOPS[COLOR_STOPS.length - 1].stop) {
    return COLOR_STOPS[COLOR_STOPS.length - 1].color;
  }

  for (let i = 1; i < COLOR_STOPS.length; i += 1) {
    const current = COLOR_STOPS[i];
    const previous = COLOR_STOPS[i - 1];
    if (clampedValue <= current.stop) {
      const range = current.stop - previous.stop || 1;
      const ratio = (clampedValue - previous.stop) / range;
      const r = Math.round(previous.color[0] + (current.color[0] - previous.color[0]) * ratio);
      const g = Math.round(previous.color[1] + (current.color[1] - previous.color[1]) * ratio);
      const b = Math.round(previous.color[2] + (current.color[2] - previous.color[2]) * ratio);
      return [r, g, b];
    }
  }
  return COLOR_STOPS[COLOR_STOPS.length - 1].color;
};

const rgbArrayToString = ([r, g, b]) => `rgb(${r}, ${g}, ${b})`;

const valueToColor = (value) => rgbArrayToString(interpolateColor(value));

const getContrastingTextColor = (rgbString) => {
  const match = rgbString.match(/rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)/);
  if (!match) return '#e2e8f0';

  const [, rRaw, gRaw, bRaw] = match.map(Number);
  const toLinear = (channel) => {
    const c = channel / 255;
    return c <= 0.03928 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  };

  const luminance = 0.2126 * toLinear(rRaw) + 0.7152 * toLinear(gRaw) + 0.0722 * toLinear(bRaw);
  return luminance > 0.42 ? '#0f172a' : '#f8fafc';
};

const HeatTooltip = ({ active, payload }) => {
  if (!active || !payload?.length) return null;
  const cell = payload[0].payload;
  const trace = cell.trace;
  const tempNet =
    (trace?.tempBonus ?? 0) - (trace?.tempPenalty ?? 0);
  const phNet = (trace?.phBonus ?? 0) - (trace?.phPenalty ?? 0);
  const clampNet = -(trace?.clampLoss ?? 0);
  const modifier = trace?.modifierBoost ?? 0;
  const scenarioBias = trace?.scenarioBias ?? 0;
  return (
    <div className="rounded-lg border border-orange-200 bg-white p-3 text-sm text-orange-700 shadow-md">
      <div className="text-xs font-medium text-orange-600/80">Cell insight</div>
      <div className="mt-1 text-sm text-orange-700">
        Stability at {cell.temp}°C, pH {cell.ph}
      </div>
      <div className="text-base font-semibold text-orange-700">Stability {cell.value.toFixed(2)}</div>
      {cell.delta !== null && (
        <div className="mt-1 text-xs text-orange-600/80">
          Reference: {cell.reference.toFixed(2)}
          {cell.delta !== null && (
            <span className={`ml-2 font-medium ${cell.delta >= 0 ? 'text-emerald-600' : 'text-rose-600'}`}>
              Δ {cell.delta >= 0 ? '+' : ''}
              {cell.delta}
            </span>
          )}
        </div>
      )}
      <svg width="120" height="40" className="mt-2">
        <polyline
          fill="none"
          stroke="url(#spark-grad)"
          strokeWidth="2"
          points={cell.sparkPoints}
        />
        <defs>
          <linearGradient id="spark-grad" x1="0" x2="1" y1="0" y2="0">
            <stop offset="0%" stopColor="#2563eb" />
            <stop offset="100%" stopColor="#10b981" />
          </linearGradient>
        </defs>
      </svg>
      {trace && (
        <div className="mt-3 space-y-1 text-xs text-slate-600">
          <div className="text-xs font-medium text-slate-500">Local breakdown</div>
          <div className="flex items-center justify-between font-mono">
            <span>Base grid</span>
            <span>{trace.base?.toFixed(2)}</span>
          </div>
          <div className="flex items-center justify-between font-mono">
            <span>Temp Δ</span>
            <span className={tempNet >= 0 ? 'text-emerald-600' : 'text-rose-600'}>
              {tempNet >= 0 ? '+' : ''}
              {tempNet.toFixed(2)}
            </span>
          </div>
          <div className="flex items-center justify-between font-mono">
            <span>pH Δ</span>
            <span className={phNet >= 0 ? 'text-emerald-600' : 'text-rose-600'}>
              {phNet >= 0 ? '+' : ''}
              {phNet.toFixed(2)}
            </span>
          </div>
          {modifier !== 0 && (
            <div className="flex items-center justify-between font-mono">
              <span>Modifier</span>
              <span className={modifier >= 0 ? 'text-emerald-600' : 'text-rose-600'}>
                {modifier >= 0 ? '+' : ''}
                {modifier.toFixed(2)}
              </span>
            </div>
          )}
          {scenarioBias !== 0 && (
            <div className="flex items-center justify-between font-mono">
              <span>Scenario focus</span>
              <span className={scenarioBias >= 0 ? 'text-emerald-600' : 'text-rose-600'}>
                {scenarioBias >= 0 ? '+' : ''}
                {scenarioBias.toFixed(2)}
              </span>
            </div>
          )}
          {clampNet !== 0 && (
            <div className="flex items-center justify-between font-mono">
              <span>Clamp</span>
              <span className={clampNet >= 0 ? 'text-emerald-600' : 'text-rose-600'}>
                {clampNet >= 0 ? '+' : ''}
                {clampNet.toFixed(2)}
              </span>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

const Heatmap = ({
  heatmap,
  candidates,
  candidateId,
  onCandidateChange,
  playbackTempIndex,
  onPlaybackChange,
  playing,
  onTogglePlay,
  tempRange,
  phRange,
  whatIf,
  valueMap,
  compareReference = false,
  reference = null,
}) => {
  const containerRef = useRef(null);
  const [dimensions, setDimensions] = useState({ width: 0, height: 0 });
  const [legendOpen, setLegendOpen] = useState(true);
  const legendGradient = useMemo(
    () =>
      `linear-gradient(to right, ${COLOR_STOPS.map(
        (stop) => `${rgbArrayToString(stop.color)} ${stop.stop * 100}%`
      ).join(', ')})`,
    []
  );
  const legendTicks = useMemo(() => [0, 0.3, 0.5, 0.7, 0.85, 1], []);

  useEffect(() => {
    if (!containerRef.current) return;
    const observer = new ResizeObserver((entries) => {
      const entry = entries[0];
      setDimensions({
        width: entry.contentRect.width,
        height: entry.contentRect.height,
      });
    });
    observer.observe(containerRef.current);
    return () => observer.disconnect();
  }, []);

  const candidate = useMemo(
    () => heatmap?.candidates?.find((c) => c.id === candidateId) ?? heatmap?.candidates?.[0],
    [heatmap, candidateId]
  );

  const cells = useMemo(() => {
    if (!candidate || !heatmap || dimensions.width === 0 || dimensions.height === 0) return [];
    const temps = heatmap.temps;
    const phs = heatmap.phs;
    const cellWidth = dimensions.width / temps.length;
    const cellHeight = dimensions.height / phs.length;

    return temps.flatMap((temp, tIndex) =>
      phs.map((ph, pIndex) => {
        const rawValue = candidate.values?.[tIndex]?.[pIndex] ?? 0;
        const mapResult = valueMap ? valueMap(candidate.id, tIndex, pIndex, rawValue) : { value: rawValue };
        const value =
          typeof mapResult === 'number' ? mapResult : mapResult?.value ?? rawValue;
        const trace = typeof mapResult === 'object' ? mapResult?.trace ?? null : null;
        const inWindow = tIndex >= tempRange[0] && tIndex <= tempRange[1] && pIndex >= phRange[0] && pIndex <= phRange[1];
        let delta = null;
        let referenceValue = null;
        if (compareReference && reference?.values?.[tIndex]?.[pIndex] !== undefined) {
          referenceValue = reference.values[tIndex][pIndex];
          delta = Number((value - referenceValue).toFixed(2));
        }
        const fill = valueToColor(value);
        const textColor = getContrastingTextColor(fill);
        const sparkPoints = temps
          .map((t, idx) => {
            const raw = candidate.values?.[idx]?.[pIndex] ?? 0;
            const v = valueMap ? valueMap(candidate.id, idx, pIndex, raw) : raw;
            const x = (idx / (temps.length - 1 || 1)) * 120;
            const y = 35 - v * 30;
            return `${x},${y}`;
          })
          .join(' ');
        const playbackStroke = playbackTempIndex === tIndex;
        const baseStroke = 'rgba(15, 23, 42, 0.35)';
        const deltaStroke = delta !== null ? (delta >= 0 ? '#22c55e' : '#f97316') : baseStroke;
        const strokeColor = playbackStroke ? '#ff3df0' : deltaStroke;
        const strokeWidth = playbackStroke ? 2.5 : delta !== null ? 2 : 1;
        return {
          x: tIndex,
          y: pIndex,
          temp,
          ph,
          value,
          fill,
          textColor,
          opacity: whatIf ? (inWindow ? 1 : 0.25) : 1,
          stroke: strokeColor,
          strokeWidth,
          cellWidth,
          cellHeight,
          trace,
          inScenario: trace?.scenarioIncluded ?? inWindow,
          sparkPoints,
          delta,
          reference: referenceValue,
        };
      })
    );
  }, [candidate, heatmap, dimensions, tempRange, phRange, whatIf, playbackTempIndex, compareReference, reference, valueMap]);

  if (!candidate || !heatmap) {
    return (
      <div className="flex h-full items-center justify-center rounded-xl border border-orange-200 bg-orange-50 text-orange-600">
        No heatmap data loaded.
      </div>
    );
  }

  const temps = heatmap.temps;
  const phs = heatmap.phs;

  const CellShape = (props) => {
    const { cx, cy, payload } = props;
    const width = payload.cellWidth * 0.88;
    const height = payload.cellHeight * 0.88;
    const minDimension = Math.min(width, height);
    const showValue = minDimension >= 28;
    const fontSize = Math.min(14, Math.max(11, minDimension * 0.4));
    return (
      <g>
        <rect
          x={cx - width / 2}
          y={cy - height / 2}
          width={width}
          height={height}
          fill={payload.fill}
          opacity={payload.opacity}
          stroke={payload.stroke}
          strokeWidth={payload.strokeWidth}
          rx={4}
          ry={4}
          className="transition-all duration-200"
        />
        {showValue && (
          <text
            x={cx}
            y={cy}
            textAnchor="middle"
            dominantBaseline="middle"
            fontSize={fontSize}
            fontWeight="600"
            fill={payload.textColor}
            opacity={payload.opacity}
          >
            {payload.value.toFixed(2)}
          </text>
        )}
      </g>
    );
  };

  return (
    <div className="flex h-full min-h-0 flex-1 flex-col gap-3">
      <div className="flex flex-wrap items-center justify-between gap-4 px-2">
        <div className="flex items-center gap-2 text-sm font-medium text-orange-700">
          <span>Candidate</span>
          <select
            value={candidateId}
            onChange={(e) => onCandidateChange(e.target.value)}
            className="rounded-md border border-orange-200 bg-white px-2 py-1 text-orange-700 focus:border-orange-300 focus:outline-none"
          >
            {candidates.map((cand) => (
              <option key={cand.id} value={cand.id}>
                {cand.id} • {cand.name}
              </option>
            ))}
          </select>
        </div>
        <div className="flex items-center gap-3 text-sm text-orange-600/70">
          <span>{temps[playbackTempIndex]}°C sweep</span>
          <div className="flex items-center gap-2">
            <button
              onClick={onTogglePlay}
              className="rounded-full border border-orange-200 bg-white p-2 text-orange-600 hover:bg-orange-50"
              aria-label={playing ? 'Pause sweep' : 'Play sweep'}
            >
              {playing ? <Pause size={16} /> : <Play size={16} />}
            </button>
            <input
              type="range"
              min={0}
              max={temps.length - 1}
              step={1}
              value={playbackTempIndex}
              onChange={(e) => onPlaybackChange(Number(e.target.value))}
              className="w-28 accent-orange-500 mr-2"
            />
          </div>
        </div>
      </div>

      <div
        ref={containerRef}
        className="relative flex-1 min-h-0 overflow-hidden rounded-2xl border border-orange-200 bg-white px-4 py-4"
      >
        {dimensions.width > 0 && dimensions.height > 0 && (
          <ResponsiveContainer width="100%" height="100%">
            <ScatterChart margin={{ top: 20, right: 20, left: 20, bottom: 20 }}>
              <XAxis
                type="number"
                dataKey="x"
                name="Temp"
                domain={[-0.5, temps.length - 0.5]}
                ticks={temps.map((_, idx) => idx)}
                tickFormatter={(idx) => `${temps[idx]}°`}
                tick={{ fill: '#64748b', fontSize: 12 }}
                axisLine={{ stroke: '#cbd5f5' }}
                tickLine={{ stroke: '#cbd5f5' }}
              />
              <YAxis
                type="number"
                dataKey="y"
                name="pH"
                domain={[-0.5, phs.length - 0.5]}
                ticks={phs.map((_, idx) => idx)}
                tickFormatter={(idx) => `${phs[idx]}`}
                tick={{ fill: '#64748b', fontSize: 12 }}
                axisLine={{ stroke: '#cbd5f5' }}
                tickLine={{ stroke: '#cbd5f5' }}
              />
              <ZAxis type="number" dataKey="value" range={[0, 1]} />
              <Tooltip content={<HeatTooltip />} cursor={false} />
              <Scatter data={cells} shape={CellShape} />
            </ScatterChart>
          </ResponsiveContainer>
        )}
        <div className="pointer-events-auto absolute right-3 top-3 flex flex-col items-end gap-2">
          <button
            type="button"
            onClick={() => setLegendOpen((prev) => !prev)}
            aria-expanded={legendOpen}
            className="flex items-center gap-2 rounded-full border border-orange-200 bg-white px-3 py-1 text-xs text-orange-700 shadow-sm transition hover:bg-orange-50 focus:outline-none focus:ring-2 focus:ring-orange-200"
          >
            <Info size={12} className="text-orange-500" />
            <span>{legendOpen ? 'Hide scale' : 'Show scale'}</span>
          </button>
          {legendOpen && (
            <div className="flex w-48 flex-col gap-2 rounded-xl border border-orange-200 bg-white p-3 shadow-lg">
              <div className="flex items-center justify-between text-xs font-medium text-orange-600/80">
                <span>Stability scale</span>
                <Info
                  size={12}
                  className="text-orange-500"
                  title="Score blends normalized assay replicates with scenario adjustments and caps at 1.0."
                />
              </div>
              <div className="h-2 rounded-full" style={{ background: legendGradient }} />
              <div className="flex justify-between text-xs text-orange-600/70">
                <span>Low</span>
                <span>High</span>
              </div>
              <div className="flex justify-between text-xs font-medium text-orange-700">
                {legendTicks.map((tick) => (
                  <span key={tick}>{tick.toFixed(1)}</span>
                ))}
              </div>
              <div className="flex items-start gap-2 text-xs leading-5 text-orange-600/80">
                <Info size={11} className="mt-[2px] flex-none text-orange-500" aria-hidden="true" />
                <p>Higher scores indicate stronger retention across the selected temperature and pH window.</p>
              </div>
            </div>
          )}
        </div>
        <div className="pointer-events-none absolute inset-0 border border-orange-200/60" />
      </div>
    </div>
  );
};

export default Heatmap;
