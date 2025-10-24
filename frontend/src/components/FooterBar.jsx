import { SlidersHorizontal, ToggleLeft, ToggleRight } from 'lucide-react';

const SliderSection = ({ label, values, range, onChange }) => {
  const [start, end] = range;
  return (
    <div className="flex flex-col gap-2 rounded-lg border border-orange-200 bg-white p-3">
      <div className="flex items-center justify-between text-xs font-medium text-orange-600/80">
        <span>{label}</span>
        <span className="text-orange-700">{`${values[start]} – ${values[end]}`}</span>
      </div>
      <div className="flex items-center gap-2 text-[11px] text-orange-500">
        <span>Min</span>
        <input
          type="range"
          min={0}
          max={values.length - 1}
          step={1}
          value={start}
          onChange={(e) => {
            const nextStart = Number(e.target.value);
            if (nextStart <= end) {
              onChange([nextStart, end]);
            }
          }}
          className="w-full accent-orange-500"
        />
      </div>
      <div className="flex items-center gap-2 text-[11px] text-orange-500">
        <span>Max</span>
        <input
          type="range"
          min={0}
          max={values.length - 1}
          step={1}
          value={end}
          onChange={(e) => {
            const nextEnd = Number(e.target.value);
            if (nextEnd >= start) {
              onChange([start, nextEnd]);
            }
          }}
          className="w-full accent-orange-500"
        />
      </div>
    </div>
  );
};

const FooterBar = ({ temps, phs, tempRange, phRange, onTempRangeChange, onPhRangeChange, whatIf, onToggleWhatIf }) => {
  return (
    <footer id="whatif-panel" className="rounded-xl border border-orange-200 bg-white p-4 shadow-sm">
      <div className="mb-3 flex items-center gap-2 text-sm font-medium text-orange-700">
        <SlidersHorizontal size={16} />
        What-if window
      </div>
      <div className="grid gap-3 md:grid-cols-[1fr_1fr_auto]">
        <SliderSection label="Temperature (°C)" values={temps} range={tempRange} onChange={onTempRangeChange} />
        <SliderSection label="pH" values={phs} range={phRange} onChange={onPhRangeChange} />
        <button
          onClick={onToggleWhatIf}
          className={`flex flex-col items-center justify-center rounded-lg border px-4 py-6 text-xs font-medium transition ${
            whatIf
              ? 'border-orange-300 bg-orange-50 text-orange-700'
              : 'border-orange-200 bg-white text-orange-500 hover:bg-orange-100'
          }`}
        >
          {whatIf ? <ToggleRight size={32} /> : <ToggleLeft size={32} />}
          <span className="mt-2">{whatIf ? 'Scenario on' : 'Scenario off'}</span>
        </button>
      </div>
    </footer>
  );
};

export default FooterBar;
