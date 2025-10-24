import { useState } from 'react';
import { ClipboardCheck, ChevronDown, ChevronUp } from 'lucide-react';

const steps = [
  'Confirm with n=3 replicates at T=45 °C, pH 7.6.',
  'Buffer: phosphate 20–50 mM.',
  'Acceptance: activity ≥0.85 vs baseline.',
  'Capture deviations in LIMS with scenario ID.',
  'Escalate anomalies to QA before release.',
];

const WetLabPanel = () => {
  const [open, setOpen] = useState(false);
  return (
    <div className="text-xs text-cyan-100">
      <button
        type="button"
        onClick={() => setOpen((prev) => !prev)}
        className="flex w-full items-center justify-between rounded-lg border border-teal-400/40 bg-teal-500/15 px-3 py-2 text-[11px] uppercase tracking-[0.2em] text-teal-100"
      >
        <span className="flex items-center gap-2">
          <ClipboardCheck size={14} />
          Wet-lab handoff checklist
        </span>
        {open ? <ChevronUp size={14} /> : <ChevronDown size={14} />}
      </button>
      {open && (
        <ol className="mt-3 space-y-1 text-cyan-100/80">
          {steps.map((step) => (
            <li key={step} className="leading-snug">
              {step}
            </li>
          ))}
        </ol>
      )}
    </div>
  );
};

export default WetLabPanel;
