import { Cpu, Activity, ThermometerSun, FlaskConical } from 'lucide-react';

const iconMap = {
  flow1_antibody: ThermometerSun,
  flow2_dnaprobe: Activity,
  flow3_enzyme: FlaskConical,
};

const FlowPicker = ({ id, flows = [], activeFlowId, onSelect }) => {
  return (
    <div id={id} className="flex gap-3">
      {flows.map((flow) => {
        const Icon = iconMap[flow.id] || Cpu;
        const active = flow.id === activeFlowId;
        return (
          <button
            key={flow.id}
            onClick={() => onSelect(flow.id)}
            className={`group flex items-center gap-2 rounded-full border px-3 py-[6px] text-xs font-medium transition ${
              active
                ? 'border-orange-300 bg-orange-50 text-orange-700 shadow-sm'
                : 'border-orange-200 bg-white text-orange-600 hover:bg-orange-50'
            }`}
          >
            <Icon size={16} className="text-orange-500" />
            <span>{flow.name}</span>
          </button>
        );
      })}
    </div>
  );
};

export default FlowPicker;
