import * as Dialog from '@radix-ui/react-dialog';
import { useEffect, useMemo, useState } from 'react';
import { getRaw } from '../lib/api.js';
import { X } from 'lucide-react';

const tabs = [
  { id: 'candidates.csv', label: 'candidates.csv' },
  { id: 'scores.json', label: 'scores.json' },
  { id: 'heatmap.json', label: 'heatmap.json' },
  { id: 'rationale.json', label: 'rationale.json' },
];

const toCsv = (rows) => {
  if (!Array.isArray(rows) || rows.length === 0) return '';
  const headers = Object.keys(rows[0]);
  const lines = [headers.join(',')];
  for (const row of rows) {
    lines.push(headers.map((h) => row[h]).join(','));
  }
  return lines.join('\n');
};

const DataModal = ({ flowId, open, onOpenChange }) => {
  const [activeTab, setActiveTab] = useState(tabs[0].id);
  const [content, setContent] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  useEffect(() => {
    if (!open || !flowId) return;
    setLoading(true);
    setError('');
    getRaw(flowId, activeTab)
      .then((value) => {
        if (activeTab.endsWith('.csv')) {
          setContent(toCsv(value));
        } else {
          setContent(JSON.stringify(value, null, 2));
        }
      })
      .catch(() => setError('Unable to load dataset in sandbox mode.'))
      .finally(() => setLoading(false));
  }, [open, flowId, activeTab]);

  useEffect(() => {
    if (open) {
      setActiveTab(tabs[0].id);
      setContent('');
      setError('');
    }
  }, [open, flowId]);

  return (
    <Dialog.Root open={open} onOpenChange={onOpenChange}>
      <Dialog.Portal>
        <Dialog.Overlay className="fixed inset-0 bg-black/30 backdrop-blur-sm" />
        <Dialog.Content className="fixed left-1/2 top-1/2 z-[11000] w-[min(90vw,900px)] -translate-x-1/2 -translate-y-1/2 rounded-2xl border border-orange-200 bg-white p-6 shadow-xl">
          <div className="mb-4 flex items-start justify-between">
            <div>
              <Dialog.Title className="text-lg font-semibold text-orange-700">
                Data payload
              </Dialog.Title>
              <Dialog.Description className="text-sm text-orange-600/70">
                Flow dataset served directly from the FastAPI mock files.
              </Dialog.Description>
            </div>
            <Dialog.Close className="rounded-full border border-orange-200 p-1 text-orange-600 hover:bg-orange-50">
              <X size={16} />
            </Dialog.Close>
          </div>

          <div className="mb-4 flex gap-2">
            {tabs.map((tab) => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`rounded-lg border px-3 py-2 text-sm font-medium ${
                  activeTab === tab.id
                    ? 'border-orange-300 bg-orange-50 text-orange-700'
                    : 'border-orange-200 bg-white text-orange-500 hover:bg-orange-50'
                }`}
              >
                {tab.label}
              </button>
            ))}
          </div>

          <div className="max-h-[420px] overflow-auto rounded-xl border border-orange-200 bg-orange-50 p-4 font-mono text-xs text-orange-700">
            {loading && <div className="text-orange-600/80">Loading…</div>}
            {!loading && error && <div className="text-rose-600">{error}</div>}
            {!loading && !error && <pre className="whitespace-pre-wrap">{content}</pre>}
          </div>
          <p className="mt-3 text-xs text-orange-600/70">
            Demo data, non-PHI; suitable for method development. Validation required for clinical use.
          </p>
        </Dialog.Content>
      </Dialog.Portal>
    </Dialog.Root>
  );
};

export default DataModal;
