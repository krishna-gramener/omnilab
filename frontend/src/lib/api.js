const envBase = typeof import.meta !== 'undefined' && import.meta.env ? import.meta.env.VITE_API_BASE : undefined;

const resolveBaseUrl = () => {
  if (envBase) {
    return envBase.replace(/\/$/, '');
  }
  if (typeof window !== 'undefined' && window.location) {
    return `${window.location.origin}/api`;
  }
  return '/api';
};

const BASE = resolveBaseUrl();

export const listFlows = () => fetch(`${BASE}/flows`).then((r) => r.json());

export const getFlowData = async (id) => {
  const [cand, scores, heat, rat] = await Promise.all([
    fetch(`${BASE}/flows/${id}/candidates`).then((r) => r.json()),
    fetch(`${BASE}/flows/${id}/scores`).then((r) => r.json()),
    fetch(`${BASE}/flows/${id}/heatmap`).then((r) => r.json()),
    fetch(`${BASE}/flows/${id}/rationale`).then((r) => r.json()),
  ]);
  return { cand, scores, heat, rat };
};

export const getRaw = (id, name) =>
  fetch(`${BASE}/flows/${id}/raw/${name}`).then((r) => r.json());

export const getPreflight = (id) =>
  fetch(`${BASE}/flows/${id}/preflight`).then((r) => r.json());

export const simulateCell = ({ flow_id, candidate_id, ph, temp }) =>
  fetch(`${BASE}/simulate`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ flow_id, candidate_id, ph, temp }),
  }).then((r) => r.json());
