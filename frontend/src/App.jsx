import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import Joyride from 'react-joyride';
import TopBar from './components/TopBar.jsx';
import HudFrame from './components/HudFrame.jsx';
import Heatmap from './components/Heatmap.jsx';
import Leaderboard from './components/Leaderboard.jsx';
import InspectorPanel from './components/InspectorPanel.jsx';
import FooterBar from './components/FooterBar.jsx';
import DataModal from './components/DataModal.jsx';
import { getFlowData, listFlows, getPreflight, simulateCell } from './lib/api.js';
import tourSteps from './tour/steps.js';

const threshold = 0.82;

const clamp01 = (value) => {
  const bounded = Math.min(0.99, Math.max(0.05, value));
  return Number(bounded.toFixed(3));
};

const TEMP_SENSITIVITY = {
  low: 0.015,
  medium: 0.045,
  high: 0.09,
};

const PH_SENSITIVITY = {
  low: 0.012,
  medium: 0.04,
  high: 0.08,
};

const HEAT_FAVOR_FLAGS = new Set(['heat_tolerant', 'shipping_friendly']);
const COLD_FAVOR_FLAGS = new Set([
  'cold_room_ready',
  'cold_boost',
  'cold_chain_ready',
  'cold_adapted',
  'low_temp_tuned',
]);
const WIDE_WINDOW_FLAGS = new Set(['wide_window']);
const HIGH_PH_FLAGS = new Set(['buffer_c_ok', 'phosphate_ok']);
const LOW_PH_FLAGS = new Set(['buffer_a_ok', 'cold_boost']);
const MID_PH_FLAGS = new Set(['buffer_b_ok', 'balanced_gc', 'locked_backbone']);

const hasAnyFlag = (flagSet, targets) => {
  for (const flag of targets) {
    if (flagSet.has(flag)) {
      return true;
    }
  }
  return false;
};

const createScenarioAwareMapper = (data, scenario, modifier, activeFlowId) => {
  if (!data) {
    return (_, __, ___, value) => ({ value });
  }

  const temps = data.heat?.temps ?? [];
  const phs = data.heat?.phs ?? [];
  const scoreMap = new Map(data.scores.map((item) => [item.id, item]));

  const maxTempIdx = Math.max(temps.length - 1, 1);
  const maxPhIdx = Math.max(phs.length - 1, 1);

  const [tempStart = 0, tempEnd = maxTempIdx] = scenario?.tempRange ?? [];
  const [phStart = 0, phEnd = maxPhIdx] = scenario?.phRange ?? [];

  const highTempFocus = tempEnd >= maxTempIdx && tempStart >= Math.max(maxTempIdx - 1, 0);
  const lowTempFocus = tempStart <= 0 && tempEnd <= Math.min(1, maxTempIdx);
  const highPhFocus = phEnd >= maxPhIdx && phStart >= Math.max(maxPhIdx - 1, 0);
  const lowPhFocus = phStart <= 0 && phEnd <= Math.min(1, maxPhIdx);

  return (candidateId, tIndex, pIndex, rawValue) => {
    const score = scoreMap.get(candidateId) ?? {};
    const flagsArray = score.flags ?? [];
    const flagSet = new Set(flagsArray);
    const tempSensitivity = TEMP_SENSITIVITY[score.sensitivity?.temp] ?? TEMP_SENSITIVITY.medium;
    const phSensitivity = PH_SENSITIVITY[score.sensitivity?.ph] ?? PH_SENSITIVITY.medium;

    const tempFraction = maxTempIdx === 0 ? 0 : tIndex / maxTempIdx;
    const coldFraction = 1 - tempFraction;
    const phFraction = maxPhIdx === 0 ? 0 : pIndex / maxPhIdx;
    const phCenterDistance = Math.abs(phFraction - 0.5);

    let adjusted = rawValue;
    const trace = {
      base: rawValue,
      tempBonus: 0,
      tempPenalty: 0,
      phBonus: 0,
      phPenalty: 0,
      scenarioBias: 0,
      modifierBoost: 0,
      jitter: 0,
      clampLoss: 0,
      final: rawValue,
      raw: rawValue,
      scenarioIncluded:
        tIndex >= tempStart && tIndex <= tempEnd && pIndex >= phStart && pIndex <= phEnd,
    };

    let tempPenalty = tempSensitivity * Math.pow(tempFraction, 1.2);
    if (hasAnyFlag(flagSet, WIDE_WINDOW_FLAGS)) {
      tempPenalty *= 0.4;
    }

    if (hasAnyFlag(flagSet, HEAT_FAVOR_FLAGS)) {
      const warmFavor = Math.pow(tempFraction, 1.25) * 0.08;
      const coldDrift = Math.pow(coldFraction, 1.2) * 0.02;
      adjusted += warmFavor;
      trace.tempBonus += warmFavor;
      adjusted -= coldDrift;
      trace.tempPenalty += coldDrift;
      tempPenalty *= 0.5;
    }

    if (hasAnyFlag(flagSet, COLD_FAVOR_FLAGS)) {
      const coldFavor = Math.pow(coldFraction, 1.25) * 0.08;
      const hotDrift = Math.pow(tempFraction, 1.1) * tempSensitivity * 0.6;
      adjusted += coldFavor;
      trace.tempBonus += coldFavor;
      tempPenalty += hotDrift;
    }

    if (highTempFocus) {
      const scenarioAdjust = hasAnyFlag(flagSet, HEAT_FAVOR_FLAGS) ? 0.05 : -0.06;
      adjusted += scenarioAdjust;
      trace.scenarioBias += scenarioAdjust;
    } else if (lowTempFocus) {
      const scenarioAdjust = hasAnyFlag(flagSet, COLD_FAVOR_FLAGS) ? 0.05 : -0.05;
      adjusted += scenarioAdjust;
      trace.scenarioBias += scenarioAdjust;
    }

    adjusted -= tempPenalty;
    trace.tempPenalty += tempPenalty;

    let phPenalty = phSensitivity * Math.pow(phCenterDistance, 1.2);
    if (hasAnyFlag(flagSet, WIDE_WINDOW_FLAGS)) {
      phPenalty *= 0.5;
    }

    if (hasAnyFlag(flagSet, HIGH_PH_FLAGS)) {
      const highBias = Math.max(0, phFraction - 0.6);
      const lowDrift = Math.max(0, 0.45 - phFraction);
      const highFavor = Math.pow(highBias, 1.2) * 0.12;
      adjusted += highFavor;
      trace.phBonus += highFavor;
      phPenalty += lowDrift * phSensitivity * 1.1;
    }

    if (hasAnyFlag(flagSet, LOW_PH_FLAGS)) {
      const lowBias = Math.max(0, 0.4 - phFraction);
      const highDrift = Math.max(0, phFraction - 0.55);
      const lowFavor = Math.pow(lowBias, 1.2) * 0.12;
      adjusted += lowFavor;
      trace.phBonus += lowFavor;
      phPenalty += highDrift * phSensitivity * 1.1;
    }

    if (hasAnyFlag(flagSet, MID_PH_FLAGS)) {
      const midBonus = Math.max(0, 0.28 - Math.abs(phFraction - 0.5));
      const edgeDrift = Math.max(0, Math.abs(phFraction - 0.5) - 0.28);
      const midFavor = Math.pow(midBonus, 1.3) * 0.08;
      adjusted += midFavor;
      trace.phBonus += midFavor;
      phPenalty += edgeDrift * phSensitivity;
    }

    if (highPhFocus) {
      const scenarioAdjust = hasAnyFlag(flagSet, HIGH_PH_FLAGS) ? 0.04 : -0.05;
      adjusted += scenarioAdjust;
      trace.scenarioBias += scenarioAdjust;
    } else if (lowPhFocus) {
      const scenarioAdjust = hasAnyFlag(flagSet, LOW_PH_FLAGS) ? 0.04 : -0.05;
      adjusted += scenarioAdjust;
      trace.scenarioBias += scenarioAdjust;
    }

    adjusted -= phPenalty;
    trace.phPenalty += phPenalty;

    if (modifier === 'mgBoost' && activeFlowId === 'flow2_dnaprobe') {
      const mgReady = flagSet.has('buffer_c_ok') || flagSet.has('mg2_plus');
      if (mgReady) {
        const boost = pIndex <= 2 ? 0.03 : 0.02;
        adjusted += boost;
        trace.modifierBoost += boost;
      }
    }

    const jitter =
      ((candidateId.charCodeAt(0) + candidateId.charCodeAt(candidateId.length - 1)) % 5) * 0.001;
    trace.jitter = jitter;

    const withJitter = adjusted + jitter;
    const clamped = clamp01(withJitter);
    trace.clampLoss = withJitter - clamped;
    trace.final = clamped;

    return { value: clamped, trace };
  };
};

const defaultTrace = (value, inScenario) => ({
  base: value,
  tempBonus: 0,
  tempPenalty: 0,
  phBonus: 0,
  phPenalty: 0,
  scenarioBias: 0,
  modifierBoost: 0,
  jitter: 0,
  clampLoss: 0,
  final: value,
  raw: value,
  scenarioIncluded: inScenario,
});

const parseScenarioCsv = (text, fallbackCandidateId) => {
  const lines = text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line.length > 0);
  if (lines.length === 0) {
    throw new Error('Scenario file contained no data.');
  }
  const headers = lines[0]
    .split(',')
    .map((h) => h.trim().toLowerCase());

  const tempIdx = headers.findIndex((header) => header === 'temp' || header === 'temperature');
  const phIdx = headers.findIndex((header) => header === 'ph' || header === 'ph_value' || header === 'phvalue');
  const valueIdx = headers.findIndex(
    (header) => header === 'value' || header === 'score' || header === 'stability'
  );
  const candidateIdx = headers.findIndex(
    (header) => header === 'candidate' || header === 'candidate_id' || header === 'id'
  );

  if (tempIdx === -1 || phIdx === -1 || valueIdx === -1) {
    throw new Error('Expected columns named temp, pH, and value.');
  }

  const overrides = new Map();
  const rows = [];
  for (let i = 1; i < lines.length; i += 1) {
    const parts = lines[i].split(',').map((part) => part.trim());
    if (parts.every((part) => part === '')) continue;

    const temp = Number(parts[tempIdx]);
    const ph = Number(parts[phIdx]);
    const rawValue = Number(parts[valueIdx]);
    const candidateId =
      candidateIdx !== -1 ? parts[candidateIdx] || fallbackCandidateId : fallbackCandidateId;

    if (!candidateId) {
      throw new Error('Candidate ID missing in file and no default candidate selected.');
    }
    if (Number.isNaN(temp) || Number.isNaN(ph) || Number.isNaN(rawValue)) {
      throw new Error(`Invalid numeric value on row ${i + 1}.`);
    }

    const entry = {
      temp,
      ph,
      value: clamp01(rawValue),
    };
    if (!overrides.has(candidateId)) {
      overrides.set(candidateId, []);
    }
    overrides.get(candidateId).push(entry);
    rows.push({ candidate: candidateId, temp, ph, value: entry.value });
  }

  if (overrides.size === 0) {
    throw new Error('No usable stability rows found in upload.');
  }

  return { overrides, rows };
};

const computeMetrics = (candidateId, data, scenario, valueMapper, breakdownAccumulator) => {
  if (!candidateId || !data?.heat) {
    return { windowAvg: 0, areaAbove: 0, maxStableTemp: '—', phWindow: 'n/a' };
  }
  const { temps, phs, candidates } = data.heat;
  const candidate = candidates.find((c) => c.id === candidateId);
  if (!candidate) {
    return { windowAvg: 0, areaAbove: 0, maxStableTemp: '—', phWindow: 'n/a' };
  }

  let breakdownEntry = null;
  if (breakdownAccumulator) {
    breakdownEntry = {
      totals: {
        count: 0,
        base: 0,
        tempBonus: 0,
        tempPenalty: 0,
        phBonus: 0,
        phPenalty: 0,
        scenarioBias: 0,
        modifierBoost: 0,
        jitter: 0,
        clampLoss: 0,
        final: 0,
      },
      cells: [],
    };
    breakdownAccumulator.set(candidateId, breakdownEntry);
  }

  let sum = 0;
  let count = 0;
  let above = 0;
  let maxStableTempIdx = -1;
  let minPhIdx = null;
  let maxPhIdx = null;

  for (let t = 0; t < temps.length; t += 1) {
    let rowStays = false;
    for (let p = 0; p < phs.length; p += 1) {
      const raw = candidate.values?.[t]?.[p] ?? 0;
      const mapped = valueMapper
        ? valueMapper(candidateId, t, p, raw)
        : { value: raw, trace: defaultTrace(raw, false) };
      const value = typeof mapped === 'number' ? mapped : mapped?.value ?? raw;
      const inScenario =
        t >= scenario.tempRange[0] &&
        t <= scenario.tempRange[1] &&
        p >= scenario.phRange[0] &&
        p <= scenario.phRange[1];
      const trace =
        typeof mapped === 'number'
          ? defaultTrace(value, inScenario)
          : mapped?.trace ?? defaultTrace(value, inScenario);
      trace.scenarioIncluded = trace.scenarioIncluded ?? inScenario;
      if (inScenario) {
        sum += value;
        count += 1;
      }
      if (value >= threshold) {
        above += 1;
        rowStays = true;
        if (minPhIdx === null || p < minPhIdx) minPhIdx = p;
        if (maxPhIdx === null || p > maxPhIdx) maxPhIdx = p;
      }
      if (breakdownEntry) {
        const cellTrace = trace ?? defaultTrace(value, inScenario);
        breakdownEntry.cells.push({
          tempIndex: t,
          phIndex: p,
          temp: temps[t],
          ph: phs[p],
          raw: raw,
          adjusted: value,
          contributions: {
            tempBonus: cellTrace.tempBonus,
            tempPenalty: cellTrace.tempPenalty,
            phBonus: cellTrace.phBonus,
            phPenalty: cellTrace.phPenalty,
            scenarioBias: cellTrace.scenarioBias,
            modifierBoost: cellTrace.modifierBoost,
            jitter: cellTrace.jitter,
            clampLoss: cellTrace.clampLoss,
          },
          inScenario: cellTrace.scenarioIncluded ?? inScenario,
        });
        if (cellTrace.scenarioIncluded ?? inScenario) {
          breakdownEntry.totals.count += 1;
          breakdownEntry.totals.base += cellTrace.base ?? raw;
          breakdownEntry.totals.tempBonus += cellTrace.tempBonus ?? 0;
          breakdownEntry.totals.tempPenalty += cellTrace.tempPenalty ?? 0;
          breakdownEntry.totals.phBonus += cellTrace.phBonus ?? 0;
          breakdownEntry.totals.phPenalty += cellTrace.phPenalty ?? 0;
          breakdownEntry.totals.scenarioBias += cellTrace.scenarioBias ?? 0;
          breakdownEntry.totals.modifierBoost += cellTrace.modifierBoost ?? 0;
          breakdownEntry.totals.jitter += cellTrace.jitter ?? 0;
          breakdownEntry.totals.clampLoss += cellTrace.clampLoss ?? 0;
          breakdownEntry.totals.final += cellTrace.final ?? value;
        }
      }
    }
    if (rowStays) {
      maxStableTempIdx = t;
    }
  }

  const windowAvg = count ? sum / count : 0;
  const areaAbove = above / (temps.length * phs.length);
  const maxStableTemp = maxStableTempIdx >= 0 ? temps[maxStableTempIdx] : '—';
  const phWindow =
    minPhIdx !== null && maxPhIdx !== null ? `${phs[minPhIdx]}–${phs[maxPhIdx]}` : 'n/a';

  return { windowAvg, areaAbove, maxStableTemp, phWindow };
};

const narratives = {
  flow1_antibody: [
    { id: 'baseline', label: 'Load Baseline', description: 'Return to default leaderboard.' },
    {
      id: 'hotday',
      label: 'Simulate Hot Day (45°C)',
      description: 'Focus on top-end temperatures and highlight heat tolerant variants.',
    },
    {
      id: 'shipping',
      label: 'Shipping Advice',
      description: 'Surface buffer recommendation for hot regions.',
    },
  ],
  flow2_dnaprobe: [
    { id: 'baseline', label: 'Load Baseline', description: 'Reset newborn screening scenario.' },
    {
      id: 'mgboost',
      label: 'Buffer Tweak (Mg²⁺ ↑)',
      description: 'Apply +0.03 stability gain to buffer-compatible probes.',
    },
    {
      id: 'coldroom',
      label: 'Low-Temp Stability',
      description: 'Focus on 4–25°C window for cold-room workflows.',
    },
  ],
  flow3_enzyme: [
    { id: 'baseline', label: 'Load Baseline', description: 'Reset enzyme assay board.' },
    {
      id: 'lostchain',
      label: 'Lost Cold Chain (10→25°C)',
      description: 'Animate temp sweep and spotlight shipping friendly variants.',
    },
    {
      id: 'specsheet',
      label: 'Spec Sheet Draft',
      description: 'Generate quick copy for datasheet output.',
    },
  ],
};

const App = () => {
  const [flows, setFlows] = useState([]);
  const [activeFlowId, setActiveFlowId] = useState(null);
  const [datasets, setDatasets] = useState({});
  const [datasetOriginals, setDatasetOriginals] = useState({});
  const [preflightMap, setPreflightMap] = useState({});
  const [selectedCandidateId, setSelectedCandidateId] = useState('');
  const [scenario, setScenario] = useState({ tempRange: [0, 0], phRange: [0, 0], whatIf: false });
  const [playing, setPlaying] = useState(false);
  const [playbackTempIndex, setPlaybackTempIndex] = useState(0);
  const [modalOpen, setModalOpen] = useState(false);
  const [loadingFlow, setLoadingFlow] = useState(false);
  const [loadError, setLoadError] = useState('');
  const [activeNarrativeId, setActiveNarrativeId] = useState('baseline');
  const [modifier, setModifier] = useState(null);
  const [narrativeNote, setNarrativeNote] = useState('');
  const [heatmapPulse, setHeatmapPulse] = useState(false);
  const [previousStats, setPreviousStats] = useState(null);
  const [simulationDetails, setSimulationDetails] = useState(null);
  const [demoMode, setDemoMode] = useState(false);
  const [toast, setToast] = useState(null);
  const [helpOpen, setHelpOpen] = useState(false);
  const [tourRunning, setTourRunning] = useState(false);
  const [hasSeenTour, setHasSeenTour] = useState(() => {
    if (typeof window === 'undefined') return true;
    return window.localStorage.getItem('rrs_seen_tour') === 'true';
  });
  const playbackTimeoutRef = useRef(null);
  const prevOrderRef = useRef(new Map());
  const toastTimeoutRef = useRef(null);
  const pulseTimeoutRef = useRef(null);
  const [macToggle, setMacToggle] = useState(false);
  const [compareReference, setCompareReference] = useState(false);
  const [tourEnabled, setTourEnabled] = useState(() => {
    if (typeof window === 'undefined') return false;
    return window.localStorage.getItem('rrs_tour_enabled') !== 'false';
  });
  const [scenarioPreview, setScenarioPreview] = useState(null);

  useEffect(() => {
    listFlows()
      .then((data) => {
        setFlows(data);
        if (data.length > 0) {
          setActiveFlowId(data[0].id);
        }
      })
      .catch(() => setLoadError('Unable to load flows. Ensure backend is running on :8080.'));
  }, []);

  useEffect(() => {
    setScenarioPreview(null);
  }, [activeFlowId]);

  useEffect(() => {
    if (!hasSeenTour && flows.length > 0 && tourEnabled) {
      setTourRunning(true);
    }
  }, [hasSeenTour, flows.length, tourEnabled]);

  useEffect(() => {
    if (!activeFlowId) return;
    if (datasets[activeFlowId]) {
      setActiveNarrativeId('baseline');
      setNarrativeNote('');
      setModifier(null);
      setPlaying(false);
      return;
    }
    setLoadingFlow(true);
    getFlowData(activeFlowId)
      .then((data) => {
        const canonical = JSON.parse(JSON.stringify(data));
        setDatasets((prev) => ({ ...prev, [activeFlowId]: canonical }));
        setDatasetOriginals((prev) => {
          if (prev[activeFlowId]) return prev;
          return { ...prev, [activeFlowId]: JSON.parse(JSON.stringify(data)) };
        });
        const firstScore = [...canonical.scores].sort((a, b) => a.rank - b.rank)[0];
        setSelectedCandidateId(firstScore?.id ?? canonical.cand[0]?.id ?? '');
        setScenario({
          tempRange: [0, canonical.heat.temps.length - 1],
          phRange: [0, canonical.heat.phs.length - 1],
          whatIf: false,
        });
        setPlaybackTempIndex(Math.max(0, Math.floor(canonical.heat.temps.length / 2)));
      })
      .catch(() => setLoadError('Failed to fetch flow dataset.'))
      .finally(() => {
        setLoadingFlow(false);
      });
  }, [activeFlowId, datasets]);

  useEffect(() => {
    if (!activeFlowId) return;
    getPreflight(activeFlowId)
      .then((info) => {
        setPreflightMap((prev) => ({ ...prev, [activeFlowId]: info }));
      })
      .catch(() => {
        setPreflightMap((prev) => ({ ...prev, [activeFlowId]: null }));
      });
  }, [activeFlowId]);

  useEffect(() => {
    if (playing && activeFlowId) {
      const data = datasets[activeFlowId];
      if (!data) return;
      const maxIndex = data.heat.temps.length;
      const id = setInterval(() => {
        setPlaybackTempIndex((prev) => (prev + 1) % maxIndex);
      }, 1400);
      return () => clearInterval(id);
    }
    return undefined;
  }, [playing, activeFlowId, datasets]);

  useEffect(() => {
    if (!activeFlowId) return;
    const data = datasets[activeFlowId];
    if (!data) return;
    setScenario({
      tempRange: [0, data.heat.temps.length - 1],
      phRange: [0, data.heat.phs.length - 1],
      whatIf: false,
    });
    const firstScore = [...data.scores].sort((a, b) => a.rank - b.rank)[0];
    setSelectedCandidateId(firstScore?.id ?? data.cand[0]?.id ?? '');
    setPlaybackTempIndex(Math.max(0, Math.floor(data.heat.temps.length / 2)));
    setActiveNarrativeId('baseline');
    setModifier(null);
    setNarrativeNote('');
    setPlaying(false);
    setSimulationDetails(null);
    setPreviousStats(null);
    setHeatmapPulse(false);
    setCompareReference(false);
    if (playbackTimeoutRef.current) {
      clearTimeout(playbackTimeoutRef.current);
      playbackTimeoutRef.current = null;
    }
  }, [activeFlowId, datasets]);

  const currentData = activeFlowId ? datasets[activeFlowId] : null;

  const valueMapper = useMemo(
    () => createScenarioAwareMapper(currentData, scenario, modifier, activeFlowId),
    [activeFlowId, currentData, modifier, scenario]
  );

  const { scenarioMetrics, scoreBreakdowns } = useMemo(() => {
    if (!currentData) return { scenarioMetrics: new Map(), scoreBreakdowns: new Map() };
    const metricsMap = new Map();
    const breakdownMap = new Map();
    for (const cand of currentData.heat.candidates) {
      const metrics = computeMetrics(cand.id, currentData, scenario, valueMapper, breakdownMap);
      metricsMap.set(cand.id, metrics);
    }
    return { scenarioMetrics: metricsMap, scoreBreakdowns: breakdownMap };
  }, [currentData, scenario, valueMapper]);

  const leaderboardBase = useMemo(() => {
    if (!currentData) return [];
    const scoreMap = new Map(currentData.scores.map((item) => [item.id, item]));
    const candidateMap = new Map(currentData.cand.map((item) => [item.id, item]));
    return currentData.rat.top5.map((entry) => {
      const score = scoreMap.get(entry.id) ?? { robustness: 0, rank: 99, flags: [] };
      const candidate = candidateMap.get(entry.id) ?? { id: entry.id, name: entry.id };
      return {
        id: entry.id,
        name: candidate.name,
        robustness: score.robustness ?? 0,
        rank: score.rank ?? 99,
        flags: score.flags ?? [],
        why: entry.why,
        watch: entry.watch,
      };
    });
  }, [currentData]);

  const leaderboardItems = useMemo(() => {
    const items = leaderboardBase.map((item) => ({
      ...item,
      scenarioScore: scenarioMetrics.get(item.id)?.windowAvg ?? item.robustness,
    }));
    if (scenario.whatIf) {
      return [...items].sort((a, b) => (b.scenarioScore ?? 0) - (a.scenarioScore ?? 0));
    }
    return items;
  }, [leaderboardBase, scenarioMetrics, scenario.whatIf]);

  const leaderboardWithDelta = useMemo(() => {
    const previousOrder = prevOrderRef.current;
    const mapped = leaderboardItems.map((item, index) => {
      const prevIndex = previousOrder.has(item.id) ? previousOrder.get(item.id) : null;
      const delta = prevIndex === null ? 0 : prevIndex - index;
      return { ...item, delta };
    });
    prevOrderRef.current = new Map(mapped.map((item, index) => [item.id, index]));
    return mapped;
  }, [leaderboardItems]);

  const selectedCandidate = currentData?.cand.find((c) => c.id === selectedCandidateId);
  const selectedRationale = currentData?.rat.top5.find((item) => item.id === selectedCandidateId);

  const preflight = preflightMap[activeFlowId] ?? null;

  const artifactStatus = useMemo(() => {
    if (!selectedCandidate) return [];
    const items = [];
    const hasSmiles = Boolean(selectedCandidate.smiles);
    items.push({
      id: 'smiles',
      label: 'SMILES',
      symbol: hasSmiles ? 'LIVE' : 'MOCK',
      variant: hasSmiles ? 'good' : 'mock',
      tooltip: hasSmiles
        ? 'SMILES parsed via RDKit ETKDG builder.'
        : 'No SMILES available; using canned conformer.',
    });

    const propkaActive = preflight?.propka_active ?? false;
    const propkaRequested = preflight?.propka_requested ?? false;
    items.push({
      id: 'propka',
      label: 'PROPKA',
      symbol: propkaActive ? 'LIVE' : propkaRequested ? 'MISS' : 'MOCK',
      variant: propkaActive ? 'good' : propkaRequested ? 'warn' : 'mock',
      tooltip: propkaActive
        ? 'Residue pKa values computed with PROPKA.'
        : propkaRequested
        ? 'PROPKA requested but falling back to mocked residue pKa map.'
        : 'PROPKA mocked for this environment.',
    });

    const dockingRequested = preflight?.docking_requested ?? false;
    const dockingActive = preflight?.docking_active ?? false;
    items.push({
      id: 'docking',
      label: 'Docking',
      symbol: dockingActive ? 'LIVE' : dockingRequested ? 'MISS' : 'MOCK',
      variant: dockingActive ? 'good' : dockingRequested ? 'warn' : 'mock',
      tooltip: dockingActive
        ? 'Docking stub executed (AutoDock Vina/DiffDock placeholder).'
        : 'Docking mocked from cached poses and scores.',
    });

    items.push({
      id: 'mmgbsa',
      label: 'MM/GBSA',
      symbol: 'MOCK',
      variant: 'mock',
      tooltip: 'Robustness proxy derived from docking score, RMSD and H-bond heuristics.',
    });
    return items;
  }, [selectedCandidate, preflight]);

  const inspectorBreakdown = scoreBreakdowns.get(selectedCandidateId) ?? null;

  const stabilityCurve = useMemo(() => {
    if (!inspectorBreakdown || !currentData?.heat) return [];
    const temps = currentData.heat.temps ?? [];
    const groups = new Map();
    inspectorBreakdown.cells.forEach((cell) => {
      const key = cell.tempIndex;
      const entry =
        groups.get(key) ?? {
          temp: cell.temp,
          rawSum: 0,
          adjustedSum: 0,
          count: 0,
          scenarioSum: 0,
          scenarioCount: 0,
        };
      entry.rawSum += cell.raw ?? 0;
      entry.adjustedSum += cell.adjusted ?? 0;
      entry.count += 1;
      if (cell.inScenario) {
        entry.scenarioSum += cell.adjusted ?? 0;
        entry.scenarioCount += 1;
      }
      groups.set(key, entry);
    });
    return temps.map((temp, index) => {
      const entry = groups.get(index);
      if (!entry || entry.count === 0) {
        return { temp, raw: 0, adjusted: 0, scenario: null };
      }
      return {
        temp,
        raw: entry.rawSum / entry.count,
        adjusted: entry.adjustedSum / entry.count,
        scenario: entry.scenarioCount ? entry.scenarioSum / entry.scenarioCount : null,
      };
    });
  }, [inspectorBreakdown, currentData]);

  const cellMatrix = useMemo(() => {
    if (!inspectorBreakdown || !currentData?.heat) return [];
    const temps = currentData.heat.temps ?? [];
    const phs = currentData.heat.phs ?? [];
    const cellMap = new Map(
      inspectorBreakdown.cells.map((cell) => [`${cell.tempIndex}:${cell.phIndex}`, cell])
    );
    return temps.map((temp, tIndex) => ({
      temp,
      cells: phs.map((ph, pIndex) => {
        const key = `${tIndex}:${pIndex}`;
        const existing = cellMap.get(key);
        if (existing) {
          return existing;
        }
        const inScenario =
          tIndex >= scenario.tempRange[0] &&
          tIndex <= scenario.tempRange[1] &&
          pIndex >= scenario.phRange[0] &&
          pIndex <= scenario.phRange[1];
        return {
          tempIndex: tIndex,
          phIndex: pIndex,
          temp,
          ph,
          raw: Number.NaN,
          adjusted: Number.NaN,
          contributions: {
            tempBonus: 0,
            tempPenalty: 0,
            phBonus: 0,
            phPenalty: 0,
            scenarioBias: 0,
            modifierBoost: 0,
            jitter: 0,
            clampLoss: 0,
          },
          inScenario,
        };
      }),
    }));
  }, [inspectorBreakdown, currentData, scenario]);

  const inspectorStats =
    scenarioMetrics.get(selectedCandidateId) ??
    computeMetrics(selectedCandidateId, currentData, scenario, valueMapper);

  const pushToast = useCallback((message) => {
    if (toastTimeoutRef.current) {
      clearTimeout(toastTimeoutRef.current);
    }
    setToast({ id: Date.now(), message });
    toastTimeoutRef.current = setTimeout(() => {
      setToast(null);
      toastTimeoutRef.current = null;
    }, 2400);
  }, []);

  const triggerHeatmapPulse = useCallback(() => {
    setHeatmapPulse(true);
    if (pulseTimeoutRef.current) {
      clearTimeout(pulseTimeoutRef.current);
    }
    pulseTimeoutRef.current = setTimeout(() => {
      setHeatmapPulse(false);
      pulseTimeoutRef.current = null;
    }, 900);
  }, []);

  const applyScenarioOverrides = useCallback(
    (overrideEntries, sourceLabel, previewRows = null) => {
      if (!activeFlowId || !currentData) {
        throw new Error('No dataset loaded. Select a flow before applying overrides.');
      }
      const temps = currentData.heat?.temps ?? [];
      const phs = currentData.heat?.phs ?? [];
      if (!temps.length || !phs.length) {
        throw new Error('Active dataset is missing heatmap grid metadata.');
      }
      const candidateIds = new Set(currentData.heat?.candidates?.map((cand) => cand.id) ?? []);
      const overrideMatrices = new Map();

      overrideEntries.forEach((entries, candidateId) => {
        if (!candidateIds.has(candidateId)) {
          throw new Error(`Candidate ${candidateId} is not part of this flow.`);
        }
        const byTemp = new Map();
        entries.forEach(({ temp, ph, value }) => {
          const tempIdx = temps.indexOf(temp);
          if (tempIdx === -1) {
            throw new Error(
              `Temperature ${temp}°C is not on this heatmap grid (allowed: ${temps.join(', ')}).`
            );
          }
          const phIdx = phs.indexOf(ph);
          if (phIdx === -1) {
            throw new Error(
              `pH ${ph} is not on this heatmap grid (allowed: ${phs.join(', ')}).`
            );
          }
          if (!byTemp.has(tempIdx)) {
            byTemp.set(tempIdx, new Map());
          }
          byTemp.get(tempIdx).set(phIdx, value);
        });
        overrideMatrices.set(candidateId, byTemp);
      });

      setDatasets((prev) => {
        const existing = prev[activeFlowId];
        if (!existing) return prev;
        const updated = JSON.parse(JSON.stringify(existing));
        updated.heat.candidates = updated.heat.candidates.map((cand) => {
          const matrix = overrideMatrices.get(cand.id);
          if (!matrix) return cand;
          const newValues = cand.values.map((row, tIndex) =>
            row.map((cellValue, pIndex) => {
              const overrideValue = matrix.get(tIndex)?.get(pIndex);
              return overrideValue !== undefined ? overrideValue : cellValue;
            })
          );
          return { ...cand, values: newValues };
        });
        return { ...prev, [activeFlowId]: updated };
      });

      setScenario((prev) => ({ ...prev, whatIf: true }));
      setNarrativeNote(`Custom stability imported from ${sourceLabel}.`);
      triggerHeatmapPulse();
      pushToast(`Scenario heatmap updated via ${sourceLabel}.`);
      if (previewRows && previewRows.length > 0) {
        setScenarioPreview({
          label: sourceLabel,
          rows: previewRows,
        });
      } else {
        setScenarioPreview(null);
      }
    },
    [activeFlowId, currentData, triggerHeatmapPulse, pushToast, setScenarioPreview]
  );

  const handleScenarioUpload = useCallback(
    async (file) => {
      if (!file || !activeFlowId || !currentData) {
        pushToast('No active dataset for scenario upload.');
        return;
      }
      try {
        const text = await file.text();
        const fallbackCandidate =
          selectedCandidateId || currentData.cand?.[0]?.id || currentData.heat?.candidates?.[0]?.id;
        if (!fallbackCandidate) {
          throw new Error('Unable to determine candidate for scenario overrides.');
        }
        const { overrides, rows } = parseScenarioCsv(text, fallbackCandidate);
        applyScenarioOverrides(overrides, file.name, rows.slice(0, 12));
      } catch (error) {
        console.error('Scenario upload failed', error);
        pushToast(error.message || 'Scenario upload failed.');
      }
    },
    [activeFlowId, currentData, selectedCandidateId, applyScenarioOverrides, pushToast]
  );

  const handleSampleScenario = useCallback(async () => {
    if (!activeFlowId || !currentData) {
      pushToast('No active dataset for scenario upload.');
      return;
    }
    try {
      const response = await fetch('/sample_scenario.csv');
      if (!response.ok) {
        throw new Error('Sample scenario file not found.');
      }
      const text = await response.text();
      const fallbackCandidate =
        selectedCandidateId || currentData.cand?.[0]?.id || currentData.heat?.candidates?.[0]?.id;
      if (!fallbackCandidate) {
        throw new Error('Unable to determine candidate for scenario overrides.');
      }
      const { overrides, rows } = parseScenarioCsv(text, fallbackCandidate);
      applyScenarioOverrides(overrides, 'Sample Scenario', rows.slice(0, 12));
    } catch (error) {
      console.error('Sample scenario failed', error);
      pushToast(error.message || 'Unable to load sample scenario.');
    }
  }, [activeFlowId, currentData, selectedCandidateId, applyScenarioOverrides, pushToast]);

  const handleTourCallback = useCallback(
    (data) => {
      const { status } = data;
      if (status === 'finished' || status === 'skipped') {
        setTourRunning(false);
        setHasSeenTour(true);
        if (typeof window !== 'undefined') {
          window.localStorage.setItem('rrs_seen_tour', 'true');
        }
      }
    },
    []
  );

  const handleNarrative = (actionId) => {
    if (!currentData) return;
    setActiveNarrativeId(actionId);
    const temps = currentData.heat.temps;
    const phs = currentData.heat.phs;
    if (playbackTimeoutRef.current) {
      clearTimeout(playbackTimeoutRef.current);
      playbackTimeoutRef.current = null;
    }
    setPreviousStats(inspectorStats);
    triggerHeatmapPulse();
    setSimulationDetails(null);

    switch (actionId) {
      case 'baseline': {
        const originalDataset = datasetOriginals[activeFlowId];
        if (originalDataset) {
          setDatasets((prev) => ({
            ...prev,
            [activeFlowId]: JSON.parse(JSON.stringify(originalDataset)),
          }));
        }
        setScenario({
          tempRange: [0, temps.length - 1],
          phRange: [0, phs.length - 1],
          whatIf: false,
        });
        const firstScore = [...currentData.scores].sort((a, b) => a.rank - b.rank)[0];
        setSelectedCandidateId(firstScore?.id ?? currentData.cand[0]?.id ?? '');
        setModifier(null);
        setNarrativeNote('');
        setPlaying(false);
        setPlaybackTempIndex(Math.max(0, Math.floor(temps.length / 2)));
        setScenarioPreview(null);
        pushToast('Baseline scenario restored; leaderboard reset.');
        break;
      }
      case 'hotday': {
        setScenario({
          tempRange: [temps.length - 1, temps.length - 1],
          phRange: [0, phs.length - 1],
          whatIf: true,
        });
        setPlaybackTempIndex(temps.length - 1);
        setSelectedCandidateId('C10');
        setNarrativeNote('C10 holds >0.85 at 45°C; other variants drop 6–9%.');
        setModifier(null);
        setPlaying(true);
        playbackTimeoutRef.current = setTimeout(() => setPlaying(false), 6000);
        pushToast('Locked on 45°C band; highlighting heat-tolerant variants.');
        break;
      }
      case 'shipping': {
        setNarrativeNote(currentData.rat.suggestions?.[0] ?? 'Apply phosphate buffer 20–50 mM.');
        setSelectedCandidateId('C12');
        setModifier(null);
        pushToast('Surfaced shipping guidance from rationale.');
        break;
      }
      case 'mgboost': {
        setModifier('mgBoost');
        setScenario((prev) => ({ ...prev, whatIf: true }));
        setNarrativeNote('Mg²⁺ tweak applied: +0.03 boost for buffer-compatible probes.');
        pushToast('Applied Mg²⁺ boost to compatible probes.');
        break;
      }
      case 'coldroom': {
        const lastColdIndex = Math.max(0, temps.findIndex((t) => t >= 25));
        setScenario({
          tempRange: [0, Math.max(lastColdIndex, 2)],
          phRange: [0, phs.length - 1],
          whatIf: true,
        });
        setPlaybackTempIndex(1);
        setModifier(null);
        setNarrativeNote('Cold-room focus (≤25°C) exposes D12 and D08 as most stable.');
        pushToast('Focused on cold-room window for neonatal workflow.');
        break;
      }
      case 'lostchain': {
        const startIdx = Math.max(0, temps.findIndex((t) => t >= 10));
        setScenario({
          tempRange: [startIdx, temps.length - 1],
          phRange: [0, phs.length - 1],
          whatIf: true,
        });
        setModifier(null);
        setSelectedCandidateId('E08');
        setNarrativeNote('Shipping alarm: highlight variants with >0.85 average up to 25°C.');
        setPlaying(true);
        playbackTimeoutRef.current = setTimeout(() => setPlaying(false), 6500);
        pushToast('Simulating lost cold-chain exposure.');
        break;
      }
      case 'specsheet': {
        setModifier(null);
        setNarrativeNote('Spec Sheet Draft: Stable 4–25°C; pH 7.2–8.0; avoid >25°C beyond 48 h.');
        pushToast('Drafted one-line spec sheet guidance.');
        break;
      }
      default:
        break;
    }
  };

  const handleTempRangeChange = (range) => {
    setPreviousStats(inspectorStats);
    setScenario((prev) => ({ ...prev, tempRange: range, whatIf: true }));
    triggerHeatmapPulse();
    pushToast('Changed temperature window → leaderboard re-ranked.');
  };

  const handlePhRangeChange = (range) => {
    setPreviousStats(inspectorStats);
    setScenario((prev) => ({ ...prev, phRange: range, whatIf: true }));
    triggerHeatmapPulse();
    pushToast('Adjusted pH window → scenario metrics updated.');
  };

  const handleWhatIfToggle = () => {
    setPreviousStats(inspectorStats);
    setScenario((prev) => ({ ...prev, whatIf: !prev.whatIf }));
    triggerHeatmapPulse();
    pushToast('What-If mode toggled.');
  };

  const handleSweep = async () => {
    if (!currentData || !activeFlowId || !selectedCandidateId) return;
    const temps = currentData.heat.temps;
    const phs = currentData.heat.phs;
    const targetTemp = temps.includes(37) ? 37 : temps[Math.min(playbackTempIndex, temps.length - 1)];
    const phIndex = Math.round((scenario.phRange[0] + scenario.phRange[1]) / 2);
    const clampedPhIndex = Math.max(0, Math.min(phs.length - 1, phIndex));
    const targetPh = phs[clampedPhIndex];

    try {
      setPreviousStats(inspectorStats);
      triggerHeatmapPulse();
      const result = await simulateCell({
        flow_id: activeFlowId,
        candidate_id: selectedCandidateId,
        ph: targetPh,
        temp: targetTemp,
      });
      setSimulationDetails(result.details);
      const msg = `Runtime chem simulate @${targetTemp}°C, pH ${targetPh}: ${result.value.toFixed(2)}`;
      setNarrativeNote('Computed using RDKit protonation and docking stub for the live cell.');
      pushToast(msg);
    } catch (error) {
      pushToast('Chem simulate failed; using heatmap baseline.');
    }
  };

  const startTour = () => {
    setTourRunning(true);
  };

  const currentFlowMeta = flows.find((f) => f.id === activeFlowId);

  const flowBanner = useMemo(() => {
    switch (activeFlowId) {
      case 'flow1_antibody':
        return 'Stability summary for antibody candidates.';
      case 'flow2_dnaprobe':
        return 'Evaluate probe robustness across screening conditions.';
      case 'flow3_enzyme':
        return 'Monitor enzyme stability under cold-chain scenarios.';
      default:
        return null;
    }
  }, [activeFlowId]);

  const successMetrics = useMemo(() => {
    if (!currentData) {
      return [];
    }
    const base = currentData.scores ?? [];
    const topScore = base?.[0]?.robustness ?? 0;
    return [
      {
        label: 'Top Robustness',
        value: topScore > 0 ? topScore.toFixed(2) : '—',
        tooltip: 'Highest baseline robustness score across candidates.',
      },
      {
        label: 'Scenario Mode',
        value: scenario.whatIf ? 'Active' : 'Baseline',
        tooltip: 'Indicates whether rankings reflect the current scenario window.',
      },
      {
        label: 'Chem Stack',
        value: macToggle ? 'Full' : 'Demo',
        tooltip: macToggle
          ? 'Full stack assumes docking and PROPKA integrations.'
          : 'Demo stack uses mocked chem operations.',
      },
    ];
  }, [currentData, scenario.whatIf, macToggle]);

  const referenceCandidate = useMemo(() => {
    if (!currentData) return null;
    const legacyId = currentData.scores.find((entry) => entry.flags?.includes('legacy_reference'))?.id
      ?? currentData.scores[currentData.scores.length - 1]?.id;
    if (!legacyId) return null;
    const heatEntry = currentData.heat.candidates.find((cand) => cand.id === legacyId);
    return heatEntry ? { id: legacyId, values: heatEntry.values } : null;
  }, [currentData]);

  useEffect(() => {
    if (!referenceCandidate && compareReference) {
      setCompareReference(false);
    }
  }, [referenceCandidate, compareReference]);

  const showTourCta = !hasSeenTour && !tourRunning;
  const demoSteps = narratives[activeFlowId] ?? [];

  useEffect(() => {
    return () => {
      if (toastTimeoutRef.current) {
        clearTimeout(toastTimeoutRef.current);
      }
      if (pulseTimeoutRef.current) {
        clearTimeout(pulseTimeoutRef.current);
      }
    };
  }, []);

  return (
    <div className="min-h-screen bg-orange-50 text-slate-900">
      {tourEnabled && (
        <Joyride
          steps={tourSteps}
          run={tourRunning}
          continuous
          showSkipButton
          callback={handleTourCallback}
          styles={{
            options: {
              primaryColor: '#2563eb',
              zIndex: 10000,
            },
          }}
        />
      )}
      {toast && (
        <div className="fixed right-6 top-6 z-40 rounded-lg border border-orange-200 bg-white px-4 py-3 text-sm text-orange-700 shadow-md">
          {toast.message}
        </div>
      )}
      <main className="mx-auto flex w-full max-w-[1280px] flex-col gap-6 px-4 py-6">
        <TopBar
          flows={flows}
          activeFlowId={activeFlowId}
          onSelectFlow={setActiveFlowId}
          onShowData={() => setModalOpen(true)}
          narrativeActions={narratives[activeFlowId] ?? []}
          activeNarrativeId={activeNarrativeId}
          onNarrative={handleNarrative}
          preflight={preflight}
          onStartTour={startTour}
          showTourCta={showTourCta}
          demoMode={demoMode}
          onToggleDemoMode={() => setDemoMode((prev) => !prev)}
          demoSteps={demoSteps}
          onDemoStep={handleNarrative}
          onShowHelp={() => setHelpOpen(true)}
          flowBanner={flowBanner}
          successMetrics={successMetrics}
          macToggle={macToggle}
          onToggleMac={() => setMacToggle((prev) => !prev)}
          showTourToggle
          tourRunning={tourRunning}
          onToggleTour={() => {
            const next = !tourEnabled;
            setTourEnabled(next);
            window.localStorage.setItem('rrs_tour_enabled', next ? 'true' : 'false');
            if (!next) {
              setTourRunning(false);
            } else {
              setTourRunning(true);
            }
          }}
          onUploadScenario={handleScenarioUpload}
          onLoadSampleScenario={handleSampleScenario}
        >
          {currentData && (
            <FooterBar
              temps={currentData.heat.temps}
              phs={currentData.heat.phs}
              tempRange={scenario.tempRange}
              phRange={scenario.phRange}
              onTempRangeChange={handleTempRangeChange}
              onPhRangeChange={handlePhRangeChange}
              whatIf={scenario.whatIf}
              onToggleWhatIf={handleWhatIfToggle}
            />
          )}
        </TopBar>

        {loadError && (
          <div className="rounded-xl border border-rose-200 bg-rose-50 p-4 text-sm text-rose-700">
            {loadError}
          </div>
        )}

        <section className="grid grid-cols-1 gap-6 lg:auto-rows-[minmax(0,520px)] lg:grid-cols-[minmax(0,2fr)_minmax(0,1fr)_minmax(0,1fr)]">
          <HudFrame
            id="heatmap-panel"
            title="pH × Temperature Map"
            subtitle="Robustness likelihood"
            className={`${heatmapPulse ? 'ring-2 ring-orange-200 border-orange-200' : ''} h-full`}
            actions={
              <div className="flex items-center gap-2">
                <button
                  type="button"
                  onClick={() => setCompareReference((prev) => !prev)}
                  className={`rounded-full border px-3 py-1 text-xs transition ${
                    compareReference
                      ? 'border-blue-400 bg-blue-50 text-blue-700'
                      : 'border-slate-200 bg-white text-slate-600 hover:bg-slate-100'
                  }`}
                  title="Overlay legacy reference for comparison"
                >
                  Compare ref
                </button>
                <button
                  type="button"
                  onClick={handleSweep}
                  className="rounded-full border border-slate-200 bg-white px-3 py-1 text-xs text-slate-600 transition hover:bg-slate-100"
                  title="Run runtime chem simulate for current candidate"
                >
                  37°C Sweep
                </button>
              </div>
            }
          >
            <div className="flex h-full min-h-[360px] flex-col overflow-hidden rounded-2xl border border-cyan-700/40 bg-[#0b1328]/80">
              {loadingFlow || !currentData ? (
                <div className="flex flex-1 items-center justify-center text-slate-500">
                  Loading heatmap…
                </div>
              ) : (
                <div className="flex-1 overflow-hidden">
                  <Heatmap
                    heatmap={currentData.heat}
                    candidates={currentData.cand}
                    candidateId={selectedCandidateId}
                    onCandidateChange={setSelectedCandidateId}
                    playbackTempIndex={playbackTempIndex}
                    onPlaybackChange={setPlaybackTempIndex}
                    playing={playing}
                    onTogglePlay={() => setPlaying((prev) => !prev)}
                    tempRange={scenario.tempRange}
                    phRange={scenario.phRange}
                    whatIf={scenario.whatIf}
                    valueMap={valueMapper}
                    compareReference={compareReference}
                    reference={referenceCandidate}
                  />
                </div>
              )}
            </div>
          </HudFrame>

          <HudFrame id="leaderboard" title="Leaderboard" subtitle="Top-5 picks" className="h-full">
            <div className="flex h-full flex-col overflow-y-auto px-1 py-1">
              <Leaderboard
                items={leaderboardWithDelta}
                activeId={selectedCandidateId}
                onSelect={setSelectedCandidateId}
                scenarioActive={scenario.whatIf}
                artifactStatus={artifactStatus}
                breakdowns={scoreBreakdowns}
              />
            </div>
          </HudFrame>

          <HudFrame id="inspector-panel" title="Inspector" subtitle="Scenario metrics" className="h-full">
            <div className="flex h-full flex-col overflow-y-auto px-1 py-1">
              <InspectorPanel
                candidate={selectedCandidate}
                stats={inspectorStats}
                previousStats={previousStats}
                rationale={selectedRationale}
                narrativeNote={narrativeNote}
                scenarioActive={scenario.whatIf}
                simulationDetails={simulationDetails}
                breakdown={inspectorBreakdown}
                curveData={stabilityCurve}
                cellMatrix={cellMatrix}
              />
            </div>
          </HudFrame>
        </section>
      </main>

      <DataModal flowId={activeFlowId ?? ''} open={modalOpen} onOpenChange={setModalOpen} />
      {helpOpen && (
        <div className="fixed inset-0 z-30 flex items-center justify-center bg-black/40 backdrop-blur-sm">
          <div className="w-[min(90vw,440px)] rounded-2xl border border-slate-200 bg-white p-6 text-sm text-slate-600 shadow-xl">
            <div className="flex items-center justify-between">
              <h2 className="text-lg font-semibold text-slate-800">How to use</h2>
              <button
                type="button"
                onClick={() => setHelpOpen(false)}
                className="rounded-full border border-slate-200 px-3 py-1 text-sm text-slate-600 hover:bg-slate-100"
              >
                Close
              </button>
            </div>
            <ol className="mt-4 space-y-2 text-slate-600">
              <li>1. Choose a flow tab (Antibody / DNA / Enzyme).</li>
              <li>2. Select a top candidate to sync the heatmap and inspector.</li>
              <li>3. Adjust pH and temperature sliders to explore scenario robustness.</li>
              <li>4. Use the narrative buttons or demo-mode steps to share insights.</li>
              <li>5. Press “Show Data” for the exact CSV/JSON, or run the 37°C sweep to call the runtime simulate endpoint.</li>
            </ol>
          </div>
        </div>
      )}
      {scenarioPreview && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 backdrop-blur-sm">
          <div className="w-[min(90vw,520px)] rounded-2xl border border-slate-200 bg-white p-6 text-sm text-slate-600 shadow-xl">
            <div className="flex items-center justify-between">
              <h2 className="text-lg font-semibold text-slate-800">
                Scenario preview • {scenarioPreview.label}
              </h2>
              <button
                type="button"
                onClick={() => setScenarioPreview(null)}
                className="rounded-full border border-slate-200 px-3 py-1 text-sm text-slate-600 hover:bg-slate-100"
              >
                Close
              </button>
            </div>
            <p className="mt-3 text-xs text-slate-500">
              Showing first {scenarioPreview.rows.length} rows. Values outside the grid are ignored.
            </p>
            <div className="mt-3 max-h-64 overflow-y-auto rounded-lg border border-slate-200 bg-slate-50">
              <table className="w-full table-auto text-left text-xs text-slate-600">
                <thead className="text-xs font-medium text-slate-500">
                  <tr>
                    <th className="px-3 py-2">Candidate</th>
                    <th className="px-3 py-2">Temp (°C)</th>
                    <th className="px-3 py-2">pH</th>
                    <th className="px-3 py-2">Value</th>
                  </tr>
                </thead>
                <tbody>
                  {scenarioPreview.rows.map((row, index) => (
                    <tr key={`${row.candidate}-${index}`} className="border-t border-slate-200">
                      <td className="px-3 py-2 font-mono text-slate-700">{row.candidate}</td>
                      <td className="px-3 py-2 font-mono text-slate-700">{row.temp}</td>
                      <td className="px-3 py-2 font-mono text-slate-700">{row.ph}</td>
                      <td className="px-3 py-2 font-mono text-slate-700">{row.value.toFixed(2)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <p className="mt-3 text-xs text-slate-500">
              The uploaded values are already applied to the heatmap and leaderboard.
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default App;
