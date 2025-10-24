const tourSteps = [
  {
    target: '#flow-tabs',
    content: 'Choose a reagent scenario to load its dataset, visuals, and mock chem pipeline.',
    disableBeacon: true,
  },
  {
    target: '#narrative-buttons',
    content: 'Story cues run preset What-Ifs, chemistry stubs, and copy updates for presenters.',
  },
  {
    target: '#heatmap-panel',
    content: 'Heatmap shows pH × temperature stability. Sweeps and runtime sims update this grid.',
  },
  {
    target: '#leaderboard',
    content: 'Top-5 ranking highlights robust picks. Click to sync Heatmap, 3D, and Inspector.',
  },
  {
    target: '#inspector-panel',
    content: 'Inspector surfaces scenario metrics, chem notes, and guidance pulled from rationale.',
  },
  {
    target: '#whatif-panel',
    content: 'Adjust pH/temperature windows to recalc scenario robustness and reorder the leaderboard.',
  },
  {
    target: '#show-data',
    content: 'Open the exact CSV/JSON powering the current view for transparency during demos.',
  },
  {
    target: '#preflight-indicator',
    content: 'Preflight badge shows which chem modules are real vs mocked in this environment.',
  },
];

export default tourSteps;
