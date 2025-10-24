# Reagent Robustness Simulator

Sci‑fi styled demo that ranks reagent candidates by robustness across pH × temperature, visualises stability heatmaps and 3D binding playback, and explains the top picks per scenario. Ships with three precomputed flows:

- **Antibody kit for hot/humid markets**
- **DNA probe for newborn screening**
- **Enzyme assay for cold-chain shipments**

No heavy compute runs at demo time—everything is backed by static mock data served from a single FastAPI service that also serves the pre-built React frontend.

---

## Project Layout

```
/backend
  ├─ app/main.py            # FastAPI entrypoint + dataset endpoints
  └─ data/…                 # flows.json + per-flow CSV/JSON mock payloads
/frontend
  ├─ index.html             # Vite boot + NGL CDN
  └─ src/
       ├─ App.jsx           # Orchestration, narratives, what-if logic
       ├─ components/…      # HUD panels, heatmap, NGL viewer, modal, etc.
       └─ lib/api.js        # REST client for FastAPI endpoints
```

---

## Prerequisites

- **Python 3.10+** with `pip`
- **Node.js 18+** with `npm`
- Optional: PROPKA binary / Python package if you want live residue pKa calls.

---

## Backend Setup (single service)

```bash
cd backend
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Build the frontend (first run or whenever the UI changes):

```bash
npm --prefix ../frontend install   # safe to re-run; installs only when needed
npm --prefix ../frontend run build
```

Then launch the unified API + UI service:

```bash
uvicorn app.main:app --reload --port 8080
```

Navigate to <http://localhost:8080> to use the app—the FastAPI process now serves both the API routes (`/api/...`) and the compiled frontend.

> NOTE: the React app now points to the backend relative to the current origin. If you ever host the API elsewhere, set `VITE_API_BASE` during the Vite build (e.g., `VITE_API_BASE=https://api.example.com npm run build`).

### ChemOps feature flags

Environment variables let you opt into heavier tooling when available:

- `USE_REAL_PROPka=true` — attempt PROPKA; falls back to mocked residue pKa map if binaries are missing.
- `USE_REAL_DOCKING=true` — placeholder for AutoDock Vina / DiffDock (currently returns deterministic stubs unless tooling is present).
- Defaults (`false`) keep everything mock-backed so the click-through demo remains instant.

### Available API Routes

- `GET /api/flows` — list flow descriptors
- `GET /api/flows/{flow_id}/candidates` — 12-row CSV as JSON records
- `GET /api/flows/{flow_id}/scores` — ranked robustness metadata
- `GET /api/flows/{flow_id}/heatmap` — pH × temperature grids
- `GET /api/flows/{flow_id}/rationale` — top-5 “why + watch” copy
- `GET /api/flows/{flow_id}/raw/{name}` — exact JSON/CSV for “Show Data” modal
- `GET /api/flows/{flow_id}/preflight` — reports RDKit / PROPKA / docking availability and active mode
- `POST /api/simulate` — recompute a single pH×Temp cell with RDKit protonation, docking stub, and MM/GBSA-like proxy (falls back to canned heuristics)

Each flow folder also contains an `ngl/` directory. Drop real PDB assets there if you want NGL to load them; otherwise the Three.js fallback renders a torus-knot target and emissive ligand sphere.

---

## Frontend Development (optional)

```bash
cd frontend
npm install
npm run dev      # serves at http://localhost:5173 with proxy to :8080
```

Use this only when you want hot module reloading during UI development. The dev server proxies `/api` requests to `http://localhost:8080`, so keep the FastAPI process running while the frontend is open. Remember to run `npm run build` afterwards so the FastAPI service serves the latest compiled assets.

### Scenario upload (CSV)

- Click **Upload Scenario** in the top bar to import a custom stability grid.
- Provide a CSV containing the columns `candidate`, `temp`, `ph`, and `value` (values are clamped to 0–1). Temperatures and pH points must match the active heatmap grid. If the `candidate` column is omitted, the currently selected candidate is overwritten.
- Only the referenced candidate(s) are overridden; select *Baseline* in the narrative panel to revert to the original dataset.
- Use the **Sample Scenario** shortcut to load `frontend/public/sample_scenario.csv`; the first few rows are previewed before/after the overrides are applied.

### Guided tour & demo mode

- First-time visitors are walked through the UI with a Joyride tour (Flow tabs → narrative buttons → heatmap → leaderboard → inspector → What-If sliders → Show Data → Preflight).
- “Start tour” stays available in the TopBar until dismissed (state persisted via `localStorage`).
- Toggle **Demo Mode** to surface presenter steps (Baseline, Hot Day, Mg²⁺ boost, etc.) with toast call-outs.
- Heatmap borders pulse, leaderboard rank delta tags render (`▲` / `▼`), and the Inspector shows a change table whenever a narrative or What-If slider updates the scenario.
- Use the **37 °C Sweep** action to hit `/api/simulate` for the active candidate and narrate “computed at runtime using RDKit protonation & docking stub.”
- Flow banner, KPI minis, ROI gear, and MAC toggle quantify impact for Revvity stakeholders; Wet-lab Handoff panel bridges simulation to bench execution.

---

## One-touch Startup

After installing prerequisites, you can launch the single FastAPI service (which also builds the frontend) with the helper script:

```bash
chmod +x start_services.sh
./start_services.sh
```

The script will:

- Create/refresh the backend virtualenv and install `requirements.txt`.
- Install frontend dependencies (first run only).
- Build the Vite frontend once per run.
- Start the FastAPI backend (serving both API and static assets) on <http://localhost:8080>.

Use `Ctrl+C` in the terminal to stop the process.

---

## Demo Playbook

1. Choose a flow from the top HUD tabs.
2. Explore the 3D binding panel (looping wobble) and heatmap sweep slider.
3. Review the neon leaderboard and inspector metrics (temp span, pH window, scenario average).
4. Use the bottom **What-if** sliders to restrict pH/temperature ranges; watch the leaderboard flip into a **Scenario Ranking** without server calls.
5. Trigger narrative buttons per flow:
   - **Antibody:** Hot Day 45 °C focus, shipping advice copy.
   - **DNA Probe:** Mg²⁺ buffer tweak (+0.03 for compatible flags), cold-room window.
   - **Enzyme:** Lost cold chain sweep, spec-sheet draft blurb.
6. Hit **37 °C Sweep** to call the runtime chem simulate and narrate the RDKit/docking provenance.
7. Click **Show Data** to pop the Radix modal with the exact CSV/JSON payloads.
8. Use **Score Breakdown** to expose the temperature/pH contributions behind each ranking, and inspect the inspector’s sparkline + grid to see per-cell deltas.

All data loads in under one second, and the UI provides hover tooltips with sparkline mini-trends for each heatmap cell.

---

## Testing & QA Notes

- No automated tests are wired (static datasets + client-side transforms).
- Manual checks:
  - Each flow loads instantly with correct metadata.
  - Leaderboard selections sync with the heatmap and 3D panel.
  - “Show Data” surfaces backend payloads verbatim.
  - What-if sliders and narrative actions update rankings client-side.
  - `/api/flows/{id}/preflight` reflects the environment flags; `/api/simulate` responds in under ~250 ms via real or mocked chem paths.

---

## Extending the Demo

- Swap in real assay data by updating the CSV/JSON files under `/backend/data/<flow>/`.
- Add more flows by registering them in `flows.json` and duplicating the dataset folder structure.
- Integrate actual NGL molecular views by placing `target.pdb` and `cand_*.pdb` files inside each flow’s `ngl/` directory.
- Layer in analytics (e.g., area under curve, stability thresholds) by expanding `computeMetrics` in `App.jsx`.

---

## Make Targets

- `make backend` — installs frontend deps, builds the UI, and starts the FastAPI service (set `USE_REAL_PROPka` / `USE_REAL_DOCKING` as needed).
- `make frontend` — runs the Vite dev server with hot reload (requires the backend to be running separately).
- `make demo-on` — shorthand for `USE_REAL_PROPka=false USE_REAL_DOCKING=false make backend`.
- `make demo-propka` — same, but enables `USE_REAL_PROPka=true`.

Enjoy the neon lab console! 🚀
