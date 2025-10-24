from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
import logging
from pathlib import Path
import json
import csv
import subprocess

from pydantic import BaseModel

from .chem import (
    load_smiles,
    embed_3d,
    protonate_for_ph,
    protein_pka,
    has_propka,
    has_rdkit,
    dock,
    mmgbsa_like,
)

app = FastAPI(title="Reagent Robustness API", version="1.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

DATA = Path(__file__).resolve().parent.parent / "data"
FRONTEND_ROOT = Path(__file__).resolve().parents[2] / "frontend"
FRONTEND_DIST = FRONTEND_ROOT / "dist"
_LOGGER = logging.getLogger("revvity.single_app")


def _maybe_build_frontend() -> None:
    if os.getenv("AUTO_BUILD_FRONTEND", "true").lower() == "false":
        _LOGGER.info("Skipping frontend build because AUTO_BUILD_FRONTEND=false")
        return

    if not FRONTEND_ROOT.exists():
        _LOGGER.warning("Frontend directory %s missing; cannot build assets.", FRONTEND_ROOT)
        return

    need_build = False
    if not FRONTEND_DIST.exists() or not any(FRONTEND_DIST.iterdir()):
        need_build = True
    else:
        try:
            dist_mtime = max((p.stat().st_mtime for p in FRONTEND_DIST.rglob("*")), default=0)
            src_mtime = max((p.stat().st_mtime for p in FRONTEND_ROOT.joinpath("src").rglob("*")), default=0)
            lock_file = FRONTEND_ROOT / "package-lock.json"
            lock_mtime = lock_file.stat().st_mtime if lock_file.exists() else 0
            if src_mtime > dist_mtime or lock_mtime > dist_mtime:
                need_build = True
        except OSError:
            need_build = True

    if not need_build:
        _LOGGER.info("Frontend assets up to date; skipping build.")
        return

    install_cmd = ["npm", "--prefix", str(FRONTEND_ROOT), "install"]
    build_cmd = ["npm", "--prefix", str(FRONTEND_ROOT), "run", "build"]

    try:
        _LOGGER.info("Installing frontend dependencies…")
        subprocess.run(install_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        _LOGGER.warning("Failed to install frontend dependencies: %s", exc)

    try:
        _LOGGER.info("Building frontend assets via npm run build…")
        subprocess.run(build_cmd, check=True)
        _LOGGER.info("Frontend assets built successfully.")
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        _LOGGER.error("Unable to build frontend assets automatically: %s", exc)


_maybe_build_frontend()


def _load_json(path: Path):
    return json.loads(path.read_text())


def _load_csv(path: Path):
    with path.open() as f:
        return list(csv.DictReader(f))


@app.get("/api/flows")
def list_flows():
    return _load_json(DATA / "flows.json")


@app.get("/api/flows/{flow_id}/candidates")
def get_candidates(flow_id: str):
    return _load_csv(DATA / f"{flow_id}/candidates.csv")


@app.get("/api/flows/{flow_id}/scores")
def get_scores(flow_id: str):
    return _load_json(DATA / f"{flow_id}/scores.json")


@app.get("/api/flows/{flow_id}/heatmap")
def get_heatmap(flow_id: str):
    return _load_json(DATA / f"{flow_id}/heatmap.json")


@app.get("/api/flows/{flow_id}/rationale")
def get_rationale(flow_id: str):
    return _load_json(DATA / f"{flow_id}/rationale.json")


@app.get("/api/flows/{flow_id}/raw/{name}")
def get_raw(flow_id: str, name: str):
    path = DATA / f"{flow_id}/{name}"
    if not path.exists():
        return {"error": "not_found"}
    if path.suffix == ".json":
        return _load_json(path)
    if path.suffix == ".csv":
        return _load_csv(path)
    return {"error": "unsupported"}


def _env_flag(key: str, default: str = "false") -> bool:
    return os.getenv(key, default).lower() == "true"


def _flow_exists(flow_id: str) -> bool:
    return (DATA / flow_id).exists()


def _find_candidate(flow_id: str, candidate_id: str) -> Optional[Dict[str, Any]]:
    rows = _load_csv(DATA / f"{flow_id}/candidates.csv")
    for row in rows:
        if row["id"] == candidate_id:
            return row
    return None


def _heatmap_lookup(flow_id: str, candidate_id: str, ph: float, temp: float) -> Dict[str, Any]:
    heat = _load_json(DATA / f"{flow_id}/heatmap.json")
    temps: List[float] = heat["temps"]
    phs: List[float] = heat["phs"]
    candidate = next((c for c in heat["candidates"] if c["id"] == candidate_id), None)
    if candidate is None:
        raise KeyError(f"Candidate {candidate_id} not present in heatmap")
    values: List[List[float]] = candidate["values"]
    try:
        t_index = temps.index(temp)
        p_index = phs.index(ph)
    except ValueError:
        # If the requested point is not on the grid, use the closest indices.
        t_index = min(range(len(temps)), key=lambda i: abs(temps[i] - temp))
        p_index = min(range(len(phs)), key=lambda i: abs(phs[i] - ph))
    base_value = values[t_index][p_index]
    return {"value": base_value, "temp_index": t_index, "ph_index": p_index, "temps": temps, "phs": phs}


def _protein_path(flow_id: str) -> Path:
    flow_target = DATA / flow_id / "ngl" / "target.pdb"
    if flow_target.exists():
        return flow_target
    placeholder = DATA / "common" / "placeholder_target.pdb"
    if placeholder.exists():
        return placeholder
    raise FileNotFoundError("No protein target available")


class SimulationRequest(BaseModel):
    flow_id: str
    candidate_id: str
    ph: float
    temp: float


@app.get("/api/flows/{flow_id}/preflight")
def get_preflight(flow_id: str):
    rdkit_active = has_rdkit()
    propka_flag = _env_flag("USE_REAL_PROPka")
    propka_active = propka_flag and has_propka()
    docking_flag = _env_flag("USE_REAL_DOCKING")
    docking_active = docking_flag and rdkit_active

    labels = []
    if rdkit_active:
        labels.append("RDKit ✔")
    else:
        labels.append("RDKit ⚠︎ mocked")
    labels.append("PROPKA ✔" if propka_active else ("PROPKA ⚠︎ mock" if propka_flag else "PROPKA ⏻ off"))
    labels.append("Docking ✔" if docking_active else ("Docking ⚠︎ stub" if docking_flag else "Docking ⏻ mock"))

    if propka_active and not propka_flag:
        propka_summary = "Mocked PROPKA"
    elif propka_active:
        propka_summary = "Real PROPKA"
    elif propka_flag and not propka_active:
        propka_summary = "Requested PROPKA but binary missing"
    else:
        propka_summary = "PROPKA mocked"

    summary = "Mocked ChemOps"
    if rdkit_active:
        summary = "RDKit enabled"
    if docking_active:
        summary += " + Docking stub"
    summary += f" | {propka_summary}"

    return {
        "flow_id": flow_id,
        "rdkit": rdkit_active,
        "propka_requested": propka_flag,
        "propka_active": propka_active,
        "docking_requested": docking_flag,
        "docking_active": docking_active,
        "labels": labels,
        "summary": summary,
    }


@app.post("/api/simulate")
def simulate_cell(payload: SimulationRequest):
    if not _flow_exists(payload.flow_id):
        raise HTTPException(status_code=404, detail="Flow not found")

    candidate = _find_candidate(payload.flow_id, payload.candidate_id)
    if candidate is None:
        raise HTTPException(status_code=404, detail="Candidate not found")

    heat_info = _heatmap_lookup(payload.flow_id, payload.candidate_id, payload.ph, payload.temp)
    base_value = float(heat_info["value"])

    smiles = candidate.get("smiles")
    ligand = None
    protonated = None
    if smiles:
        ligand = load_smiles(smiles)
        protonated = protonate_for_ph(embed_3d(ligand), payload.ph)

    protein_path = _protein_path(payload.flow_id)
    try:
        docking_result = dock(
            protonated,
            protein_path,
            center=(0.0, 0.0, 0.0),
            size=(20.0, 20.0, 20.0),
            flow_id=payload.flow_id,
            candidate_id=payload.candidate_id,
        )
    except RuntimeError:
        docking_result = {"pose_pdbqt": "REMARK fallback", "score": -7.2, "source": "mock"}
    dock_score = float(docking_result.get("score", -7.2))

    # Quick heuristic features based on temperature / pH drift.
    temp_penalty = abs(payload.temp - heat_info["temps"][heat_info["temp_index"]]) / max(heat_info["temps"])
    ph_penalty = abs(payload.ph - heat_info["phs"][heat_info["ph_index"]]) / max(heat_info["phs"])
    rmsd = round(0.6 + temp_penalty * 1.4 + ph_penalty * 0.6, 2)
    hbonds = max(1, int(base_value * 5))

    proxy_value = mmgbsa_like({"dock_score": dock_score, "rmsd": rmsd, "hbonds": hbonds})
    blended = round(min(0.99, 0.6 * base_value + 0.4 * proxy_value), 3)

    return {
        "value": blended,
        "details": {
            "dock_score": dock_score,
            "rmsd": rmsd,
            "hbonds": hbonds,
            "source": docking_result.get("source", "mock"),
            "base_grid": base_value,
        },
    }


if FRONTEND_DIST.exists():
    app.mount("/", StaticFiles(directory=FRONTEND_DIST, html=True), name="frontend")
else:
    _LOGGER.warning(
        "Frontend build directory %s missing. Run `npm install` and `npm run build` inside frontend/.",
        FRONTEND_DIST,
    )
