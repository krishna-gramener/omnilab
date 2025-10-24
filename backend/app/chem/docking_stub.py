from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

from .smiles import has_rdkit

_DATA_ROOT = Path(__file__).resolve().parent.parent.parent / "data"


def _mock_score(candidate_id: str | None, center: Iterable[float]) -> float:
    token = f"{candidate_id}:{','.join(map(lambda x: f'{x:.2f}', center))}"
    digest = hashlib.sha256(token.encode()).digest()[0]
    return round(-6.5 - (digest / 255.0) * 2.5, 2)


def _load_pose(flow_id: str | None, candidate_id: str | None) -> str:
    if not flow_id or not candidate_id:
        return "REMARK MOCKED POSE\n"
    pose_path = _DATA_ROOT / flow_id / "poses" / f"{candidate_id}.pdbqt"
    if pose_path.exists():
        return pose_path.read_text()
    return f"REMARK GENERATED PLACEHOLDER FOR {candidate_id}\n"


def dock(
    ligand_mol3d: Any,
    protein_pdb: str | Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    *,
    flow_id: str | None = None,
    candidate_id: str | None = None,
) -> Dict[str, Any]:
    """Run a docking stub or return mocked data depending on feature flags."""

    want_real = os.getenv("USE_REAL_DOCKING", "false").lower() == "true"

    if not want_real:
        score = _mock_score(candidate_id, center)
        pose = _load_pose(flow_id, candidate_id)
        return {"pose_pdbqt": pose, "score": score, "source": "mock"}

    if not has_rdkit():
        raise RuntimeError("RDKit required for real docking but unavailable")

    # Placeholder for real docking integration; in demo mode we return a synthetic
    # value but tag the source as 'stub'.
    score = _mock_score(candidate_id, center)
    pose = _load_pose(flow_id, candidate_id)
    return {"pose_pdbqt": pose, "score": score, "source": "stub"}
