from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict

try:  # pragma: no cover - optional dependency
    from propka import run as propka_run
    _HAS_PROPKA = True
except ImportError:  # pragma: no cover
    propka_run = None  # type: ignore
    _HAS_PROPKA = False


_MOCK_VALUES = {
    "ASP45": 3.9,
    "GLU67": 4.2,
    "HIS112": 6.8,
    "TYR155": 9.8,
}


def has_propka() -> bool:
    return _HAS_PROPKA


def protein_pka(pdb_path: str) -> Dict[str, float]:
    """Return a residue→pKa map using PROPKA when enabled.

    If USE_REAL_PROPka is not truthy or PROPKA is unavailable, returns a canned dict
    so the application remains functional.
    """

    want_real = os.getenv("USE_REAL_PROPka", "false").lower() == "true"
    pdb_file = Path(pdb_path)

    if not want_real or not _HAS_PROPKA or not pdb_file.exists():
        return _MOCK_VALUES.copy()

    try:
        # PROPKA's Python API expects argv-style parameters. We capture output JSON
        # into a temporary directory and parse the result.
        output_dir = pdb_file.parent / "_propka_cache"
        output_dir.mkdir(exist_ok=True)
        args = ["--json", str(pdb_file), "--out", str(output_dir)]
        propka_run.main(args)  # type: ignore[arg-type]
        json_path = next(output_dir.glob("*.json"))
        data = json.loads(json_path.read_text())
        residues = data.get("residues", [])
        return {f"{res['residue_name']}{res['residue_number']}": float(res["pka"]) for res in residues}
    except Exception:  # pragma: no cover - fallback, PROPKA optional
        return _MOCK_VALUES.copy()
