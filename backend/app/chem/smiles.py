from __future__ import annotations

import os
from typing import Optional

try:  # pragma: no cover - runtime availability guard
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdMolStandardize

    _HAS_RDKIT = True
except ImportError:  # pragma: no cover - graceful fallback
    Chem = None  # type: ignore
    AllChem = None  # type: ignore
    rdMolStandardize = None  # type: ignore
    _HAS_RDKIT = False


def has_rdkit() -> bool:
    """Return True when RDKit is importable."""

    return _HAS_RDKIT


def load_smiles(smiles: str):
    """Parse a SMILES string into an RDKit Mol, raising if unavailable.

    Returns None if RDKit is missing and `USE_REAL_DOCKING` flag is false to let the
    caller fall back to canned data quietly.
    """

    if not _HAS_RDKIT:
        if os.getenv("USE_REAL_DOCKING", "false").lower() == "true":
            raise RuntimeError("RDKit not available but real docking requested")
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Unable to parse SMILES: {smiles}")
    return mol


def embed_3d(mol):
    """Generate a 3D conformer using ETKDG; returns a copy.

    When RDKit is missing, returns None to let upstream use canned coordinates.
    """

    if not _HAS_RDKIT or mol is None:
        return None

    mol3d = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    AllChem.EmbedMolecule(mol3d, params)
    AllChem.UFFOptimizeMolecule(mol3d, maxIters=200)
    return mol3d


def protonate_for_ph(mol3d, ph: float):
    """Perform crude protonation adjustments for the requested pH.

    Uses RDKit's rdMolStandardize module. Falls back to the input molecule if RDKit
    (or standardizer) is missing. This is intentionally lightweight for demo use.
    """

    if not _HAS_RDKIT or mol3d is None or rdMolStandardize is None:
        return mol3d

    adjuster = rdMolStandardize.Uncharger()
    clean = rdMolStandardize.Cleanup(mol3d)
    reionized = rdMolStandardize.Reionize(clean)
    neutral = adjuster.uncharge(reionized)

    # Basic heuristic: for acidic pH, add hydrogens; for basic, remove some.
    ph = float(ph)
    if ph < 6.5:
        neutral = Chem.AddHs(neutral, explicitOnly=True)
    elif ph > 8.5:
        neutral = Chem.RemoveHs(neutral)

    return neutral
