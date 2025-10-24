"""Chemistry utilities with graceful fallbacks for the demo pipeline."""

from .smiles import load_smiles, embed_3d, protonate_for_ph, has_rdkit
from .propka_bridge import protein_pka, has_propka
from .docking_stub import dock
from .scoring import mmgbsa_like

__all__ = [
    "load_smiles",
    "embed_3d",
    "protonate_for_ph",
    "has_rdkit",
    "has_propka",
    "protein_pka",
    "dock",
    "mmgbsa_like",
]
