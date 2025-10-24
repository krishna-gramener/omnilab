from __future__ import annotations

from math import exp
from typing import Mapping


def mmgbsa_like(features: Mapping[str, float | int]) -> float:
    """Compute a bounded robustness proxy from docking and structural features."""

    dock_score = float(features.get("dock_score", -7.0))
    rmsd = float(features.get("rmsd", 1.0))
    hbonds = float(features.get("hbonds", 3))

    # Convert docking score (negative) into positive contribution.
    dock_term = max(0.0, min(1.0, 1 / (1 + exp(0.8 * (dock_score + 7.0)))))
    rmsd_term = max(0.0, min(1.0, exp(-rmsd / 2.5)))
    hb_term = max(0.0, min(1.0, hbonds / 6.0))

    value = 0.55 + 0.25 * dock_term + 0.15 * rmsd_term + 0.1 * hb_term
    return round(min(0.99, max(0.0, value)), 3)
