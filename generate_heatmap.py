#!/usr/bin/env python3
"""
Generate heatmap data for antibody candidates using RDKit.
This script loads candidates from a CSV file, processes their SMILES strings,
and generates stability heatmap values across a range of pH and temperature values.
"""

import os
import csv
import json
import random
import numpy as np
from pathlib import Path
from typing import Dict, List, Any

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski
    from rdkit.Chem import rdMolStandardize
    HAS_RDKIT = True
except ImportError:
    print("Warning: RDKit not available. Using mock data generation.")
    HAS_RDKIT = False

# Constants
PH_VALUES = [6.8, 7.0, 7.4, 7.8, 8.2]
TEMP_VALUES = [25, 30, 37, 40, 45]
CANDIDATES_PATH = Path("backend/data/flow1_antibody/candidates.csv")
OUTPUT_PATH = Path("backend/data/flow1_antibody/generated_heatmap.json")

def load_candidates() -> List[Dict[str, str]]:
    """Load candidates from CSV file."""
    candidates = []
    with open(CANDIDATES_PATH, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append(row)
    return candidates

def calculate_stability(mol, ph: float, temp: float) -> float:
    """
    Calculate stability score (0-1) based on molecular properties and conditions.
    
    This uses a combination of RDKit descriptors and environmental factors to 
    estimate stability. For a real application, this would be replaced with
    actual experimental data or more sophisticated models.
    """
    # Base stability factors
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    # pH effects - molecules are generally more stable near neutral pH
    ph_optimum = 7.4
    ph_factor = 1.0 - (abs(ph - ph_optimum) / 2.0)  # Penalize deviation from optimum
    
    # Temperature effects - higher temps generally decrease stability
    temp_factor = 1.0 - ((temp - 25) / 40.0)  # Penalize higher temperatures
    
    # Molecular property effects
    property_factor = 0.0
    
    # LogP contribution: moderate logP (1-3) is generally good for stability
    if 1.0 <= logp <= 3.0:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, abs(logp - 2.0) / 3.0))
    
    # TPSA contribution: moderate TPSA often correlates with stability
    if 40 <= tpsa <= 120:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, abs(tpsa - 80) / 80))
    
    # H-bond donors/acceptors: balanced numbers are often better
    hb_balance = 1.0 - min(1.0, abs(hbd - hba) / max(1, hbd + hba))
    property_factor += 0.1 * hb_balance
    
    # Rotatable bonds: fewer is generally better for stability
    if rotatable_bonds <= 5:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, (rotatable_bonds - 5) / 10))
    
    # Combine factors with weights
    stability = (0.4 * ph_factor) + (0.3 * temp_factor) + (0.3 * property_factor)
    
    # Add a small random component to create variation
    random_factor = random.uniform(-0.05, 0.05)
    stability += random_factor
    
    # Ensure result is between 0 and 1
    return max(0.0, min(1.0, stability))

def generate_mock_stability(candidate_id: str, ph: float, temp: float) -> float:
    """Generate mock stability data when RDKit is not available."""
    # Use candidate ID as a seed for reproducibility
    seed = sum(ord(c) for c in candidate_id)
    random.seed(seed + int(ph * 10) + temp)
    
    # Base value between 0.5 and 0.9
    base = 0.5 + (int(candidate_id[1:]) % 5) / 10.0
    
    # pH effect - optimal around 7.4
    ph_effect = 1.0 - abs(ph - 7.4) / 2.0
    
    # Temperature effect - lower is better
    temp_effect = 1.0 - (temp - 25) / 40.0
    
    # Combine with some randomness
    stability = base * 0.4 + ph_effect * 0.3 + temp_effect * 0.3
    stability += random.uniform(-0.05, 0.05)
    
    return max(0.1, min(0.99, stability))

def process_candidate(candidate: Dict[str, str]) -> Dict[str, Any]:
    """Process a single candidate and generate its heatmap values."""
    candidate_id = candidate['id']
    smiles = candidate['smiles']
    result = {"id": candidate_id, "values": []}
    
    if HAS_RDKIT:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Warning: Could not parse SMILES for {candidate_id}: {smiles}")
                raise ValueError("Invalid SMILES")
                
            # Generate values for each temperature row
            for temp in TEMP_VALUES:
                row = []
                for ph in PH_VALUES:
                    stability = calculate_stability(mol, ph, temp)
                    row.append(round(stability, 2))
                result["values"].append(row)
                
        except Exception as e:
            print(f"Error processing {candidate_id}: {e}")
            # Fall back to mock data
            result["values"] = generate_mock_values(candidate_id)
    else:
        # Generate mock data if RDKit is not available
        result["values"] = generate_mock_values(candidate_id)
        
    return result

def generate_mock_values(candidate_id: str) -> List[List[float]]:
    """Generate mock stability values for a candidate."""
    values = []
    for temp in TEMP_VALUES:
        row = []
        for ph in PH_VALUES:
            stability = generate_mock_stability(candidate_id, ph, temp)
            row.append(round(stability, 2))
        values.append(row)
    return values

def main():
    """Main function to generate heatmap data."""
    candidates = load_candidates()
    print(f"Loaded {len(candidates)} candidates from {CANDIDATES_PATH}")
    
    results = {
        "temps": TEMP_VALUES,
        "phs": PH_VALUES,
        "candidates": []
    }
    
    for candidate in candidates:
        print(f"Processing {candidate['id']}: {candidate['name']}")
        result = process_candidate(candidate)
        results["candidates"].append(result)
    
    # Write results to JSON file
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Heatmap data written to {OUTPUT_PATH}")

if __name__ == "__main__":
    main()
