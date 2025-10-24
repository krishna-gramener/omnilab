#!/usr/bin/env python3
"""
Generate heatmap data for DNA probe and enzyme candidates using RDKit.
This script loads candidates from CSV files, processes their SMILES strings,
and generates stability heatmap values across specified pH and temperature ranges.
"""

import os
import csv
import json
import random
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Tuple

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
FLOW_CONFIGS = {
    "flow2_dnaprobe": {
        "candidates_path": Path("backend/data/flow2_dnaprobe/candidates.csv"),
        "output_path": Path("backend/data/flow2_dnaprobe/generated_heatmap.json"),
        "ph_values": [6.5, 6.8, 7.2, 7.6],
        "temp_values": [4, 10, 25, 37],
        "stability_weights": {
            "ph_optimum": 7.0,  # DNA probes prefer slightly acidic conditions
            "temp_optimum": 15,  # DNA probes often work best at lower temps
            "ph_weight": 0.35,
            "temp_weight": 0.35,
            "property_weight": 0.3,
            "random_range": (-0.04, 0.04)
        }
    },
    "flow3_enzyme": {
        "candidates_path": Path("backend/data/flow3_enzyme/candidates.csv"),
        "output_path": Path("backend/data/flow3_enzyme/generated_heatmap.json"),
        "ph_values": [7.0, 7.4, 8.0, 8.5],
        "temp_values": [4, 10, 18, 25],
        "stability_weights": {
            "ph_optimum": 7.4,  # Most enzymes prefer neutral to slightly basic pH
            "temp_optimum": 20,  # Enzymes often have optimal activity around room temp
            "ph_weight": 0.3,
            "temp_weight": 0.4,
            "property_weight": 0.3,
            "random_range": (-0.05, 0.05)
        }
    }
}

def load_candidates(path: Path) -> List[Dict[str, str]]:
    """Load candidates from CSV file."""
    candidates = []
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append(row)
    return candidates

def calculate_stability_dna_probe(mol, ph: float, temp: float, weights: Dict) -> float:
    """
    Calculate stability score (0-1) for DNA probes based on molecular properties and conditions.
    
    DNA probes are generally more stable at lower temperatures and neutral pH.
    """
    if not HAS_RDKIT or mol is None:
        return 0.5  # Fallback if RDKit is not available
    
    # Base stability factors
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    # pH effects - DNA probes are generally more stable near neutral pH
    ph_optimum = weights["ph_optimum"]
    ph_factor = 1.0 - (abs(ph - ph_optimum) / 2.0)  # Penalize deviation from optimum
    
    # Temperature effects - DNA probes generally prefer lower temperatures
    temp_optimum = weights["temp_optimum"]
    temp_factor = 1.0 - (abs(temp - temp_optimum) / 40.0)  # Penalize deviation from optimum
    
    # Molecular property effects
    property_factor = 0.0
    
    # LogP contribution: moderate logP is generally good for stability
    if 0.5 <= logp <= 2.5:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, abs(logp - 1.5) / 3.0))
    
    # TPSA contribution: higher TPSA often correlates with DNA probe stability
    if 60 <= tpsa <= 140:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, abs(tpsa - 100) / 80))
    
    # H-bond donors/acceptors: balanced numbers are often better
    hb_balance = 1.0 - min(1.0, abs(hbd - hba) / max(1, hbd + hba))
    property_factor += 0.1 * hb_balance
    
    # Rotatable bonds: fewer is generally better for stability
    if rotatable_bonds <= 4:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, (rotatable_bonds - 4) / 8))
    
    # Combine factors with weights
    stability = (
        weights["ph_weight"] * ph_factor + 
        weights["temp_weight"] * temp_factor + 
        weights["property_weight"] * property_factor
    )
    
    # Add a small random component to create variation
    random_factor = random.uniform(*weights["random_range"])
    stability += random_factor
    
    # Ensure result is between 0 and 1
    return max(0.1, min(0.99, stability))

def calculate_stability_enzyme(mol, ph: float, temp: float, weights: Dict) -> float:
    """
    Calculate stability score (0-1) for enzymes based on molecular properties and conditions.
    
    Enzymes are generally more sensitive to temperature and have specific pH optima.
    """
    if not HAS_RDKIT or mol is None:
        return 0.5  # Fallback if RDKit is not available
    
    # Base stability factors
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    # pH effects - enzymes have specific pH optima
    ph_optimum = weights["ph_optimum"]
    ph_factor = 1.0 - (abs(ph - ph_optimum) / 2.0)  # Penalize deviation from optimum
    
    # Temperature effects - enzymes are often sensitive to higher temperatures
    temp_optimum = weights["temp_optimum"]
    temp_factor = 1.0 - (abs(temp - temp_optimum) / 30.0)  # Penalize deviation from optimum
    
    # Molecular property effects
    property_factor = 0.0
    
    # LogP contribution: moderate logP is generally good for stability
    if 1.0 <= logp <= 3.5:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, abs(logp - 2.25) / 3.0))
    
    # TPSA contribution: moderate TPSA often correlates with enzyme stability
    if 50 <= tpsa <= 150:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, abs(tpsa - 100) / 100))
    
    # H-bond donors/acceptors: balanced numbers are often better
    hb_balance = 1.0 - min(1.0, abs(hbd - hba) / max(1, hbd + hba))
    property_factor += 0.1 * hb_balance
    
    # Rotatable bonds: fewer is generally better for stability
    if rotatable_bonds <= 6:
        property_factor += 0.2
    else:
        property_factor += 0.1 * (1.0 - min(1.0, (rotatable_bonds - 6) / 10))
    
    # Combine factors with weights
    stability = (
        weights["ph_weight"] * ph_factor + 
        weights["temp_weight"] * temp_factor + 
        weights["property_weight"] * property_factor
    )
    
    # Add a small random component to create variation
    random_factor = random.uniform(*weights["random_range"])
    stability += random_factor
    
    # Ensure result is between 0 and 1
    return max(0.1, min(0.99, stability))

def generate_mock_stability(candidate_id: str, ph: float, temp: float, flow_type: str, weights: Dict) -> float:
    """Generate mock stability data when RDKit is not available."""
    # Use candidate ID as a seed for reproducibility
    seed = sum(ord(c) for c in candidate_id)
    random.seed(seed + int(ph * 10) + temp)
    
    # Base value between 0.5 and 0.9
    base = 0.5 + (int(candidate_id[1:]) % 5) / 10.0
    
    # pH effect - optimal around specified optimum
    ph_optimum = weights["ph_optimum"]
    ph_effect = 1.0 - abs(ph - ph_optimum) / 2.0
    
    # Temperature effect - optimal around specified optimum
    temp_optimum = weights["temp_optimum"]
    temp_effect = 1.0 - abs(temp - temp_optimum) / 40.0
    
    # Combine with some randomness
    stability = (
        base * 0.4 + 
        ph_effect * weights["ph_weight"] + 
        temp_effect * weights["temp_weight"]
    )
    stability += random.uniform(*weights["random_range"])
    
    return max(0.1, min(0.99, stability))

def process_candidate(candidate: Dict[str, str], flow_type: str, ph_values: List[float], 
                     temp_values: List[float], weights: Dict) -> Dict[str, Any]:
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
            for temp in temp_values:
                row = []
                for ph in ph_values:
                    if flow_type == "flow2_dnaprobe":
                        stability = calculate_stability_dna_probe(mol, ph, temp, weights)
                    else:  # flow3_enzyme
                        stability = calculate_stability_enzyme(mol, ph, temp, weights)
                    row.append(round(stability, 2))
                result["values"].append(row)
                
        except Exception as e:
            print(f"Error processing {candidate_id}: {e}")
            # Fall back to mock data
            result["values"] = generate_mock_values(candidate_id, flow_type, ph_values, temp_values, weights)
    else:
        # Generate mock data if RDKit is not available
        result["values"] = generate_mock_values(candidate_id, flow_type, ph_values, temp_values, weights)
        
    return result

def generate_mock_values(candidate_id: str, flow_type: str, ph_values: List[float], 
                        temp_values: List[float], weights: Dict) -> List[List[float]]:
    """Generate mock stability values for a candidate."""
    values = []
    for temp in temp_values:
        row = []
        for ph in ph_values:
            stability = generate_mock_stability(candidate_id, ph, temp, flow_type, weights)
            row.append(round(stability, 2))
        values.append(row)
    return values

def process_flow(flow_type: str, config: Dict) -> None:
    """Process all candidates for a specific flow."""
    candidates = load_candidates(config["candidates_path"])
    print(f"Loaded {len(candidates)} candidates from {config['candidates_path']}")
    
    results = {
        "temps": config["temp_values"],
        "phs": config["ph_values"],
        "candidates": []
    }
    
    for candidate in candidates:
        print(f"Processing {candidate['id']}: {candidate['name']}")
        result = process_candidate(
            candidate, 
            flow_type, 
            config["ph_values"], 
            config["temp_values"],
            config["stability_weights"]
        )
        results["candidates"].append(result)
    
    # Write results to JSON file
    os.makedirs(os.path.dirname(config["output_path"]), exist_ok=True)
    with open(config["output_path"], 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Heatmap data written to {config['output_path']}")

def main():
    """Main function to generate heatmap data for multiple flows."""
    # Set a fixed seed for reproducibility
    random.seed(42)
    
    for flow_type, config in FLOW_CONFIGS.items():
        print(f"\nProcessing {flow_type}...")
        process_flow(flow_type, config)

if __name__ == "__main__":
    main()
