import sys
import os
from rdkit import Chem
from typing import List, Dict, Any

# Add project root to path
sys.path.append(os.getcwd())

from chemtools.taxonomy.dynamic_classifier import DynamicClassifier
from examples.sample_compounds import ARYL_HALIDES

def test_samples():
    data_path = "chemtools/taxonomy/data"
    classifier = DynamicClassifier(data_path)
    
    # We'll test a few representative samples from ARYL_HALIDES
    # and some manual ones to cover different cases.
    samples = [
        ARYL_HALIDES[0],  # Bromobenzene
        ARYL_HALIDES[1],  # Chlorobenzene
        ARYL_HALIDES[4],  # 4-Bromobenzonitrile
        {
            "name": "Benzyl bromide",
            "smiles": "BrCc1ccccc1"
        },
        {
            "name": "Phenylboronic acid",
            "smiles": "c1ccc(B(O)O)cc1"
        },
        {
            "name": "Phenyl triflate",
            "smiles": "c1ccc(OS(=O)(=O)C(F)(F)F)cc1"
        }
    ]

    print(f"{'Compound Name':<25} | {'SMILES':<20} | {'Best Classification':<20} | {'Score'}")
    print("-" * 90)

    for s in samples:
        name = s["name"]
        smiles = s["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
            
        # Simulate detection logic
        # We define potential matches for the primary functional group of each sample.
        # The goal is to show that the classifier picks the most specific one.
        
        matches = []
        
        if name == "Bromobenzene":
            matches = [
                {"id": "Any-X", "A": "Any_Scaffold", "B": "X"},
                {"id": "Csp2-X", "A": "sp2_Carbon", "B": "X"},
                {"id": "Ar-Br", "A": "Ar", "B": "Br"},
                {"id": "Ar-X", "A": "Ar", "B": "X"},
            ]
        elif name == "Chlorobenzene":
            matches = [
                {"id": "Any-X", "A": "Any_Scaffold", "B": "X"},
                {"id": "Csp2-X", "A": "sp2_Carbon", "B": "X"},
                {"id": "Ar-Cl", "A": "Ar", "B": "Cl"},
                {"id": "Ar-X", "A": "Ar", "B": "X"},
            ]
        elif name == "4-Bromobenzonitrile":
            matches = [
                {"id": "Any-X", "A": "Any_Scaffold", "B": "X"},
                {"id": "Ar-Br", "A": "Ar", "B": "Br"},
                {"id": "Ar-X", "A": "Ar", "B": "X"},
            ]
        elif name == "Benzyl bromide":
            matches = [
                {"id": "Any-X", "A": "Any_Scaffold", "B": "X"},
                {"id": "R-Br", "A": "R", "B": "Br"},
                {"id": "Bn-Br", "A": "Bn", "B": "Br"},
            ]
        elif name == "Phenylboronic acid":
            matches = [
                {"id": "Csp2-B(OH)2", "A": "sp2_Carbon", "B": "B(OH)2"},
                {"id": "Ar-B(OH)2", "A": "Ar", "B": "B(OH)2"},
            ]
        elif name == "Phenyl triflate":
            matches = [
                {"id": "Any-Sulfonate", "A": "Any_Scaffold", "B": "Sulfonate"},
                {"id": "Csp2-Sulfonate", "A": "sp2_Carbon", "B": "Sulfonate"},
                {"id": "Ar-OTf", "A": "Ar", "B": "OTf"},
            ]

        # Score and find the best
        best_match = None
        max_score = -1
        
        for m in matches:
            score = classifier.get_compound_score(m)
            if score > max_score:
                max_score = score
                best_match = m
        
        if best_match:
            print(f"{name:<25} | {smiles:<20} | {best_match['id']:<20} | {max_score}")

if __name__ == "__main__":
    test_samples()
