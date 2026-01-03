import sys
import os
from rdkit import Chem

# Add project root to path
sys.path.append(os.getcwd())

from chemtools.taxonomy.dynamic_classifier import DynamicClassifier

def run_poc():
    data_path = "chemtools/taxonomy/data"
    classifier = DynamicClassifier(data_path)
    
    print("=== Chemical Motif Specificity PoC ===")
    print(f"Loading taxonomy from: {data_path}")
    
    # Example 1: Chlorobenzene
    # It matches: Any-X, Csp2-X, Ar-X, Ar-Cl
    print("\nScenario: Classifying Chlorobenzene (c1ccccc1Cl)")
    
    # Simulated matches for the C-Cl bond
    bond_atoms = (0, 6) # Carbon 0 and Chlorine 6
    
    matches = [
        {"id": "Any-X",   "A": "Any_Scaffold", "B": "X"},
        {"id": "Csp2-X",  "A": "sp2_Carbon",   "B": "X"},
        {"id": "Ar-X",    "A": "Ar",           "B": "X"},
        {"id": "Ar-Cl",   "A": "Ar",           "B": "Cl"},
    ]
    
    scored_matches = []
    for m in matches:
        score = classifier.get_compound_score(m)
        scored_matches.append({**m, "score": score})
        print(f"  Detected: {m['id']:10} | Score: {score}")
        
    # Filtering logic
    max_score = max(m["score"] for m in scored_matches)
    winners = [m["id"] for m in scored_matches if m["score"] == max_score]
    
    print(f"\nResult: The most specific classification is: {winners}")
    print("-" * 40)

    # Example 2: Benzyl Bromide
    print("\nScenario: Classifying Benzyl Bromide (c1ccccc1CH2Br)")
    matches_2 = [
        {"id": "Any-X",   "A": "Any_Scaffold", "B": "X"},
        {"id": "R-X",     "A": "R",            "B": "X"},
        {"id": "R-Br",    "A": "R",            "B": "Br"},
        {"id": "Bn-Br",   "A": "Bn",           "B": "Br"},
    ]
    
    scored_matches_2 = []
    for m in matches_2:
        score = classifier.get_compound_score(m)
        scored_matches_2.append({**m, "score": score})
        print(f"  Detected: {m['id']:10} | Score: {score}")
        
    max_score_2 = max(m["score"] for m in scored_matches_2)
    winners_2 = [m["id"] for m in scored_matches_2 if m["score"] == max_score_2]
    
    print(f"\nResult: The most specific classification is: {winners_2}")
    print("Note: Bn-Br wins because 'Bn' (priority 2) is more specific than 'R' (priority 1).")

if __name__ == "__main__":
    run_poc()
