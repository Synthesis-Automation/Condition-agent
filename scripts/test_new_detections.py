import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path.cwd()))

from chemtools import detect_reaction

def test_new_detections():
    test_cases = [
        {
            "name": "Alcohol Oxidation",
            "smiles": "OCC1CCCCC1>>O=CC1CCCCC1",
            "expected": "oxidation_alcohol_to_carbonyl"
        },
        {
            "name": "Carbonyl Reduction",
            "smiles": "O=CC1CCCCC1>>OCC1CCCCC1",
            "expected": "reduction_carbonyl_to_alcohol"
        },
        {
            "name": "Nitro Reduction",
            "smiles": "O=[N+]([O-])c1ccccc1>>Nc1ccccc1",
            "expected": "reduction_nitro_to_amine"
        },
        {
            "name": "Aromatic Halogenation",
            "smiles": "c1ccccc1.BrBr>>Brc1ccccc1",
            "expected": "halogenation_aromatic"
        }
    ]
    
    print(f"{'Test Case':<30} | {'Detected ID':<30} | {'Status'}")
    print("-" * 75)
    
    for case in test_cases:
        result = detect_reaction(case["smiles"])
        detected_id = result.get("family", "None")
        status = "PASS" if detected_id == case["expected"] else "FAIL"
        print(f"{case['name']:<30} | {detected_id:<30} | {status}")

if __name__ == "__main__":
    test_new_detections()
