
import sys
from pathlib import Path
from chemtools.rule import RuleEngine
from chemtools.taxonomy.rule_db import resolve_rule_db_v2
from chemtools import detect_reaction

TEST_REACTIONS = [
    # Sonogashira
    ("Ic1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1", "Sonogashira"),
    ("Brc1ccc(C(F)(F)F)cc1.C#CCCCCC>>C(#CCCCCC)c1ccc(C(F)(F)F)cc1", "Sonogashira"),
    
    # Amide Formation
    ("O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1", "Amide Formation"),
    ("O=C(Cl)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1", "Amide Formation"),
    
    # Reductive Amination
    ("O=Cc1ccccc1.Nc1ccccc1>>c1ccc(NCc2ccccc2)cc1", "Reductive Amination"),
    ("O=C1CCCCC1.N1CCCCC1>>c1ccc(N2CCCCC2)cc1", "Reductive Amination"),
    
    # RCM
    ("C=CCCC=C>>C1=CCCC1", "RCM"),
    ("C=CCCCCC=C>>C1=CCCCC1", "RCM"),
    
    # Suzuki
    ("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1", "Suzuki"),
    ("Clc1ccc(C(F)(F)F)cc1.OB(O)c1ccncc1>>c1ccc(-c2ccncc2)cc1", "Suzuki"),
    
    # Cu C-N
    ("Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1", "Cu C-N"),
    ("Brc1ccccc1.n1cc[nH]c1>>c1ccc(-n2ccnc2)cc1", "Cu C-N"),
    
    # SNAr
    ("Fc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>O=[N+]([O-])c1ccc(Nc2ccccc2)cc1", "SNAr"),
    ("Clc1ccnc(Cl)c1.NCCCC>>c1cc(NCCCC)nc(Cl)c1", "SNAr"),
    
    # C-O Coupling
    ("Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1", "C-O Coupling")
]

def test_families():
    print("==================================================")
    print("CROSS-FAMILY RULE DATABASE TEST")
    print("==================================================\n")
    
    success_count = 0
    total = len(TEST_REACTIONS)
    
    for smiles, expected_label in TEST_REACTIONS:
        print(f"Testing {expected_label}: {smiles}")
        try:
            # Detect reaction type
            detection = detect_reaction(smiles)
            family = detection.get("family", "unknown")
            print(f"  Detected Family: {family}")
            
            # Resolve rule database
            db_stem = resolve_rule_db_v2(family)
            if not db_stem:
                print(f"  [FAIL] No rule database found for family: {family}")
                continue
            
            db_path = Path("data/rule_db_v2") / f"{db_stem}.json"
            print(f"  Using DB: {db_path.name}")
            
            # Run recommendation
            engine = RuleEngine.from_file(db_path)
            rec = engine.recommend(smiles)
            
            print(f"  [SUCCESS] Rule: {rec.base_rule.name}")
            success_count += 1
            
        except Exception as e:
            print(f"  [ERROR] {e}")
        print("-" * 30)
        
    print(f"\nSummary: {success_count}/{total} successful recommendations.")

if __name__ == "__main__":
    test_families()
