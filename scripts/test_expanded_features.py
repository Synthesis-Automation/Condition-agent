"""
Test the expanded calculable_features.json with various molecules
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.featurizers import calculable

# Test molecules covering various new features
test_molecules = [
    # Protecting groups
    ("CC(=O)Oc1ccccc1", "Acetate ester (phenyl acetate)"),
    ("CC(C)(C)OC(=O)Nc1ccccc1", "Boc-protected aniline"),
    ("c1ccc(COc2ccccc2)cc1", "Benzyl ether"),
    ("O[Si](C)(C)C(C)(C)C", "TBS ether"),
    
    # Heterocycles
    ("c1ccoc1", "Furan"),
    ("c1ccsc1", "Thiophene"),
    ("c1c[nH]cc1", "Pyrrole"),
    ("c1cnc[nH]1", "Imidazole"),
    ("c1ccc2ncccc2c1", "Quinoline"),
    ("C1COCCN1", "Morpholine"),
    ("C1CNCNC1", "Piperazine"),
    
    # Reactive intermediates
    ("C1CO1", "Epoxide (ethylene oxide)"),
    ("CC(=N=[N+]=[N-])C", "Diazo compound"),
    ("CCN=[N+]=[N-]", "Azide"),
    ("C=CC(=O)C", "α,β-Unsaturated ketone (MVK)"),
    ("BrCC(=O)C", "α-Bromo ketone"),
    
    # Medicinal chemistry
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Ibuprofen"),
    ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
    ("FC(F)(F)c1ccccc1", "Trifluorotoluene"),
    ("O=S(=O)(N)c1ccccc1", "Benzenesulfonamide"),
    
    # Organometallics
    ("c1ccc(B(O)O)cc1", "Phenylboronic acid"),
    ("c1ccc(B2OC(C)(C)C(C)(C)O2)cc1", "Phenyl Bpin"),
    ("Brc1ccccc1", "Bromobenzene (electrophile)"),
    
    # Drug-like molecules
    ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Caffeine"),
    ("CC(C)NCC(COc1ccccc1)O", "Propranolol (simplified)"),
]

print("=" * 80)
print("EXPANDED CALCULABLE FEATURES TEST")
print("=" * 80)

for smiles, name in test_molecules:
    print(f"\n{name}")
    print(f"SMILES: {smiles}")
    print("-" * 80)
    
    try:
        features = calculable.detect_all_features(smiles)
        
        # Filter to only True boolean values and non-zero counts
        active_features = {}
        for k, v in features.items():
            if isinstance(v, bool) and v:
                active_features[k] = v
            elif isinstance(v, (int, float)) and v != 0:
                active_features[k] = v
        
        if active_features:
            print(f"Active features ({len(active_features)}):")
            for feat, val in sorted(active_features.items())[:15]:  # Show first 15
                print(f"  ✓ {feat}: {val}")
            
            if len(active_features) > 15:
                print(f"  ... and {len(active_features) - 15} more")
        else:
            print("  No features detected (may need RDKit)")
            
    except Exception as e:
        print(f"  ERROR: {e}")

print("\n" + "=" * 80)
print("TEST COMPLETE")
print("=" * 80)

# Test derived features specifically
print("\n" + "=" * 80)
print("DERIVED FEATURE TESTS")
print("=" * 80)

derived_tests = [
    ("CC(C)(C)OC(=O)Nc1ccccc1", "protected_amine_present", "Boc-aniline"),
    ("c1ccc(COc2ccccc2)cc1", "protected_alcohol_present", "Benzyl ether"),
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "drug_like", "Ibuprofen"),
    ("c1ccc(B(O)O)cc1", "cross_coupling_nucleophile", "Phenylboronic acid"),
    ("Brc1ccccc1", "cross_coupling_electrophile", "Bromobenzene"),
    ("C=CC(=O)C", "michael_acceptor_present", "MVK"),
    ("CCN=[N+]=[N-]", "explosive_risk", "Ethyl azide"),
    ("FC(F)(F)c1ccccc1", "fluorinated_motif_present", "Trifluorotoluene"),
]

for smiles, feature, name in derived_tests:
    try:
        features = calculable.detect_all_features(smiles)
        result = features.get(feature, False)
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status} | {name:25s} | {feature:30s} = {result}")
    except Exception as e:
        print(f"✗ ERROR | {name:25s} | {feature:30s} = {e}")

print("\n" + "=" * 80)
