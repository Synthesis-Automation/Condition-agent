#!/usr/bin/env python3
"""
Simple validation of calculable features using sample compounds.
Shows what features are being detected correctly.
"""

import sys
from pathlib import Path

# Add parent directory to path
repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))
sys.path.insert(0, str(repo_root / "tests"))

from sample_compounds import ALL_SAMPLE_COMPOUNDS
from chemtools.featurizers.calculable import detect_all_features
from chemtools.util.rdkit_helpers import rdkit_available

def main():
    """Run simple validation."""
    if not rdkit_available():
        print("❌ RDKit is not available. Cannot run tests.")
        return 1
    
    print("=" * 80)
    print("🧪 Calculable Features Detection Summary")
    print("=" * 80)
    
    print(f"\nTesting {len(ALL_SAMPLE_COMPOUNDS)} compounds...")
    
    # Track all detected features
    all_features_detected = set()
    feature_counts = {}
    successful_parses = 0
    failed_parses = 0
    
    for compound in ALL_SAMPLE_COMPOUNDS:
        smiles = compound["smiles"]
        name = compound["name"]
        
        features = detect_all_features(smiles)
        
        if features is None:
            failed_parses += 1
            print(f"  ❌ Failed to parse: {name}")
            continue
        
        successful_parses += 1
        
        # Count detected features
        for feature_name, feature_value in features.items():
            if feature_value is True:
                all_features_detected.add(feature_name)
                feature_counts[feature_name] = feature_counts.get(feature_name, 0) + 1
    
    print(f"\n📊 Results:")
    print(f"  Successfully parsed: {successful_parses}/{len(ALL_SAMPLE_COMPOUNDS)}")
    print(f"  Failed to parse: {failed_parses}")
    print(f"  Unique features detected: {len(all_features_detected)}")
    
    print(f"\n🔍 Top 20 Most Common Features:")
    sorted_features = sorted(feature_counts.items(), key=lambda x: x[1], reverse=True)
    for i, (feature, count) in enumerate(sorted_features[:20], 1):
        pct = count / successful_parses * 100
        print(f"  {i:2d}. {feature:35s} {count:3d} ({pct:5.1f}%)")
    
    # Show some specific examples
    print(f"\n📝 Example Detections:")
    
    examples = [
        ("Bromobenzene", ["ArBr_present", "sp2_bromide_present", "aryl_halide_present"]),
        ("Phenylboronic acid", ["boronic_acid_present", "sp2_boron_present"]),
        ("Aniline", ["primary_amine_present"]),
        ("4-Bromopyridine", ["heteroaryl_halide_present", "pyridine_present"])
    ]
    
    for name, expected in examples:
        compound = next((c for c in ALL_SAMPLE_COMPOUNDS if c["name"] == name), None)
        if compound:
            features = detect_all_features(compound["smiles"])
            if features:
                detected = [k for k, v in features.items() if v is True]
                print(f"\n  {name}:")
                print(f"    SMILES: {compound['smiles']}")
                print(f"    Expected: {', '.join(expected)}")
                print(f"    Detected ({len(detected)}): {', '.join(sorted(detected))}")
                
                # Check matches
                found = [f for f in expected if features.get(f)]
                missing = [f for f in expected if not features.get(f)]
                if found:
                    print(f"    ✓ Found: {', '.join(found)}")
                if missing:
                    print(f"    ✗ Missing: {', '.join(missing)}")
    
    print("\n" + "=" * 80)
    print("✅ Validation complete!")
    print("=" * 80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
