#!/usr/bin/env python3
"""
Test the calculable features system using the sample compounds library.

This script validates that the feature detection system correctly identifies
molecular features across a diverse set of real-world compounds.
"""

import sys
from pathlib import Path

# Add parent directory to path
repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(repo_root))
sys.path.insert(0, str(repo_root / "tests"))

from sample_compounds import (
    ALL_SAMPLE_COMPOUNDS,
    get_compounds_by_role,
    get_compounds_by_reaction,
    get_compounds_by_feature
)
from chemtools.featurizers.calculable import detect_all_features, detect_feature
from chemtools.util.rdkit_helpers import rdkit_available

def test_feature_detection_accuracy():
    """Test that expected features are correctly detected for each compound."""
    print("=" * 80)
    print("Testing Feature Detection Accuracy")
    print("=" * 80)
    
    total_compounds = len(ALL_SAMPLE_COMPOUNDS)
    total_features_expected = sum(len(c["features"]) for c in ALL_SAMPLE_COMPOUNDS)
    correct_detections = 0
    incorrect_detections = 0
    
    results = {
        "correct": [],
        "missing": [],
        "unexpected": []
    }
    
    for compound in ALL_SAMPLE_COMPOUNDS:
        smiles = compound["smiles"]
        name = compound["name"]
        expected_features = set(compound["features"])
        
        # Detect all features
        detected = detect_all_features(smiles)
        
        if detected is None:
            print(f"❌ Failed to parse: {name} ({smiles})")
            incorrect_detections += len(expected_features)
            continue
        
        # Get detected features (those that are True)
        detected_features = {k for k, v in detected.items() if v is True}
        
        # Check matches
        correct = expected_features & detected_features
        missing = expected_features - detected_features
        unexpected = detected_features - expected_features
        
        correct_detections += len(correct)
        incorrect_detections += len(missing)
        
        if missing or unexpected:
            results["missing"].extend([
                {"compound": name, "smiles": smiles, "feature": f}
                for f in missing
            ])
            results["unexpected"].extend([
                {"compound": name, "smiles": smiles, "feature": f}
                for f in unexpected
            ])
        else:
            results["correct"].append({
                "compound": name,
                "smiles": smiles,
                "features": list(correct)
            })
    
    # Print summary
    print(f"\n📊 Summary:")
    print(f"  Total compounds tested: {total_compounds}")
    print(f"  Total expected features: {total_features_expected}")
    print(f"  Correctly detected: {correct_detections}")
    print(f"  Missed features: {incorrect_detections}")
    print(f"  Accuracy: {correct_detections / total_features_expected * 100:.1f}%")
    
    # Show some examples
    if results["correct"]:
        print(f"\n✅ Perfect matches: {len(results['correct'])}")
        for example in results["correct"][:3]:
            print(f"  - {example['compound']}: {', '.join(example['features'])}")
    
    if results["missing"]:
        print(f"\n⚠️  Missing features: {len(results['missing'])}")
        for example in results["missing"][:5]:
            print(f"  - {example['compound']}: expected '{example['feature']}'")
    
    if results["unexpected"]:
        print(f"\n🔍 Unexpected features: {len(results['unexpected'])}")
        for example in results["unexpected"][:5]:
            print(f"  - {example['compound']}: detected '{example['feature']}'")
    
    return results


def test_reaction_specific_features():
    """Test feature detection for specific reaction types."""
    print("\n" + "=" * 80)
    print("Testing Reaction-Specific Feature Detection")
    print("=" * 80)
    
    # Test Suzuki-Miyaura coupling compounds
    print("\n🔬 Suzuki-Miyaura Coupling Compounds:")
    suzuki_compounds = get_compounds_by_reaction("Suzuki-Miyaura")
    print(f"  Total: {len(suzuki_compounds)}")
    
    # Count key features
    ar_x_count = 0
    boronic_count = 0
    for compound in suzuki_compounds:
        features = detect_all_features(compound["smiles"])
        if features:
            if features.get("ArX"):
                ar_x_count += 1
            if features.get("boronic_acid_present") or features.get("boronic_ester_present"):
                boronic_count += 1
    
    print(f"  Aryl halides (ArX): {ar_x_count}")
    print(f"  Boronic acids/esters: {boronic_count}")
    
    # Test Buchwald-Hartwig coupling compounds
    print("\n🔬 Buchwald-Hartwig Coupling Compounds:")
    bh_compounds = get_compounds_by_reaction("Buchwald-Hartwig")
    print(f"  Total: {len(bh_compounds)}")
    
    # Count key features
    ar_x_count = 0
    amine_count = 0
    for compound in bh_compounds:
        features = detect_all_features(compound["smiles"])
        if features:
            if features.get("ArX"):
                ar_x_count += 1
            if features.get("primary_amine_present") or features.get("secondary_amine_present"):
                amine_count += 1
    
    print(f"  Aryl halides (ArX): {ar_x_count}")
    print(f"  Amines: {amine_count}")


def test_electrophile_features():
    """Test feature detection for electrophiles."""
    print("\n" + "=" * 80)
    print("Testing Electrophile Feature Detection")
    print("=" * 80)
    
    electrophiles = get_compounds_by_role("electrophile")
    print(f"\nTotal electrophiles: {len(electrophiles)}")
    
    # Count different electrophile types
    feature_counts = {
        "ArCl": 0,
        "ArBr": 0,
        "ArI": 0,
        "triflate_present": 0,
        "tosylate_present": 0,
        "vinyl_halide_present": 0,
        "alkyl_halide_present": 0
    }
    
    for compound in electrophiles:
        features = detect_all_features(compound["smiles"])
        if features:
            for feature_name in feature_counts.keys():
                if features.get(feature_name):
                    feature_counts[feature_name] += 1
    
    print("\nElectrophile breakdown:")
    for feature, count in sorted(feature_counts.items(), key=lambda x: x[1], reverse=True):
        if count > 0:
            print(f"  {feature}: {count}")


def test_nucleophile_features():
    """Test feature detection for nucleophiles."""
    print("\n" + "=" * 80)
    print("Testing Nucleophile Feature Detection")
    print("=" * 80)
    
    nucleophiles = get_compounds_by_role("nucleophile")
    print(f"\nTotal nucleophiles: {len(nucleophiles)}")
    
    # Count different nucleophile types
    feature_counts = {
        "boronic_acid_present": 0,
        "boronic_ester_present": 0,
        "primary_amine_present": 0,
        "secondary_amine_present": 0,
        "alcohol_present": 0,
        "thiol_present": 0,
        "terminal_alkyne_present": 0,
        "grignard_reagent": 0,
        "organozinc_present": 0
    }
    
    for compound in nucleophiles:
        features = detect_all_features(compound["smiles"])
        if features:
            for feature_name in feature_counts.keys():
                if features.get(feature_name):
                    feature_counts[feature_name] += 1
    
    print("\nNucleophile breakdown:")
    for feature, count in sorted(feature_counts.items(), key=lambda x: x[1], reverse=True):
        if count > 0:
            print(f"  {feature}: {count}")


def test_polarity_features():
    """Test polarity feature detection across compounds."""
    print("\n" + "=" * 80)
    print("Testing Polarity Feature Detection")
    print("=" * 80)
    
    polarity_counts = {
        "low_polarity": 0,
        "moderate_polarity": 0,
        "high_polarity": 0
    }
    
    examples = {
        "low_polarity": [],
        "moderate_polarity": [],
        "high_polarity": []
    }
    
    for compound in ALL_SAMPLE_COMPOUNDS:
        features = detect_all_features(compound["smiles"])
        if features:
            for polarity in polarity_counts.keys():
                if features.get(polarity):
                    polarity_counts[polarity] += 1
                    if len(examples[polarity]) < 3:
                        examples[polarity].append(compound["name"])
    
    print("\nPolarity distribution:")
    for polarity, count in polarity_counts.items():
        print(f"  {polarity}: {count}")
        if examples[polarity]:
            print(f"    Examples: {', '.join(examples[polarity])}")


def test_beta_hydride_elimination():
    """Test β-hydride elimination risk detection."""
    print("\n" + "=" * 80)
    print("Testing β-Hydride Elimination Risk Detection")
    print("=" * 80)
    
    beta_hydride_compounds = []
    no_beta_hydride_compounds = []
    
    for compound in ALL_SAMPLE_COMPOUNDS:
        features = detect_all_features(compound["smiles"])
        if features:
            if features.get("beta_hydride"):
                beta_hydride_compounds.append(compound["name"])
            else:
                no_beta_hydride_compounds.append(compound["name"])
    
    print(f"\nCompounds with β-hydride risk: {len(beta_hydride_compounds)}")
    if beta_hydride_compounds[:5]:
        print(f"  Examples: {', '.join(beta_hydride_compounds[:5])}")
    
    print(f"\nCompounds without β-hydride risk: {len(no_beta_hydride_compounds)}")
    if no_beta_hydride_compounds[:5]:
        print(f"  Examples: {', '.join(no_beta_hydride_compounds[:5])}")


def test_specific_compound_details():
    """Show detailed feature detection for specific interesting compounds."""
    print("\n" + "=" * 80)
    print("Detailed Feature Detection Examples")
    print("=" * 80)
    
    # Select some interesting compounds to showcase
    interesting_names = [
        "4-Bromoanisole",
        "2-Bromopyridine",
        "Phenylboronic acid",
        "tert-Butyl carbamate",
        "4-Chlorobenzyl alcohol"
    ]
    
    for name in interesting_names:
        compound = next((c for c in ALL_SAMPLE_COMPOUNDS if c["name"] == name), None)
        if compound:
            print(f"\n🔍 {name}")
            print(f"   SMILES: {compound['smiles']}")
            print(f"   Role: {compound['role']}")
            print(f"   Reactions: {', '.join(compound['reaction_types'])}")
            
            features = detect_all_features(compound["smiles"])
            if features:
                detected = [k for k, v in features.items() if v is True]
                print(f"   Detected features ({len(detected)}):")
                for feature in sorted(detected):
                    expected = "✓" if feature in compound["features"] else "⚠️"
                    print(f"     {expected} {feature}")
            
            if compound["notes"]:
                print(f"   Notes: {compound['notes']}")


def main():
    """Run all tests."""
    if not rdkit_available():
        print("❌ RDKit is not available. Cannot run tests.")
        return 1
    
    print("\n" + "=" * 80)
    print("🧪 CALCULABLE FEATURES VALIDATION TEST SUITE")
    print("   Using Sample Compounds Library")
    print("=" * 80)
    
    # Run all test functions
    results = test_feature_detection_accuracy()
    test_reaction_specific_features()
    test_electrophile_features()
    test_nucleophile_features()
    test_polarity_features()
    test_beta_hydride_elimination()
    test_specific_compound_details()
    
    print("\n" + "=" * 80)
    print("✅ Test suite completed!")
    print("=" * 80)
    
    # Return exit code based on accuracy
    accuracy = len(results["correct"]) / len(ALL_SAMPLE_COMPOUNDS) * 100
    if accuracy >= 80:
        print(f"\n✅ PASS: {accuracy:.1f}% of compounds had perfect feature detection")
        return 0
    else:
        print(f"\n⚠️  WARN: Only {accuracy:.1f}% of compounds had perfect feature detection")
        return 1


if __name__ == "__main__":
    sys.exit(main())
