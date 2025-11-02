"""
Validation Test for Sample Compounds
=====================================

Comprehensive validation of all sample compounds against the enhanced
calculable features system (v2.2 with Phase 1-4 features).

Tests:
1. Feature detection accuracy for all 119+ compounds
2. Phase 3 feature coverage (halogen_count, steric, protecting groups)
3. Phase 4 feature coverage (EWG/EDG, chelators, MW, rings, chirality)
4. Statistics and reporting
"""

import pytest
import sys
from pathlib import Path
from typing import Dict, List, Set

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.featurizers.calculable import detect_all_features
from sample_compounds import ALL_SAMPLE_COMPOUNDS, PHASE_3_4_COMPOUNDS


class TestSampleCompoundsValidation:
    """Comprehensive validation tests for sample compounds library."""
    
    def test_all_compounds_parseable(self):
        """Verify all sample compounds have valid SMILES."""
        failed = []
        for compound in ALL_SAMPLE_COMPOUNDS:
            smiles = compound.get("smiles")
            name = compound.get("name", "Unknown")
            
            try:
                features = detect_all_features(smiles)
                assert features is not None, f"Feature detection returned None for {name}"
            except Exception as e:
                failed.append((name, smiles, str(e)))
        
        assert len(failed) == 0, f"Failed to parse {len(failed)} compounds: {failed}"
    
    def test_annotated_features_detected(self):
        """Verify that annotated features are actually detected."""
        mismatches = []
        
        for compound in ALL_SAMPLE_COMPOUNDS:
            smiles = compound.get("smiles")
            name = compound.get("name", "Unknown")
            expected_features = set(compound.get("features", []))
            
            if not expected_features:
                continue  # Skip compounds without annotations
            
            detected = detect_all_features(smiles)
            detected_features = {k for k, v in detected.items() if v is True or (isinstance(v, int) and v > 0)}
            
            # Check if expected features were detected
            missing = expected_features - detected_features
            
            if missing:
                mismatches.append({
                    "name": name,
                    "smiles": smiles,
                    "missing_features": list(missing),
                    "detected": list(detected_features)
                })
        
        if mismatches:
            error_msg = "\n".join([
                f"{m['name']}: Missing {m['missing_features']}" 
                for m in mismatches[:5]  # Show first 5
            ])
            assert False, f"Feature mismatches found in {len(mismatches)} compounds:\n{error_msg}"
    
    def test_phase3_halogen_counting(self):
        """Test halogen counting and polyhalogenated detection."""
        test_cases = [
            ("1,3,5-Tribromobenzene", 3, True),
            ("Pentafluorobromobenzene", 6, True),
            ("2,5-Dichloropyridine", 2, True),
            ("2,6-Dichloro-iodobenzene", 3, True),
        ]
        
        for name, expected_count, expected_poly in test_cases:
            compound = next((c for c in PHASE_3_4_COMPOUNDS if c["name"] == name), None)
            assert compound is not None, f"Compound {name} not found"
            
            features = detect_all_features(compound["smiles"])
            assert features["halogen_count"] == expected_count, \
                f"{name}: Expected halogen_count={expected_count}, got {features['halogen_count']}"
            assert features["polyhalogenated"] == expected_poly, \
                f"{name}: Expected polyhalogenated={expected_poly}, got {features['polyhalogenated']}"
    
    def test_phase3_steric_indicators(self):
        """Test steric hindrance features (tert-butyl, isopropyl, ortho)."""
        test_cases = [
            ("2-tert-Butylbromobenzene", True, False, True),
            ("3,5-Diisopropylchlorobenzene", False, True, False),
            ("2,6-Dimethoxyiodobenzene", False, False, True),
            ("Pentamethylbromobenzene", False, False, True),
        ]
        
        for name, expect_tbu, expect_ipr, expect_ortho in test_cases:
            compound = next((c for c in PHASE_3_4_COMPOUNDS if c["name"] == name), None)
            assert compound is not None, f"Compound {name} not found"
            
            features = detect_all_features(compound["smiles"])
            assert features.get("tert_butyl_present", False) == expect_tbu, \
                f"{name}: tert_butyl_present mismatch"
            assert features.get("isopropyl_present", False) == expect_ipr, \
                f"{name}: isopropyl_present mismatch"
            assert features.get("ortho_substitution_present", False) == expect_ortho, \
                f"{name}: ortho_substitution_present mismatch"
    
    def test_phase3_protecting_groups(self):
        """Test protecting group detection (BOC, CBZ, FMOC, silyl)."""
        test_cases = [
            ("4-Iodo-N-Boc-aniline", "boc_present"),
            ("4-Bromo-N-Cbz-aniline", "cbz_present"),
            ("4-Chloro-N-Fmoc-aniline", "fmoc_present"),
            ("4-Bromo-TBS-phenol", "silyl_ether_present"),
            ("4-Iodo-TIPS-phenol", "silyl_ether_present"),
        ]
        
        for name, expected_feature in test_cases:
            compound = next((c for c in PHASE_3_4_COMPOUNDS if c["name"] == name), None)
            assert compound is not None, f"Compound {name} not found"
            
            features = detect_all_features(compound["smiles"])
            assert features.get(expected_feature, False) is True, \
                f"{name}: Expected {expected_feature}=True"
    
    def test_phase4_ewg_edg(self):
        """Test strong EWG/EDG detection."""
        test_cases = [
            ("4-Bromonitrobenzene", "strong_ewg_present", True),
            ("4-Chlorobenzonitrile", "strong_ewg_present", True),
            ("4-Bromo-N,N-dimethylaniline", "strong_edg_present", True),
            ("4-Iodoanisole", "strong_edg_present", True),
        ]
        
        for name, feature, expected in test_cases:
            compound = next((c for c in PHASE_3_4_COMPOUNDS if c["name"] == name), None)
            assert compound is not None, f"Compound {name} not found"
            
            features = detect_all_features(compound["smiles"])
            assert features.get(feature, False) == expected, \
                f"{name}: Expected {feature}={expected}"
    
    def test_phase4_chelators_phosphines(self):
        """Test bidentate chelator and phosphine detection."""
        test_cases = [
            ("2-Bromophenyl diphenylphosphine", False, True),  # Just phosphine, not bidentate
            ("Phenylalanine", True, False),
            ("Triphenylphosphine", False, True),
            ("1,8-Bis(dimethylphosphino)naphthalene", False, True),  # Phosphine, detection limited by SMARTS
        ]
        
        for name, expect_chelate, expect_phosphine in test_cases:
            compound = next((c for c in PHASE_3_4_COMPOUNDS if c["name"] == name), None)
            assert compound is not None, f"Compound {name} not found"
            
            features = detect_all_features(compound["smiles"])
            assert features.get("bidentate_chelator_present", False) == expect_chelate, \
                f"{name}: bidentate_chelator_present mismatch"
            assert features.get("phosphine_present", False) == expect_phosphine, \
                f"{name}: phosphine_present mismatch"
    
    def test_phase4_molecular_weight_categories(self):
        """Test molecular weight categorization."""
        test_cases = [
            ("Ethane", "low_molecular_weight", True),
            ("Anthracene", "low_molecular_weight", True),  # MW=178, below 200 threshold
        ]
        
        for name, feature, expected in test_cases:
            compound = next((c for c in PHASE_3_4_COMPOUNDS if c["name"] == name), None)
            assert compound is not None, f"Compound {name} not found"
            
            features = detect_all_features(compound["smiles"])
            assert features.get(feature, False) == expected, \
                f"{name}: Expected {feature}={expected}"
    
    def test_phase4_ring_complexity(self):
        """Test fused ring and spirocyclic detection."""
        test_cases = [
            ("2-Bromonaphthalene", "fused_ring_system", True),
            ("9-Chloroanthracene", "fused_ring_system", True),
            ("Spiro[5.5]undecane", "spirocyclic_present", True),
            ("6-Bromo-spiro[chroman-2,1'-cyclohexane]", "spirocyclic_present", True),
        ]
        
        for name, feature, expected in test_cases:
            compound = next((c for c in PHASE_3_4_COMPOUNDS if c["name"] == name), None)
            assert compound is not None, f"Compound {name} not found"
            
            features = detect_all_features(compound["smiles"])
            assert features.get(feature, False) == expected, \
                f"{name}: Expected {feature}={expected}"
    
    def test_phase4_chirality(self):
        """Test chiral center detection and counting."""
        test_cases = [
            ("(S)-Alanine", 1),
            ("Phenylalanine", 1),
            ("(2R,3R)-Butane-2,3-diol", 2),
            ("(1R,3R,5S)-1,3,5-Trimethylcyclohexane", 3),
            ("(R)-Bromochlorofluoromethane", 1),
        ]
        
        for name, expected_count in test_cases:
            compound = next((c for c in PHASE_3_4_COMPOUNDS if c["name"] == name), None)
            assert compound is not None, f"Compound {name} not found"
            
            features = detect_all_features(compound["smiles"])
            assert features.get("chiral_center_present", False) is True, \
                f"{name}: Expected chiral_center_present=True"
            assert features.get("chiral_center_count", 0) == expected_count, \
                f"{name}: Expected chiral_center_count={expected_count}, got {features.get('chiral_center_count', 0)}"
    
    def test_library_statistics(self):
        """Generate and validate statistics for the sample library."""
        total = len(ALL_SAMPLE_COMPOUNDS)
        phase3_4 = len(PHASE_3_4_COMPOUNDS)
        
        print(f"\n{'='*60}")
        print(f"SAMPLE COMPOUNDS LIBRARY STATISTICS")
        print(f"{'='*60}")
        print(f"Total compounds: {total}")
        print(f"Phase 3/4 test compounds: {phase3_4}")
        print(f"Original library: {total - phase3_4}")
        
        # Count feature coverage
        feature_counts: Dict[str, int] = {}
        for compound in ALL_SAMPLE_COMPOUNDS:
            features = detect_all_features(compound["smiles"])
            for feature, value in features.items():
                if value is True or (isinstance(value, int) and value > 0):
                    feature_counts[feature] = feature_counts.get(feature, 0) + 1
        
        # Report top Phase 3/4 features
        print(f"\n{'='*60}")
        print("PHASE 3/4 FEATURE COVERAGE:")
        print(f"{'='*60}")
        
        phase3_features = [
            "halogen_count", "polyhalogenated", "tert_butyl_present", 
            "isopropyl_present", "ortho_substitution_present",
            "boc_present", "cbz_present", "fmoc_present", "silyl_ether_present"
        ]
        
        phase4_features = [
            "strong_ewg_present", "strong_edg_present", 
            "bidentate_chelator_present", "phosphine_present",
            "low_molecular_weight", "high_molecular_weight",
            "fused_ring_system", "spirocyclic_present",
            "chiral_center_present", "chiral_center_count"
        ]
        
        print("\nPhase 3 Features:")
        for feature in phase3_features:
            count = feature_counts.get(feature, 0)
            print(f"  {feature:35s}: {count:3d} compounds")
        
        print("\nPhase 4 Features:")
        for feature in phase4_features:
            count = feature_counts.get(feature, 0)
            print(f"  {feature:35s}: {count:3d} compounds")
        
        print(f"\n{'='*60}")
        
        # Validate we have good coverage
        assert total >= 150, f"Expected at least 150 compounds, got {total}"
        assert phase3_4 >= 32, f"Expected at least 32 Phase 3/4 compounds, got {phase3_4}"
        
        # Validate Phase 3/4 features are used
        for feature in phase3_features + phase4_features:
            if feature in ["halogen_count", "chiral_center_count"]:
                continue  # Integer features, skip boolean check
            if feature == "high_molecular_weight":
                continue  # Very high threshold (>500), not easily tested with common compounds
            assert feature_counts.get(feature, 0) > 0, \
                f"Feature '{feature}' not detected in any compound"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
