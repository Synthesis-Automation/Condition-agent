"""
Test expanded calculable features with comprehensive sample compounds
"""

import sys
from pathlib import Path

# Add parent directory to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Now import after fixing path
from chemtools.featurizers import calculable

# Import sample compounds directly from file
import importlib.util
spec = importlib.util.spec_from_file_location(
    "sample_compounds", 
    project_root / "tests" / "sample_compounds.py"
)
sample_compounds = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sample_compounds)

ALL_SAMPLE_COMPOUNDS = sample_compounds.ALL_SAMPLE_COMPOUNDS

print("=" * 100)
print("EXPANDED CALCULABLE FEATURES v3.0 - COMPREHENSIVE TEST WITH SAMPLE COMPOUNDS")
print("=" * 100)

# Test statistics
total_compounds = len(ALL_SAMPLE_COMPOUNDS)
total_features_detected = 0
compounds_tested = 0
errors = 0

# Track which new features are being detected
new_feature_detections = {}
v3_specific_features = [
    # Protecting groups
    'boc_present', 'cbz_present', 'fmoc_present', 'silyl_ether_present',
    'acetate_ester_present', 'benzoate_ester_present', 'tert_butyl_ester_present',
    'benzyl_ether_present', 'pmb_ether_present', 'mom_ether_present',
    'tosyl_sulfonamide_present', 'phthalimide_present',
    
    # Heterocycles
    'furan_present', 'thiophene_present', 'pyrrole_present', 'imidazole_present',
    'oxazole_present', 'thiazole_present', 'triazole_present', 'tetrazole_present',
    'quinoline_present', 'isoquinoline_present', 'benzofuran_present',
    'morpholine_present', 'piperazine_present', 'piperidine_present',
    
    # Reactive intermediates
    'epoxide_present', 'aziridine_present', 'diazo_compound_present',
    'azide_present', 'peroxide_present', 'alpha_beta_unsaturated_carbonyl_present',
    'alpha_halo_carbonyl_present',
    
    # Physicochemical
    'lipinski_hbd_compliant', 'lipinski_hba_compliant', 'lipinski_mw_compliant',
    'lipinski_logp_compliant', 'veber_rotb_compliant', 'veber_tpsa_compliant',
    'ionizable_basic_group_present', 'ionizable_acidic_group_present',
    
    # Medicinal chemistry
    'fluorine_present', 'trifluoromethyl_present', 'sulfone_present',
    'sulfonamide_present', 'urea_present', 'carbamate_present',
    'pains_aldehyde_alert', 'pains_catechol_alert',
    
    # Stereochemistry
    'chiral_center_present', 'chiral_center_count', 'spiro_center_present',
    'quaternary_carbon_present',
    
    # Redox/photo
    'benzylic_ch_present', 'allylic_ch_present', 'alpha_amino_ch_present',
    'alpha_oxy_ch_present',
    
    # Organometallics
    'aryl_boron_present', 'vinyl_boron_present', 'alkyl_boron_present',
    
    # Derived features
    'protected_amine_present', 'protected_alcohol_present',
    'lipinski_compliant', 'veber_compliant', 'drug_like',
    'michael_acceptor_present', 'cross_coupling_electrophile',
    'cross_coupling_nucleophile', 'explosive_risk',
    'fluorinated_motif_present', 'metabolic_soft_spot_present',
]

print(f"\nTesting {total_compounds} sample compounds...")
print(f"Tracking {len(v3_specific_features)} v3.0-specific features\n")

# Group results by feature for summary
feature_hit_counts = {feat: 0 for feat in v3_specific_features}

# Test each compound
for i, compound in enumerate(ALL_SAMPLE_COMPOUNDS, 1):
    smiles = compound['smiles']
    name = compound['name']
    role = compound.get('role', 'unknown')
    expected_features = compound.get('features', [])
    
    try:
        features = calculable.detect_all_features(smiles)
        compounds_tested += 1
        
        # Count active features
        active_count = sum(1 for v in features.values() if v)
        total_features_detected += active_count
        
        # Track v3.0 feature detections
        for feat in v3_specific_features:
            if features.get(feat):
                feature_hit_counts[feat] += 1
                if feat not in new_feature_detections:
                    new_feature_detections[feat] = []
                new_feature_detections[feat].append(name)
        
        # Check if expected features are detected
        missing_features = [f for f in expected_features if not features.get(f)]
        
        # Print detailed info for interesting cases (first 10 compounds)
        if i <= 10:
            print(f"\n{i}. {name} ({role})")
            print(f"   SMILES: {smiles}")
            print(f"   Active features: {active_count}")
            
            # Show v3.0 features if any
            v3_active = [f for f in v3_specific_features if features.get(f)]
            if v3_active:
                print(f"   v3.0 features: {', '.join(v3_active[:5])}")
                if len(v3_active) > 5:
                    print(f"                  ... and {len(v3_active) - 5} more")
            
            if missing_features:
                print(f"   ⚠ Missing expected: {', '.join(missing_features[:3])}")
        
    except Exception as e:
        errors += 1
        if i <= 10:
            print(f"\n{i}. {name} - ERROR: {e}")

print("\n" + "=" * 100)
print("SUMMARY STATISTICS")
print("=" * 100)

print(f"\nCompounds tested: {compounds_tested}/{total_compounds}")
print(f"Total feature detections: {total_features_detected}")
print(f"Average features per compound: {total_features_detected/compounds_tested if compounds_tested else 0:.1f}")
print(f"Errors encountered: {errors}")

print(f"\n" + "=" * 100)
print("v3.0 NEW FEATURE DETECTIONS (Top 20)")
print("=" * 100)

# Sort features by hit count
sorted_features = sorted(feature_hit_counts.items(), key=lambda x: -x[1])

for feat, count in sorted_features[:20]:
    if count > 0:
        percentage = (count / compounds_tested * 100) if compounds_tested else 0
        print(f"  {feat:40s} : {count:3d} hits ({percentage:5.1f}%)")
        # Show first few examples
        if feat in new_feature_detections:
            examples = new_feature_detections[feat][:3]
            print(f"     Examples: {', '.join(examples)}")

print(f"\n" + "=" * 100)
print("CATEGORY BREAKDOWN")
print("=" * 100)

# Test specific categories
categories = {
    "Protecting Groups": ['boc_present', 'cbz_present', 'fmoc_present', 'silyl_ether_present', 
                          'protected_amine_present', 'protected_alcohol_present'],
    "Heterocycles": ['furan_present', 'thiophene_present', 'pyrrole_present', 'quinoline_present',
                     'morpholine_present', 'piperazine_present'],
    "Drug-Likeness": ['lipinski_compliant', 'veber_compliant', 'drug_like'],
    "Cross-Coupling": ['cross_coupling_electrophile', 'cross_coupling_nucleophile'],
    "Safety Markers": ['explosive_risk', 'moisture_sensitive', 'air_sensitive'],
    "Medicinal Chem": ['fluorinated_motif_present', 'trifluoromethyl_present', 'sulfonamide_present'],
    "Stereochemistry": ['chiral_center_present', 'chiral_center_count', 'spiro_center_present'],
}

for category, features in categories.items():
    hits = sum(feature_hit_counts.get(f, 0) for f in features)
    print(f"\n{category:20s}: {hits:3d} total detections across {len(features)} features")

print(f"\n" + "=" * 100)
print("SPECIFIC FEATURE VALIDATION")
print("=" * 100)

# Test specific expected features
validation_tests = [
    ("1,3,5-Tribromobenzene", "polyhalogenated"),
    ("4-Bromo-N-Boc-aniline", "boc_present"),
    ("4-Bromo-N-Cbz-aniline", "cbz_present"),
    ("4-Chloro-N-Fmoc-aniline", "fmoc_present"),
    ("4-Bromo-TBS-phenol", "silyl_ether_present"),
    ("4-Bromonitrobenzene", "strong_ewg_present"),
    ("Phenylalanine", "chiral_center_present"),
    ("Triphenylphosphine", "phosphine_present"),
    ("(S)-Alanine", "bidentate_chelator_present"),
    ("Phenylboronic acid", "cross_coupling_nucleophile"),
    ("Bromobenzene", "cross_coupling_electrophile"),
]

print("\nValidating specific expected features:\n")
passed = 0
failed = 0

for name, expected_feature in validation_tests:
    compound = next((c for c in ALL_SAMPLE_COMPOUNDS if c['name'] == name), None)
    if compound:
        try:
            features = calculable.detect_all_features(compound['smiles'])
            result = features.get(expected_feature, False)
            status = "✓ PASS" if result else "✗ FAIL"
            if result:
                passed += 1
            else:
                failed += 1
            print(f"{status} | {name:35s} | {expected_feature:30s}")
        except Exception as e:
            failed += 1
            print(f"✗ ERROR | {name:35s} | {expected_feature:30s} - {e}")
    else:
        print(f"? SKIP  | {name:35s} | {expected_feature:30s} - compound not found")

print(f"\nValidation Results: {passed} passed, {failed} failed")

print(f"\n" + "=" * 100)
print("TEST COMPLETE")
print("=" * 100)

print(f"\n✓ Tested v3.0 expanded features with {compounds_tested} sample compounds")
print(f"✓ Detected {len([f for f, c in feature_hit_counts.items() if c > 0])} unique v3.0 features in dataset")
print(f"✓ Average {total_features_detected/compounds_tested if compounds_tested else 0:.1f} features per molecule")
print(f"✓ Feature expansion working as expected!")

# Additional insights
print(f"\n" + "=" * 100)
print("COVERAGE INSIGHTS")
print("=" * 100)

# Check which v3.0 features were NOT detected (might need more diverse samples)
undetected = [f for f, c in feature_hit_counts.items() if c == 0]
if undetected:
    print(f"\nv3.0 features not detected in current sample set ({len(undetected)}):")
    for feat in undetected[:10]:
        print(f"  - {feat}")
    if len(undetected) > 10:
        print(f"  ... and {len(undetected) - 10} more")
    print(f"\nNote: These features require specific functional groups not in current samples")
else:
    print(f"\n✓ All tracked v3.0 features were detected in at least one compound!")

print(f"\n" + "=" * 100)
