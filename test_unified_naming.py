"""
Quick test to verify unified naming convention.

Tests that:
1. Rule DB files have been renamed correctly
2. Dataset files match rule DB naming pattern
3. Family name mapping works correctly
"""

import os
from pathlib import Path

# Project root
ROOT = Path(__file__).parent

def test_file_naming_alignment():
    """Test that rule DB and dataset files align in naming."""
    
    # Check conditionDB files
    condition_db_dir = ROOT / "data" / "conditionDB"
    expected_rule_dbs = {
        "C_N_Coupling_Cu_db.json",
        "C_N_Coupling_Pd_db.json",
        "C_N_Coupling_Ni_db.json",
        "Suzuki_db.json",
        "Amide_formation_db.json",
    }
    
    actual_rule_dbs = set(f.name for f in condition_db_dir.glob("*.json"))
    
    print("=" * 70)
    print("Rule DB Files (data/conditionDB/)")
    print("=" * 70)
    for db in sorted(actual_rule_dbs):
        status = "✅" if db in expected_rule_dbs else "❌"
        print(f"  {status} {db}")
    
    # Check reaction dataset files
    dataset_dir = ROOT / "data" / "reaction_dataset"
    expected_datasets = {
        "C_N_Coupling_Cu.jsonl",
        "C_N_Coupling_Pd.jsonl",
        "C_N_Coupling_Ni.jsonl",
        "Suzuki.jsonl",
        "Amide_formation.jsonl",
    }
    
    actual_datasets = set(f.name for f in dataset_dir.glob("*.jsonl"))
    
    print("\n" + "=" * 70)
    print("Dataset Files (data/reaction_dataset/)")
    print("=" * 70)
    for ds in sorted(actual_datasets):
        status = "✅" if ds in expected_datasets else "⚠️ "
        print(f"  {status} {ds}")
    
    # Check alignment
    print("\n" + "=" * 70)
    print("Naming Alignment Check")
    print("=" * 70)
    
    # Extract base names
    rule_bases = {f.replace("_db.json", "") for f in expected_rule_dbs}
    dataset_bases = {f.replace(".jsonl", "") for f in expected_datasets}
    
    for rule_base in sorted(rule_bases):
        dataset_base = rule_base
        has_dataset = dataset_base in dataset_bases
        status = "✅" if has_dataset else "❌"
        print(f"  {status} {rule_base:25s} -> {dataset_base}.jsonl {'(found)' if has_dataset else '(MISSING)'}")
    
    # Test missing files
    missing_rule = expected_rule_dbs - actual_rule_dbs
    missing_dataset = expected_datasets - actual_datasets
    
    all_good = len(missing_rule) == 0 and len(missing_dataset) == 0
    
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    if missing_rule:
        print(f"  ❌ Missing rule DBs: {', '.join(sorted(missing_rule))}")
    else:
        print(f"  ✅ All expected rule DBs present ({len(expected_rule_dbs)} files)")
    
    if missing_dataset:
        print(f"  ⚠️  Missing datasets: {', '.join(sorted(missing_dataset))}")
    else:
        print(f"  ✅ All expected datasets present ({len(expected_datasets)} files)")
    
    if all_good:
        print("\n  🎉 All files properly aligned!")
    
    return all_good


def test_family_mapping():
    """Test that family name mapping works correctly."""
    print("\n" + "=" * 70)
    print("Family Name Mapping Test")
    print("=" * 70)
    
    # Import after path setup
    import sys
    sys.path.insert(0, str(ROOT))
    
    from chemtools.recommend.utils import canonical_family
    
    test_cases = [
        # Detection name -> Expected dataset name
        ("Ullmann_CN", "C_N_Coupling_Cu"),
        ("Buchwald_CN", "C_N_Coupling_Pd"),
        ("C_N_Coupling_Cu", "C_N_Coupling_Cu"),
        ("C_N_Coupling_Pd", "C_N_Coupling_Pd"),
        ("C_N_Coupling_Ni", "C_N_Coupling_Ni"),
        ("Suzuki_CC", "Suzuki"),
        ("Suzuki", "Suzuki"),
        ("Amide_Coupling", "Amide_formation"),
        ("Amide_formation", "Amide_formation"),
    ]
    
    all_pass = True
    for input_name, expected_output in test_cases:
        actual_output = canonical_family(input_name)
        passed = actual_output == expected_output
        all_pass = all_pass and passed
        status = "✅" if passed else "❌"
        print(f"  {status} {input_name:20s} -> {actual_output:20s} (expected: {expected_output})")
    
    print("\n" + "=" * 70)
    if all_pass:
        print("  🎉 All family mappings correct!")
    else:
        print("  ❌ Some family mappings failed")
    
    return all_pass


if __name__ == "__main__":
    files_ok = test_file_naming_alignment()
    mapping_ok = test_family_mapping()
    
    print("\n" + "=" * 70)
    print("FINAL RESULT")
    print("=" * 70)
    if files_ok and mapping_ok:
        print("  ✅ All tests passed! Unified naming is working correctly.")
        exit(0)
    else:
        print("  ❌ Some tests failed. Please review the output above.")
        exit(1)
