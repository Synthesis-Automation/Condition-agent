"""
Diagnostic script to check for issues in chemtools/analysis and /taxonomy
"""
import json
from pathlib import Path
from typing import List, Dict, Any

def check_taxonomy_data_integrity():
    """Check taxonomy JSON files for data integrity issues"""
    print("=" * 60)
    print("CHECKING TAXONOMY DATA INTEGRITY")
    print("=" * 60)
    
    data_dir = Path("chemtools/taxonomy/data")
    issues: List[str] = []
    
    # Load all JSON files
    try:
        with open(data_dir / "manifest.json", encoding="utf-8") as f:
            manifest = json.load(f)
        with open(data_dir / "reaction_categories.json", encoding="utf-8") as f:
            reaction_categories = json.load(f)
        with open(data_dir / "reaction_types.json", encoding="utf-8") as f:
            reaction_types = json.load(f)
        with open(data_dir / "reactant_types.json", encoding="utf-8") as f:
            reactant_types = json.load(f)
        with open(data_dir / "reagent_roles.json", encoding="utf-8") as f:
            reagent_roles = json.load(f)
        with open(data_dir / "reagent_families.json", encoding="utf-8") as f:
            reagent_families = json.load(f)
        with open(data_dir / "aliases.json", encoding="utf-8") as f:
            aliases = json.load(f)
    except Exception as e:
        print(f"❌ ERROR loading JSON files: {e}")
        import traceback
        traceback.print_exc()
        return []
    
    # Check 1: Duplicate IDs
    print("\n1. Checking for duplicate IDs...")
    
    reaction_cat_ids = [r["id"] for r in reaction_categories]
    if len(reaction_cat_ids) != len(set(reaction_cat_ids)):
        issues.append("Duplicate reaction category IDs found")
        
    reaction_type_ids = [r["id"] for r in reaction_types]
    if len(reaction_type_ids) != len(set(reaction_type_ids)):
        issues.append("Duplicate reaction type IDs found")
        
    reactant_type_ids = [r["id"] for r in reactant_types]
    if len(reactant_type_ids) != len(set(reactant_type_ids)):
        issues.append("Duplicate reactant type IDs found")
        
    reagent_role_ids = [r["id"] for r in reagent_roles]
    if len(reagent_role_ids) != len(set(reagent_role_ids)):
        issues.append("Duplicate reagent role IDs found")
        
    reagent_family_ids = [r["id"] for r in reagent_families]
    if len(reagent_family_ids) != len(set(reagent_family_ids)):
        issues.append("Duplicate reagent family IDs found")
    
    if not issues:
        print("   ✓ No duplicate IDs found")
    
    # Check 2: Missing SMARTS patterns
    print("\n2. Checking for missing SMARTS patterns in reactant types...")
    reactants_without_smarts = []
    for rt in reactant_types:
        if not rt.get("smarts"):
            members_with_smarts = [m for m in rt.get("members", []) if m.get("smarts")]
            if not members_with_smarts:
                reactants_without_smarts.append(rt["id"])
    
    if reactants_without_smarts:
        issues.append(f"Reactant types without SMARTS: {', '.join(reactants_without_smarts[:5])}...")
        print(f"   ⚠ {len(reactants_without_smarts)} reactant types have no SMARTS patterns")
    else:
        print("   ✓ All reactant types have SMARTS patterns")
    
    # Check 3: Reactant member SMARTS validation
    print("\n3. Checking reactant member SMARTS...")
    member_issues = []
    for rt in reactant_types:
        for member in rt.get("members", []):
            if not member.get("id"):
                member_issues.append(f"Member in {rt['id']} has no ID")
            if not member.get("smarts"):
                member_issues.append(f"Member {member.get('id', 'unknown')} in {rt['id']} has no SMARTS")
    
    if member_issues:
        issues.append(f"{len(member_issues)} reactant member issues")
        print(f"   ⚠ {len(member_issues)} member issues found")
        for issue in member_issues[:5]:
            print(f"      - {issue}")
    else:
        print("   ✓ All reactant members valid")
    
    # Check 4: Reference integrity
    print("\n4. Checking reference integrity...")
    reaction_cat_set = {r["id"] for r in reaction_categories}
    reactant_type_set = {r["id"] for r in reactant_types}
    reagent_role_set = {r["id"] for r in reagent_roles}
    reagent_family_set = {r["id"] for r in reagent_families}
    
    ref_issues = []
    for rt in reaction_types:
        # Check category_id
        if rt.get("category_id") and rt["category_id"] not in reaction_cat_set:
            ref_issues.append(f"ReactionType {rt['id']} references unknown category {rt['category_id']}")
        
        # Check reactant requirements
        for req in rt.get("reactants", []):
            reactant_id = req.get("reactant_type_id")
            if reactant_id and reactant_id not in reactant_type_set:
                ref_issues.append(f"ReactionType {rt['id']} references unknown reactant {reactant_id}")
        
        # Check required roles
        for role_req in rt.get("required_roles", []):
            role_id = role_req.get("role_id")
            if role_id and role_id not in reagent_role_set:
                ref_issues.append(f"ReactionType {rt['id']} references unknown role {role_id}")
            
            default_family = role_req.get("default_family_id")
            if default_family and default_family not in reagent_family_set:
                ref_issues.append(f"ReactionType {rt['id']} references unknown family {default_family}")
    
    for family in reagent_families:
        if family.get("role_id") and family["role_id"] not in reagent_role_set:
            ref_issues.append(f"ReagentFamily {family['id']} references unknown role {family['role_id']}")
    
    if ref_issues:
        issues.append(f"{len(ref_issues)} reference integrity issues")
        print(f"   ⚠ {len(ref_issues)} reference issues found")
        for issue in ref_issues[:5]:
            print(f"      - {issue}")
        if len(ref_issues) > 5:
            print(f"      ... and {len(ref_issues) - 5} more")
    else:
        print("   ✓ All references valid")
    
    # Check 5: Alias integrity
    print("\n5. Checking alias integrity...")
    alias_issues = []
    all_entity_ids = {
        "reaction_type": {r["id"] for r in reaction_types},
        "reactant_type": reactant_type_set,
        "reagent_role": reagent_role_set,
        "reagent_family": reagent_family_set,
        "reaction_category": reaction_cat_set,
    }
    
    for alias in aliases:
        entity_type = alias.get("entity_type")
        entity_id = alias.get("entity_id")
        
        if entity_type not in all_entity_ids:
            alias_issues.append(f"Alias '{alias.get('alias')}' has unknown entity_type '{entity_type}'")
        elif entity_id and entity_id not in all_entity_ids.get(entity_type, set()):
            alias_issues.append(f"Alias '{alias.get('alias')}' references unknown {entity_type} '{entity_id}'")
    
    if alias_issues:
        issues.append(f"{len(alias_issues)} alias issues")
        print(f"   ⚠ {len(alias_issues)} alias issues found")
        for issue in alias_issues[:5]:
            print(f"      - {issue}")
    else:
        print("   ✓ All aliases valid")
    
    # Summary
    print("\n" + "=" * 60)
    if issues:
        print(f"❌ FOUND {len(issues)} ISSUE CATEGORIES:")
        for i, issue in enumerate(issues, 1):
            print(f"  {i}. {issue}")
    else:
        print("✓ NO ISSUES FOUND - Taxonomy data is valid!")
    print("=" * 60)
    
    return issues


def check_analysis_module_functionality():
    """Test analysis module functions"""
    print("\n" + "=" * 60)
    print("CHECKING ANALYSIS MODULE FUNCTIONALITY")
    print("=" * 60)
    
    issues = []
    
    try:
        from chemtools.analysis import (
            classify_reactant_smiles,
            classify_reactant_category,
            get_all_reactant_matches,
            analyze_reaction,
            normalize_reactant_identifier,
        )
        print("\n✓ All imports successful")
    except Exception as e:
        print(f"\n❌ Import error: {e}")
        return [f"Import error: {e}"]
    
    # Test 1: Classify simple reactants
    print("\n1. Testing reactant classification...")
    test_smiles = [
        ("c1ccccc1Br", "ArBr"),  # aryl bromide
        ("CCBr", "Alkyl-Br"),  # alkyl bromide
        ("c1ccccc1B(O)O", "ArB(OH)2"),  # aryl boronic acid
        ("CCN", "RNH2"),  # primary amine
    ]
    
    for smiles, expected_category in test_smiles:
        try:
            result = classify_reactant_category(smiles)
            if result:
                print(f"   ✓ {smiles} -> {result}")
            else:
                issues.append(f"No classification for {smiles} (expected {expected_category})")
                print(f"   ⚠ {smiles} -> None (expected {expected_category})")
        except Exception as e:
            issues.append(f"Error classifying {smiles}: {e}")
            print(f"   ❌ {smiles} -> Error: {e}")
    
    # Test 2: Analyze a simple reaction
    print("\n2. Testing reaction analysis...")
    test_reactions = [
        "c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1",  # Suzuki
        "c1ccccc1Br.CCN>>c1ccccc1NCC",  # Buchwald-Hartwig
    ]
    
    for rxn_smiles in test_reactions:
        try:
            result = analyze_reaction(rxn_smiles)
            family = result.get("family", {}).get("canonical_id")
            print(f"   ✓ Reaction analyzed -> family: {family}")
            
            # Check structure
            if "reactants" not in result:
                issues.append(f"Missing 'reactants' in reaction analysis")
            if "family" not in result:
                issues.append(f"Missing 'family' in reaction analysis")
                
        except Exception as e:
            issues.append(f"Error analyzing reaction: {e}")
            print(f"   ❌ Error: {e}")
    
    # Test 3: Normalize identifiers
    print("\n3. Testing identifier normalization...")
    test_ids = [
        ("ArBr", "aryl bromide alias"),
        ("Buchwald-Hartwig", "reaction type alias"),
        ("Suzuki-Miyaura", "reaction type alias"),
    ]
    
    for test_id, desc in test_ids:
        try:
            result = normalize_reactant_identifier(test_id)
            print(f"   ✓ '{test_id}' -> {result}")
        except Exception as e:
            issues.append(f"Error normalizing {test_id}: {e}")
            print(f"   ❌ '{test_id}' -> Error: {e}")
    
    # Summary
    print("\n" + "=" * 60)
    if issues:
        print(f"❌ FOUND {len(issues)} FUNCTIONAL ISSUES:")
        for i, issue in enumerate(issues, 1):
            print(f"  {i}. {issue}")
    else:
        print("✓ ALL FUNCTIONALITY TESTS PASSED!")
    print("=" * 60)
    
    return issues


def check_missing_tests():
    """Check for missing test coverage"""
    print("\n" + "=" * 60)
    print("CHECKING TEST COVERAGE")
    print("=" * 60)
    
    issues = []
    
    # Expected test files
    expected_tests = {
        "test_analysis_reactants.py": "Tests for chemtools/analysis/reactants.py",
        "test_analysis_reactions.py": "Tests for chemtools/analysis/reactions.py",
        "test_analysis_integration.py": "Integration tests for analysis module",
        "test_taxonomy_models.py": "Tests for taxonomy models",
        "test_taxonomy_validation.py": "Tests for taxonomy validation",
    }
    
    tests_dir = Path("tests")
    existing_tests = {f.name for f in tests_dir.glob("test_*.py")} if tests_dir.exists() else set()
    
    print("\nExpected test files:")
    for test_file, desc in expected_tests.items():
        if test_file in existing_tests:
            print(f"   ✓ {test_file} - exists")
        else:
            print(f"   ❌ {test_file} - MISSING ({desc})")
            issues.append(f"Missing test file: {test_file}")
    
    # Summary
    print("\n" + "=" * 60)
    if issues:
        print(f"⚠ {len(issues)} test files missing")
    else:
        print("✓ All expected test files present")
    print("=" * 60)
    
    return issues


def main():
    print("\n" + "█" * 60)
    print("DIAGNOSTIC CHECK FOR CHEMTOOLS/ANALYSIS & /TAXONOMY")
    print("█" * 60)
    
    all_issues = []
    
    # Run all checks
    all_issues.extend(check_taxonomy_data_integrity())
    all_issues.extend(check_analysis_module_functionality())
    all_issues.extend(check_missing_tests())
    
    # Final summary
    print("\n" + "█" * 60)
    print("FINAL SUMMARY")
    print("█" * 60)
    
    if all_issues:
        print(f"\n❌ TOTAL: {len(all_issues)} issues found")
        print("\nRECOMMENDATIONS:")
        print("1. Fix reference integrity issues in taxonomy JSON files")
        print("2. Add missing SMARTS patterns for reactant types")
        print("3. Create comprehensive test suite for analysis module")
        print("4. Validate all aliases point to existing entities")
    else:
        print("\n✓✓✓ ALL CHECKS PASSED - System is healthy! ✓✓✓")
    
    print("█" * 60 + "\n")
    
    return len(all_issues)


if __name__ == "__main__":
    exit(main())
