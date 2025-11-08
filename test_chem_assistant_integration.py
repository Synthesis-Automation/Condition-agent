"""
Quick verification test for chem_assistant integration with new system.

Verifies:
1. Rule database paths are updated to v2
2. Sonogashira file is correctly named and accessible
3. Unified recommender with validation is working
4. chem_assistant can access the new features
"""

from pathlib import Path
from chem_assistant.chemtools_wrapper import (
    _RULE_DB_SEARCH_PATHS,
    _FAMILY_TO_RULE_DB,
    unified_recommender_tool,
    UNIFIED_RECOMMENDER_AVAILABLE,
)


def test_rule_db_paths():
    """Test that rule DB paths are updated to v2."""
    print("=" * 80)
    print("TEST 1: Rule Database Paths")
    print("=" * 80)
    
    print("\nConfigured search paths:")
    for i, path in enumerate(_RULE_DB_SEARCH_PATHS, 1):
        print(f"  {i}. {path}")
        if "rule_db_v2" in str(path):
            print(f"     ✅ Uses v2 path")
        elif "rule_db" in str(path) and "rule_db_v2" not in str(path):
            print(f"     ❌ Uses old v1 path")
    
    # Check if any v2 paths exist
    v2_paths = [p for p in _RULE_DB_SEARCH_PATHS if "rule_db_v2" in str(p)]
    if v2_paths:
        print(f"\n✅ v2 paths configured: {len(v2_paths)}")
    else:
        print(f"\n❌ No v2 paths found!")
    
    # Check if paths exist
    existing_paths = [p for p in _RULE_DB_SEARCH_PATHS if p.exists()]
    print(f"\nExisting paths: {len(existing_paths)}/{len(_RULE_DB_SEARCH_PATHS)}")
    for path in existing_paths:
        if path.name == "rule_db_v2":
            rule_count = len(list(path.glob("*.json")))
            print(f"  ✅ {path} ({rule_count} rule files)")


def test_sonogashira_mapping():
    """Test that sonogashira is correctly mapped."""
    print("\n" + "=" * 80)
    print("TEST 2: Sonogashira File Mapping")
    print("=" * 80)
    
    sonogashira_aliases = [
        k for k, v in _FAMILY_TO_RULE_DB.items()
        if "sonogashira" in k.lower()
    ]
    
    print(f"\nSonogashira aliases found: {len(sonogashira_aliases)}")
    for alias in sonogashira_aliases:
        target = _FAMILY_TO_RULE_DB[alias]
        print(f"  '{alias}' → {target}")
    
    if sonogashira_aliases:
        target_file = _FAMILY_TO_RULE_DB[sonogashira_aliases[0]]
        print(f"\nTarget file: {target_file}")
        
        # Try to find the file
        for search_path in _RULE_DB_SEARCH_PATHS:
            if not search_path.exists():
                continue
            rule_file = search_path / f"{target_file}.json"
            if rule_file.exists():
                print(f"✅ Found: {rule_file}")
                
                # Check file content
                import json
                with open(rule_file, 'r') as f:
                    data = json.load(f)
                    schema_version = data.get("schema_version", "unknown")
                    source_type = data.get("source_type", "unknown")
                    print(f"   Schema: {schema_version}, Type: {source_type}")
                break
        else:
            print(f"❌ File not found: {target_file}.json")
    else:
        print("❌ No sonogashira mappings found!")


def test_unified_recommender_validation():
    """Test unified recommender with validation enabled."""
    print("\n" + "=" * 80)
    print("TEST 3: Unified Recommender with Validation")
    print("=" * 80)
    
    if not UNIFIED_RECOMMENDER_AVAILABLE:
        print("⚠️  Unified recommender not available (DRFP not installed)")
        return
    
    # Test Sonogashira reaction
    rxn = "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"
    print(f"\nQuery: {rxn}")
    print("(Sonogashira coupling: aryl bromide + terminal alkyne)")
    
    try:
        result = unified_recommender_tool(
            reaction_smiles=rxn,
            top_k=3,
            validate_rules=True
        )
        
        if result.get("success"):
            recs = result.get("recommendations", [])
            validation_enabled = result.get("validation_enabled", False)
            
            print(f"\n✅ Got {len(recs)} recommendations")
            print(f"   Validation enabled: {validation_enabled}")
            
            for rec in recs:
                print(f"\n   {rec['rank']}. {rec['name']}")
                print(f"      Type: {rec['source_type']}, Similarity: {rec['similarity']}")
                print(f"      Family: {rec['family']}")
                if "validated" in rec:
                    print(f"      Validated: {rec['validated']}")
            
            # Check if sonogashira is in results
            sonogashira_found = any(
                "sonogashira" in rec.get("name", "").lower() or
                "sonogashira" in rec.get("family", "").lower()
                for rec in recs
            )
            
            if sonogashira_found:
                print("\n✅ Sonogashira rule/protocol found in results!")
            else:
                print("\n⚠️  Sonogashira not in top results (may have low similarity)")
        else:
            print(f"❌ Request failed: {result.get('error')}")
    
    except Exception as e:
        print(f"❌ Error: {e}")


def test_validation_parameter():
    """Test that validation parameter is properly passed."""
    print("\n" + "=" * 80)
    print("TEST 4: Validation Parameter")
    print("=" * 80)
    
    if not UNIFIED_RECOMMENDER_AVAILABLE:
        print("⚠️  Skipping - DRFP not available")
        return
    
    rxn = "Ic=Cc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)C=Cc1ccccc1"
    print(f"\nQuery: {rxn}")
    print("(Alkenyl iodide Sonogashira - should match if validation works)")
    
    try:
        # Test with validation enabled
        result_validated = unified_recommender_tool(
            reaction_smiles=rxn,
            top_k=10,
            validate_rules=True
        )
        
        # Test with validation disabled
        result_unvalidated = unified_recommender_tool(
            reaction_smiles=rxn,
            top_k=10,
            validate_rules=False
        )
        
        count_validated = len(result_validated.get("recommendations", []))
        count_unvalidated = len(result_unvalidated.get("recommendations", []))
        
        print(f"\nResults:")
        print(f"  With validation: {count_validated} recommendations")
        print(f"  Without validation: {count_unvalidated} recommendations")
        
        if count_validated <= count_unvalidated:
            print("\n✅ Validation is filtering results (fewer or equal with validation)")
        else:
            print("\n⚠️  Unexpected: More results with validation")
        
        # Check validation flag
        val_flag_validated = result_validated.get("validation_enabled", None)
        val_flag_unvalidated = result_unvalidated.get("validation_enabled", None)
        
        print(f"\nValidation flags:")
        print(f"  With validate_rules=True: {val_flag_validated}")
        print(f"  With validate_rules=False: {val_flag_unvalidated}")
        
        if val_flag_validated and not val_flag_unvalidated:
            print("\n✅ Validation flag correctly set!")
        else:
            print("\n⚠️  Validation flag issue")
    
    except Exception as e:
        print(f"❌ Error: {e}")


if __name__ == "__main__":
    print("\n" + "🧪 " * 30)
    print("CHEM_ASSISTANT INTEGRATION VERIFICATION")
    print("🧪 " * 30)
    
    test_rule_db_paths()
    test_sonogashira_mapping()
    test_unified_recommender_validation()
    test_validation_parameter()
    
    print("\n" + "=" * 80)
    print("VERIFICATION COMPLETE")
    print("=" * 80)
    print()
