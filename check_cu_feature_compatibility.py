#!/usr/bin/env python
"""Check if all feature tokens in C_N_Coupling_Cu_db.json are defined in calculable_features.json."""

import json
from pathlib import Path
from typing import Set

def extract_features_from_rule_db(db: dict) -> Set[str]:
    """Extract all feature tokens referenced in a rule database."""
    features = set()
    
    # From applies_if
    if "applies_if" in db:
        if "all" in db["applies_if"]:
            features.update(db["applies_if"]["all"])
        if "any" in db["applies_if"]:
            features.update(db["applies_if"]["any"])
    
    # From base_rules
    for rule in db.get("base_rules", []):
        if "reactant_features" in rule:
            if "any" in rule["reactant_features"]:
                features.update(rule["reactant_features"]["any"])
            if "all" in rule["reactant_features"]:
                features.update(rule["reactant_features"]["all"])
        
        if "electrophile_features" in rule:
            if "any" in rule["electrophile_features"]:
                features.update(rule["electrophile_features"]["any"])
            if "all" in rule["electrophile_features"]:
                features.update(rule["electrophile_features"]["all"])
    
    # From modifiers - extract feature tokens (not symptom: tokens)
    for mod in db.get("modifiers", []):
        if "when" in mod:
            for condition in mod["when"]:
                # Skip symptom: and any generic tokens
                if not condition.startswith("symptom:") and condition != "any":
                    features.add(condition)
    
    return features

def main():
    print("="*70)
    print("FEATURE TOKEN COMPATIBILITY CHECK")
    print("="*70)
    
    # Load Cu C-N coupling database
    cu_db_path = Path("data/rule_db/C_N_Coupling_Cu_db.json")
    with open(cu_db_path, encoding='utf-8') as f:
        cu_db = json.load(f)
    
    # Load calculable features
    calc_features_path = Path("chemtools/featurizers/calculable_features.json")
    with open(calc_features_path, encoding='utf-8') as f:
        calc_features_data = json.load(f)
    
    # Extract available feature tokens (from both features and derived_shortcuts)
    available_features = set()
    if "features" in calc_features_data and isinstance(calc_features_data["features"], list):
        for feat in calc_features_data["features"]:
            if "token" in feat:
                available_features.add(feat["token"])
    if "derived_shortcuts" in calc_features_data and isinstance(calc_features_data["derived_shortcuts"], list):
        for feat in calc_features_data["derived_shortcuts"]:
            if "token" in feat:
                available_features.add(feat["token"])
    
    # Extract feature tokens used in rule db
    used_features = extract_features_from_rule_db(cu_db)
    
    print(f"\nFeature tokens used in C_N_Coupling_Cu_db.json: {len(used_features)}")
    print(f"Available features in calculable_features.json: {len(available_features)}")
    
    # Check compatibility
    missing_features = used_features - available_features
    
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    
    if missing_features:
        print(f"\n❌ MISSING FEATURES ({len(missing_features)}):")
        print("The following features are used in the rule database but NOT defined:")
        for feat in sorted(missing_features):
            print(f"   - {feat}")
        print("\n⚠️  These features need to be added to calculable_features.json")
        return False
    else:
        print("\n✅ ALL FEATURES COMPATIBLE!")
        print("All feature tokens used in the rule database are properly defined.")
        
        # Show what features are used
        print(f"\n📋 Features used from applies_if:")
        if "applies_if" in cu_db:
            if "all" in cu_db["applies_if"]:
                for feat in cu_db["applies_if"]["all"]:
                    print(f"   ✓ {feat}")
            if "any" in cu_db["applies_if"]:
                for feat in cu_db["applies_if"]["any"]:
                    print(f"   ✓ {feat}")
        
        print(f"\n📋 Additional features used in base_rules and modifiers:")
        rule_features = used_features - set(cu_db.get("applies_if", {}).get("all", [])) - set(cu_db.get("applies_if", {}).get("any", []))
        for feat in sorted(rule_features):
            if feat in available_features:
                print(f"   ✓ {feat}")
        
        return True

if __name__ == "__main__":
    success = main()
    print("\n" + "="*70)
    exit(0 if success else 1)
