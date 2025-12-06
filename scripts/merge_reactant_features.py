#!/usr/bin/env python3
"""
Merge reactant type features into calculable_features.json
"""

import json
from pathlib import Path

# Paths
REPO_ROOT = Path(__file__).parent.parent
CALCULABLE_FILE = REPO_ROOT / "chemtools" / "taxonomy" / "data" / "calculable_features.json"
GENERATED_FILE = REPO_ROOT / "scripts" / "reactant_features_generated.json"
BACKUP_FILE = REPO_ROOT / "chemtools" / "taxonomy" / "data" / "calculable_features_backup.json"

def main():
    print("=" * 80)
    print("Merging Reactant Type Features into calculable_features.json")
    print("=" * 80)
    
    # Load existing calculable features
    print(f"\nLoading existing features from: {CALCULABLE_FILE}")
    with open(CALCULABLE_FILE) as f:
        calculable_data = json.load(f)
    
    current_features = len(calculable_data.get("features", []))
    current_derived = len(calculable_data.get("derived_shortcuts", []))
    print(f"Current features: {current_features}")
    print(f"Current derived shortcuts: {current_derived}")
    
    # Load generated reactant features
    print(f"\nLoading generated features from: {GENERATED_FILE}")
    with open(GENERATED_FILE) as f:
        generated_data = json.load(f)
    
    new_features = len(generated_data["features"])
    new_derived = len(generated_data["derived_shortcuts"])
    print(f"New reactant features: {new_features}")
    print(f"New reactant derived: {new_derived}")
    
    # Create backup
    print(f"\nCreating backup: {BACKUP_FILE}")
    with open(BACKUP_FILE, 'w') as f:
        json.dump(calculable_data, f, indent=2)
    
    # Merge features
    print("\nMerging features...")
    calculable_data["features"].extend(generated_data["features"])
    
    # Add derived shortcuts section if it doesn't exist
    if "derived_shortcuts" not in calculable_data:
        calculable_data["derived_shortcuts"] = []
    
    calculable_data["derived_shortcuts"].extend(generated_data["derived_shortcuts"])
    
    # Update version and description
    old_version = calculable_data.get("version", "unknown")
    calculable_data["version"] = "2025-11-09.v4.0-reactant-types"
    calculable_data["description"] = (
        calculable_data.get("description", "") + 
        " | EXTENDED: Added 225 reactant type features from taxonomy system (176 members + 49 categories)."
    )
    
    # Update schema notes
    if "categories" in calculable_data.get("schema_notes", {}):
        if "reactant_types" not in calculable_data["schema_notes"]["categories"]:
            calculable_data["schema_notes"]["categories"].append("reactant_types")
    
    # Add expansion note
    if "expansion_notes" in calculable_data.get("schema_notes", {}):
        calculable_data["schema_notes"]["expansion_notes"].append(
            "Added reactant type features from taxonomy (ArBr_reactant, ArX_reactant, etc.)"
        )
    
    # Save merged file
    print(f"\nSaving merged features to: {CALCULABLE_FILE}")
    with open(CALCULABLE_FILE, 'w') as f:
        json.dump(calculable_data, f, indent=2, ensure_ascii=False)
    
    print("\n✓ Successfully merged features!")
    print(f"\nFinal counts:")
    print(f"  - Total features: {len(calculable_data['features'])} (+{new_features})")
    print(f"  - Total derived: {len(calculable_data['derived_shortcuts'])} (+{new_derived})")
    print(f"  - Version: {old_version} → {calculable_data['version']}")
    
    print(f"\n✓ Backup saved to: {BACKUP_FILE}")
    print("\n" + "=" * 80)
    print("Merge complete! Next steps:")
    print("1. Test the updated calculable.py with new features")
    print("2. Add utility functions for reactant type queries")
    print("3. Create backward-compatible wrappers")
    print("=" * 80)


if __name__ == "__main__":
    main()
