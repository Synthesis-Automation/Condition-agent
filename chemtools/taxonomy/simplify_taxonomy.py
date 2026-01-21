#!/usr/bin/env python3
"""
Simplify taxonomy by:
1. Adding "-" prefix to all substituent group IDs
2. Removing all "name" fields (ID is the display name)
3. Simplifying compound ID generation to A+B (no separator)
"""

import json
from pathlib import Path
from typing import Dict, List, Any


def simplify_groups(groups_file: Path) -> Dict[str, str]:
    """
    Update groups:
    - Add "-" prefix to all substituent IDs (if not already present)
    - Remove all "name" fields
    - Return mapping of old_id -> new_id
    """
    with open(groups_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    groups = data.get('groups', [])
    id_mapping = {}
    
    print(f"Processing {len(groups)} groups...")
    
    for group in groups:
        old_id = group['id']
        kind = group.get('kind', '')
        
        # Determine new ID
        if kind == 'substituent':
            # Add "-" prefix if not already present
            if not old_id.startswith('-'):
                new_id = f"-{old_id}"
                id_mapping[old_id] = new_id
                group['id'] = new_id
                print(f"  Substituent: {old_id:20s} → {new_id}")
            else:
                id_mapping[old_id] = old_id  # Already has prefix
        else:
            # Scaffold: keep as-is
            id_mapping[old_id] = old_id
        
        # Remove "name" field
        if 'name' in group:
            del group['name']
    
    # Save updated groups
    with open(groups_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write('\n')
    
    print(f"\n✓ Updated {groups_file.name}")
    print(f"  - Renamed {len([v for k,v in id_mapping.items() if k!=v])} substituent IDs")
    print(f"  - Removed 'name' field from all {len(groups)} groups")
    
    return id_mapping


def update_compounds(compounds_file: Path, id_mapping: Dict[str, str]):
    """
    Update compound A/B references to use new group IDs
    """
    with open(compounds_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    compounds = data.get('compounds', [])
    update_count = 0
    
    print(f"\nProcessing {len(compounds)} compounds...")
    
    for compound in compounds:
        updated = False
        
        # Update A reference
        if 'A' in compound:
            old_a = compound['A']
            new_a = id_mapping.get(old_a, old_a)
            if new_a != old_a:
                compound['A'] = new_a
                updated = True
        
        # Update B reference
        if 'B' in compound:
            old_b = compound['B']
            new_b = id_mapping.get(old_b, old_b)
            if new_b != old_b:
                compound['B'] = new_b
                updated = True
        
        if updated:
            update_count += 1
    
    # Save updated compounds
    with open(compounds_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write('\n')
    
    print(f"✓ Updated {compounds_file.name}")
    print(f"  - Updated {update_count} compound references")


def update_motif_registry():
    """
    Update motif_registry.py to use simple A+B concatenation
    """
    registry_file = Path("chemtools/featurizers/motif_registry.py")
    
    with open(registry_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find and replace the ID generation logic
    old_code = '''        # Auto-generate compound ID from A-B pattern if not explicitly provided
        if "id" not in entry:
            group_a = str(entry.get("A") or "")
            group_b = str(entry.get("B") or "")
            if group_a and group_b:
                # Normalize display names: strip _Subst suffix for cleaner IDs
                display_a = group_a.replace("_Subst", "")
                display_b = group_b.replace("_Subst", "")
                entry["id"] = f"{display_a}-{display_b}"'''
    
    new_code = '''        # Auto-generate compound ID from A+B (no separator needed)
        if "id" not in entry:
            group_a = str(entry.get("A") or "")
            group_b = str(entry.get("B") or "")
            if group_a and group_b:
                # Simple concatenation: A+B (substituents already have "-" prefix)
                entry["id"] = f"{group_a}{group_b}"'''
    
    if old_code in content:
        content = content.replace(old_code, new_code)
        
        with open(registry_file, 'w', encoding='utf-8') as f:
            f.write(content)
        
        print(f"\n✓ Updated {registry_file.name}")
        print(f"  - Simplified ID generation to A+B concatenation")
    else:
        print(f"\n⚠ Could not find target code in {registry_file.name}")


def main():
    print("=" * 70)
    print("Simplifying Organic Taxonomy")
    print("=" * 70)
    print()
    
    taxonomy_dir = Path("chemtools/taxonomy/data")
    groups_file = taxonomy_dir / "organic_groups.v1.3.json"
    compounds_file = taxonomy_dir / "organic_compounds.v1.3.json"
    
    # Step 1: Update groups and get ID mapping
    id_mapping = simplify_groups(groups_file)
    
    # Step 2: Update compound references
    update_compounds(compounds_file, id_mapping)
    
    # Step 3: Update motif registry code
    update_motif_registry()
    
    print("\n" + "=" * 70)
    print("✓ Taxonomy Simplification Complete!")
    print("=" * 70)
    print("\nKey improvements:")
    print("  • All substituent IDs have '-' prefix (e.g., -Cl, -Br, -NH2)")
    print("  • Removed 'name' field (ID is the display name)")
    print("  • Compound IDs are simple A+B concatenation:")
    print("    - Ar + -Cl = 'Ar-Cl'")
    print("    - Ar + -Ar = 'Ar-Ar'")
    print("    - Alkyl + -OH = 'Alkyl-OH'")
    print("\nNext: Run validation to verify all changes")


if __name__ == "__main__":
    main()
