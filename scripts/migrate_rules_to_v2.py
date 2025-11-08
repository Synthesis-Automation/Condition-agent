"""
Migrate rule databases to v2.0 schema format.

This script converts rule databases from v1 to v2.0 format,
adding required metadata and schema fields.
"""

import json
from pathlib import Path
from datetime import datetime


def migrate_rule_to_v2(old_data: dict, filename: str) -> dict:
    """Convert rule database from v1 to v2 format."""
    
    # Extract family from reaction_type or name
    family = old_data.get('reaction_type', old_data.get('name', 'unknown'))
    family = family.replace('-', '_').replace(' ', '_')
    
    # Generate ID from filename
    rule_id = filename.replace('_db.json', '').lower() + '_v2'
    
    # Create v2 structure
    v2_data = {
        "schema_version": "2.0",
        "source_type": "rule",
        "metadata": {
            "id": rule_id,
            "name": old_data.get('name', family),
            "version": "v2.0",
            "created_date": datetime.now().strftime('%Y-%m-%d'),
            "status": "active",
            "tags": extract_tags(old_data, family)
        },
        "reaction": {
            "family": family,
            "reference_reactions": old_data.get('reaction', {}).get('reference_reactions', []),
            "scope": {
                "scope_type": "broad",
                "compatible_functional_groups": [],
                "incompatible_functional_groups": []
            },
            "notes": old_data.get('evaluation', '')
        },
        "applies_if": old_data.get('applies_if', {}),
        "default_rule": old_data.get('default_rule', {}),
        "base_rules": old_data.get('base_rules', []),
        "modifiers": old_data.get('modifiers', [])
    }
    
    return v2_data


def extract_tags(data: dict, family: str) -> list:
    """Extract tags from various fields."""
    tags = []
    
    # Add family as base tag
    if 'suzuki' in family.lower():
        tags.extend(['suzuki', 'palladium', 'cross-coupling', 'boronic-acid'])
    elif 'c_n' in family.lower() or 'buchwald' in family.lower():
        tags.extend(['C-N-coupling', 'amination', 'cross-coupling'])
        if 'pd' in family.lower():
            tags.append('palladium')
        if 'cu' in family.lower():
            tags.append('copper')
    elif 'c_o' in family.lower():
        tags.extend(['C-O-coupling', 'etherification', 'cross-coupling'])
    elif 'rcm' in family.lower():
        tags.extend(['RCM', 'metathesis', 'ruthenium', 'ring-closing'])
    elif 'sonogashira' in family.lower():
        tags.extend(['sonogashira', 'palladium', 'copper', 'alkyne'])
    elif 'snar' in family.lower():
        tags.extend(['SNAr', 'nucleophilic-substitution', 'aromatic'])
    elif 'amide' in family.lower():
        tags.extend(['amide-formation', 'coupling', 'peptide'])
    elif 'reductive' in family.lower() and 'amination' in family.lower():
        tags.extend(['reductive-amination', 'amine-synthesis'])
    
    return tags


def main():
    rule_db_dir = Path("data/rule_db")
    output_dir = Path("data/rule_db_v2")
    output_dir.mkdir(exist_ok=True)
    
    print("="*80)
    print("Migrating rule databases to v2.0 format")
    print("="*80)
    
    migrated_count = 0
    
    for json_file in rule_db_dir.glob("*.json"):
        if json_file.name == 'sonogashira_v2.json':
            print(f"\nSkipping {json_file.name} (already v2)")
            continue
        
        print(f"\nMigrating {json_file.name}...")
        
        # Read v1 file
        with open(json_file, 'r', encoding='utf-8') as f:
            old_data = json.load(f)
        
        # Check if it has reference_reactions
        if 'reaction' not in old_data or 'reference_reactions' not in old_data.get('reaction', {}):
            print(f"  ⚠️  No reference_reactions found - skipping (run add_reference_reactions.py first)")
            continue
        
        # Convert to v2
        v2_data = migrate_rule_to_v2(old_data, json_file.name)
        
        # Save v2 file
        output_path = output_dir / json_file.name
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(v2_data, f, indent=2, ensure_ascii=False)
        
        print(f"  ✅ Migrated to {output_path}")
        migrated_count += 1
    
    print("\n" + "="*80)
    print(f"Summary: {migrated_count} files migrated to v2.0")
    print("="*80)


if __name__ == '__main__':
    main()
