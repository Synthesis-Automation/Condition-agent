"""
Fix common validation issues in migrated v2 rule databases.
"""

import json
from pathlib import Path


def fix_additives_field(data: dict):
    """Convert additives from string to array if needed."""
    modified = False
    
    # Fix in default_rule
    if 'default_rule' in data and 'conditions' in data['default_rule']:
        conditions = data['default_rule']['conditions']
        if 'additives' in conditions and isinstance(conditions['additives'], str):
            conditions['additives'] = [conditions['additives']]
            modified = True
    
    # Fix in base_rules
    if 'base_rules' in data:
        for rule in data['base_rules']:
            if 'conditions' in rule and 'additives' in rule['conditions']:
                if isinstance(rule['conditions']['additives'], str):
                    rule['conditions']['additives'] = [rule['conditions']['additives']]
                    modified = True
    
    return modified


def main():
    rule_db_v2_dir = Path("data/rule_db_v2")
    
    print("="*80)
    print("Fixing validation issues in v2 rule databases")
    print("="*80)
    
    fixed_count = 0
    
    for json_file in rule_db_v2_dir.glob("*.json"):
        # Read file
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Fix issues
        modified = fix_additives_field(data)
        
        if modified:
            # Save back
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            
            print(f"✅ Fixed {json_file.name}")
            fixed_count += 1
    
    print(f"\nSummary: {fixed_count} files fixed")


if __name__ == '__main__':
    main()
