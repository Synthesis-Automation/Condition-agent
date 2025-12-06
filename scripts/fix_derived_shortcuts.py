"""
Update derived_shortcuts to use _present suffix instead of _reactant
"""

import json
from pathlib import Path

def update_derived_shortcuts(json_path: Path):
    """Update derived expressions to use _present suffix"""
    
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    updated_count = 0
    
    # Update derived_shortcuts
    for shortcut in data.get('derived_shortcuts', []):
        token = shortcut.get('token', '')
        derive_expr = shortcut.get('derive', '')
        
        # Update token if it ends with _reactant
        if token.endswith('_reactant'):
            new_token = token[:-9] + '_present'
            shortcut['token'] = new_token
            updated_count += 1
            print(f"  Token: {token} → {new_token}")
        
        # Update derive expression - replace all _reactant with _present
        if '_reactant' in derive_expr:
            new_expr = derive_expr.replace('_reactant', '_present')
            shortcut['derive'] = new_expr
            updated_count += 1
            if updated_count <= 15:  # Show first few
                print(f"  Expr: Updated expression for {shortcut['token']}")
    
    # Save
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"\n✓ Updated {updated_count} references in derived_shortcuts")
    return updated_count

def main():
    json_path = Path('chemtools/taxonomy/data/calculable_features.json')
    print(f"Updating derived_shortcuts in: {json_path}\n")
    
    count = update_derived_shortcuts(json_path)
    
    print(f"\n{'='*70}")
    print(f"COMPLETE: Updated {count} _reactant → _present references")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
