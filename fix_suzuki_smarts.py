"""
Fix SMARTS patterns in Suzuki_protocols.json to match the actual reaction SMILES

The issue: SMARTS patterns use B(O[H])O[H] (explicit H) but reaction SMILES use B(O)O (implicit H)

Solution: Update SMARTS patterns to be more flexible and match both forms
"""

import json
from pathlib import Path

def fix_suzuki_smarts():
    """Fix SMARTS patterns in Suzuki_protocols.json"""
    
    # Path to the file
    protocol_file = Path("data/protocol_db/Suzuki_protocols.json")
    
    print("=" * 80)
    print("Fixing SMARTS Patterns in Suzuki_protocols.json")
    print("=" * 80)
    print()
    
    # Load the file
    with open(protocol_file, 'r', encoding='utf-8') as f:
        protocols = json.load(f)
    
    print(f"Loaded {len(protocols)} protocols")
    print()
    
    # Fix each protocol
    fixes_made = 0
    for idx, protocol in enumerate(protocols):
        reaction = protocol.get('reaction', {})
        reaction_smiles = reaction.get('reaction_smiles', '')
        old_smarts = reaction.get('reaction_SMARTS', [])
        
        if not old_smarts:
            continue
        
        print(f"Protocol [{idx}]: {reaction.get('family', 'Unknown')}")
        print(f"  Reaction: {reaction_smiles[:80]}...")
        print(f"  Old SMARTS: {old_smarts}")
        
        # Determine the correct SMARTS pattern based on reaction type
        new_smarts = []
        fixed = False
        
        for smarts in old_smarts:
            # Fix boronic acid pattern: B(O[H])O[H] -> B(O)O
            # Make it flexible to accept both implicit and explicit H
            if 'B(O[H])O[H]' in smarts:
                # Replace with implicit H version (matches the actual SMILES)
                new_pattern = smarts.replace('B(O[H])O[H]', 'B(O)O')
                new_smarts.append(new_pattern)
                fixed = True
                print(f"  ✓ Fixed: {new_pattern}")
            else:
                new_smarts.append(smarts)
        
        if fixed:
            protocol['reaction']['reaction_SMARTS'] = new_smarts
            fixes_made += 1
        
        print()
    
    # Save the updated file
    if fixes_made > 0:
        with open(protocol_file, 'w', encoding='utf-8') as f:
            json.dump(protocols, f, indent=2, ensure_ascii=False)
        
        print("=" * 80)
        print(f"✅ Fixed {fixes_made} protocol(s)")
        print(f"✅ Updated file: {protocol_file}")
        print("=" * 80)
    else:
        print("=" * 80)
        print("No fixes needed")
        print("=" * 80)

if __name__ == "__main__":
    fix_suzuki_smarts()
