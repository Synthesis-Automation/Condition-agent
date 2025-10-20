import json
from pathlib import Path

# List all protocols
db = Path('data/protocol_db')
files = sorted(db.glob('*.json'))

print(f"Found {len(files)} protocols in database:")
print("=" * 80)
print()

for i, f in enumerate(files, 1):
    with open(f, encoding='utf-8') as fp:
        protocol = json.load(fp)
    
    family = protocol.get('reaction', {}).get('family', 'N/A')
    smarts = protocol.get('reaction', {}).get('reaction_SMARTS', [])
    reaction_smiles = protocol.get('reaction', {}).get('reaction_smiles', 'N/A')
    
    print(f"{i}. {f.stem}")
    print(f"   Family: {family}")
    print(f"   Reaction: {reaction_smiles}")
    print(f"   SMARTS: {smarts}")
    print()

# Check if there are any Sonogashira protocols
sonogashira = []
for f in files:
    if 'sonogashira' in f.stem.lower():
        sonogashira.append(f)
    else:
        try:
            with open(f, encoding='utf-8') as fp:
                protocol = json.load(fp)
                family = protocol.get('reaction', {}).get('family', '').lower()
                if 'sonogashira' in family:
                    sonogashira.append(f)
        except:
            pass

print("=" * 80)
print(f"Sonogashira protocols found: {len(sonogashira)}")
if sonogashira:
    for f in sonogashira:
        print(f"  - {f.stem}")
else:
    print("  None - This explains why Sonogashira reactions didn't match!")
