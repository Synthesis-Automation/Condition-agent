import json

# Read ligand file
with open('data/reagents/ligand.json', 'r') as f:
    data = json.load(f)

# Get entry at index 16
entry = data[16]

print("=" * 60)
print("LIGAND AT INDEX 16 (Missing SMILES)")
print("=" * 60)
print(f"Name:         {entry.get('name', 'UNKNOWN')}")
print(f"ID:           {entry.get('id', 'UNKNOWN')}")
print(f"Abbreviation: {entry.get('abbreviation', [])}")
print(f"Aliases:      {entry.get('aliases', [])}")
print(f"SMILES:       {entry.get('smiles', 'MISSING ❌')}")
print("=" * 60)
