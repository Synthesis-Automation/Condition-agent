"""Check taxonomy for motif definitions."""
import json
from pathlib import Path

# Load organic compounds taxonomy
path = Path('chemtools/taxonomy/data/organic_compounds.v1.3.json')
with open(path) as f:
    data = json.load(f)

compounds = data.get('compounds', [])
print(f'Total compounds in taxonomy: {len(compounds)}')

print()
print('=== Ar-Br type motifs ===')
for c in compounds:
    a = c.get('A', '')
    b = c.get('B', '')
    if 'Br' in b:
        print(f"  {a}{b}")

print()
print('=== B(OH)2 type motifs ===')
for c in compounds:
    a = c.get('A', '')
    b = c.get('B', '')
    if 'B(' in b:
        print(f"  {a}{b}")

print()
print('=== HeteroAr scaffolds ===')
for c in compounds:
    a = c.get('A', '')
    b = c.get('B', '')
    if 'HeteroAr' in a:
        print(f"  {a}{b}")

print()
print('=== Furan/Pyrimidine scaffolds ===')
for c in compounds:
    a = c.get('A', '')
    b = c.get('B', '')
    if 'Furan' in a or 'Pyrimidine' in a:
        print(f"  {a}{b}")
