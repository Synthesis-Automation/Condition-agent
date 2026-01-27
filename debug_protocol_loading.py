"""Debug protocol CSV loading."""
import csv

with open('data/HTE_db/protocols/protocols_all_v2_hte.csv', 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    print(f'Total rows in protocol CSV: {len(rows)}')
    if rows:
        r = rows[0]
        print(f'\nFirst row details:')
        print(f'  reaction_id: {r.get("reaction_id", "MISSING")}')
        print(f'  detected_reaction_type: {r.get("detected_reaction_type", "MISSING")}')
        print(f'  reaction_smiles: {r.get("reaction_smiles", "MISSING")[:50]}...')
        print(f'  catalyst: {r.get("catalyst", "MISSING")}')

# Now test loader
from chemtools.precedent.loader import _make_row_from_csv

print('\n\nTesting loader on first protocol row:')
result = _make_row_from_csv(rows[0], row_index=0, file_family="protocols_all_v2_hte")
if result:
    print('✓ Successfully parsed row')
    print(f'  rxn_type: {result.get("rxn_type")}')
    print(f'  conditions core: {result.get("conditions", {}).get("catalyst")}')
else:
    print('✗ Failed to parse row')
