"""Check file_family values in loaded data."""
from chemtools.precedent.loader import _load

_load.cache_clear()
rows = _load()

print(f'Total rows: {len(rows):,}\n')

# Get unique file_family values
families = set(r.get('file_family', 'MISSING') for r in rows)
print(f'Unique file_family values ({len(families)}):')
for fam in sorted(families):
    count = sum(1 for r in rows if r.get('file_family') == fam)
    print(f'  {fam}: {count:,}')

# Check a protocol row specifically
print('\n\nChecking first few rows for protocols source:')
from chemtools.precedent.loader import _iter_literature_files
import os

for path in _iter_literature_files():
    if 'protocols' in path.lower():
        print(f'\nFound protocol file: {os.path.basename(path)}')
        # Find rows from this file
        proto_family = path.replace('.csv', '').replace('\\', '/').split('/')[-1]
        print(f'  Expected file_family pattern: {proto_family}')
        
        matching = [r for r in rows if proto_family in str(r.get('file_family', ''))]
        print(f'  Rows with this file_family: {len(matching)}')
        if matching:
            print(f'  Sample rxn_type: {matching[0].get("rxn_type")}')
            print(f'  Sample file_family: {matching[0].get("file_family")}')
