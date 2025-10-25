"""Check all matches for propene"""
from classify_reactant import get_all_matches, load_reactant_types

rt = load_reactant_types()
matches = get_all_matches('C=CC', rt)

print('All matches for C=CC (propene):')
for m in matches:
    print(f'  - {m["member_type"]:25s} {m["name"]}')
