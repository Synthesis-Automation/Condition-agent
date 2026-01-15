import pandas as pd
import sys
from app.A_convert_to_hte_format import _split_csv_list, _dedupe_list, CAS_PATTERN
from chemtools.reagent.lookup import find_reagent

df = pd.read_csv(r'c:\Git-softwares\Condition-agent\data\reaction_dataset\C-O Coupling.csv', nrows=1)
row = df.iloc[0].to_dict()

print('Testing CAS extraction and lookup:')
print(f'\nsolvent_cas column value: {row.get("solvent_cas")}')

all_cas = []
for cas in _split_csv_list(row.get("solvent_cas", "")):
    all_cas.append(cas)
    
print(f'Split CAS list: {all_cas}')

for cas in _dedupe_list(all_cas):
    print(f'\n  Processing CAS: {cas}')
    print(f'    Pattern match: {CAS_PATTERN.match(cas)}')
    if not cas or not CAS_PATTERN.match(cas):
        print(f'    SKIPPED - no match')
        continue
        
    hit = None
    for r_type in ['metal_catalyst', 'ligand', 'base', 'additive', 'condensation_agent', 'other_reagent', 'solvent', 'acid', 'oxidant']:
        hit = find_reagent(cas, r_type)
        if hit:
            print(f'    Found as {r_type}')
            print(f'    Name: {hit.get("name")}')
            print(f'    Role: {hit.get("role")}')
            print(f'    Abbrev: {hit.get("abbreviation")}')
            break
    if not hit:
        print(f'    NOT FOUND')
