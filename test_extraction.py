import pandas as pd
import sys
from app.A_convert_to_hte_format import extract_reagents

df = pd.read_csv(r'c:\Git-softwares\Condition-agent\data\reaction_dataset\C-O Coupling.csv', nrows=1)
row = df.iloc[0].to_dict()

print('Source data:')
print(f'  catalyst_cas: {row.get("catalyst_cas")}')
print(f'  solvent_cas: {row.get("solvent_cas")}')
print(f'  reagent_cas: {row.get("reagent_cas")}')

result = extract_reagents({}, row)
print('\nExtracted reagents:')
for k, v in result.items():
    if v:
        print(f'  {k}: {v}')
