#!/usr/bin/env python3
"""Check Cu-catalyzed C-N coupling data"""
from chemtools.HTE import HTEAnalytics

a = HTEAnalytics()
df = a.list_reactant_pairs('C-N', 'Cu', min_experiments=1)

print(f'Total pairs with min_exp=1: {len(df)}')
print(f'Total experiments: {df["Num_Experiments"].sum()}')
print('\nTop 10 pairs:')
print(df.head(10)[['Reactant_A_Type', 'Reactant_B_Type', 'Num_Experiments', 'Success_Rate']].to_string(index=False))

print('\n\nWith min_experiments=5:')
df5 = a.list_reactant_pairs('C-N', 'Cu', min_experiments=5)
print(f'Total pairs: {len(df5)}')

print('\n\nWith min_experiments=10:')
df10 = a.list_reactant_pairs('C-N', 'Cu', min_experiments=10)
print(f'Total pairs: {len(df10)}')
