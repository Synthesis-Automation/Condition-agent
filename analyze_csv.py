import pandas as pd

# Read the CSV
df = pd.read_csv(r'c:\Git-softwares\Condition-agent\data-processor\other_data\z-Score Peaks with FG.csv')

print(f'Total rows: {len(df)}')
print(f'\nColumns: {list(df.columns)}')
print(f'\nReaction Types:')
print(df['Reaction Type'].value_counts())

# Check reagent columns
reagent_cols = ['Additive', 'Base', 'Catalyst', 'Coupling Reagent', 'Solvent', 'Ligand', 'Secondary Solvent', 'Tertiary Solvent']
print(f'\n\nReagent columns sample:')
for col in reagent_cols:
    if col in df.columns:
        unique_count = df[col].nunique()
        print(f'{col}: {unique_count} unique values')
