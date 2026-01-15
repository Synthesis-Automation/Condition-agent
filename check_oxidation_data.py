import pandas as pd

df = pd.read_csv(r'data\reaction_dataset\Oxidation_of_Aromatic_Side_Chains.csv')
print(f'Total rows: {len(df)}')
print(f'Rows with catalyst_cas: {df["catalyst_cas"].notna().sum()}')
print(f'Rows with catalyst_amd: {df["catalyst_amd"].notna().sum()}')
print(f'Rows with reaction_smiles: {df["reaction_smiles"].notna().sum()}')

print('\nFirst 3 rows:')
for i in range(min(3, len(df))):
    r = df.iloc[i]
    print(f'\nRow {i}:')
    print(f'  catalyst_cas: {r.get("catalyst_cas")}')
    print(f'  catalyst_amd: {r.get("catalyst_amd")}')
    print(f'  reaction_smiles: {str(r.get("reaction_smiles"))[:80]}')
    print(f'  Has >>: {">> " in str(r.get("reaction_smiles"))}')
