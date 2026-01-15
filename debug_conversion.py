import pandas as pd
import sys
sys.path.insert(0, '.')

from app.A_convert_to_hte_format import _csv_row_to_record, extract_reagents

df = pd.read_csv(r'data\reaction_dataset\Oxidation_of_Aromatic_Side_Chains.csv')

print("Testing first few rows:\n")
for i in range(min(5, len(df))):
    row = df.iloc[i].to_dict()
    record = _csv_row_to_record(row)
    
    print(f"Row {i}:")
    print(f"  SMILES: {record.get('smiles')[:50]}...")
    print(f"  Has '>>': {'>>' in record.get('smiles', '')}")
    
    if ">>" not in record.get('smiles', ''):
        print(f"  SKIPPED: No >> in SMILES")
        continue
    
    reagents = extract_reagents(record, csv_row=row)
    print(f"  Catalyst: {reagents.get('catalyst', 'NONE')}")
    print(f"  Would drop (no catalyst): {not reagents.get('catalyst')}")
    print()
