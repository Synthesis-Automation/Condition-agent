"""Extract and analyze reaction types from z-Score CSV and reaction dataset."""

import pandas as pd
from pathlib import Path
import json

# Load z-Score CSV
csv_path = Path(__file__).parent / "z-Score Peaks with FG.csv"
df = pd.read_csv(csv_path, encoding='utf-8-sig')

print("=" * 70)
print("REACTION TYPES FROM Z-SCORE CSV")
print("=" * 70)
print(df['Reaction Type'].value_counts())
print(f"\nTotal unique reaction types: {df['Reaction Type'].nunique()}")
print(f"Total reactions: {len(df)}")

# Get unique reaction types
zscore_types = set(df['Reaction Type'].dropna().unique())

# Try to load reaction dataset types
dataset_types = set()
reaction_data_dir = Path(__file__).parent.parent.parent / "data" / "reaction_dataset"

if reaction_data_dir.exists():
    print("\n" + "=" * 70)
    print("REACTION TYPES FROM REACTION DATASET")
    print("=" * 70)
    
    # Look for JSONL files
    jsonl_files = list(reaction_data_dir.glob("*.jsonl"))
    
    for jsonl_file in jsonl_files:
        print(f"\nAnalyzing: {jsonl_file.name}")
        try:
            with open(jsonl_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    if line.strip():
                        try:
                            data = json.loads(line)
                            if 'family' in data:
                                dataset_types.add(data['family'])
                            elif 'reaction_type' in data:
                                dataset_types.add(data['reaction_type'])
                        except json.JSONDecodeError:
                            continue
        except Exception as e:
            print(f"  Error reading {jsonl_file.name}: {e}")
    
    print(f"\nUnique reaction types from dataset: {sorted(dataset_types)}")
    print(f"Total: {len(dataset_types)}")

# Combine all types
print("\n" + "=" * 70)
print("COMBINED REACTION TYPES")
print("=" * 70)
all_types = sorted(zscore_types | dataset_types)
print(f"\nTotal unique reaction types: {len(all_types)}\n")

for idx, rtype in enumerate(all_types, 1):
    source = []
    if rtype in zscore_types:
        count = len(df[df['Reaction Type'] == rtype])
        source.append(f"z-Score: {count}")
    if rtype in dataset_types:
        source.append("Dataset")
    print(f"{idx:2}. {rtype:<35} ({', '.join(source)})")
