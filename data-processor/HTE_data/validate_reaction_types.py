"""Validate and analyze the canonical reaction type taxonomy coverage."""

from collections import defaultdict
from pathlib import Path

import pandas as pd

from chemtools.reagent import (
    get_reaction_type_definitions,
)

def validate_reaction_types():
    # Load reaction types JSON
    reaction_data = get_reaction_type_definitions()
    
    # Load z-Score CSV
    csv_path = Path(__file__).parent / "z-Score Peaks with FG.csv"
    df = pd.read_csv(csv_path, encoding='utf-8-sig')
    zscore_types = set(df['Reaction Type'].dropna().unique())
    
    print("=" * 80)
    print("REACTION TYPES JSON VALIDATION")
    print("=" * 80)
    
    # Count reactions by category
    total_reactions = 0
    category_counts = defaultdict(int)
    all_reaction_ids = set()
    all_aliases = set()
    
    print(f"\n{'Category':<40} {'Reactions':<10}")
    print("-" * 80)
    
    for category_key, category_info in reaction_data.items():
        category_name = category_info['category']
        reaction_count = len(category_info['reactions'])
        category_counts[category_name] = reaction_count
        total_reactions += reaction_count
        
        print(f"{category_name:<40} {reaction_count:<10}")
        
        # Collect all IDs and aliases
        for reaction in category_info['reactions']:
            all_reaction_ids.add(reaction['id'])
            all_aliases.update(reaction.get('aliases', []))
    
    print("-" * 80)
    print(f"{'TOTAL':<40} {total_reactions:<10}")
    
    print(f"\n\nTotal unique reaction IDs: {len(all_reaction_ids)}")
    print(f"Total aliases: {len(all_aliases)}")
    
    # Check coverage of z-Score dataset
    print("\n" + "=" * 80)
    print("COVERAGE ANALYSIS: z-Score Dataset Reaction Types")
    print("=" * 80)
    
    # Create mapping from z-Score types to JSON reactions
    zscore_to_json = {}
    
    for zscore_type in zscore_types:
        matched = False
        
        # Direct ID match
        for category_info in reaction_data.values():
            for reaction in category_info['reactions']:
                # Check ID match (with normalization)
                if reaction['id'].replace('-', ' ').replace('_', ' ').lower() == zscore_type.replace('-', ' ').replace('_', ' ').lower():
                    zscore_to_json[zscore_type] = reaction['id']
                    matched = True
                    break
                
                # Check alias match
                for alias in reaction.get('aliases', []):
                    if alias.lower() == zscore_type.lower():
                        zscore_to_json[zscore_type] = reaction['id']
                        matched = True
                        break
                
                # Check if zscore_type contains the reaction name
                if zscore_type.lower() in reaction['name'].lower() or reaction['name'].lower() in zscore_type.lower():
                    zscore_to_json[zscore_type] = reaction['id']
                    matched = True
                    break
            
            if matched:
                break
    
    # Report coverage
    covered = len(zscore_to_json)
    total_zscore = len(zscore_types)
    coverage_pct = (covered / total_zscore * 100) if total_zscore > 0 else 0
    
    print(f"\nCovered: {covered}/{total_zscore} ({coverage_pct:.1f}%)")
    
    print("\n" + "-" * 80)
    print(f"{'z-Score Type':<40} {'Count':<10} {'Mapped To':<30}")
    print("-" * 80)
    
    # Get counts
    type_counts = df['Reaction Type'].value_counts()
    
    for zscore_type in sorted(zscore_types):
        count = type_counts.get(zscore_type, 0)
        mapped = zscore_to_json.get(zscore_type, "NOT MAPPED")
        status = "✓" if zscore_type in zscore_to_json else "✗"
        print(f"{status} {zscore_type:<38} {count:<10} {mapped:<30}")
    
    # Show unmapped types
    unmapped = [t for t in zscore_types if t not in zscore_to_json]
    if unmapped:
        print("\n" + "=" * 80)
        print(f"UNMAPPED TYPES ({len(unmapped)}):")
        print("=" * 80)
        for utype in sorted(unmapped):
            count = type_counts.get(utype, 0)
            print(f"  - {utype} ({count} reactions)")
    
    # Show reaction dataset compatibility
    print("\n" + "=" * 80)
    print("REACTION DATASET COMPATIBILITY")
    print("=" * 80)
    
    dataset_families = {
        'Amide formation': 'Amide-coupling',
        'Suzuki': 'Suzuki-Miyaura',
        'C_N_Coupling': 'CN-Coupling',
        'C_O_Coupling': 'CO-Coupling',
        'C_S_Coupling': 'CS-Coupling'
    }
    
    for family, expected_id in dataset_families.items():
        exists = expected_id in all_reaction_ids
        status = "✓" if exists else "✗"
        print(f"{status} {family:<30} → {expected_id:<30} {'EXISTS' if exists else 'MISSING'}")
    
    return reaction_data

if __name__ == "__main__":
    validate_reaction_types()
