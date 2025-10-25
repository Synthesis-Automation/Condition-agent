"""Process z-Score Peaks CSV to use standardized reactant and reaction types."""

from collections import defaultdict
from pathlib import Path

import pandas as pd

from chemtools.reagent import (
    build_reaction_lookup,
    build_reactant_lookup,
)


def split_mixed_values(value):
    """Split mixed reactant values (e.g., 'ArBr, ArCl' -> ['ArBr', 'ArCl'])."""
    if pd.isna(value) or value == '':
        return []
    
    # Split by comma
    parts = [p.strip() for p in str(value).split(',')]
    return parts

def process_csv(csv_path, reactant_alias_map, id_to_category, reaction_alias_map):
    """Process the CSV file with standardized types."""
    print(f"Reading CSV from: {csv_path}")
    df = pd.read_csv(csv_path, encoding='utf-8-sig')
    
    print(f"Original shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")
    
    # Add new standardized columns
    df['Reactant_Type_Electrophile'] = ''
    df['Reactant_Type_Nucleophile'] = ''
    df['Reactant_Category_Electrophile'] = ''
    df['Reactant_Category_Nucleophile'] = ''
    df['Reaction_Type_Standardized'] = ''
    
    # Process each row
    stats = {
        'total_rows': len(df),
        'electrophile_mapped': 0,
        'nucleophile_mapped': 0,
        'reaction_mapped': 0,
        'electrophile_unmapped': defaultdict(int),
        'nucleophile_unmapped': defaultdict(int),
        'reaction_unmapped': defaultdict(int),
    }
    
    for idx, row in df.iterrows():
        # Process Aryl Halide (Electrophile)
        aryl_halide = row.get('Aryl Halide', '')
        if pd.notna(aryl_halide) and aryl_halide != '':
            parts = split_mixed_values(aryl_halide)
            mapped_ids = []
            mapped_categories = []
            
            for part in parts:
                key = part.strip()
                canonical = reactant_alias_map.get(key.lower())
                if canonical:
                    mapped_ids.append(canonical)
                    mapped_categories.append(id_to_category.get(canonical, 'Unknown'))
                    stats['electrophile_mapped'] += 1
                else:
                    stats['electrophile_unmapped'][part] += 1
            
            unique_ids = list(dict.fromkeys(mapped_ids))
            unique_categories = list(dict.fromkeys(mapped_categories))
            df.at[idx, 'Reactant_Type_Electrophile'] = ', '.join(unique_ids)
            df.at[idx, 'Reactant_Category_Electrophile'] = ', '.join(unique_categories)
        
        # Process N-Nucleophile/Boronate Type
        nucleophile = row.get('N-Nucleophile/Boronate Type', '')
        if pd.notna(nucleophile) and nucleophile != '':
            parts = split_mixed_values(nucleophile)
            mapped_ids = []
            mapped_categories = []
            
            for part in parts:
                # Handle special cases
                key = part.strip()
                if 'RCO2H or M' in key:
                    key = 'RCO2H or M'

                canonical = reactant_alias_map.get(key.lower())
                if canonical:
                    mapped_ids.append(canonical)
                    mapped_categories.append(id_to_category.get(canonical, 'Unknown'))
                    stats['nucleophile_mapped'] += 1
                else:
                    stats['nucleophile_unmapped'][part] += 1
            
            unique_ids = list(dict.fromkeys(mapped_ids))
            unique_categories = list(dict.fromkeys(mapped_categories))
            df.at[idx, 'Reactant_Type_Nucleophile'] = ', '.join(unique_ids)
            df.at[idx, 'Reactant_Category_Nucleophile'] = ', '.join(unique_categories)
        
        # Process Reaction Type
        reaction_type = row.get('Reaction Type', '')
        if pd.notna(reaction_type) and reaction_type != '':
            canonical_rt = reaction_alias_map.get(str(reaction_type).strip().lower())
            if canonical_rt:
                df.at[idx, 'Reaction_Type_Standardized'] = canonical_rt
                stats['reaction_mapped'] += 1
            else:
                stats['reaction_unmapped'][reaction_type] += 1
    
    return df, stats

def save_processed_csv(df, output_path):
    """Save the processed DataFrame."""
    df.to_csv(output_path, index=False, encoding='utf-8-sig')
    print(f"\nProcessed CSV saved to: {output_path}")

def print_statistics(stats):
    """Print processing statistics."""
    print("\n" + "=" * 80)
    print("PROCESSING STATISTICS")
    print("=" * 80)
    
    print(f"\nTotal rows processed: {stats['total_rows']:,}")
    print(f"Electrophile reactants mapped: {stats['electrophile_mapped']:,}")
    print(f"Nucleophile reactants mapped: {stats['nucleophile_mapped']:,}")
    print(f"Reactions mapped: {stats['reaction_mapped']:,}")
    
    if stats['electrophile_unmapped']:
        print(f"\nUnmapped electrophiles ({len(stats['electrophile_unmapped'])}):")
        for key, count in sorted(stats['electrophile_unmapped'].items(), key=lambda x: -x[1])[:20]:
            print(f"  - {key}: {count} occurrences")
    
    if stats['nucleophile_unmapped']:
        print(f"\nUnmapped nucleophiles ({len(stats['nucleophile_unmapped'])}):")
        for key, count in sorted(stats['nucleophile_unmapped'].items(), key=lambda x: -x[1])[:20]:
            print(f"  - {key}: {count} occurrences")
    
    if stats['reaction_unmapped']:
        print(f"\nUnmapped reactions ({len(stats['reaction_unmapped'])}):")
        for key, count in sorted(stats['reaction_unmapped'].items(), key=lambda x: -x[1]):
            print(f"  - {key}: {count} occurrences")

if __name__ == "__main__":
    # Paths
    base_path = Path(__file__).parent
    input_csv = base_path / "z-Score Peaks with FG.csv"
    output_csv = base_path / "z-Score Peaks with FG_STANDARDIZED.csv"
    
    # Build lookup tables from the shared ChemTools taxonomy
    print("Loading reactant and reaction type definitions...")
    reactant_alias_map, id_to_category = build_reactant_lookup()
    reaction_meta_index, reaction_alias_map = build_reaction_lookup()
    
    print(f"Reactant aliases available: {len(reactant_alias_map)}")
    print(f"Reaction aliases available: {len(reaction_alias_map)} (covering {len(reaction_meta_index)} reactions)")
    
    # Process CSV
    print("\nProcessing CSV...")
    df_processed, stats = process_csv(input_csv, reactant_alias_map, id_to_category, reaction_alias_map)
    
    # Save processed CSV
    save_processed_csv(df_processed, output_csv)
    
    # Print statistics
    print_statistics(stats)
    
    # Sample output
    print("\n" + "=" * 80)
    print("SAMPLE PROCESSED DATA (first 5 rows)")
    print("=" * 80)
    cols_to_show = [
        'Aryl Halide', 'Reactant_Type_Electrophile', 'Reactant_Category_Electrophile',
        'N-Nucleophile/Boronate Type', 'Reactant_Type_Nucleophile', 'Reactant_Category_Nucleophile',
        'Reaction Type', 'Reaction_Type_Standardized'
    ]
    print(df_processed[cols_to_show].head().to_string())
