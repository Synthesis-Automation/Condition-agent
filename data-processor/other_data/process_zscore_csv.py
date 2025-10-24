"""Process z-Score Peaks CSV to use standardized reactant and reaction types."""

import pandas as pd
import json
from pathlib import Path
from collections import defaultdict

def load_type_definitions():
    """Load reactant and reaction type definitions."""
    base_path = Path(__file__).parent
    
    with open(base_path / "reactant_types.json", 'r', encoding='utf-8') as f:
        reactant_types = json.load(f)
    
    with open(base_path / "reaction_types.json", 'r', encoding='utf-8') as f:
        reaction_types = json.load(f)
    
    return reactant_types, reaction_types

def create_reactant_mapping(reactant_types):
    """Create mapping from CSV values to standardized reactant type IDs."""
    # Build all possible IDs and names
    all_ids = set()
    id_to_category = {}
    
    for category_key, category_data in reactant_types.items():
        all_ids.add(category_key)
        id_to_category[category_key] = category_key
        
        for member in category_data.get("members", []):
            member_id = member["id"]
            all_ids.add(member_id)
            id_to_category[member_id] = category_key
    
    # Manual mapping for CSV-specific values
    csv_to_standard = {
        # From z-Score CSV "Aryl Halide" column
        "ArBr": "ArBr",
        "ArCl": "ArCl",
        "ArI": "ArI",
        "ArF": "ArF",
        "ArOH": "ArOH",
        "ArOSO2R": "ArOSO2R",
        "Alkyl-Br": "Alkyl-Br",
        "Alkyl-Cl": "Alkyl-Cl",
        "Alkyl-I": "Alkyl-I",
        "Alkyl-OSO2R": "Alkyl-OSO2R",
        "alkene-Br": "alkene-Br",
        "alkene-I": "alkene-I",
        
        # From "N-Nucleophile/Boronate Type" column
        "ArB(OR)2": "ArB(OR)2",
        "ArB(OH)2": "ArB(OH)2",
        "ArBF3K": "ArBF3K",
        "RNH2 a-branch": "RNH2-a-branch",
        "RNH2": "RNH2",
        "R2NH a-branch": "R2NH-a-branch",
        "R2NH": "R2NH",
        "ArNH2": "ArNH2",
        "ArNHR": "ArNHR",
        "Ar2NH": "Ar2NH",
        "arom. NH": "arom-NH",
        "alkeneB(OR)2": "alkeneB(OR)2",
        "Alkyl-B(OH)2": "Alkyl-B(OH)2",
        "Alkyl-B(OR)2": "Alkyl-B(OR)2",
        "Alkyl-BF3K": "Alkyl-BF3K",
        "Alkyl-OH a-branch": "ROH-a-branch",
        "Alkyl-OH primary": "ROH-primary",
        "Lactam": "Lactam",
        "RCONH2": "RCONH2",
        "ROC(O)NR2": "Carbamate",
        "Urea": "Urea",
        "Alkyl-M": "Alkyl-M",
        "Alkyl-H acidic": "Alkyl-H-acidic",
        "Alkyl-H": "Alkyl-H",
        "alkene": "alkene",
        "alkyne": "alkyne",
        "Ar-H": "Ar-H",
        "enolether": "enol-ether",
        "RSH": "RSH",
        
        # Mixed entries (split these)
        "RCO2H or M": "RCO2H",  # Carboxylic acid
    }
    
    return csv_to_standard, id_to_category

def create_reaction_mapping(reaction_types):
    """Create mapping from CSV reaction names to standardized reaction type IDs."""
    # Build reverse mapping from aliases and names to IDs
    csv_to_standard = {}
    
    for category_data in reaction_types.values():
        for reaction in category_data["reactions"]:
            reaction_id = reaction["id"]
            
            # Map from reaction name
            csv_to_standard[reaction["name"]] = reaction_id
            
            # Map from aliases
            for alias in reaction.get("aliases", []):
                csv_to_standard[alias] = reaction_id
            
            # Map from ID itself
            csv_to_standard[reaction_id] = reaction_id
    
    # Manual mappings for CSV-specific names
    csv_specific = {
        "Buchwald-Hartwig": "Buchwald-Hartwig",
        "Suzuki-Miyaura": "Suzuki-Miyaura",
        "Suzuki-Miyaura, in situ": "Suzuki-Miyaura-in-situ",
        "Arylation, acidic C-H": "Arylation-acidic-C-H",
        "Amide coupling": "Amide-coupling",
        "CN-Coupling": "CN-Coupling",
        "CO-Coupling": "CO-Coupling",
        "Condensation": "Condensation",
        "CH-Activation": "CH-Activation",
        "Negishi, in-situ": "Negishi-in-situ",
        "Cyclization": "Cyclization",
        "Borylation, Miyaura": "Borylation-Miyaura",
        "Alkylation": "Alkylation",
        "Deprotection": "Deprotection",
        "Negishi": "Negishi",
        "Heck": "Heck",
        "CC-Coupling": "CC-Coupling",
        "SNAr": "SNAr",
        "Hydrolysis": "Hydrolysis",
        "Salt formation": "Salt-formation",
        "Sonogashira": "Sonogashira",
        "Stetter": "Stetter",
        "Cyanation": "Cyanation",
        "Oxidation": "Oxidation",
        "Activation": "Activation",
        "Hydrodehalogenation": "Hydrodehalogenation",
        "Glycosidation": "Glycosidation",
        "Stannylation": "Stannylation",
        "Dehydration": "Dehydration",
        "Dimerization, reductive": "Dimerization-reductive",
        "Mitsunobu": "Mitsunobu",
        "CS-Coupling": "CS-Coupling",
        "Borylation, C-H": "Borylation-C-H",
        "Wittig": "Wittig",
        "Sandmeyer": "Sandmeyer",
        "Addition": "Addition",
        "Hydration": "Hydration",
        "Reduction": "Reduction",
        "Deoxyfluorination": "Deoxyfluorination",
        "Chlorination": "Halogenation-Chlorination",
        "Fluorination, oxidative": "Fluorination-oxidative",
        "Protection": "Protection",
    }
    
    csv_to_standard.update(csv_specific)
    
    return csv_to_standard

def split_mixed_values(value):
    """Split mixed reactant values (e.g., 'ArBr, ArCl' -> ['ArBr', 'ArCl'])."""
    if pd.isna(value) or value == '':
        return []
    
    # Split by comma
    parts = [p.strip() for p in str(value).split(',')]
    return parts

def process_csv(csv_path, reactant_mapping, reaction_mapping, id_to_category):
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
                if part in reactant_mapping:
                    std_id = reactant_mapping[part]
                    mapped_ids.append(std_id)
                    mapped_categories.append(id_to_category.get(std_id, 'Unknown'))
                    stats['electrophile_mapped'] += 1
                else:
                    stats['electrophile_unmapped'][part] += 1
            
            df.at[idx, 'Reactant_Type_Electrophile'] = ', '.join(mapped_ids)
            df.at[idx, 'Reactant_Category_Electrophile'] = ', '.join(set(mapped_categories))
        
        # Process N-Nucleophile/Boronate Type
        nucleophile = row.get('N-Nucleophile/Boronate Type', '')
        if pd.notna(nucleophile) and nucleophile != '':
            parts = split_mixed_values(nucleophile)
            mapped_ids = []
            mapped_categories = []
            
            for part in parts:
                # Handle special cases
                if 'RCO2H or M' in part:
                    part = 'RCO2H or M'
                
                if part in reactant_mapping:
                    std_id = reactant_mapping[part]
                    mapped_ids.append(std_id)
                    mapped_categories.append(id_to_category.get(std_id, 'Unknown'))
                    stats['nucleophile_mapped'] += 1
                else:
                    stats['nucleophile_unmapped'][part] += 1
            
            df.at[idx, 'Reactant_Type_Nucleophile'] = ', '.join(mapped_ids)
            df.at[idx, 'Reactant_Category_Nucleophile'] = ', '.join(set(mapped_categories))
        
        # Process Reaction Type
        reaction_type = row.get('Reaction Type', '')
        if pd.notna(reaction_type) and reaction_type != '':
            if reaction_type in reaction_mapping:
                df.at[idx, 'Reaction_Type_Standardized'] = reaction_mapping[reaction_type]
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
    
    # Load type definitions
    print("Loading reactant and reaction type definitions...")
    reactant_types, reaction_types = load_type_definitions()
    
    # Create mappings
    print("Creating mappings...")
    reactant_mapping, id_to_category = create_reactant_mapping(reactant_types)
    reaction_mapping = create_reaction_mapping(reaction_types)
    
    print(f"Reactant mappings created: {len(reactant_mapping)}")
    print(f"Reaction mappings created: {len(reaction_mapping)}")
    
    # Process CSV
    print("\nProcessing CSV...")
    df_processed, stats = process_csv(input_csv, reactant_mapping, reaction_mapping, id_to_category)
    
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
