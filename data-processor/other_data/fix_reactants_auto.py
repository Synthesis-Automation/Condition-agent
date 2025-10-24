"""Automatically fix reactants field in reaction_types.json to reference reactant_types.json properly."""

import json
from pathlib import Path
import re

def load_json_files():
    """Load both JSON files."""
    base_path = Path(__file__).parent
    
    with open(base_path / "reaction_types.json", 'r', encoding='utf-8') as f:
        reaction_types = json.load(f)
    
    with open(base_path / "reactant_types.json", 'r', encoding='utf-8') as f:
        reactant_types = json.load(f)
    
    return reaction_types, reactant_types

def get_reactant_mapping():
    """Define mappings from free text to reactant type IDs."""
    return {
        # Exact matches
        "ArX": "ArX*",
        "Ar-X": "ArX*",
        "aryl halide": "ArX*",
        "aryl halides": "ArX*",
        
        "R-X": "Alkyl-X",
        "RX": "Alkyl-X",
        "alkyl halide": "Alkyl-X",
        
        "vinyl halide": "VinylX*",
        
        # Amines
        "R2NH or RNH2": ["Aliphatic-amine"],
        "RNH2": ["Aliphatic-amine"],
        "R2NH": ["Aliphatic-amine"],
        "amine": ["Aliphatic-amine"],
        "aniline": ["Aniline-type"],
        "ArNH2": ["Aniline-type"],
        "N-nucleophile": ["Aliphatic-amine", "Aniline-type"],
        
        # Alcohols
        "ROH or ArOH": ["ROH", "ArOH"],
        "ROH": "ROH",
        "ArOH": "ArOH",
        "alcohol": "ROH",
        "phenol": "ArOH",
        
        # Boronic acids
        "ArB(OH)2 or ArB(OR)2": "ArB*",
        "ArB(OH)2": "ArB*",
        "ArB(OR)2": "ArB*",
        "boronic acid": "ArB*",
        
        # Organometallics
        "RZnX": "Organometallic",
        "RMgX (Grignard reagent)": "Organometallic",
        "Grignard": "Organometallic",
        "R-SnR3 (organostannane)": "Organometallic",
        
        # Carbonyls
        "Carbonyl compound": ["Aldehyde", "Ketone"],
        "Carbonyl compounds": ["Aldehyde", "Ketone"],
        "carbonyl": ["Aldehyde", "Ketone"],
        "R2C=O": ["Aldehyde", "Ketone"],
        
        # Unsaturated
        "alkene": "Alkene",
        "terminal alkyne": "Alkyne",
        "alkyne": "Alkyne",
        
        # Acids/esters
        "Carboxylic acid": "Acyl-source",
        "carboxylic acid": "Acyl-source",
        "Acyl source": "Acyl-source",
        "α-halo ester": "Acyl-source",
        "Phosphonate ester": "Acyl-source",
        "ester": "Acyl-source",
        
        # Thiols
        "RSH": "RSH",
        "thiol": "RSH",
        "metal thiolate": "RSH",
        "SR": "RSH",
        
        # C-H donors
        "C-H substrate": ["Alkyl-C-H", "ArH"],
        "Acidic C-H substrate": "Alkyl-C-H",
        
        # General/unclear
        "Nucleophile": ["Aliphatic-amine", "Aniline-type", "ROH", "ArOH", "RSH"],
        "Nucleophile (NH2, OH, SR, etc.)": ["Aliphatic-amine", "Aniline-type", "ROH", "ArOH", "RSH"],
        "Nucleophile or electrophile": ["Aliphatic-amine", "ArX*", "Alkyl-X"],
        "Coupling partner": ["ArX*", "Organometallic"],
        
        # Substrates
        "Substrate": [],  # Too generic, skip
        "Substrate with functional group": [],
        "Protected substrate": [],
        "Bifunctional substrate": [],
        "Unsaturated substrate": ["Alkene", "Alkyne"],
        "Substrate with OH or two OH groups": "ROH",
        
        # Reagents
        "B2pin2 or similar": "RB*",
        "B2pin2 or HBpin": "RB*",
        "Zn": "Organometallic",
        "Oxidant": [],  # Not a reactant type
        "H2O": [],  # Not a reactant type
        "H source (H2, formate, etc.)": [],
        "CN source (NaCN, KCN, Zn(CN)2, etc.)": [],
        "Chlorinating agent": [],
        "Electrophilic F source (Selectfluor, NFSI, etc.)": [],
        "Fluorinating reagent (DAST, Deoxo-Fluor, etc.)": [],
        "NaNO2/HX": [],
        "CuX": [],
        "Protecting reagent": [],
        "Deprotecting reagent": [],
        "Activating reagent": [],
        "Ph3P=CR2": [],
        "α,β-unsaturated carbonyl": ["Aldehyde", "Ketone"],
        "Glycosyl donor": [],
        "Acid": [],
        "Base": [],
        "Various": [],
    }

def normalize_reactants(reaction_types, mapping):
    """Normalize all reactants fields."""
    changes = []
    
    for category_key, category_data in reaction_types.items():
        for reaction in category_data["reactions"]:
            old_reactants = reaction["reactants"].copy()
            new_reactants = []
            
            for reactant_str in old_reactants:
                # Try exact match first
                if reactant_str in mapping:
                    mapped = mapping[reactant_str]
                    if isinstance(mapped, list):
                        new_reactants.extend(mapped)
                    elif mapped:  # Not empty string
                        new_reactants.append(mapped)
                    else:
                        # Empty mapping means skip (reagent, not reactant)
                        pass
                else:
                    # Try to find partial matches
                    found = False
                    for key, value in mapping.items():
                        if key.lower() in reactant_str.lower() or reactant_str.lower() in key.lower():
                            if isinstance(value, list):
                                new_reactants.extend(value)
                            elif value:
                                new_reactants.append(value)
                            found = True
                            break
                    
                    if not found:
                        # Keep original if no mapping found
                        new_reactants.append(reactant_str)
            
            # Remove duplicates while preserving order
            seen = set()
            unique_reactants = []
            for r in new_reactants:
                if r not in seen:
                    seen.add(r)
                    unique_reactants.append(r)
            
            if unique_reactants != old_reactants:
                changes.append({
                    "reaction": reaction["id"],
                    "old": old_reactants,
                    "new": unique_reactants
                })
                reaction["reactants"] = unique_reactants
    
    return changes

def save_reaction_types(reaction_types):
    """Save the updated reaction types JSON."""
    output_path = Path(__file__).parent / "reaction_types_UPDATED.json"
    
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(reaction_types, f, indent=2, ensure_ascii=False)
    
    return output_path

if __name__ == "__main__":
    print("Loading JSON files...")
    reaction_types, reactant_types = load_json_files()
    
    print("Getting reactant mapping...")
    mapping = get_reactant_mapping()
    
    print("Normalizing reactants fields...")
    changes = normalize_reactants(reaction_types, mapping)
    
    print(f"\nMade {len(changes)} changes:")
    print("=" * 80)
    
    for change in changes:
        print(f"\nReaction: {change['reaction']}")
        print(f"  Old: {change['old']}")
        print(f"  New: {change['new']}")
    
    output_path = save_reaction_types(reaction_types)
    print(f"\n\nUpdated reaction_types saved to: {output_path}")
    print("\nPlease review the changes, then rename to reaction_types.json if satisfied.")
