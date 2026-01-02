"""
Add category-level SMARTS patterns to reactant_types.json

This enables hierarchical matching:
1. First check if molecule matches category (broad pattern)
2. Then find best specific member match within that category
"""

import json
from pathlib import Path

_REACTANT_TYPES_PATH = Path(__file__).resolve().parents[2] / "chemtools" / "taxonomy" / "data" / "reactant_types.json"
if not _REACTANT_TYPES_PATH.exists():
    raise SystemExit(
        "Legacy script: reactant_types.json has been removed; use organic_compounds.v1.3.json instead."
    )

from chemtools.reagent import (
    clear_reactant_type_cache,
    get_reactant_types_file,
)

# Category-level SMARTS patterns (broader than member patterns)
CATEGORY_SMARTS = {
    "ArX*": "c[Br,Cl,I,F,$([OX2][SX4](=O)(=O))]",  # Aryl + any halogen or sulfonate
    "Heterocyclic-halide": "[n,o,s]1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Br,Cl,I,F,$([OX2][SX4](=O)(=O))]",  # Heteroaryl + leaving group
    "VinylX*": "[CX3]=[CX3][Br,Cl,I,F,$([OX2][SX4](=O)(=O))]",  # Vinyl + leaving group
    "Alkyl-X": "[CX4][Br,Cl,I,F,$([OX2][SX4](=O)(=O))]",  # Alkyl + leaving group
    "Allylic-halide": "[CX4][CX3]=[CX3]",  # Allylic position (general)
    "Benzylic-halide": "[CX4]c",  # Benzylic position (general)
    "RSO2Cl": "[#6][SX4](=O)(=O)[Cl]",  # Sulfonyl chloride
    "Acyl-source-electrophile": "[#6][CX3](=O)[$([Cl]),$([OX2][CX3](=O))]",  # Acyl chloride or anhydride
    "ArB*": "c[B]",  # Aryl-boron
    "RB*": "[CX4,CX3][B]",  # Alkyl/alkenyl-boron
    "Organometallic": "[#6][Mg,Zn,Li,Al,Cu]",  # Any organometallic
    "Aliphatic-amine": "[CX4][NX3;H1,H2;!$(NC=O)]",  # Aliphatic amine (primary or secondary)
    "Aniline-type": "c[NX3;H1,H2;!$(NC=O)]",  # Aromatic amine
    "Enamine": "[NX3][#6]=[CX3]",  # Enamine structure
    "Imines": "[NX2]=[CX3]",  # Imine structure
    "ROH": "[CX4][OX2H]",  # Aliphatic alcohol
    "ArOH": "c[OX2H]",  # Phenol
    "RSH": "[#6][SX2H]",  # Thiol
    "Alkyl-C-H": "[CX4;H1,H2,H3]",  # Aliphatic C-H
    "ArH": "[c,n,o,s]:[c,n,o,s]",  # Aromatic C-H (aryl or heteroaryl)
    "Amide-type": "[#6,SX4][CX3,SX4](=O)[NX3]",  # Amide or sulfonamide
    "Acyl-source": "[#6][CX3](=O)[OX2,O-]",  # Carboxylic acid/ester/salt
    "Alkene": "[CX3]=[CX3]",  # Alkene
    "Azide": "[#6][NX2]=[NX2]=[NX1]",  # Azide
    "Nitrile": "[#6][CX2]#[NX1]",  # Nitrile
    "Alkyne": "[CX2]#[CX2]",  # Alkyne
    "Aldehyde": "[#6][CX3H1](=O)",  # Aldehyde
    "Ketone": "[#6][CX3](=O)[#6]",  # Ketone
}


def add_category_smarts(input_file, output_file):
    """Add category-level SMARTS patterns to reactant_types.json"""
    
    with open(input_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    categories_updated = 0
    categories_missing = []
    
    for category_name, category_data in data.items():
        if category_name in CATEGORY_SMARTS:
            # Add smarts field at category level (before members)
            category_data['smarts'] = CATEGORY_SMARTS[category_name]
            categories_updated += 1
            print(f"✓ Added SMARTS to {category_name}: {CATEGORY_SMARTS[category_name]}")
        else:
            categories_missing.append(category_name)
            print(f"⚠ Missing SMARTS for {category_name}")
    
    # Save with proper formatting
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"\n{'='*60}")
    print(f"Updated {categories_updated} categories")
    if categories_missing:
        print(f"Missing SMARTS for: {', '.join(categories_missing)}")
    print(f"Saved to: {output_file}")
    print(f"{'='*60}")


if __name__ == "__main__":
    target_path = str(get_reactant_types_file())
    add_category_smarts(target_path, target_path)
    clear_reactant_type_cache()
