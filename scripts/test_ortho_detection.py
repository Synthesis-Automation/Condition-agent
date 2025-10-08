"""
Test ortho substitution detection for SCDB-SUZ-ARBR-ORTHO-XPhos matching.

This script demonstrates how the SCDB matcher detects ortho-substituted
aryl halides and routes them to the appropriate XPhos conditions.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rdkit import Chem
from chemtools.scdb_matcher.loader import load_db
from chemtools.scdb_matcher.matcher import match

def analyze_ortho_substitution(smiles: str, name: str):
    """
    Analyze a reaction SMILES for ortho substitution detection.
    """
    print(f"\n{'='*80}")
    print(f"Testing: {name}")
    print(f"{'='*80}")
    print(f"SMILES: {smiles}")
    print()
    
    # Load the Suzuki database
    suzuki_db = load_db("data/conditionDB/suzuki_db.json")
    
    # Match the reaction
    result = match(suzuki_db, smiles)
    
    if result:
        print(f"âœ?Match Found!")
        print(f"  Entry ID: {result.entry_id}")
        print(f"  Entry Name: {result.entry_name or 'N/A'}")
        print(f"  Match Type: {result.match_type}")
        print(f"  Priority: {result.priority}")
        print()
        
        # Display detected features from trace
        if 'observed_features' in result.trace:
            print("Detected Features:")
            features = result.trace['observed_features']
            if 'numeric_features' in features:
                numeric = features['numeric_features']
                if 'electrophile.ortho_sub_count' in numeric:
                    ortho_count = numeric['electrophile.ortho_sub_count']
                    print(f"  - Ortho substitution count: {ortho_count}")
                if 'electrophile.ring_hetero_count' in numeric:
                    print(f"  - Ring heteroatom count: {numeric['electrophile.ring_hetero_count']}")
            if 'set_features' in features:
                sets = features['set_features']
                if 'electrophile.lg_class' in sets:
                    print(f"  - Leaving group(s): {', '.join(sets['electrophile.lg_class'])}")
        
        # Display recommended conditions
        print()
        print("Recommended Conditions:")
        conds = result.conditions
        if 'pd_source' in conds:
            pd_sources = conds['pd_source']
            if isinstance(pd_sources, list):
                print(f"  - Pd Source: {', '.join(pd_sources[:2])}")
            else:
                print(f"  - Pd Source: {pd_sources}")
        if 'ligands' in conds or 'ligand' in conds:
            ligands = conds.get('ligands', conds.get('ligand'))
            if isinstance(ligands, list):
                print(f"  - Ligand: {', '.join(ligands[:2])}")
            else:
                print(f"  - Ligand: {ligands}")
        if 'base' in conds:
            bases = conds['base']
            if isinstance(bases, list):
                print(f"  - Base: {', '.join(bases[:2])}")
            else:
                print(f"  - Base: {bases}")
        if 'temperature_C' in conds:
            temps = conds['temperature_C']
            if isinstance(temps, list):
                print(f"  - Temperature: {temps[0]}-{temps[1]}Â°C")
            elif isinstance(temps, (int, float)):
                print(f"  - Temperature: {temps}Â°C")
    else:
        print("âœ?No match found")

def main():
    """
    Test various aryl halides with different ortho substitution patterns.
    """
    
    print("="*80)
    print(" ORTHO SUBSTITUTION DETECTION TEST")
    print(" Testing SCDB-SUZ-ARBR-ORTHO-XPhos Entry")
    print("="*80)
    
    test_cases = [
        {
            "smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
            "name": "Simple bromobenzene (NO ortho substitution)"
        },
        {
            "smiles": "Brc1ccccc1C.c1ccc(B(O)O)cc1>>Cc1ccccc1-c1ccccc1",
            "name": "2-Bromotoluene (1 ortho substituent)"
        },
        {
            "smiles": "Brc1cc(C)ccc1C.c1ccc(B(O)O)cc1>>Cc1ccc(C)c(-c2ccccc2)c1",
            "name": "2-Bromo-4,5-dimethylbenzene (2 ortho substituents)"
        },
        {
            "smiles": "Brc1c(C)cccc1C.COc1ccc(B(O)O)cc1>>Cc1cccc(C)c1-c1ccc(OC)cc1",
            "name": "2,6-Dimethylbromobenzene (2 ortho substituents)"
        },
        {
            "smiles": "Brc1ccccc1C(C)C.c1ccc(B(O)O)cc1>>CC(C)c1ccccc1-c1ccccc1",
            "name": "2-Isopropylbromobenzene (1 bulky ortho substituent)"
        },
        {
            "smiles": "Brc1ccccc1OCC.c1ccc(B(O)O)cc1C>>CCOc1ccccc1-c1ccccc1C",
            "name": "2-Ethoxybromobenzene (1 ortho + nucleophile has 1 ortho)"
        },
    ]
    
    for test in test_cases:
        analyze_ortho_substitution(test["smiles"], test["name"])
    
    print()
    print("="*80)
    print(" SUMMARY")
    print("="*80)
    print()
    print("The SCDB matcher detects ortho substitution by:")
    print()
    print("1. Identifying the aromatic carbon attached to the leaving group (ipso position)")
    print("2. Finding the 6-membered aromatic ring containing this carbon")
    print("3. Checking the ortho positions (adjacent carbons in the ring)")
    print("4. Counting non-hydrogen substituents at ortho positions")
    print()
    print("Feature: 'electrophile.ortho_sub_count'")
    print("  - Value 0: No ortho substituents (uses standard conditions)")
    print("  - Value 1: One ortho substituent (may need XPhos)")
    print("  - Value 2+: Two ortho substituents (definitely needs XPhos)")
    print()
    print("SCDB-SUZ-ARBR-ORTHO-XPhos requirements:")
    print("  - electrophile.lg_class: Br")
    print("  - electrophile.ortho_sub_count: >= 1")
    print()
    print("When matched, uses bulky biaryl phosphine ligands (XPhos, Ligand 95)")
    print("to accommodate steric hindrance at higher temperatures (75-90Â°C).")
    print()
    print("="*80)

if __name__ == "__main__":
    main()
