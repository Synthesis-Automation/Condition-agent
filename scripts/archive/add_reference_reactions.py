"""
Add reference_reactions to all rule database files.

This script adds representative reaction SMILES to each rule database
for DRFP-based similarity matching.
"""

import json
from pathlib import Path

# Define reference reactions for each rule database
REFERENCE_REACTIONS = {
    "Suzuki_db.json": [
        # Aryl bromide + aryl boronic acid (standard)
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        # Aryl chloride + heteroaryl boronic acid (harder substrate)
        "Clc1ccc(C(F)(F)F)cc1.OB(O)c1ccncc1>>c1ccc(-c2ccncc2)cc1",
        # Aryl iodide + aryl boronic ester
        "Ic1ccc(C#N)cc1.c1ccc(B2OC(C)(C)C(C)(C)O2)cc1>>N#Cc1ccc(-c2ccccc2)cc1",
        # Heteroaryl bromide (pyridine)
        "Brc1cccnc1.OB(O)c1ccc(OC)cc1>>COc1ccc(-c2cccnc2)cc1",
        # Vinyl bromide
        "BrC=Cc1ccccc1.OB(O)c1ccc(C)cc1>>Cc1ccc(/C=C/c2ccccc2)cc1",
        # Aryl triflate
        "O=S(=O)(Oc1ccccc1)C(F)(F)F.OB(O)c1ccc(F)cc1>>Fc1ccc(-c2ccccc2)cc1",
        # Ortho-substituted aryl bromide
        "Brc1ccccc1C.OB(O)c1ccccc1>>Cc1ccccc1-c1ccccc1",
        # Electron-rich aryl chloride
        "Clc1ccc(OC)cc1.OB(O)c1ccc(C(C)(C)C)cc1>>COc1ccc(-c2ccc(C(C)(C)C)cc2)cc1"
    ],
    
    "C_N_Coupling_Pd_db.json": [
        # Aryl bromide + aniline (Buchwald-Hartwig)
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        # Aryl chloride + secondary amine
        "Clc1ccc(C(F)(F)F)cc1.CN(C)C>>CN(C)c1ccc(C(F)(F)F)cc1",
        # Heteroaryl bromide + primary amine
        "Brc1ccncc1.NCc1ccccc1>>c1ccc(CNc2ccncc2)cc1",
        # Aryl iodide + morpholine
        "Ic1ccc(C#N)cc1.C1COCCN1>>N#Cc1ccc(N2CCOCC2)cc1",
        # Aryl bromide + ammonia surrogate
        "Brc1ccccc1.CC(C)(C)OC(=O)N>>CC(C)(C)OC(=O)Nc1ccccc1",
        # Ortho-substituted aryl bromide + aniline
        "Brc1ccccc1C.Nc1ccc(OC)cc1>>COc1ccc(Nc2ccccc2C)cc1"
    ],
    
    "C_N_Coupling_Cu_db.json": [
        # Aryl iodide + aniline (Ullmann)
        "Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        # Aryl bromide + aliphatic amine
        "Brc1ccc(C(F)(F)F)cc1.NCCCC>>c1ccc(NCCCC)cc1",
        # Heteroaryl iodide + primary amine
        "Ic1ccncc1.NCc1ccccc1>>c1ccc(CNc2ccncc2)cc1",
        # Aryl iodide + secondary amine
        "Ic1ccc(C#N)cc1.CN(C)C>>CN(C)c1ccc(C#N)cc1",
        # Aryl bromide + heterocyclic amine
        "Brc1ccccc1.c1cnccn1>>n1ccnc(-c2ccccc2)c1"
    ],
    
    "C_O_coupling_db.json": [
        # Aryl bromide + phenol
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1",
        # Aryl chloride + aliphatic alcohol
        "Clc1ccc(C(F)(F)F)cc1.OCCCC>>c1ccc(OCCCC)cc1",
        # Heteroaryl bromide + phenol
        "Brc1ccncc1.Oc1ccc(OC)cc1>>COc1ccc(Oc2ccncc2)cc1",
        # Aryl iodide + benzyl alcohol
        "Ic1ccc(C#N)cc1.OCc1ccccc1>>N#Cc1ccc(OCc2ccccc2)cc1",
        # Aryl triflate + methanol
        "O=S(=O)(Oc1ccccc1)C(F)(F)F.CO>>COc1ccccc1"
    ],
    
    "RCM_db.json": [
        # Simple diene (5-membered ring)
        "C=CCCC=C>>C1=CCCC1",
        # Substituted diene (6-membered ring)
        "C=CCCCCC=C>>C1=CCCCC1",
        # Diene with functional groups
        "C=CCC(=O)CC=C>>C1=CCC(=O)C1",
        # Styrene diene
        "C=CCc1ccc(CC=C)cc1>>C=CCc1ccc(CC1)cc1",
        # Trisubstituted alkene formation
        "C=C(C)CCC=C>>CC1=CCCC1",
        # Large ring (macrocycle)
        "C=CCCCCCCCC=C>>C1=CCCCCCCC1"
    ],
    
    "sonogashira_db.json": [
        # Already done - skip
    ],
    
    "SNAr_db.json": [
        # Aryl fluoride + aniline
        "Fc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>O=[N+]([O-])c1ccc(Nc2ccccc2)cc1",
        # Heteroaryl chloride + amine
        "Clc1ccnc(Cl)c1.NCCCC>>c1cc(NCCCC)nc(Cl)c1",
        # Activated aryl chloride + phenoxide
        "Clc1cc([N+](=O)[O-])cc([N+](=O)[O-])c1.Oc1ccccc1>>O=[N+]([O-])c1cc(Oc2ccccc2)cc([N+](=O)[O-])c1",
        # Pyridine N-oxide SNAr
        "Clc1cccc[n+]1[O-].Nc1ccccc1>>c1ccc(Nc2cccc[n+]2[O-])cc1",
        # Dinitrochlorobenzene + thiol
        "Clc1c([N+](=O)[O-])cccc1[N+](=O)[O-].SCc1ccccc1>>O=[N+]([O-])c1cccc([N+](=O)[O-])c1SCc1ccccc1"
    ],
    
    "amide_formation_db.json": [
        # Carboxylic acid + amine (EDC coupling)
        "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1",
        # Acid chloride + amine
        "O=C(Cl)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1",
        # Activated ester + amine
        "O=C(On1nnc2ccccc21)c1ccccc1.NCCCC>>O=C(NCCCC)c1ccccc1",
        # Carboxylic acid + aniline (HATU)
        "O=C(O)c1ccc(OC)cc1.Nc1ccc(F)cc1>>COc1ccc(C(=O)Nc2ccc(F)cc2)cc1",
        # Aliphatic acid + secondary amine
        "O=C(O)CCCC.CN(C)C>>O=C(CCCC)N(C)C"
    ],
    
    "reductive_amination_db.json": [
        # Benzaldehyde + aniline
        "O=Cc1ccccc1.Nc1ccccc1>>c1ccc(NCc2ccccc2)cc1",
        # Aliphatic aldehyde + primary amine
        "O=CCCCC.NCc1ccccc1>>c1ccc(CNCCCCC)cc1",
        # Ketone + secondary amine
        "O=C(C)c1ccccc1.CN(C)C>>CN(C)C(C)c1ccccc1",
        # Cyclic ketone + aniline
        "O=C1CCCCC1.Nc1ccccc1>>c1ccc(NC2CCCCC2)cc1",
        # Aromatic aldehyde + heterocyclic amine
        "O=Cc1ccc(OC)cc1.C1CCNCC1>>COc1ccc(CN2CCCCC2)cc1"
    ]
}


def add_reference_reactions_to_file(file_path: Path, reference_reactions: list):
    """Add reference_reactions to a rule database file."""
    print(f"\nProcessing {file_path.name}...")
    
    # Read existing file
    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Check if reference_reactions already exists
    if 'reaction' not in data:
        data['reaction'] = {}
    
    if 'reference_reactions' in data['reaction']:
        print(f"  ⚠️  reference_reactions already exists - skipping")
        return False
    
    # Add reference_reactions
    data['reaction']['reference_reactions'] = reference_reactions
    
    # Save back
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"  ✅ Added {len(reference_reactions)} reference reactions")
    return True


def main():
    rule_db_dir = Path("data/rule_db")
    
    print("="*80)
    print("Adding reference_reactions to rule databases")
    print("="*80)
    
    updated_count = 0
    skipped_count = 0
    
    for filename, reactions in REFERENCE_REACTIONS.items():
        if filename == "sonogashira_db.json":
            print(f"\nSkipping {filename} (already has v2 version)")
            skipped_count += 1
            continue
        
        file_path = rule_db_dir / filename
        
        if not file_path.exists():
            print(f"\n⚠️  {filename} not found - skipping")
            continue
        
        if reactions and add_reference_reactions_to_file(file_path, reactions):
            updated_count += 1
        else:
            skipped_count += 1
    
    print("\n" + "="*80)
    print(f"Summary: {updated_count} updated, {skipped_count} skipped")
    print("="*80)


if __name__ == '__main__':
    main()
