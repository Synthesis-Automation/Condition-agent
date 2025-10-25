"""
Add SMARTS patterns to reactant_types.json for automatic detection.

This script reads the current reactant_types.json and adds SMARTS patterns
to each member type for automatic substrate classification.
"""

import json

# Load existing reactant types
with open('data-processor/other_data/reactant_types.json', 'r', encoding='utf-8') as f:
    reactant_types = json.load(f)

# SMARTS patterns for each member type
SMARTS_PATTERNS = {
    # ArX* members
    "ArBr": "c[Br]",
    "ArCl": "c[Cl]",
    "ArI": "c[I]",
    "ArF": "c[F]",
    "ArOSO2R": "c[OX2][SX4](=O)(=O)[#6]",
    "ArOMs": "c[OX2][SX4](=O)(=O)[CH3]",
    "ArOTf": "c[OX2][SX4](=O)(=O)[CF3]",
    "ArOTs": "c[OX2][SX4](=O)(=O)c1ccc([CH3])cc1",
    
    # Heterocyclic-halide members
    "HetArBr": "[n,o,s]1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Br]",
    "HetArCl": "[n,o,s]1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[Cl]",
    "HetArI": "[n,o,s]1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[I]",
    "HetArOTf": "[n,o,s]1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[OX2][SX4](=O)(=O)[CF3]",
    "PyridineBr": "n1ccccc1[Br]",
    "PyrimidineCl": "n1cnccn1[Cl]",
    "ThiazoleBr": "s1ccnc1[Br]",
    "IndoleBr": "c1ccc2[nH]ccc2c1[Br]",
    
    # VinylX* members
    "alkene-Br": "[CX3]=[CX3][Br]",
    "alkene-I": "[CX3]=[CX3][I]",
    "alkene-Cl": "[CX3]=[CX3][Cl]",
    "VinylOTf": "[CX3]=[CX3][OX2][SX4](=O)(=O)[CF3]",
    
    # Alkyl-X members
    "Alkyl-Br": "[CX4][Br]",
    "Alkyl-Cl": "[CX4][Cl]",
    "Alkyl-I": "[CX4][I]",
    "Alkyl-OSO2R": "[CX4][OX2][SX4](=O)(=O)[#6]",
    
    # Allylic-halide members  
    "Allyl-Br": "[CX4;H1,H2,H3][CX3]=[CX3].[Br]",
    "Allyl-Cl": "[CX4;H1,H2,H3][CX3]=[CX3].[Cl]",
    "Allyl-I": "[CX4;H1,H2,H3][CX3]=[CX3].[I]",
    "Allyl-OAc": "[CX4;H1,H2,H3][CX3]=[CX3].[OX2][CX3](=O)[CH3]",
    "Allyl-OCO2R": "[CX4;H1,H2,H3][CX3]=[CX3].[OX2][CX3](=O)[OX2][#6]",
    
    # Benzylic-halide members
    "Bn-Br": "[CX4;H1,H2,H3]c.[Br]",
    "Bn-Cl": "[CX4;H1,H2,H3]c.[Cl]",
    "Bn-I": "[CX4;H1,H2,H3]c.[I]",
    
    # RSO2Cl members
    "RSO2Cl": "[#6][SX4](=O)(=O)[Cl]",
    "MsCl": "[CH3][SX4](=O)(=O)[Cl]",
    
    # Acyl-source-electrophile members
    "RCOCl": "[#6][CX3](=O)[Cl]",
    "PivCl": "[CX4]([CH3])([CH3])([CH3])[CX3](=O)[Cl]",
    "Anhydride": "[#6][CX3](=O)[OX2][CX3](=O)[#6]",
    
    # ArB* members
    "ArB(OH)2": "c[B]([OH])[OH]",
    "ArB(OR)2": "c[B]([OX2][#6])[OX2][#6]",
    "ArBF3K": "c[B-]([F])([F])[F]",
    
    # RB* members
    "Alkyl-B(OH)2": "[CX4][B]([OH])[OH]",
    "Alkyl-B(OR)2": "[CX4][B]([OX2][#6])[OX2][#6]",
    "Alkyl-BF3K": "[CX4][B-]([F])([F])[F]",
    "alkeneB(OR)2": "[CX3]=[CX3][B]([OX2][#6])[OX2][#6]",
    "alkeneB(OH)2": "[CX3]=[CX3][B]([OH])[OH]",
    
    # Organometallic members
    "Alkyl-M": "[CX4][Mg,Zn,Li]",
    "Ar-M": "c[Mg,Zn,Li]",
    
    # Aliphatic-amine members
    "RNH2": "[CX4][NX3;H2;!$(NC=O)]",
    "RNH2-a-branch": "[CX4;H1]([#6])([#6])[NX3;H2;!$(NC=O)]",
    "R2NH": "[CX4][NX3;H1;!$(NC=O)][CX4]",
    "R2NH-a-branch": "[CX4;H1]([#6])([#6])[NX3;H1;!$(NC=O)]",
    "R3N": "[NX3;H0;!$(NC=O);!$(N=*)]([CX4])([CX4])[CX4]",
    
    # Aniline-type members
    "ArNH2": "c[NX3;H2;!$(NC=O)]",
    "ArNHR": "c[NX3;H1;!$(NC=O)][CX4]",
    "Ar2NH": "c[NX3;H1;!$(NC=O)]c",
    "arom-NH": "[n,o,s]1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[NH2,NH1]",
    
    # Enamine members
    "enamine": "[NX3]([#6])[#6]=[CX3]",
    "morpholine-enamine": "O1CCN([#6]=[CX3])CC1",
    "pyrrolidine-enamine": "C1CCN([#6]=[CX3])C1",
    
    # Imines members
    "imine": "[NX2]=[CX3]",
    "Ar-imine": "c[CX3]=[NX2]",
    "N-protected-imine": "[NX2]([SX4](=O)(=O)[#6])=[CX3]",
    
    # ROH members
    "ROH-primary": "[CX4;H2][OX2H]",
    "ROH-secondary": "[CX4;H1]([#6])[OX2H]",
    "ROH-tertiary": "[CX4;H0]([#6])([#6])[OX2H]",
    "ROH-a-branch": "[CX4;H1]([#6])([#6])[OX2H]",
    
    # ArOH members
    "ArOH": "c[OX2H]",
    
    # RSH members
    "RSH": "[#6][SX2H]",
    "ArSH": "c[SX2H]",
    
    # Azide members
    "R-N3": "[CX4][NX2]=[NX2]=[NX1]",
    "Ar-N3": "c[NX2]=[NX2]=[NX1]",
    "NaN3": "[Na+].[N-]=[N+]=[N-]",
    
    # Alkyl-C-H members
    "Alkyl-H": "[CX4;H1,H2,H3]",
    "Alkyl-H-acidic": "[CX4;H1,H2,H3][CX3]=O",
    
    # ArH members
    "Ar-H": "c[H]",
    "Hetero-ArH": "[n,o,s]1:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:[c,n,o,s]:1[H]",
    
    # Amide-type members
    "RCONH2": "[#6][CX3](=O)[NX3;H2]",
    "RCONHR": "[#6][CX3](=O)[NX3;H1][#6]",
    "Lactam": "[CX3r](=O)[NX3r]",
    "Carbamate": "[OX2][CX3](=O)[NX3]",
    "Sulfonyl-amide": "[SX4](=O)(=O)[NX3]",
    "Urea": "[NX3][CX3](=O)[NX3]",
    
    # Acyl-source members
    "RCO2H": "[#6][CX3](=O)[OX2H]",
    "RCO2M": "[#6][CX3](=O)[O-].[Na+,K+,Li+]",
    "RCO2OR": "[#6][CX3](=O)[OX2][#6]",
    
    # Alkene members
    "alkene": "[CX3]=[CX3]",
    "R-alkene": "[CX4][CX3]=[CX3][CX4]",
    "Ar-alkene": "c[CX3]=[CX3]",
    "enol-ether": "[OX2][#6]=[CX3]",
    
    # Nitrile members
    "R-CN": "[CX4][CX2]#[NX1]",
    "Ar-CN": "c[CX2]#[NX1]",
    
    # Alkyne members
    "alkyne": "[CX2]#[CX2]",
    "R-alkyne": "[CX4][CX2]#[CX2]",
    "Ar-alkyne": "c[CX2]#[CX2]",
    
    # Aldehyde members
    "RCHO": "[CX4][CX3H1](=O)",
    "ArCHO": "c[CX3H1](=O)",
    
    # Ketone members
    "RCOR": "[#6][CX3](=O)[#6]",
    "ArCOR": "c[CX3](=O)[#6]",
}

# Add SMARTS patterns to each member
updated_count = 0
missing_patterns = []

for category, data in reactant_types.items():
    for member in data['members']:
        member_id = member['id']
        if member_id in SMARTS_PATTERNS:
            member['smarts'] = SMARTS_PATTERNS[member_id]
            updated_count += 1
        else:
            missing_patterns.append(f"{category}/{member_id}")

# Save updated reactant_types.json
with open('data-processor/other_data/reactant_types.json', 'w', encoding='utf-8') as f:
    json.dump(reactant_types, f, indent=2, ensure_ascii=False)

print(f"✅ Updated {updated_count} members with SMARTS patterns")
print(f"📊 Total categories: {len(reactant_types)}")
print(f"📊 Total members: {sum(len(v['members']) for v in reactant_types.values())}")

if missing_patterns:
    print(f"\n⚠️  Missing SMARTS patterns for {len(missing_patterns)} members:")
    for mp in missing_patterns:
        print(f"   - {mp}")
else:
    print("\n✅ All members have SMARTS patterns!")

print("\n📁 Updated: data-processor/other_data/reactant_types.json")
