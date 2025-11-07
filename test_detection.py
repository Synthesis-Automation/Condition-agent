"""Quick test to see what families are detected"""
from chemtools import detect_reaction

reactions = {
    "Suzuki": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    "Amide": "CC(=O)Cl.NCC>>CC(=O)NCC",
    "SN2": "CBr.NCC>>CCNCC",
    "SNAr": "[O-][N+](=O)c1ccc([N+](=O)[O-])cc1F.NH2C2H4OH>>[O-][N+](=O)c1ccc([N+](=O)[O-])cc1NCCO",
    "Radical": "CC.CCl4>>CCl+CCCl3",
}

for name, smiles in reactions.items():
    try:
        result = detect_reaction(smiles)
        print(f"{name:15} -> {result.get('family', 'Unknown'):30} (conf: {result.get('confidence', 0):.2f})")
    except Exception as e:
        print(f"{name:15} -> ERROR: {e}")
