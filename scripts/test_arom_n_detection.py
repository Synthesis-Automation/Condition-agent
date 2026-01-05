from rdkit import Chem
from chemtools.featurizers.molecule import featurize_molecule

def test_arom_n():
    # N-phenylpyrrole
    smiles = "c1ccc(N2C=CC=C2)cc1"
    
    result = featurize_molecule(smiles)
    
    print(f"SMILES: {smiles}")
    print("Detected Motifs:")
    for motif in result["motifs"]:
        print(f"  - {motif}")

if __name__ == "__main__":
    test_arom_n()
