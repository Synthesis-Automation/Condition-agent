from chemtools.featurizers.unified import featurize_molecule
import json

def test_motifs(smiles):
    payload = featurize_molecule(smiles)
    motifs = payload['molecule']['motifs']
    print(f"SMILES: {smiles}")
    print(f"Motifs: {[m['compound_id'] for m in motifs]}")

if __name__ == "__main__":
    test_motifs("C=C")
    test_motifs("Brc1ccccc1")
    test_motifs("c1ccc([Zn]Cl)cc1")
