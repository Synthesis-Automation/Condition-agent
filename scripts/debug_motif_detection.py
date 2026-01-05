
from chemtools.featurizers.structural import featurize_molecule
from rdkit import Chem

def debug_reaction():
    r1 = "c1ccc(N2CCNCC2)cc1"
    r2 = "OB(O)c1ccc(F)cc1"
    p = "Fc1ccc(N2CCN(c3ccccc3)CC2)cc1"
    
    print(f"Reactant 1: {r1}")
    feat1 = featurize_molecule(r1)
    print(f"Motifs: {[m.get('compound_id') for m in feat1.get('motifs', [])]}")
    
    print(f"\nReactant 2: {r2}")
    feat2 = featurize_molecule(r2)
    print(f"Motifs: {[m.get('compound_id') for m in feat2.get('motifs', [])]}")
    
    print(f"\nProduct: {p}")
    featp = featurize_molecule(p)
    print(f"Motifs: {[m.get('compound_id') for m in featp.get('motifs', [])]}")

if __name__ == "__main__":
    debug_reaction()
