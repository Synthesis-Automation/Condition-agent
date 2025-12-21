
from rdkit import Chem
mol = Chem.MolFromSmiles("c1cc[nH]c1")
smarts = "[NX3H1;!$([N][C,c]=O);!$([N]S(=O)(=O))]"
pattern = Chem.MolFromSmarts(smarts)
print(f"Pyrrole match: {mol.HasSubstructMatch(pattern)}")

mol2 = Chem.MolFromSmiles("CCNCC")
print(f"Diethylamine match: {mol2.HasSubstructMatch(pattern)}")
