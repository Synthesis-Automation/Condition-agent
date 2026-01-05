from rdkit import Chem

def test_smarts():
    smiles = "c1ccc(N2C=CC=C2)cc1"
    mol = Chem.MolFromSmiles(smiles)
    
    # Test [c:1]-[n:2]
    smarts = "[c:1]-[n:2]"
    query = Chem.MolFromSmarts(smarts)
    
    matches = mol.GetSubstructMatches(query)
    print(f"SMILES: {smiles}")
    print(f"SMARTS: {smarts}")
    print(f"Matches: {matches}")

    # Test [c:1]-[NX3;H0:2]
    smarts2 = "[c:1]-[NX3;H0:2]"
    query2 = Chem.MolFromSmarts(smarts2)
    matches2 = mol.GetSubstructMatches(query2)
    print(f"SMARTS: {smarts2}")
    print(f"Matches: {matches2}")

    # Test [c:1]-[#7X3;H0:2]
    smarts3 = "[c:1]-[#7X3;H0:2]"
    query3 = Chem.MolFromSmarts(smarts3)
    matches3 = mol.GetSubstructMatches(query3)
    print(f"SMARTS: {smarts3}")
    print(f"Matches: {matches3}")

if __name__ == "__main__":
    test_smarts()
