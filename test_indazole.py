from rdkit import Chem
from rdkit.Chem import AllChem

rxn_smiles = 'IC1=NNC2=C1C=CC=C2>>N#CC3=NNC4=CC=CC=C43'
pattern = 'IC1=NNC2=C1C=CC=C2>>N#CC3=NNC4=CC=CC=C43'

rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
pat_rxn = AllChem.ReactionFromSmarts(pattern)

print(f'Reaction parsed: {rxn is not None}')
print(f'Pattern parsed: {pat_rxn is not None}')

rxn_r = rxn.GetReactants()
rxn_p = rxn.GetProducts()
pat_r = pat_rxn.GetReactants()
pat_p = pat_rxn.GetProducts()

print(f'\nRxn reactants: {len(rxn_r)}, products: {len(rxn_p)}')
print(f'Pat reactants: {len(pat_r)}, products: {len(pat_p)}')

# Sanitize
for mol in rxn_r:
    Chem.SanitizeMol(mol)
for mol in rxn_p:
    Chem.SanitizeMol(mol)

print(f'\nRxn R0 SMILES: {Chem.MolToSmiles(rxn_r[0])}')
print(f'Pat R0 SMILES: {Chem.MolToSmiles(pat_r[0])}')
print(f'\nRxn P0 SMILES: {Chem.MolToSmiles(rxn_p[0])}')
print(f'Pat P0 SMILES: {Chem.MolToSmiles(pat_p[0])}')

print(f'\nReactant match: {rxn_r[0].HasSubstructMatch(pat_r[0])}')
print(f'Product match: {rxn_p[0].HasSubstructMatch(pat_p[0])}')

# Check if they're exactly equal
print(f'\nReactant exact match: {Chem.MolToSmiles(rxn_r[0]) == Chem.MolToSmiles(pat_r[0])}')
print(f'Product exact match: {Chem.MolToSmiles(rxn_p[0]) == Chem.MolToSmiles(pat_p[0])}')
