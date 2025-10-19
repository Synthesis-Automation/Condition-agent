"""Test if RDKit can automatically identify reaction centers via atom mapping"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions

# Test reaction: Sonogashira coupling
rxn_smiles = "NC(=O)c1ccccc1I.C#CCCCC>>NC(=O)c1ccccc1C#CCCCC"

print("="*80)
print("REACTION CENTER IDENTIFICATION TEST")
print("="*80)
print(f"\nReaction: {rxn_smiles}")

# Try to parse reaction
try:
    rxn = AllChem.ReactionFromSmarts(rxn_smiles)
    if rxn:
        print(f"✓ Reaction parsed successfully")
        
        # Try to get reactants and products
        print(f"\nReactants: {rxn.GetNumReactantTemplates()}")
        print(f"Products: {rxn.GetNumProductTemplates()}")
        
        # Check if RDKit can automatically identify changed atoms
        # This requires comparing reactants to products
        reactant_mols = [Chem.MolFromSmiles(s) for s in rxn_smiles.split('>>')[0].split('.')]
        product_mol = Chem.MolFromSmiles(rxn_smiles.split('>>')[1])
        
        print(f"\nReactant 1 atoms: {reactant_mols[0].GetNumAtoms()}")
        print(f"Reactant 2 atoms: {reactant_mols[1].GetNumAtoms()}")
        print(f"Product atoms: {product_mol.GetNumAtoms()}")
        
except Exception as e:
    print(f"✗ Error: {e}")

print("\n" + "="*80)
print("ATOM MAPPING APPROACHES")
print("="*80)

print("""
1. **Manual Atom Mapping** (Most Reliable)
   - Provide mapped SMILES in protocol JSON
   - Example: [C:1]#[C:2][C:3][C:4].[c:5][I:6]>>[C:1]#[C:2]-[c:5]
   - PROS: Explicit, unambiguous
   - CONS: Requires manual work or automated mapper

2. **RDKit Reaction Mapping** (Complex)
   - Use RDKit's ReactionFromSmarts with mapping
   - Or use external tools like RXNMapper, ChemAxon
   - PROS: Can be automated
   - CONS: Not always accurate, requires setup

3. **Heuristic + Validation** (Current Approach)
   - Identify likely reaction centers by functional groups
   - Compare reactants → products to validate
   - PROS: No extra input needed
   - CONS: Can fail with complex molecules

4. **Hybrid Approach** (Recommended)
   - Use mapped SMILES if available in protocol
   - Fall back to heuristics if not mapped
   - Add validation step to detect mismatches
""")

print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)
print("""
Best approach for your use case:

1. **Update Protocol Schema** to include optional atom-mapped SMILES:
   {
     "reaction_smiles": "NC(=O)c1ccccc1I.C#CCCCC>>NC(=O)c1ccccc1C#CCCCC",
     "reaction_smiles_mapped": "[NH2:7][C:8](=[O:9])[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[I:10].[C:11]#[C:12][CH2:13][CH2:14][CH2:15][CH3:16]>>[NH2:7][C:8](=[O:9])[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[C:11]#[C:12][CH2:13][CH2:14][CH2:15][CH3:16]"
   }

2. **Use RXNMapper** or similar tool to generate mappings automatically
   - GitHub: rxn4chemistry/rxnmapper
   - Can process batch of reactions

3. **SMARTS Generator** checks for mapped SMILES first:
   - If mapped → use mapping to identify reaction center
   - If not mapped → fall back to heuristics + warn user
""")
