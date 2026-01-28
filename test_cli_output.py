"""
Test script to verify CLI output with new extended format.
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[0]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction

# Test molecule
print("=" * 80)
print("TEST 1: Simple Molecule (Bromobenzene)")
print("=" * 80)
mol_result = featurize_molecule("c1ccccc1Br", options={"detailed": True})
print(f"Schema: {mol_result.get('schema_version')}")
print(f"SMILES: {mol_result.get('smiles')}")
print(f"Motifs: {[m['id'] for m in mol_result.get('motifs', [])]}")
print(f"Has extended: {'extended' in mol_result}")
if 'extended' in mol_result:
    extended = mol_result['extended']
    print(f"Extended keys: {list(extended.keys())}")
    if 'analyses' in extended:
        print(f"  - analyses: {len(extended['analyses'])} entries")

print("\n" + "=" * 80)
print("TEST 2: Suzuki Coupling Reaction")
print("=" * 80)
rxn_result = featurize_reaction(
    "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1",
    options={"detailed": True, "confirm_coupling_products": True}
)
print(f"Schema: {rxn_result.get('schema_version')}")
print(f"Reaction SMILES: {rxn_result.get('reaction_smiles')}")
print(f"Reaction Key: {rxn_result.get('reaction_key')}")
print(f"Reaction Type: {rxn_result.get('reaction_type')}")
print(f"Confidence: {rxn_result.get('confidence')}")
print(f"Has extended: {'extended' in rxn_result}")
if 'extended' in rxn_result:
    extended = rxn_result['extended']
    print(f"Extended keys: {list(extended.keys())}")
    if 'roles' in extended:
        roles = extended['roles']
        print(f"  - roles.reaction_type: {roles.get('reaction_type')}")
        print(f"  - roles.num_reactants: {roles.get('num_reactants')}")
    if 'agent_roles' in extended:
        agent_roles = extended['agent_roles']
        print(f"  - agent_roles.agent_count: {agent_roles.get('agent_count')}")
        print(f"  - agent_roles.role_counts: {agent_roles.get('role_counts')}")

# Test reactants
print("\nReactants:")
for idx, reactant in enumerate(rxn_result.get('reactants', []), 1):
    print(f"  Reactant {idx}:")
    print(f"    SMILES: {reactant.get('smiles')}")
    print(f"    Motifs: {[m['id'] for m in reactant.get('motifs', [])]}")
    print(f"    Has extended: {'extended' in reactant}")

print("\n✅ All tests completed!")
