"""
Quick reference for using the updated CLI with extended output.
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction

# Example 1: Molecule with detailed analysis
print("=" * 60)
print("MOLECULE EXAMPLE: 4-Bromobenzonitrile")
print("=" * 60)

mol = featurize_molecule("N#Cc1ccc(Br)cc1", options={"detailed": True})

print(f"SMILES: {mol['smiles']}")
print(f"Schema: {mol['schema_version']}")
print(f"\nCore Fields ({len([k for k in mol.keys() if k != 'extended'])} fields):")
for key in mol.keys():
    if key != 'extended':
        print(f"  - {key}")

print(f"\nMotifs Found:")
for m in mol['motifs']:
    print(f"  - {m['id']} (rank: {m['rank']})")

extended = mol.get('extended', {})
print(f"\nExtended Analysis:")
print(f"  Per-motif analysis entries: {len(extended.get('per_motif_analysis', []))}")
print(f"  SNAr feasibility entries: {len(extended.get('snar_feasibility', []))}")

for analysis in extended.get('per_motif_analysis', []):
    print(f"\n  {analysis['motif_id']}:")
    if 'steric' in analysis:
        print(f"    Steric: {analysis['steric']}")
    if 'electronic' in analysis:
        print(f"    Electronic: {analysis['electronic']}")
    if 'nearby_groups' in analysis:
        print(f"    Nearby: {analysis['nearby_groups']}")

# Example 2: Reaction with detailed analysis
print("\n" + "=" * 60)
print("REACTION EXAMPLE: Buchwald-Hartwig Amination")
print("=" * 60)

rxn = featurize_reaction(
    "c1ccccc1Br.c1ccccc1N>>c1ccccc1Nc1ccccc1",
    options={
        "detailed": True,
        "confirm_coupling_products": True,
    }
)

print(f"Reaction SMILES: {rxn['reaction_smiles']}")
print(f"Reaction Key: {rxn['reaction_key']}")
print(f"Reaction Type: {rxn['reaction_type']} (confidence: {rxn['confidence']})")
print(f"Feasibility: {rxn['feasibility']}")

extended = rxn.get('extended', {})
print(f"\nExtended Analysis:")
print(f"  Detection matches: {len(extended.get('detection', {}).get('matches', []))}")

aggregates = extended.get('aggregates', {})
print(f"\n  Aggregates:")
for key, value in sorted(aggregates.items()):
    print(f"    {key}: {value}")

role_class = extended.get('role_classification', {})
print(f"\n  Role Classification:")
reactant_roles = role_class.get('reactants', {})
print(f"    Reaction type: {reactant_roles.get('reaction_type')}")
print(f"    Num reactants: {reactant_roles.get('num_reactants')}")
print(f"    Reactant details:")
for r in reactant_roles.get('reactants', []):
    print(f"      - {r.get('name')} ({r.get('role')})")

agent_roles = role_class.get('agents', {})
print(f"\n  Agent Roles:")
print(f"    Agent count: {agent_roles.get('agent_count', 0)}")
print(f"    Role flags: {list(agent_roles.get('flags', {}).keys())}")

print("\n" + "=" * 60)
print("KEY TAKEAWAYS FOR CONDITIONS RECOMMENDATION")
print("=" * 60)
print("""
1. Motifs with ranks identify functional groups
2. Steric/electronic scores quantify reactivity
3. Reaction type + confidence classify the transformation
4. Aggregates provide reaction-wide statistics
5. Role classification identifies reactant/agent roles
6. All information is in a structured, queryable format

Use this data to:
- Match reactions to condition databases
- Filter by reaction type and substrate features
- Assess compatibility with known protocols
- Generate condition recommendations
""")
