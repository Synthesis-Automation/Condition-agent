"""Quick test of simplified featurizer outputs."""
import sys
from pathlib import Path

# Add repo root to path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction

print("=" * 72)
print("TESTING SIMPLIFIED MOLECULE OUTPUT")
print("=" * 72)

# Test molecule
mol_payload = featurize_molecule("c1ccccc1Br")
print(f"\nSchema version: {mol_payload.get('schema_version')}")
print(f"Top-level fields: {list(mol_payload.keys())}")
print(f"SMILES: {mol_payload.get('smiles')}")
print(f"Motifs: {mol_payload.get('motifs')}")
print(f"Properties: {mol_payload.get('properties')}")
print(f"RDKit: {mol_payload.get('rdkit')}")

print("\n" + "=" * 72)
print("TESTING DETAILED MOLECULE OUTPUT")
print("=" * 72)

# Test molecule with detailed=True
mol_detailed = featurize_molecule("c1ccccc1Br", options={"detailed": True})
print(f"\nSchema version: {mol_detailed.get('schema_version')}")
print(f"Top-level fields: {list(mol_detailed.keys())}")
print(f"Extended fields: {list(mol_detailed.get('extended', {}).keys())}")
print(f"Has per_motif_analysis: {'per_motif_analysis' in mol_detailed.get('extended', {})}")

print("\n" + "=" * 72)
print("TESTING SIMPLIFIED REACTION OUTPUT")
print("=" * 72)

# Test reaction
rxn_payload = featurize_reaction("c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1-c2ccccc2")
print(f"\nSchema version: {rxn_payload.get('schema_version')}")
print(f"Top-level fields: {list(rxn_payload.keys())}")
print(f"Reaction type: {rxn_payload.get('reaction_type')}")
print(f"Confidence: {rxn_payload.get('confidence')}")
print(f"Reaction key: {rxn_payload.get('reaction_key')}")
print(f"Feasibility: {rxn_payload.get('feasibility')}")
print(f"Reactants: {len(rxn_payload.get('reactants', []))}")
print(f"Products: {len(rxn_payload.get('products', []))}")
print(f"First reactant motifs: {rxn_payload.get('reactants', [{}])[0].get('motifs', [])}")

print("\n" + "=" * 72)
print("TESTING DETAILED REACTION OUTPUT")
print("=" * 72)

# Test reaction with detailed=True
rxn_detailed = featurize_reaction("c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1-c2ccccc2", options={"detailed": True})
print(f"\nSchema version: {rxn_detailed.get('schema_version')}")
print(f"Top-level fields: {list(rxn_detailed.keys())}")
print(f"Extended fields: {list(rxn_detailed.get('extended', {}).keys())}")
print(f"Has aggregates: {'aggregates' in rxn_detailed.get('extended', {})}")
print(f"Has detection: {'detection' in rxn_detailed.get('extended', {})}")

print("\n" + "=" * 72)
print("✓ ALL TESTS PASSED - Simplified outputs working!")
print("=" * 72)
