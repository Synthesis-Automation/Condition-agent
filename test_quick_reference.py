"""
Quick Reference: Clean Featurizer System
=========================================

Two-Tier Architecture (v2 only - No Legacy)
"""

from chemtools.featurizers.formatters.molecule import featurize_molecule
from chemtools.featurizers.formatters.reaction import featurize_reaction

# ============================================================================
# CORE FORMAT (DEFAULT) - Simplified for common use
# ============================================================================

# Molecule: 6 fields (smiles, motifs, properties, rdkit, kind, schema_version)
mol_core = featurize_molecule("c1ccccc1Br")
print("Core Molecule:")
print(f"  Fields: {list(mol_core.keys())}")
print(f"  Motifs: {mol_core['motifs']}")
print(f"  Properties: {mol_core['properties']}")
print()

# Reaction: 9 fields (reaction_smiles, reaction_type, confidence, reaction_key,
#                     reactants, products, feasibility, kind, schema_version)
rxn_core = featurize_reaction("c1ccc(Br)cc1.CC(=O)O>>CC(=O)c1ccccc1")
print("Core Reaction:")
print(f"  Fields: {list(rxn_core.keys())}")
print(f"  Type: {rxn_core['reaction_type']}")
print(f"  Key: {rxn_core['reaction_key']}")
print()

# ============================================================================
# EXTENDED FORMAT - Detailed analysis when needed
# ============================================================================

# Molecule: 7 fields (core + extended section)
mol_ext = featurize_molecule("c1ccccc1Br", options={"detailed": True})
print("Extended Molecule:")
print(f"  Fields: {list(mol_ext.keys())}")
print(f"  Has extended: {'extended' in mol_ext}")
if 'extended' in mol_ext:
    print(f"  Extended keys: {list(mol_ext['extended'].keys())}")
print()

# Reaction: 10 fields (core + extended section)
rxn_ext = featurize_reaction("c1ccc(Br)cc1.CC(=O)O>>CC(=O)c1ccccc1", options={"detailed": True})
print("Extended Reaction:")
print(f"  Fields: {list(rxn_ext.keys())}")
print(f"  Has extended: {'extended' in rxn_ext}")
if 'extended' in rxn_ext:
    print(f"  Extended keys: {list(rxn_ext['extended'].keys())}")
print()

# ============================================================================
# KEY DIFFERENCES FROM OLD SYSTEM
# ============================================================================
print("=" * 72)
print("What Changed:")
print("  ❌ No more nested .molecule or .reaction wrappers")
print("  ❌ No more legacy=True option")
print("  ❌ No more compound_id/rank_score fields (now id/rank)")
print("  ✅ Direct access: result['smiles'] instead of result['molecule']['smiles']")
print("  ✅ Smaller outputs: 54-56% fewer fields by default")
print("  ✅ Single schema: All outputs are v2")
print("=" * 72)
