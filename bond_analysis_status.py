#!/usr/bin/env python3
"""
Bond Breaking and Formation Analysis - Feature Status
=======================================================

Summary of what exists and what's missing for analyzing bond changes in reactions.
"""

print("="*80)
print("BOND ANALYSIS CAPABILITIES IN CONDITION-AGENT")
print("="*80)

print("\n📋 WHAT EXISTS:")
print("-" * 80)

print("\n1. ✅ Reaction Center Detector (chemtools/util/reaction_center_detector.py)")
print("   Location: chemtools.util.reaction_center_detector")
print("   Functions:")
print("   • identify_changed_atoms_from_mapped_smiles()")
print("     - Input: Atom-mapped reaction SMILES")
print("     - Returns: changed_atoms, broken_bonds, formed_bonds, spectator_atoms")
print("   • generate_smarts_from_mapped_reaction()")
print("   • compare_unmapped_reaction_to_find_changes() [stub only]")
print()
print("   Capabilities:")
print("   ✓ Identifies which atoms changed in reaction")
print("   ✓ Detects broken bonds (reactant → product)")
print("   ✓ Detects formed bonds (reactant → product)")
print("   ✓ Identifies spectator atoms (unchanged)")
print()
print("   Requirements:")
print("   ⚠ Needs ATOM-MAPPED SMILES (e.g., [C:1]#[C:2].[c:3][I:4]>>[C:1]#[C:2]-[c:3])")
print("   ⚠ Requires RDKit to be installed")

print("\n2. ✅ ML-Based Reaction Detection (rxn-insight)")
print("   Location: chemtools._ml_helpers, chemtools.detection")
print("   Library: rxn-insight (installed)")
print("   Functions:")
print("   • detect_reaction() - Unified API")
print("   • Reaction type classification")
print("   • Reaction name prediction")
print()
print("   Capabilities:")
print("   ✓ Classifies reaction type without atom mapping")
print("   ✓ Provides reaction names and classes")
print("   ✓ Confidence scores")
print()
print("   Limitations:")
print("   ✗ Does NOT provide atom mapping")
print("   ✗ Does NOT identify specific bonds broken/formed")

print("\n3. ✅ Reactant Classification (Two-Pass Approach)")
print("   Location: chemtools.analysis")
print("   Functions:")
print("   • classify_reactants_with_context()")
print("   • classify_reactant_smiles()")
print()
print("   Capabilities:")
print("   ✓ Identifies functional groups in reactants")
print("   ✓ Categorizes reactants by role")
print("   ✓ Context-aware classification")
print()
print("   Limitations:")
print("   ✗ Does NOT identify which bonds break/form")

print("\n" + "="*80)
print("❌ WHAT'S MISSING:")
print("="*80)

print("\n1. ❌ Automatic Atom Mapping")
print("   Current Status: NOT IMPLEMENTED")
print("   What's needed:")
print("   • RXNMapper integration (not currently used)")
print("   • Automatic atom mapping for unmapped reactions")
print("   • Alternative: graphormer/LocalMapper/other tools")

print("\n2. ❌ Bond Change Analysis for Unmapped Reactions")
print("   Current Status: STUB ONLY")
print("   Function exists but returns:")
print("   'Automatic reaction center detection without atom mapping is unreliable.'")
print("   ")
print("   What's needed:")
print("   • Maximum Common Substructure (MCS) implementation")
print("   • Heuristic bond change detection")
print("   • Integration with atom mapping tools")

print("\n3. ❌ Public API Exposure")
print("   Current Status: NOT EXPOSED")
print("   • reaction_center_detector is in util/ but not in chemtools.__init__")
print("   • Not available through 'from chemtools import ...'")
print("   • Not documented in public API")

print("\n" + "="*80)
print("💡 RECOMMENDATIONS:")
print("="*80)

print("\n1. For Atom-Mapped Reactions (WORKS NOW):")
print("   >>> from chemtools.util.reaction_center_detector import identify_changed_atoms_from_mapped_smiles")
print("   >>> mapped = '[c:1]1[c:2][I:3].[C:4]#[C:5]>>[c:1]1[c:2][C:4]#[C:5]'")
print("   >>> result = identify_changed_atoms_from_mapped_smiles(mapped)")
print("   >>> print(result['broken_bonds'])  # [(2, 3)] - C-I bond")
print("   >>> print(result['formed_bonds'])  # [(2, 4)] - new C-C bond")

print("\n2. For Unmapped Reactions (NEED TO ADD):")
print("   Option A: Add RXNMapper integration")
print("   >>> from rxnmapper import RXNMapper")
print("   >>> mapper = RXNMapper()")
print("   >>> mapped = mapper.get_attention_guided_atom_maps(['Br.C#C>>C#C'])")
print("   >>> # Then use reaction_center_detector")
print()
print("   Option B: Use rxn-insight's capabilities (if it has mapping)")
print("   >>> # Check if rxn-insight provides atom mapping")

print("\n3. Expose in Public API:")
print("   Add to chemtools/__init__.py:")
print("   >>> from chemtools.util.reaction_center_detector import (")
print("   >>>     identify_changed_atoms_from_mapped_smiles,")
print("   >>>     analyze_bond_changes  # New wrapper function")
print("   >>> )")

print("\n" + "="*80)
print("🎯 ANSWER TO YOUR QUESTION:")
print("="*80)
print("""
Q: Do we use RXNMapper to do this?
A: NO - Currently RXNMapper is NOT being used.

What IS being used:
  • rxn-insight: For reaction type classification (ML-based)
  • reaction_center_detector: For bond analysis (REQUIRES atom-mapped input)

What COULD be added:
  • RXNMapper: For automatic atom mapping (then feed to reaction_center_detector)
  • This would enable bond analysis for unmapped reactions

Current workflow for bond analysis:
  1. User provides atom-mapped SMILES manually
  2. Use reaction_center_detector.identify_changed_atoms_from_mapped_smiles()
  3. Get broken_bonds and formed_bonds

Recommended workflow (NOT YET IMPLEMENTED):
  1. User provides unmapped SMILES
  2. Use RXNMapper to add atom mapping automatically
  3. Use reaction_center_detector.identify_changed_atoms_from_mapped_smiles()
  4. Get broken_bonds and formed_bonds
""")

print("="*80)
