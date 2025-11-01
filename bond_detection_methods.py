#!/usr/bin/env python3
"""
Bond Breaking/Formation Detection Methods - Complete Overview
==============================================================

This shows ALL available methods for detecting bond changes in reactions,
not just RXNMapper.
"""

print("="*80)
print("BOND BREAKING/FORMATION DETECTION METHODS")
print("="*80)

print("""
RXNMapper is NOT the only way! Here are all available methods:

┌─────────────────────────────────────────────────────────────────────────────┐
│ METHOD 1: Manual Atom Mapping (WORKS NOW - NO DEPENDENCIES)                │
└─────────────────────────────────────────────────────────────────────────────┘

If you have atom-mapped SMILES (manually created or from literature):

from chemtools import identify_changed_atoms_from_mapped_smiles

# Manually mapped Sonogashira
mapped = "[c:1][Br:2].[C:3]#[C:4]>>[c:1][C:3]#[C:4]"

result = identify_changed_atoms_from_mapped_smiles(mapped)
# Returns: broken_bonds, formed_bonds, changed_atoms

✅ Pros:
  • Works immediately (no extra dependencies)
  • 100% accurate if mapping is correct
  • Uses only RDKit (already installed)

❌ Cons:
  • Requires manual atom mapping
  • Time-consuming for many reactions


┌─────────────────────────────────────────────────────────────────────────────┐
│ METHOD 2: RXNMapper Auto-Mapping (AUTOMATED - RECOMMENDED)                 │
└─────────────────────────────────────────────────────────────────────────────┘

Automatic atom mapping using ML attention-guided approach:

from chemtools import analyze_bond_changes

# Unmapped reaction
unmapped = "Brc1ccccc1.C#C>>c1ccccc1C#C"

result = analyze_bond_changes(unmapped, auto_map=True)
# Automatically maps, then analyzes

✅ Pros:
  • Fully automatic
  • High accuracy (transformer-based)
  • Handles complex reactions
  • Confidence scores provided

❌ Cons:
  • Requires: pip install rxnmapper (optional dependency)
  • Slower than pre-mapped SMILES
  • ~100MB model download on first use


┌─────────────────────────────────────────────────────────────────────────────┐
│ METHOD 3: RDKit MCS (Maximum Common Substructure) - NOT YET IMPLEMENTED    │
└─────────────────────────────────────────────────────────────────────────────┘

Uses RDKit's built-in MCS algorithm to find unchanged parts:

# COULD BE IMPLEMENTED:
from rdkit.Chem import rdFMCS

def find_bond_changes_with_mcs(reaction_smiles):
    reactants, products = reaction_smiles.split('>>')
    r_mol = Chem.MolFromSmiles(reactants)
    p_mol = Chem.MolFromSmiles(products)
    
    # Find maximum common substructure
    mcs = rdFMCS.FindMCS([r_mol, p_mol])
    # Atoms NOT in MCS are the reaction center
    # Compare bonds to find changes

✅ Pros:
  • No additional dependencies (uses RDKit)
  • Works without atom mapping
  • Deterministic algorithm

❌ Cons:
  • Less accurate than RXNMapper
  • Struggles with rearrangements
  • Can miss reaction center in complex cases
  • NOT YET IMPLEMENTED (would need ~100 lines of code)


┌─────────────────────────────────────────────────────────────────────────────┐
│ METHOD 4: ReactionFromSmarts + SMARTS Patterns (PARTIAL - WORKS NOW)       │
└─────────────────────────────────────────────────────────────────────────────┘

Use RDKit's ChemicalReaction to analyze reactions:

# ALREADY PARTIALLY USED in chemtools/router.py
from rdkit.Chem import AllChem

rxn = AllChem.ReactionFromSmarts('[C:1]Br.[C:2]#C>>[C:1][C:2]#C')
# Can extract which bonds appear in products but not reactants

✅ Pros:
  • RDKit native (no extra dependencies)
  • Fast pattern matching
  • Good for known reaction types

❌ Cons:
  • Requires predefined SMARTS for each reaction type
  • Only works for known patterns
  • Can't handle novel reactions
  • Currently only used for detection, not bond analysis


┌─────────────────────────────────────────────────────────────────────────────┐
│ METHOD 5: Graph-Based Difference (COULD BE IMPLEMENTED)                    │
└─────────────────────────────────────────────────────────────────────────────┘

Compare molecular graphs directly:

# COULD BE IMPLEMENTED:
def graph_based_bond_diff(reaction_smiles):
    r_mol, p_mol = parse_reaction(reaction_smiles)
    
    # Build adjacency matrices
    r_bonds = set((i,j) for i,j in r_mol.GetBonds())
    p_bonds = set((i,j) for i,j in p_mol.GetBonds())
    
    broken = r_bonds - p_bonds
    formed = p_bonds - r_bonds

✅ Pros:
  • Simple algorithm
  • No ML required
  • Fast

❌ Cons:
  • Requires some way to align atoms (back to mapping problem)
  • Less reliable than MCS
  • NOT YET IMPLEMENTED


┌─────────────────────────────────────────────────────────────────────────────┐
│ METHOD 6: Other Atom Mapping Tools (ALTERNATIVES TO RXNMAPPER)             │
└─────────────────────────────────────────────────────────────────────────────┘

Alternative tools for atom mapping:

1. GraphMapper (IBM RXN4Chemistry)
   pip install rxn-graph-mapper
   • Graph-based approach
   • Fast but less accurate than RXNMapper

2. LocalMapper 
   pip install localmapper
   • Local alignment approach
   • Good for similar reactants/products

3. Graphormer (Microsoft)
   • Transformer-based
   • State-of-the-art accuracy
   • Heavier than RXNMapper

4. ChemAxon Reactor
   • Commercial tool
   • Very accurate
   • Requires license

5. NameRxn (Name-based)
   • Uses reaction names
   • Limited to known reactions

All of these would work with our existing infrastructure - just replace
the RXNMapper call in _atom_mapping.py with the alternative tool.


┌─────────────────────────────────────────────────────────────────────────────┐
│ CURRENT IMPLEMENTATION SUMMARY                                              │
└─────────────────────────────────────────────────────────────────────────────┘

What's IMPLEMENTED NOW:
  ✅ Method 1: Manual atom mapping → bond analysis (WORKS)
  ✅ Method 2: RXNMapper auto-mapping → bond analysis (WORKS)
  ✅ Method 4: SMARTS patterns for detection (PARTIAL)

What's POSSIBLE but NOT implemented:
  ⏳ Method 3: RDKit MCS approach (~100 lines to add)
  ⏳ Method 5: Graph-based difference (~50 lines to add)
  ⏳ Method 6: Alternative mapping tools (easy to swap in)

What's RECOMMENDED:
  🏆 For unmapped reactions: RXNMapper (Method 2)
  🏆 For mapped reactions: Direct analysis (Method 1)
  🏆 For no dependencies: Implement MCS (Method 3)


┌─────────────────────────────────────────────────────────────────────────────┐
│ DECISION MATRIX                                                             │
└─────────────────────────────────────────────────────────────────────────────┘

Choose your method based on:

┌──────────────────┬─────────┬──────────┬────────────┬──────────────┐
│ Method           │ Deps    │ Accuracy │ Speed      │ Use Case     │
├──────────────────┼─────────┼──────────┼────────────┼──────────────┤
│ Manual Mapping   │ RDKit   │ ⭐⭐⭐⭐⭐ │ ⚡⚡⚡⚡⚡   │ Pre-mapped   │
│ RXNMapper        │ +ML     │ ⭐⭐⭐⭐   │ ⚡⚡        │ Auto (best)  │
│ MCS              │ RDKit   │ ⭐⭐⭐     │ ⚡⚡⚡⚡     │ Simple rxns  │
│ SMARTS           │ RDKit   │ ⭐⭐       │ ⚡⚡⚡⚡⚡   │ Known types  │
│ Graph Diff       │ RDKit   │ ⭐⭐       │ ⚡⚡⚡⚡⚡   │ Quick/dirty  │
└──────────────────┴─────────┴──────────┴────────────┴──────────────┘

""")

print("="*80)
print("ANSWER: NO, RXNMapper is NOT the only way!")
print("="*80)
print("""
You have multiple options:

1. ✅ ALREADY WORKING: Manual atom-mapped SMILES
   → Uses identify_changed_atoms_from_mapped_smiles()
   → No extra dependencies

2. ✅ ALREADY WORKING: RXNMapper (recommended for automation)
   → Uses analyze_bond_changes() with auto_map=True
   → Requires: pip install rxnmapper

3. ⏳ CAN ADD: RDKit MCS approach
   → Would work without RXNMapper
   → Just needs ~100 lines of implementation

4. ⏳ CAN ADD: Alternative mapping tools
   → GraphMapper, LocalMapper, Graphormer
   → Easy to swap in place of RXNMapper

5. ⏳ CAN ADD: Graph-based simple difference
   → Very simple algorithm
   → Lower accuracy but no ML

💡 RECOMMENDATION:
  • Use RXNMapper for production (best accuracy)
  • Fall back to MCS if RXNMapper not installed
  • We could implement MCS as a fallback!
""")

print("="*80)
