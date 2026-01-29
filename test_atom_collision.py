"""Test per-reactant motif selection to fix atom index collision bug."""

from chemtools.featurizers.formatters.aggregation import select_primary_motifs_by_atom

# Simulate TWO reactants with atom index collision
# Reactant 1 (piperazine): has RCH2-NHR at atom 3
# Reactant 2 (aniline): has Ar-Cl at atom 3
reactant1_motifs = [
    {"id": "RCH2-NHR", "a_idx": 3, "reactivity_weight": 0},
    {"id": "RCH2-NR2", "a_idx": 2, "reactivity_weight": 0},
]

reactant2_motifs = [
    {"id": "Ar-Cl", "a_idx": 3, "reactivity_weight": 100},  # Higher weight
    {"id": "Ar-NH2", "a_idx": 1, "reactivity_weight": 0},
]

print("=== WRONG WAY: Combine first, then select (causes atom collision) ===")
combined_wrong = reactant1_motifs + reactant2_motifs
selected_wrong = select_primary_motifs_by_atom(combined_wrong)
print(f"Input: {[m['id'] for m in combined_wrong]}")
print(f"Output: {[m['id'] for m in selected_wrong]}")
print(f"RCH2-NHR preserved? {'RCH2-NHR' in [m['id'] for m in selected_wrong]}")
print()

print("=== RIGHT WAY: Select per reactant, then combine ===")
selected_r1 = select_primary_motifs_by_atom(reactant1_motifs)
selected_r2 = select_primary_motifs_by_atom(reactant2_motifs)
combined_right = selected_r1 + selected_r2
print(f"Reactant 1 output: {[m['id'] for m in selected_r1]}")
print(f"Reactant 2 output: {[m['id'] for m in selected_r2]}")
print(f"Combined: {[m['id'] for m in combined_right]}")
print(f"RCH2-NHR preserved? {'RCH2-NHR' in [m['id'] for m in combined_right]}")
