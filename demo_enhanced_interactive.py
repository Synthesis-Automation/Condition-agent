"""
Quick Demo: Enhanced Interactive App with Reaction-Specific Reactant Analysis

This shows the new feature added to reaction_analysis_interactive.py
"""

from chemtools import detect_reaction
from app.reaction_analysis_interactive import display_detection_only

print("=" * 80)
print("ENHANCED INTERACTIVE APP - NEW REACTANT CLASSIFICATION FEATURE")
print("=" * 80)
print()
print("The interactive app now includes reaction-specific reactant analysis!")
print()

# Example 1: Suzuki
print("\n" + "=" * 80)
print("EXAMPLE 1: SUZUKI-MIYAURA COUPLING")
print("=" * 80)
suzuki = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
display_detection_only(suzuki, use_ml=True)

# Example 2: Buchwald-Hartwig
print("\n" + "=" * 80)
print("EXAMPLE 2: BUCHWALD-HARTWIG C-N COUPLING")
print("=" * 80)
bh = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
display_detection_only(bh, use_ml=True)

print("\n" + "=" * 80)
print("WHAT'S NEW:")
print("=" * 80)
print("""
The interactive app now shows:

1. ✓ Reaction type detection (as before)
2. ✓ ML vs Rule-based comparison (as before)
3. ✓ Functional groups detected (as before)
4. ✓ Catalyst detection (as before)
5. ✓✓ NEW: Reaction-Specific Reactant Classification ✓✓

For each detected reaction type, it now:
- Shows the EXPECTED reactant types (e.g., ArX* + ArB* for Suzuki)
- Classifies each actual reactant (e.g., ArBr = aryl bromide)
- Indicates whether reactants match expectations (✓ or ⚠)
- Shows alternative classifications if applicable

This makes it easy to verify that:
- For Suzuki: You have an aryl halide + boron reagent
- For Buchwald-Hartwig: You have an aryl halide + amine
- For other reactions: Appropriate reactant pairs
""")
print("=" * 80)
