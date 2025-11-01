"""Quick test of interactive bond analysis."""

import sys
from io import StringIO
from app.reaction_analysis_interactive import display_bond_analysis

# Test with Suzuki-Miyaura coupling
print("Testing Suzuki-Miyaura Coupling:")
print("=" * 80)
smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
print(f"SMILES: {smiles}")
print()

display_bond_analysis(smiles)
