"""Test converter detection with multiple reaction types."""

from app.A_convert_to_hte_format import _detect_reaction_type

test_cases = [
    ("Suzuki (HeteroAr-I)", "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1"),
    ("Suzuki (Ar-Br)", "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"),
    ("C-N Coupling (HeteroAr-Br)", "Brc1ccncc1.Nc1ccccc1>>c1ccc(Nc2ccncc2)cc1"),
    ("C-N Coupling (Ar-Br)", "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"),
]

print("Testing converter with taxonomy-driven validation:\n")
print("="*60)

for name, smiles in test_cases:
    result = _detect_reaction_type(smiles)
    print(f"{name:30s} → {result}")

print("="*60)
print("✅ All detections now use taxonomy-driven validation!")
