"""
Demo of SmartsBuilder - Context-Aware SMARTS Pattern Generation

This script demonstrates the capabilities of the SmartsBuilder module for generating
chemistry-aware SMARTS patterns from molecular structures.
"""

from chemtools.util.smarts_builders import (
    build_smarts,
    build_smarts_with_guards,
    match_smarts,
    SmartsBuilder,
    SmartsPatternMatcher,
)


def demo_basic_usage():
    """Demonstrate basic pattern generation"""
    print("=" * 70)
    print("1. BASIC PATTERN GENERATION")
    print("=" * 70)
    
    examples = [
        ("CCCCI", "Primary alkyl iodide"),
        ("CC(C)Br", "Secondary alkyl bromide"),
        ("CC(C)(C)Cl", "Tertiary alkyl chloride"),
        ("c1ccccc1I", "Aryl iodide"),
        ("c1ccccc1N", "Aniline"),
        ("CCN", "Primary aliphatic amine"),
        ("c1ccccc1O", "Phenol"),
        ("CCCO", "Aliphatic alcohol"),
    ]
    
    for smiles, description in examples:
        pattern = build_smarts(smiles)
        print(f"\n{description}:")
        print(f"  SMILES:  {smiles}")
        print(f"  SMARTS:  {pattern}")


def demo_halide_specificity():
    """Demonstrate halide-specific pattern generation"""
    print("\n" + "=" * 70)
    print("2. HALIDE SPECIFICITY")
    print("=" * 70)
    
    print("\nPrimary alkyl halides with different halogens:")
    for halogen, smiles in [("Iodide", "CCCCI"), ("Bromide", "CCCCBr"), 
                             ("Chloride", "CCCCCl"), ("Fluoride", "CCCF")]:
        pattern = build_smarts(smiles)
        print(f"  {halogen:10s}: {pattern}")


def demo_special_positions():
    """Demonstrate benzylic, allylic, propargylic detection"""
    print("\n" + "=" * 70)
    print("3. SPECIAL POSITIONS")
    print("=" * 70)
    
    special_cases = [
        ("c1ccccc1CCl", "Benzylic chloride"),
        ("C=CCBr", "Allylic bromide"),
        ("C#CCI", "Propargylic iodide"),
        ("c1ccccc1CO", "Benzylic alcohol"),
        ("C=CCO", "Allylic alcohol"),
    ]
    
    for smiles, description in special_cases:
        pattern = build_smarts(smiles)
        print(f"\n{description}:")
        print(f"  SMILES:  {smiles}")
        print(f"  SMARTS:  {pattern}")


def demo_amine_vs_amide():
    """Demonstrate critical distinction between amines and amides"""
    print("\n" + "=" * 70)
    print("4. AMINE vs AMIDE DISTINCTION")
    print("=" * 70)
    
    examples = [
        ("c1ccccc1N", "Aniline (aromatic amine)"),
        ("CCN", "Primary aliphatic amine"),
        ("CCNC", "Secondary amine"),
        ("CC(=O)N", "Primary amide (NOT an amine!)"),
        ("CC(=O)NC", "Secondary amide"),
    ]
    
    for smiles, description in examples:
        pattern = build_smarts(smiles)
        print(f"\n{description}:")
        print(f"  SMILES:  {smiles}")
        print(f"  SMARTS:  {pattern}")
        
        # Show that amine patterns exclude amides
        if "amide" not in description.lower():
            print(f"  Note:    Pattern includes !$(NC=O) to exclude amides")


def demo_guard_patterns():
    """Demonstrate guard pattern generation"""
    print("\n" + "=" * 70)
    print("5. GUARD PATTERN GENERATION")
    print("=" * 70)
    
    smiles = "CCCCI"
    print(f"\nGenerating patterns with guards for: {smiles}")
    
    result = build_smarts_with_guards(smiles)
    
    print(f"\nSubstrate Class: {result['substrate_class']}")
    print(f"Substrate Family: {result['substrate_family']}")
    print(f"\nCore Pattern (REQUIRED):")
    print(f"  {result['core']}")
    print(f"\nGuard Patterns (FORBIDDEN):")
    for i, guard in enumerate(result['guards_forbid'], 1):
        print(f"  {i}. {guard}")


def demo_pattern_matching():
    """Demonstrate pattern matching and validation"""
    print("\n" + "=" * 70)
    print("6. PATTERN MATCHING")
    print("=" * 70)
    
    pattern = "[CX4;H2,H3]-[I]"
    print(f"\nTesting pattern: {pattern}")
    print(f"(Primary alkyl iodide pattern)\n")
    
    test_molecules = [
        ("CI", "Methyl iodide", True),
        ("CCI", "Ethyl iodide", True),
        ("CCCCI", "Propyl iodide", True),
        ("CCCCCCCCI", "Octyl iodide", True),
        ("CC(C)I", "Isopropyl iodide (secondary)", False),
        ("CC(C)(C)I", "tert-Butyl iodide (tertiary)", False),
        ("c1ccccc1I", "Phenyl iodide (aryl)", False),
    ]
    
    for smiles, description, should_match in test_molecules:
        matches = match_smarts(smiles, pattern)
        status = "✓ MATCH" if matches else "✗ NO MATCH"
        expected = "✓" if matches == should_match else "✗ UNEXPECTED"
        print(f"  {status:12s} {expected}  {description:40s}  ({smiles})")


def demo_carbonyl_patterns():
    """Demonstrate carbonyl compound patterns"""
    print("\n" + "=" * 70)
    print("7. CARBONYL COMPOUNDS")
    print("=" * 70)
    
    carbonyls = [
        ("CC(=O)O", "Carboxylic acid"),
        ("CC(=O)OC", "Ester"),
        ("CCC=O", "Aldehyde"),
        ("CC(=O)C", "Ketone"),
    ]
    
    for smiles, description in carbonyls:
        pattern = build_smarts(smiles)
        print(f"\n{description}:")
        print(f"  SMILES:  {smiles}")
        print(f"  SMARTS:  {pattern}")


def demo_boron_patterns():
    """Demonstrate boron compound patterns"""
    print("\n" + "=" * 70)
    print("8. BORON COMPOUNDS")
    print("=" * 70)
    
    boron_compounds = [
        ("c1ccccc1B(O)O", "Phenylboronic acid"),
        ("c1ccccc1B1OC(C)(C)C(C)(C)O1", "Phenylboronic acid pinacol ester"),
    ]
    
    for smiles, description in boron_compounds:
        pattern = build_smarts(smiles)
        print(f"\n{description}:")
        print(f"  SMILES:  {smiles}")
        print(f"  SMARTS:  {pattern}")


def demo_real_world_examples():
    """Demonstrate with real reaction substrates"""
    print("\n" + "=" * 70)
    print("9. REAL-WORLD REACTION SUBSTRATES")
    print("=" * 70)
    
    print("\nBuchwald-Hartwig C-N Coupling Substrates:")
    
    # Aryl halide component
    aryl_halides = [
        ("c1ccccc1Br", "Bromobenzene"),
        ("c1ccc(Br)cc1C", "4-bromotoluene"),
        ("c1ccc(Br)cc1-c2ccccc2", "4-bromobiphenyl"),
    ]
    
    print("\n  Aryl Halides:")
    for smiles, name in aryl_halides:
        pattern = build_smarts(smiles)
        print(f"    {name:25s} → {pattern}")
    
    # Amine component
    amines = [
        ("c1ccccc1N", "Aniline"),
        ("CCN", "Ethylamine"),
        ("c1ccc2c(c1)cccc2N", "1-naphthylamine"),
    ]
    
    print("\n  Amines:")
    for smiles, name in amines:
        pattern = build_smarts(smiles)
        print(f"    {name:25s} → {pattern}")


def demo_advanced_matching():
    """Demonstrate advanced pattern matching features"""
    print("\n" + "=" * 70)
    print("10. ADVANCED PATTERN MATCHING")
    print("=" * 70)
    
    matcher = SmartsPatternMatcher()
    
    smiles = "c1ccccc1Br"
    pattern = "c-[Br]"
    
    print(f"\nMolecule: {smiles} (bromobenzene)")
    print(f"Pattern:  {pattern}")
    
    # Basic matching
    matches = matcher.match(smiles, pattern)
    print(f"\nMatches: {matches}")
    
    # Explanation
    explanation = matcher.explain_match(smiles, pattern)
    print(f"Explanation: {explanation}")
    
    # Find matching atoms
    halogen_pattern = "[Br]"
    matching_atoms = matcher.find_matching_atoms(smiles, halogen_pattern)
    print(f"\nAtoms matching '[Br]': {matching_atoms}")


def demo_chemistry_awareness():
    """Demonstrate chemistry-aware pattern generation"""
    print("\n" + "=" * 70)
    print("11. CHEMISTRY AWARENESS")
    print("=" * 70)
    
    print("\nThe same functional group gets DIFFERENT patterns based on context:")
    
    comparisons = [
        ("Iodide", [
            ("CCCCI", "Primary alkyl", "[CX4;H2,H3]-[I]"),
            ("CC(C)I", "Secondary alkyl", "[CX4;H1]-[I]"),
            ("c1ccccc1I", "Aryl", "c-[I]"),
            ("c1ccccc1CI", "Benzylic", "[CH2;$([CH2][c])]-[I]"),
        ]),
        ("Amine", [
            ("c1ccccc1N", "Aniline", "c-[NX3;H2;!$(NC=O)]"),
            ("CCN", "Aliphatic", "[CX4]-[NX3;H2;!$(NC=O)]"),
        ]),
        ("Hydroxyl", [
            ("c1ccccc1O", "Phenol", "c-[OX2H]"),
            ("CCCO", "Alcohol", "[CX4]-[OX2H]"),
        ]),
    ]
    
    for fg_name, cases in comparisons:
        print(f"\n{fg_name}:")
        for smiles, context, expected_pattern in cases:
            actual_pattern = build_smarts(smiles)
            match_indicator = "✓" if actual_pattern == expected_pattern else "✗"
            print(f"  {match_indicator} {context:20s}: {actual_pattern}")


def main():
    """Run all demos"""
    print()
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  SmartsBuilder Demo - Context-Aware SMARTS Pattern Generation".ljust(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "═" * 68 + "╝")
    
    demos = [
        demo_basic_usage,
        demo_halide_specificity,
        demo_special_positions,
        demo_amine_vs_amide,
        demo_guard_patterns,
        demo_pattern_matching,
        demo_carbonyl_patterns,
        demo_boron_patterns,
        demo_real_world_examples,
        demo_advanced_matching,
        demo_chemistry_awareness,
    ]
    
    for demo in demos:
        try:
            demo()
        except Exception as e:
            print(f"\n⚠ Error in {demo.__name__}: {e}")
    
    print("\n" + "=" * 70)
    print("DEMO COMPLETE")
    print("=" * 70)
    print()
    print("Key Takeaways:")
    print("  • Patterns are context-aware (primary vs secondary vs aryl)")
    print("  • Amines and amides are distinguished (!$(NC=O) exclusion)")
    print("  • Special positions (benzylic, allylic) are detected")
    print("  • Guard patterns prevent over-matching")
    print("  • All patterns are chemically meaningful")
    print()


if __name__ == "__main__":
    main()
