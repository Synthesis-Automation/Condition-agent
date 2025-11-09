"""
Test script for the image generator with sample reactions.

This script tests the rendering.py module's capability to generate
images for both molecules and reactions using examples from sample_reactions.py.
"""

import sys
from pathlib import Path

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent))

from chemtools.visualization.rendering import (
    render_molecule_image,
    render_reaction_image,
    rdkit_available,
)

# Import from tests directory
try:
    from tests.sample_reactions import SAMPLE_REACTIONS
except ImportError:
    # Fallback: define a small subset for testing
    SAMPLE_REACTIONS = [
        "Select a sample reaction...",
        "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 (Suzuki - Simple Ph-Ph)",
        "Clc1ccc(C#N)cc1.c1ccc(B(O)O)cc1>>N#Cc1ccc(-c2ccccc2)cc1 (Suzuki - Electron-poor ArCl)",
        "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (Buchwald-Hartwig - Diphenylamine)",
        "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1 (Sonogashira - Diphenylacetylene)",
        "Brc1ccccc1.C=C>>C(=Cc1ccccc1) (Heck - Simple styrene)",
    ]


def test_molecule_rendering():
    """Test molecule rendering with various SMILES."""
    print("\n" + "="*70)
    print("TESTING MOLECULE RENDERING")
    print("="*70)
    
    # Test molecules
    test_molecules = [
        ("c1ccccc1", "benzene"),
        ("CCO", "ethanol"),
        ("c1ccc(Br)cc1", "bromobenzene"),
        ("Nc1ccccc1", "aniline"),
        ("c1ccc(C(F)(F)F)cc1", "trifluorotoluene"),
        ("c1ccncc1", "pyridine"),
        ("c1ccc2[nH]ccc2c1", "indole"),
    ]
    
    output_dir = Path("test_images/molecules")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    success_count = 0
    fail_count = 0
    
    for smiles, name in test_molecules:
        try:
            output_path = output_dir / f"{name}.png"
            result = render_molecule_image(
                smiles=smiles,
                output_path=output_path,
                size=(400, 300),
                legend=name
            )
            print(f"✓ {name:20s} → {result}")
            success_count += 1
        except Exception as e:
            print(f"✗ {name:20s} → Error: {e}")
            fail_count += 1
    
    print(f"\nMolecule rendering: {success_count} passed, {fail_count} failed")
    return success_count, fail_count


def test_reaction_rendering():
    """Test reaction rendering with sample reactions."""
    print("\n" + "="*70)
    print("TESTING REACTION RENDERING")
    print("="*70)
    
    output_dir = Path("test_images/reactions")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Select diverse test reactions from different categories
    test_reactions = []
    
    # Extract SMILES from sample reactions (first 20 actual reactions, skip the label)
    for i, rxn_str in enumerate(SAMPLE_REACTIONS[1:21], 1):
        # Parse reaction string: "SMILES (Description)"
        if ">>" in rxn_str:
            smiles_part = rxn_str.split("(")[0].strip()
            desc_part = rxn_str.split("(")[1].rstrip(")") if "(" in rxn_str else f"reaction_{i}"
            test_reactions.append((smiles_part, desc_part, i))
    
    success_count = 0
    fail_count = 0
    
    for smiles, description, idx in test_reactions:
        try:
            # Create safe filename from description
            safe_name = "".join(c if c.isalnum() or c in (' ', '-', '_') else '_' 
                               for c in description[:50])
            safe_name = safe_name.replace(' ', '_')
            output_path = output_dir / f"rxn_{idx:02d}_{safe_name}.png"
            
            result = render_reaction_image(
                reaction_smiles=smiles,
                output_path=output_path,
                size=(960, 320)
            )
            print(f"✓ Reaction {idx:2d}: {description[:60]:60s} → {result.name}")
            success_count += 1
        except Exception as e:
            print(f"✗ Reaction {idx:2d}: {description[:60]:60s} → Error: {str(e)[:50]}")
            fail_count += 1
    
    print(f"\nReaction rendering: {success_count} passed, {fail_count} failed")
    return success_count, fail_count


def test_svg_rendering():
    """Test SVG output format."""
    print("\n" + "="*70)
    print("TESTING SVG FORMAT")
    print("="*70)
    
    output_dir = Path("test_images/svg")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Test molecule SVG
    try:
        mol_path = output_dir / "benzene.svg"
        result = render_molecule_image(
            smiles="c1ccccc1",
            output_path=mol_path,
            size=(400, 300),
            image_format="svg",
            legend="Benzene (SVG)"
        )
        print(f"✓ Molecule SVG: {result}")
        mol_success = True
    except Exception as e:
        print(f"✗ Molecule SVG: Error: {e}")
        mol_success = False
    
    # Test reaction SVG
    try:
        rxn_path = output_dir / "suzuki_reaction.svg"
        result = render_reaction_image(
            reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
            output_path=rxn_path,
            size=(960, 320),
            image_format="svg"
        )
        print(f"✓ Reaction SVG: {result}")
        rxn_success = True
    except Exception as e:
        print(f"✗ Reaction SVG: Error: {e}")
        rxn_success = False
    
    return 2 if (mol_success and rxn_success) else (1 if (mol_success or rxn_success) else 0), \
           0 if (mol_success and rxn_success) else (1 if (mol_success or rxn_success) else 2)


def test_error_handling():
    """Test error handling with invalid inputs."""
    print("\n" + "="*70)
    print("TESTING ERROR HANDLING")
    print("="*70)
    
    output_dir = Path("test_images/errors")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    test_cases = [
        ("invalid_smiles", "InvalidSMILES", "molecule", "Invalid SMILES should raise error"),
        ("", "", "molecule", "Empty SMILES should raise error"),
        ("c1ccccc1>>invalid", "c1ccccc1>>invalid", "reaction", "Invalid product SMILES should raise error"),
        (">>c1ccccc1", ">>c1ccccc1", "reaction", "Missing reactants should raise error"),
    ]
    
    expected_errors = 0
    caught_errors = 0
    
    for name, smiles, render_type, desc in test_cases:
        expected_errors += 1
        try:
            output_path = output_dir / f"{name}.png"
            if render_type == "molecule":
                render_molecule_image(smiles=smiles, output_path=output_path)
            else:
                render_reaction_image(reaction_smiles=smiles, output_path=output_path)
            print(f"✗ {desc}: Did not raise error (unexpected)")
        except Exception as e:
            print(f"✓ {desc}: Correctly raised {type(e).__name__}")
            caught_errors += 1
    
    print(f"\nError handling: {caught_errors}/{expected_errors} errors correctly caught")
    return caught_errors, expected_errors - caught_errors


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("IMAGE GENERATOR TEST SUITE")
    print("="*70)
    
    # Check RDKit availability
    if not rdkit_available():
        print("\n❌ RDKit is not available. Cannot run tests.")
        print("Install RDKit to use image generation features:")
        print("  conda install -c conda-forge rdkit")
        return 1
    
    print("\n✓ RDKit is available")
    
    # Run all test suites
    results = []
    
    results.append(("Molecules", *test_molecule_rendering()))
    results.append(("Reactions", *test_reaction_rendering()))
    results.append(("SVG Format", *test_svg_rendering()))
    results.append(("Error Handling", *test_error_handling()))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    total_pass = 0
    total_fail = 0
    
    for name, passed, failed in results:
        total = passed + failed
        status = "✓ PASS" if failed == 0 else "⚠ PARTIAL" if passed > 0 else "✗ FAIL"
        print(f"{status:10s} {name:20s}: {passed}/{total} tests passed")
        total_pass += passed
        total_fail += failed
    
    print("-" * 70)
    print(f"{'OVERALL':10s} {'':20s}: {total_pass}/{total_pass + total_fail} tests passed")
    
    if total_fail == 0:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n⚠️  {total_fail} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
