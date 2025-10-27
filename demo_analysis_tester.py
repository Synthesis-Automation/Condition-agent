#!/usr/bin/env python3
"""
Demo script showing the analysis module interactive tester in action.

This demonstrates the output without requiring interactive input.
"""

from test_analysis_interactive import (
    run_analysis,
    show_examples,
    list_reaction_types,
    print_header,
)


def demo():
    """Run a demonstration of the analysis module."""
    
    print_header("ANALYSIS MODULE TESTER - DEMONSTRATION", "=")
    
    print("This demo shows the interactive tester output for several reactions.")
    print("To use the interactive version, run: python test_analysis_interactive.py")
    print()
    
    # Show examples first
    show_examples()
    
    print("=" * 80)
    print()
    
    # Demo 1: Auto-detect Suzuki
    print_header("DEMO 1: Auto-Detect Suzuki-Miyaura Coupling")
    smiles1 = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    run_analysis(smiles1, reaction_type=None)
    
    print("\n" + "=" * 80 + "\n")
    input("Press Enter to continue to next demo...")
    
    # Demo 2: User-specified Buchwald-Hartwig
    print_header("DEMO 2: User-Specified Buchwald-Hartwig C-N")
    smiles2 = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    run_analysis(smiles2, reaction_type="buchwald_hartwig_c_n")
    
    print("\n" + "=" * 80 + "\n")
    input("Press Enter to continue to next demo...")
    
    # Demo 3: Auto-detect Sonogashira
    print_header("DEMO 3: Auto-Detect Sonogashira")
    smiles3 = "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"
    run_analysis(smiles3, reaction_type=None)
    
    print("\n" + "=" * 80 + "\n")
    input("Press Enter to continue to next demo...")
    
    # Demo 4: Multi-functional substrate
    print_header("DEMO 4: Multi-Functional Substrate (ArBr with ArNH2)")
    smiles4 = "Brc1ccc(N)cc1.c1ccc(B(O)O)cc1>>Nc1ccc(-c2ccccc2)cc1"
    print("This demonstrates context-aware classification:")
    print("The substrate has BOTH ArBr and ArNH2 functional groups.")
    print("In Suzuki coupling, ArBr is reactive (electrophile), ArNH2 is not.")
    print()
    run_analysis(smiles4, reaction_type="suzuki_miyaura")
    
    print("\n" + "=" * 80 + "\n")
    
    # Show how to list reaction types
    print_header("AVAILABLE REACTION TYPES")
    print("The interactive tester can show all available reaction types:")
    print()
    list_reaction_types()
    
    print("\n" + "=" * 80)
    print()
    print("✓ Demo complete!")
    print()
    print("To use the interactive version:")
    print("  python test_analysis_interactive.py")
    print()


if __name__ == "__main__":
    try:
        demo()
    except KeyboardInterrupt:
        print("\n\nDemo interrupted. Goodbye!")
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
