#!/usr/bin/env python3
"""
Simple demo of the analysis module showing multiple reaction types.

Run this to see automated analysis without interactive prompts.
"""

from test_analysis_interactive import run_analysis, print_header


def simple_demo():
    """Run automated demo of several reactions."""
    
    print_header("REACTION ANALYSIS MODULE - AUTOMATED DEMO", "=")
    
    print("Demonstrating reaction type detection and reactant classification")
    print("for common organic reactions using the Two-Pass Approach.")
    print()
    
    reactions = [
        ("Suzuki-Miyaura (auto-detect)", 
         "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1", 
         None),
        
        ("Buchwald-Hartwig (user-specified)", 
         "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1", 
         "buchwald_hartwig_c_n"),
        
        ("Sonogashira (auto-detect)", 
         "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1", 
         None),
        
        ("Multi-functional substrate (Suzuki)", 
         "Brc1ccc(N)cc1.c1ccc(B(O)O)cc1>>Nc1ccc(-c2ccccc2)cc1", 
         "suzuki_miyaura"),
    ]
    
    for i, (name, smiles, rxn_type) in enumerate(reactions, 1):
        print("=" * 80)
        print(f"\nDEMO {i}: {name}")
        print("=" * 80)
        
        run_analysis(smiles, rxn_type)
        
        if i < len(reactions):
            print("\n" + "─" * 80 + "\n")
    
    print("\n" + "=" * 80)
    print("\n✓ Demo complete!")
    print("\nFor interactive testing, run: python test_analysis_interactive.py")
    print()


if __name__ == "__main__":
    try:
        simple_demo()
    except KeyboardInterrupt:
        print("\n\nDemo interrupted. Goodbye!")
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
