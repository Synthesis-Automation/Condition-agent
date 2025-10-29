#!/usr/bin/env python3
"""
Interactive CLI tester for the analysis module.

Tests reaction type detection and reactant classification with the Two-Pass Approach.
Users can input reaction SMILES and optionally specify the reaction type.
"""

import sys
from pathlib import Path
from typing import Optional

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.analysis.reaction_context import (
    classify_reactants_with_context,
    ReactionClassification,
)
from chemtools import detect_reaction
from chemtools.analysis._registry import get_registry


def print_header(text: str, char: str = "="):
    """Print a formatted header."""
    width = 80
    print()
    print(char * width)
    print(text.center(width))
    print(char * width)
    print()


def print_section(text: str):
    """Print a section divider."""
    print()
    print(f"--- {text} " + "-" * (70 - len(text)))
    print()


def get_reaction_smiles() -> Optional[str]:
    """Prompt user for reaction SMILES."""
    print("Enter reaction SMILES (or 'quit' to exit):")
    print("Format: reactants>>products")
    print("Example: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1")
    print()
    
    smiles = input("> ").strip()
    
    if smiles.lower() in ['quit', 'exit', 'q']:
        return None
    
    if not smiles or '>>' not in smiles:
        print("❌ Invalid SMILES format. Must contain '>>' separator.")
        return None
    
    return smiles


def list_reaction_types() -> None:
    """Display available reaction types from registry."""
    registry = get_registry()
    if not registry:
        print("⚠ Registry not loaded - cannot list reaction types")
        return
    
    reaction_types = sorted(registry.reaction_types.keys())
    print("\nAvailable Reaction Types:")
    print("-" * 80)
    
    # Define explicit category mappings for clear organization
    categories = {
        'C-C Coupling': [
            'suzuki_miyaura', 'suzuki_miyaura_in_situ', 'negishi', 'negishi_in_situ',
            'heck', 'sonogashira', 'stille', 'kumada', 'grignard_addition',
            'organolithium_addition', 'reformatsky', 'cc_coupling'
        ],
        'C-N Coupling': [
            'buchwald_hartwig_c_n', 'chan_lam', 'ullmann_cn', 'cn_coupling',
            'reductive_amination', 'cyanation'
        ],
        'C-O Coupling': [
            'mitsunobu', 'ullmann_ether', 'williamson_ether', 'co_coupling',
            'esterification', 'glycosidation'
        ],
        'C-S Coupling': [
            'cs_coupling'
        ],
        'Amide Formation': [
            'amide_coupling'
        ],
        'C-H Activation': [
            'ch_activation', 'ch_alkylation', 'borylation_c_h', 'arylation_acidic_c_h'
        ],
        'Reduction': [
            'hydrogenation', 'carbonyl_reduction', 'reduction', 'birch_reduction',
            'hydrodehalogenation'
        ],
        'Oxidation': [
            'oxidation', 'alcohol_oxidation', 'baeyer_villiger', 'epoxidation',
            'fluorination_oxidative'
        ],
        'Metathesis': [
            'cross_metathesis', 'ring_closing_metathesis'
        ],
        'Condensation': [
            'aldol_condensation', 'claisen_condensation', 'condensation',
            'horner_wadsworth_emmons', 'wittig'
        ],
        'Halogenation/Substitution': [
            'halogenation_chlorination', 'finkelstein', 'deoxyfluorination',
            'sn2_substitution', 'snar', 'sandmeyer'
        ],
        'Addition': [
            'addition', 'michael_addition', 'diels_alder', 'hydroboration',
            'hydration'
        ],
        'Other': []  # Catch-all for reactions not explicitly categorized
    }
    
    # Build reverse mapping for quick lookup
    categorized = set()
    for cat_reactions in categories.values():
        categorized.update(cat_reactions)
    
    # Add uncategorized reactions to 'Other'
    categories['Other'] = [rt for rt in reaction_types if rt not in categorized]
    
    # Display by category
    for category, reactions in categories.items():
        if reactions:
            reactions_sorted = sorted([rt for rt in reactions if rt in reaction_types])
            if reactions_sorted:
                print(f"\n  {category}:")
                for i, rt in enumerate(reactions_sorted, 1):
                    print(f"    {i:2d}. {rt}")
    
    print()


def get_reaction_type_choice() -> Optional[str]:
    """Prompt user for reaction type selection."""
    print("\nReaction Type Selection:")
    print("  1. Auto-detect (recommended)")
    print("  2. Specify reaction type manually")
    print("  3. List all available reaction types")
    print()
    
    choice = input("Select option (1-3): ").strip()
    
    if choice == '1' or not choice:
        return None  # Auto-detect
    
    elif choice == '2':
        reaction_type = input("Enter reaction type (e.g., suzuki_miyaura, buchwald_hartwig_c_n): ").strip()
        if reaction_type:
            return reaction_type
        return None
    
    elif choice == '3':
        list_reaction_types()
        reaction_type = input("\nEnter reaction type from the list above: ").strip()
        if reaction_type:
            return reaction_type
        return None
    
    else:
        print("⚠ Invalid choice, using auto-detect")
        return None


def format_confidence(confidence: float) -> str:
    """Format confidence with color indicator."""
    if confidence >= 0.8:
        return f"{confidence:.2f} (High ✓)"
    elif confidence >= 0.5:
        return f"{confidence:.2f} (Medium ~)"
    else:
        return f"{confidence:.2f} (Low ⚠)"


def display_reaction_analysis(result: ReactionClassification, smiles: str) -> None:
    """Display formatted reaction analysis results."""
    
    print_header("REACTION ANALYSIS RESULTS")
    
    # Input
    print(f"Input SMILES: {smiles}")
    print()
    
    # Reaction Type Detection
    print_section("Reaction Type Detection")
    
    print(f"Detected Family:    {result.reaction_type}")
    print(f"Confidence:         {format_confidence(result.reaction_confidence)}")
    print(f"Detection Method:   {result.detection_method}")
    
    if result.detection_method == 'user_provided':
        print("                    (User specified - confidence set to 1.0)")
    elif result.detection_method == 'ml_detected':
        print("                    (ML model: rxn-insight)")
    elif result.detection_method == 'rule_based':
        print("                    (SMARTS pattern matching)")
    
    # Reactant Classification
    print_section("Reactant Classification (Two-Pass Approach)")
    
    if not result.reactants:
        print("⚠ No reactants classified")
        print()
        print("Possible reasons:")
        print("  - SMILES parsing failed")
        print("  - No SMARTS patterns matched")
        print("  - Reactants don't match expected types for this reaction")
    else:
        print(f"Total Reactants Found: {len(result.reactants)}")
        print(f"Multi-functional Groups: {'Yes' if result.has_multi_functional else 'No'}")
        print()
        
        for i, reactant in enumerate(result.reactants, 1):
            print(f"Reactant {i}:")
            print(f"  Category:      {reactant.category}")
            print(f"  Member Type:   {reactant.member_type}")
            print(f"  Name:          {reactant.name}")
            
            if reactant.role:
                print(f"  Role:          {reactant.role}")
            
            print(f"  Position:      {reactant.position}")
            print(f"  Expected:      {'Yes ✓' if reactant.is_expected else 'No'}")
            
            if reactant.alternative_matches:
                print(f"  Alternatives:  {len(reactant.alternative_matches)} other functional groups found")
                for alt in reactant.alternative_matches[:3]:  # Show first 3
                    print(f"                 - {alt.member_type} ({alt.category})")
            
            print()


def display_detection_only(smiles: str) -> None:
    """Display only reaction type detection (no user-provided type)."""
    
    print_section("Detecting Reaction Type...")
    
    try:
        detection_result = detect_reaction(smiles, use_ml=False)
        
        detected_family = detection_result.get('family', 'Unknown')
        confidence = detection_result.get('confidence', 0.0)
        
        print(f"Detected Family:  {detected_family}")
        print(f"Confidence:       {format_confidence(confidence)}")
        
        if detected_family == 'Unknown':
            print()
            print("⚠ Could not detect reaction type")
            print("  Consider specifying the reaction type manually")
        
    except Exception as e:
        print(f"❌ Error during detection: {e}")


def run_analysis(smiles: str, reaction_type: Optional[str] = None) -> None:
    """Run the complete analysis workflow."""
    
    try:
        # If no reaction type specified, show detection first
        if not reaction_type:
            display_detection_only(smiles)
            print()
        
        # Run Two-Pass classification
        print_section("Running Two-Pass Classification...")
        
        result = classify_reactants_with_context(
            reaction_smiles=smiles,
            reaction_type=reaction_type,
            auto_detect=True
        )
        
        # Display results
        display_reaction_analysis(result, smiles)
        
    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()


def show_examples() -> None:
    """Display example reactions."""
    print_header("EXAMPLE REACTIONS")
    
    examples = [
        ("Suzuki-Miyaura Coupling",
         "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
         "suzuki_miyaura"),
        
        ("Buchwald-Hartwig C-N",
         "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
         "buchwald_hartwig_c_n"),
        
        ("Sonogashira",
         "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",
         "sonogashira"),
        
        ("Heck Reaction",
         "Brc1ccccc1.C=C>>C(=Cc1ccccc1)",
         "heck"),
        
        ("Amide Coupling",
         "c1ccc(C(=O)O)cc1.Nc1ccccc1>>c1ccc(C(=O)Nc2ccccc2)cc1",
         "amide_coupling"),
    ]
    
    for i, (name, smiles, rxn_type) in enumerate(examples, 1):
        print(f"{i}. {name}")
        print(f"   SMILES: {smiles}")
        print(f"   Type: {rxn_type}")
        print()
    
    print("Copy any SMILES above to test!")
    print()


def main():
    """Main interactive loop."""
    
    print_header("REACTION ANALYSIS MODULE - INTERACTIVE TESTER", "=")
    
    print("This tool tests:")
    print("  1. Reaction type detection (auto or manual)")
    print("  2. Reactant classification using the Two-Pass Approach")
    print("  3. Context-aware functional group analysis")
    print()
    print("Commands:")
    print("  - Type 'examples' to see example reactions")
    print("  - Type 'types' to list available reaction types")
    print("  - Type 'quit' or 'exit' to quit")
    print()
    
    # Check registry
    registry = get_registry()
    if registry:
        print(f"✓ Registry loaded: {len(registry.reaction_types)} reaction types, "
              f"{len(registry.reactant_types)} reactant types")
    else:
        print("⚠ Warning: Registry not loaded - analysis may fail")
    
    print()
    print("=" * 80)
    
    while True:
        print()
        
        # Get reaction SMILES
        smiles = get_reaction_smiles()
        
        if smiles is None:
            print("\nGoodbye!")
            break
        
        # Handle special commands
        if smiles.lower() == 'examples':
            show_examples()
            continue
        
        if smiles.lower() in ['types', 'list']:
            list_reaction_types()
            continue
        
        # Get reaction type preference
        reaction_type = get_reaction_type_choice()
        
        # Run analysis
        run_analysis(smiles, reaction_type)
        
        # Ask to continue
        print()
        print("=" * 80)
        continue_choice = input("\nAnalyze another reaction? (y/n): ").strip().lower()
        if continue_choice not in ['y', 'yes', '']:
            print("\nGoodbye!")
            break


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Goodbye!")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
