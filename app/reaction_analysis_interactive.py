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
from chemtools.analysis.reactants import (
    classify_reactant_smiles,
    required_reactant_categories,
    iter_reactant_matches,
)


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
            'sn2_substitution', 'snar_cn', 'snar_co', 'snar_cs', 'sandmeyer'
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


def get_reaction_type_choice() -> tuple[Optional[str], bool]:
    """
    Prompt user for reaction type selection.
    
    Returns:
        (reaction_type, use_ml) tuple
        - reaction_type: Specified type or None for auto-detect
        - use_ml: Whether to use ML-enhanced detection
    """
    print("\nReaction Type Selection:")
    print("  1. Auto-detect with ML (recommended)")
    print("  2. Auto-detect (rule-based only)")
    print("  3. Specify reaction type manually")
    print("  4. List all available reaction types")
    print()
    
    choice = input("Select option (1-4): ").strip()
    
    if choice == '1' or not choice:
        return None, True  # Auto-detect with ML
    
    elif choice == '2':
        return None, False  # Auto-detect rule-based only
    
    elif choice == '3':
        reaction_type = input("Enter reaction type (e.g., suzuki_miyaura, buchwald_hartwig_c_n): ").strip()
        if reaction_type:
            return reaction_type, True
        return None, True
    
    elif choice == '4':
        list_reaction_types()
        reaction_type = input("\nEnter reaction type from the list above: ").strip()
        if reaction_type:
            return reaction_type, True
        return None, True
    
    else:
        print("⚠ Invalid choice, using auto-detect with ML")
        return None, True


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


def _display_reactant_classification(family: str, details: dict) -> None:
    """
    Display reaction-specific reactant classification.
    
    Args:
        family: Detected reaction family
        details: Detection details containing reactants
    """
    if family == 'Unknown':
        return
    
    # Map internal family names to taxonomy reaction IDs
    reaction_id_map = {
        'suzuki_miyaura': 'Suzuki-Miyaura',
        'buchwald_hartwig_c_n': 'Buchwald-Hartwig-C-N',
        'ullmann_cn': 'Ullmann-CN',
        'ullmann_co': 'Ullmann-CO',
        'sonogashira': 'Sonogashira',
        'negishi': 'Negishi',
        'heck': 'Heck',
        'stille': 'Stille',
        'cn_coupling': 'Buchwald-Hartwig-C-N',  # Generic C-N coupling
        'co_coupling': 'Ullmann-CO',  # Generic C-O coupling
        'aldol_condensation': 'Aldol',
        'wittig': 'Wittig',
        'grignard': 'Grignard',
    }
    
    reaction_id = reaction_id_map.get(family)
    if not reaction_id:
        return
    
    # Get expected reactant types for this reaction
    expected_types = required_reactant_categories(reaction_id)
    if not expected_types:
        return
    
    reactants = details.get('reactants', [])
    if not reactants:
        return
    
    print()
    print(f"  Reactant Classification for {reaction_id}:")
    print(f"  {'-' * 38}")
    
    # Show expected types
    print(f"  Expected: {' + '.join(['/'.join(group) for group in expected_types])}")
    print()
    
    # Classify each reactant
    for i, reactant_smiles in enumerate(reactants, 1):
        match = classify_reactant_smiles(reactant_smiles)
        if match:
            # Check if it matches expected types
            is_expected = any(match.category in group for group in expected_types)
            status = "✓" if is_expected else "⚠"
            
            print(f"  Reactant {i}: {reactant_smiles}")
            print(f"    {status} Type: {match.category} ({match.member_type})")
            print(f"      Name: {match.name}")
            
            # Show if it matches other categories too
            all_matches = iter_reactant_matches(reactant_smiles)
            if len(all_matches) > 1:
                other_categories = [m.category for m in all_matches[1:4] if m.category != match.category]
                if other_categories:
                    print(f"      Also: {', '.join(other_categories)}")
        else:
            print(f"  Reactant {i}: {reactant_smiles}")
            print(f"    ⚠ No classification match")
        print()


def display_detection_only(smiles: str, use_ml: bool = True) -> Optional[dict]:
    """
    Display reaction type detection with detailed breakdown.
    
    Args:
        smiles: Reaction SMILES
        use_ml: Whether to use ML-enhanced detection (default: True)
    
    Returns:
        Detection result dict or None if error
    """
    
    print_section("Detecting Reaction Type...")
    
    try:
        # Run unified detection
        detection_result = detect_reaction(smiles, use_ml=use_ml)
        
        detected_family = detection_result.get('family', 'Unknown')
        confidence = detection_result.get('confidence', 0.0)
        method = detection_result.get('method', 'unknown')
        details = detection_result.get('details', {})
        
        # Main detection result
        print(f"Detected Family:  {detected_family}")
        print(f"Confidence:       {format_confidence(confidence)}")
        print(f"Detection Method: {method}")
        
        # Show detection breakdown if available
        if details:
            print()
            print("Detection Breakdown:")
            print("-" * 40)
            
            # Rule-based prediction
            rule_pred = details.get('rule_prediction', {})
            if rule_pred:
                print(f"  Rule-based:  {rule_pred.get('family', 'N/A')} "
                      f"(conf: {rule_pred.get('confidence', 0):.2f})")
            
            # ML prediction (if ML was used)
            ml_pred = details.get('ml_prediction', {})
            if ml_pred and use_ml:
                ml_family = ml_pred.get('family', 'N/A')
                ml_conf = ml_pred.get('confidence')
                ml_name = ml_pred.get('rxn_name', '')
                ml_class = ml_pred.get('rxn_class', '')
                
                if ml_family and ml_family != 'N/A':
                    conf_str = f"{ml_conf:.2f}" if isinstance(ml_conf, (int, float)) else "N/A"
                    print(f"  ML-based:    {ml_family} (conf: {conf_str})")
                    if ml_name:
                        print(f"               ML name: \"{ml_name}\"")
                    if ml_class:
                        print(f"               ML class: \"{ml_class}\"")
                else:
                    print(f"  ML-based:    Not available or failed")
            
            # Functional groups detected
            fg = details.get('functional_groups', {})
            if fg:
                detected_groups = [k for k, v in fg.items() if v]
                if detected_groups:
                    print()
                    print(f"  Functional Groups Detected ({len(detected_groups)}):")
                    for group in sorted(detected_groups)[:8]:  # Show first 8
                        print(f"    • {group}")
                    if len(detected_groups) > 8:
                        print(f"    ... and {len(detected_groups) - 8} more")
            
            # Catalysts detected
            catalysts = details.get('catalysts', [])
            if catalysts:
                print()
                print(f"  Catalysts Detected: {', '.join(catalysts)}")
            
            # Reaction-specific reactant analysis
            _display_reactant_classification(detected_family, details)
        
        # Agreement status
        agreement = detection_result.get('agreement')
        status = detection_result.get('status', '')
        if agreement is not None and use_ml:
            print()
            if agreement:
                print("  ✓ Rule-based and ML predictions agree")
            else:
                if status == 'conflict':
                    print("  ⚠ Conflict: Rule-based and ML predictions differ")
                    print("    Using ML prediction (higher confidence)")
                elif status == 'rule_only':
                    print("  ℹ ML detection not available - using rule-based only")
        
        if detected_family == 'Unknown':
            print()
            print("⚠ Could not detect reaction type")
            print("  Suggestions:")
            print("    • Check if reaction SMILES is valid")
            print("    • Try specifying the reaction type manually")
            print("    • Reaction may not be in supported families")
        
        return detection_result
        
    except Exception as e:
        print(f"❌ Error during detection: {e}")
        import traceback
        traceback.print_exc()
        return None


def run_analysis(smiles: str, reaction_type: Optional[str] = None, use_ml: bool = True) -> None:
    """
    Run the complete analysis workflow.
    
    Args:
        smiles: Reaction SMILES
        reaction_type: Optional user-specified reaction type
        use_ml: Whether to use ML-enhanced detection (default: True)
    """
    
    try:
        # If no reaction type specified, show detailed detection first
        detection_result = None
        if not reaction_type:
            detection_result = display_detection_only(smiles, use_ml=use_ml)
            print()
            
            # Offer to use detected type or specify manually
            if detection_result and detection_result.get('family') != 'Unknown':
                print("Options:")
                print(f"  1. Use detected type: {detection_result['family']}")
                print("  2. Specify different reaction type")
                print("  3. Skip classification (detection only)")
                print()
                
                opt = input("Select option (1-3, default=1): ").strip()
                
                if opt == '2':
                    reaction_type = input("Enter reaction type: ").strip()
                    if not reaction_type:
                        print("⚠ No type entered, using detected type")
                elif opt == '3':
                    print("\n✓ Detection complete. Skipping classification.")
                    return
                # opt == '1' or default: use detected type (reaction_type remains None)
        
        # Run Two-Pass classification
        print_section("Running Two-Pass Classification...")
        
        result = classify_reactants_with_context(
            reaction_smiles=smiles,
            reaction_type=reaction_type,
            auto_detect=True
        )
        
        # Display results
        display_reaction_analysis(result, smiles)
        
        # Show additional insights
        if detection_result and result.reactants:
            print_section("Analysis Summary")
            
            total_reactants = len(result.reactants)
            expected_reactants = sum(1 for r in result.reactants if r.is_expected)
            
            print(f"Reactants classified:     {total_reactants}")
            print(f"Expected reactants:       {expected_reactants}")
            print(f"Unexpected reactants:     {total_reactants - expected_reactants}")
            print(f"Detection confidence:     {format_confidence(result.reaction_confidence)}")
            print(f"Multi-functional groups:  {'Yes' if result.has_multi_functional else 'No'}")
            
            if result.reaction_confidence < 0.5:
                print()
                print("⚠ Low confidence detection - results may be unreliable")
                print("  Consider verifying the reaction type manually")
        
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
    
    print("This tool tests the unified detection system:")
    print("  1. Reaction type detection (unified detect_reaction API)")
    print("     • ML-enhanced detection (rxn-insight integration)")
    print("     • Rule-based detection (SMARTS pattern matching)")
    print("     • Catalyst-aware refinements (Pd/Cu/Ni detection)")
    print()
    print("  2. Reactant classification using the Two-Pass Approach")
    print("     • Context-aware functional group analysis")
    print("     • Multi-functional group detection")
    print("     • Role-based reactant categorization")
    print()
    print("Commands:")
    print("  - Type 'examples' to see example reactions")
    print("  - Type 'types' to list available reaction types")
    print("  - Type 'quit' or 'exit' to quit")
    print()
    
    # Check ML availability
    try:
        from chemtools._ml_helpers import is_available as ml_available
        ml_status = "✓ Available" if ml_available() else "✗ Not installed"
    except:
        ml_status = "✗ Not available"
    
    print(f"ML Detection (rxn-insight): {ml_status}")
    
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
        reaction_type, use_ml = get_reaction_type_choice()
        
        # Run analysis
        run_analysis(smiles, reaction_type, use_ml=use_ml)
        
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
