"""
Interactive CLI Tester for Condition Recommendation System

This script provides a detailed, step-by-step view of how the condition
recommendation system works. It's useful for:
- Understanding the internal workflow
- Debugging issues
- Testing different reactions
- Learning how the taxonomy system works

Usage:
    python test_recommendation_interactive.py
"""

import sys
import time
from typing import Dict, Any, Optional
from pathlib import Path

# Color codes for terminal output
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

    @classmethod
    def disable(cls):
        cls.HEADER = ''
        cls.BLUE = ''
        cls.CYAN = ''
        cls.GREEN = ''
        cls.YELLOW = ''
        cls.RED = ''
        cls.BOLD = ''
        cls.UNDERLINE = ''
        cls.END = ''


def print_step(step_num: int, title: str):
    """Print a step header."""
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}STEP {step_num}: {title}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}")


def print_substep(text: str):
    """Print a substep."""
    print(f"{Colors.BLUE}  → {text}{Colors.END}")


def print_success(text: str):
    """Print success message."""
    print(f"{Colors.GREEN}  ✓ {text}{Colors.END}")


def print_error(text: str):
    """Print error message."""
    print(f"{Colors.RED}  ✗ {text}{Colors.END}")


def print_info(key: str, value: Any):
    """Print key-value info."""
    print(f"{Colors.BOLD}  {key}:{Colors.END} {value}")


def print_warning(text: str):
    """Print warning message."""
    print(f"{Colors.YELLOW}  ⚠ {text}{Colors.END}")


def show_menu():
    """Display the main menu."""
    print(f"\n{Colors.BOLD}{Colors.HEADER}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.HEADER}Condition Recommendation System - Interactive Tester{Colors.END}")
    print(f"{Colors.BOLD}{Colors.HEADER}{'='*80}{Colors.END}")
    print("\nThis tool shows step-by-step how condition recommendation works.")
    print("\nSelect a test reaction:")
    print(f"\n{Colors.BOLD}1.{Colors.END} Suzuki-Miyaura Coupling (should find ~1900 experiments)")
    print("   Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
    print(f"\n{Colors.BOLD}2.{Colors.END} Amide Formation (should find experiments)")
    print("   CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1")
    print(f"\n{Colors.BOLD}3.{Colors.END} C-N Coupling (limited data)")
    print("   Brc1ccccc1.NCc1ccccc1>>c1ccc(NCc2ccccc2)cc1")
    print(f"\n{Colors.BOLD}4.{Colors.END} Complex Suzuki (with protecting group)")
    print("   CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1")
    print(f"\n{Colors.BOLD}5.{Colors.END} Simple Substitution (edge case)")
    print("   Brc1ccccc1.O>>Oc1ccccc1")
    print(f"\n{Colors.BOLD}6.{Colors.END} Custom input (enter your own SMILES)")
    print(f"\n{Colors.BOLD}0.{Colors.END} Exit")


def get_test_reaction(choice: str) -> Optional[str]:
    """Get reaction SMILES for the selected test."""
    reactions = {
        '1': "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        '2': "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1",
        '3': "Brc1ccccc1.NCc1ccccc1>>c1ccc(NCc2ccccc2)cc1",
        '4': "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1",
        '5': "Brc1ccccc1.O>>Oc1ccccc1",
    }
    return reactions.get(choice)


def test_recommendation_step_by_step(rxn_smiles: str, top_k: int = 5):
    """
    Run complete recommendation workflow with detailed step-by-step output.

    Args:
        rxn_smiles: Reaction SMILES string
        top_k: Number of recommendations to show
    """

    # ==========================================================================
    # STEP 1: Initialize Components
    # ==========================================================================
    print_step(1, "Initialize LLM Client and Analyzer")

    try:
        print_substep("Importing modules...")
        from llmtools.clients import LLMClient
        from reaction_agent import ReactionSMILESAnalyzer
        from reaction_agent.condition_bridge import ConditionBridge
        print_success("Modules imported")

        print_substep("Creating LLM client (gpt-4o-mini)...")
        start_time = time.time()
        client = LLMClient(provider="openai", model="gpt-4o-mini")
        elapsed = time.time() - start_time
        print_success(f"LLM client created ({elapsed:.3f}s)")
        print_info("Model", "gpt-4o-mini")
        print_info("Provider", "openai")

        print_substep("Creating ReactionSMILESAnalyzer...")
        analyzer = ReactionSMILESAnalyzer(client)
        print_success("Analyzer created")

        print_substep("Creating ConditionBridge...")
        bridge = ConditionBridge()
        print_success("Bridge created")
        print_info("HTE Database Path", bridge.hte_db_path)

    except Exception as e:
        print_error(f"Initialization failed: {e}")
        import traceback
        traceback.print_exc()
        return

    # ==========================================================================
    # STEP 2: Analyze Reaction (Tiers 1-4)
    # ==========================================================================
    print_step(2, "Analyze Reaction (Tiers 1-4)")

    print_info("Input SMILES", rxn_smiles)

    try:
        print_substep("Running three-tier analysis...")
        print("  (This may take 15-25 seconds - calling LLM APIs)")

        analysis_start = time.time()
        analysis_result = analyzer.analyze(
            rxn_smiles,
            mode="fast",  # Fast mode for testing
            validate=False  # Skip validation for speed
        )
        analysis_time = time.time() - analysis_start

        print_success(f"Analysis complete ({analysis_time:.2f}s)")

        # Show analysis results
        input_data = analysis_result.get('input', {})
        print_info("Clean SMILES", input_data.get('rxn_smiles_clean', 'N/A'))

        # Tier 1 (automatic interpretation)
        if 'auto_interpretation' in analysis_result:
            auto = analysis_result['auto_interpretation']
            interp = auto.get('interpretation', {})
            print_info("Tier 1 Complexity", interp.get('complexity', 'N/A'))
            patterns = interp.get('patterns_detected', [])
            if patterns:
                print_info("Tier 1 Patterns", ', '.join(patterns))

        # Tier 2 (quick glance)
        if 'quick_glance' in analysis_result:
            qg = analysis_result['quick_glance']
            if qg.get('success'):
                rxn_types = qg.get('reaction_types', [])
                print_info("Tier 2 Reaction Types", ', '.join(rxn_types[:3]) if rxn_types else 'None')
                print_info("Tier 2 Confidence", f"{qg.get('confidence', 0):.2f}")
                print_info("Tier 2 Complexity", qg.get('complexity', 'N/A'))

        # Tier 3 (interpretation)
        if 'interpretation' in analysis_result:
            interp = analysis_result['interpretation']
            print_info("Tier 3 Class", interp.get('overall_class', 'N/A'))
            print_info("Tier 3 Confidence", f"{interp.get('confidence', 0):.2f}")

        # Metadata
        meta = analysis_result.get('metadata', {})
        print_info("Model Used", meta.get('model', 'N/A'))
        print_info("Total Tokens", meta.get('total_tokens', 'N/A'))
        print_info("Latency", f"{meta.get('latency_ms', 0)/1000:.2f}s")

    except Exception as e:
        print_error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return

    # ==========================================================================
    # STEP 3: Extract SMILES Components
    # ==========================================================================
    print_step(3, "Extract SMILES Components")

    try:
        print_substep("Parsing reaction SMILES...")
        reactant_a, reactant_b, product = bridge.extract_smiles(analysis_result)

        if not reactant_a:
            print_error("Failed to extract reactants")
            return

        print_success("SMILES components extracted")
        print_info("Reactant A", reactant_a)
        print_info("Reactant B", reactant_b or "(none)")
        print_info("Product", product or "(none)")

        if not product:
            print_warning("No product SMILES found!")
            print_warning("HTERecommender requires product for taxonomy detection")
            print_warning("Recommendation will likely fail")

    except Exception as e:
        print_error(f"SMILES extraction failed: {e}")
        import traceback
        traceback.print_exc()
        return

    # ==========================================================================
    # STEP 4: Call HTERecommender (Taxonomy Detection)
    # ==========================================================================
    print_step(4, "Query HTE Database (Taxonomy Detection)")

    try:
        print_substep("Calling HTERecommender.recommend()...")
        print("  This will internally:")
        print("    1. Construct full reaction SMILES")
        print("    2. Call featurize_reaction() for taxonomy detection")
        print("    3. Detect reaction type (e.g., Suzuki_miyaura)")
        print("    4. Extract reacted/formed motifs (e.g., Ar-Br, Ar-B(OH)2)")
        print("    5. Match against HTE database")
        print("    6. Rank by Z-score")

        rec_start = time.time()
        recommendations = bridge.recommend_conditions(
            analysis_result,
            top_k=top_k,
            min_experiments=2
        )
        rec_time = time.time() - rec_start

        print_success(f"Recommendation complete ({rec_time:.3f}s)")

        # Show taxonomy detection results
        print(f"\n{Colors.BOLD}  Taxonomy Detection Results:{Colors.END}")
        print_info("  Detected Reaction Type", recommendations.predicted_reaction_type or "Unknown")
        print_info("  Confidence", f"{recommendations.reaction_type_confidence:.2f}")

        if recommendations.reacted_motifs:
            print_info("  Reacted Motifs", ', '.join(recommendations.reacted_motifs))
        if recommendations.formed_motifs:
            print_info("  Formed Motifs", ', '.join(recommendations.formed_motifs))
        if recommendations.spectator_motifs:
            print_info("  Spectator Motifs", ', '.join(recommendations.spectator_motifs))

        if recommendations.query_reaction_key:
            print_info("  Reaction Key", recommendations.query_reaction_key[:80] + "...")

        print_info("  Total Experiments", recommendations.total_matching_experiments)
        print_info("  Recommendations Found", len(recommendations.recommendations))

    except Exception as e:
        print_error(f"Recommendation query failed: {e}")
        import traceback
        traceback.print_exc()
        return

    # ==========================================================================
    # STEP 5: Display Recommendations
    # ==========================================================================
    print_step(5, f"Display Top {top_k} Recommendations")

    if not recommendations.recommendations:
        print_warning("No recommendations found")
        print("\n  Possible reasons:")
        print("    - Reaction type not in HTE experiments database")
        print("    - Reactant motifs not recognized by taxonomy")
        print("    - No experimental data for this reaction class")
        return

    print(f"\n{Colors.BOLD}  Showing top {min(top_k, len(recommendations.recommendations))} conditions:{Colors.END}\n")

    for i, rec in enumerate(recommendations.recommendations[:top_k], 1):
        # Color code by Z-score
        if rec.avg_z_score > 1.0:
            score_color = Colors.GREEN
            score_label = "Excellent"
        elif rec.avg_z_score > 0.0:
            score_color = Colors.YELLOW
            score_label = "Good"
        else:
            score_color = Colors.RED
            score_label = "Poor"

        print(f"  {Colors.BOLD}Rank {i}:{Colors.END}")
        print(f"    {score_color}Z-Score: {rec.avg_z_score:.2f} ({score_label}){Colors.END}")
        print(f"    Confidence: {rec.confidence_score:.0f}%")
        print(f"    Experiments: {rec.num_experiments}")

        print(f"\n    {Colors.BOLD}Conditions:{Colors.END}")
        if rec.catalyst:
            print(f"      Catalyst:  {rec.catalyst}")
        if rec.ligand:
            print(f"      Ligand:    {rec.ligand}")
        if rec.base:
            print(f"      Base:      {rec.base}")
        if rec.solvent:
            print(f"      Solvent:   {rec.solvent}")
        if rec.additive:
            print(f"      Additive:  {rec.additive}")

        print(f"\n    {Colors.BOLD}Performance:{Colors.END}")
        print(f"      Success rate: {rec.success_rate:.1f}%")
        print(f"      Avg yield:    {rec.avg_yield:.1f}%")
        print(f"      Median yield: {rec.median_yield:.1f}%")

        print()  # Blank line between recommendations

    # ==========================================================================
    # STEP 6: Summary
    # ==========================================================================
    print_step(6, "Summary")

    total_time = analysis_time + rec_time
    print_info("Total Time", f"{total_time:.2f}s")
    print_info("  Analysis Time", f"{analysis_time:.2f}s ({analysis_time/total_time*100:.1f}%)")
    print_info("  Recommendation Time", f"{rec_time:.3f}s ({rec_time/total_time*100:.1f}%)")

    if recommendations.recommendations:
        print_success(f"Successfully found {len(recommendations.recommendations)} condition sets")
        print(f"\n{Colors.BOLD}Key Insights:{Colors.END}")
        print(f"  • Reaction type auto-detected as '{recommendations.predicted_reaction_type}'")
        print(f"  • Confidence: {recommendations.reaction_type_confidence:.0%}")
        print(f"  • Matched {recommendations.total_matching_experiments} experiments from database")
        print(f"  • Recommendation adds <1% overhead to analysis time")
    else:
        print_warning("No recommendations found")
        print("\nThis is expected for:")
        print("  • Unusual reaction types not in HTE database")
        print("  • Novel transformations without experimental precedent")
        print("  • Reactions with unrecognized substrate motifs")


def main():
    """Main interactive loop."""

    print(f"\n{Colors.CYAN}Starting interactive tester...{Colors.END}")

    # Check for API key
    import os
    if not os.getenv('OPENAI_API_KEY'):
        print(f"\n{Colors.RED}ERROR: OPENAI_API_KEY environment variable not set{Colors.END}")
        print("\nPlease set it:")
        print("  export OPENAI_API_KEY='your-key-here'")
        sys.exit(1)

    while True:
        show_menu()

        choice = input(f"\n{Colors.BOLD}Enter your choice (0-6): {Colors.END}").strip()

        if choice == '0':
            print(f"\n{Colors.CYAN}Goodbye!{Colors.END}")
            break

        elif choice == '6':
            # Custom input
            print(f"\n{Colors.BOLD}Enter custom reaction SMILES:{Colors.END}")
            print("Format: reactants>>products")
            print("Example: Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
            rxn_smiles = input(f"{Colors.BOLD}SMILES: {Colors.END}").strip()

            if not rxn_smiles:
                print(f"{Colors.YELLOW}No input provided, returning to menu{Colors.END}")
                continue

            if '>>' not in rxn_smiles:
                print(f"{Colors.RED}Invalid SMILES: missing '>>' separator{Colors.END}")
                continue

        elif choice in ['1', '2', '3', '4', '5']:
            rxn_smiles = get_test_reaction(choice)
        else:
            print(f"{Colors.RED}Invalid choice. Please enter 0-6.{Colors.END}")
            continue

        # Ask for number of recommendations
        top_k_input = input(f"\n{Colors.BOLD}Number of recommendations to show (default: 5): {Colors.END}").strip()
        try:
            top_k = int(top_k_input) if top_k_input else 5
            top_k = max(1, min(top_k, 20))  # Clamp to 1-20
        except ValueError:
            top_k = 5

        # Run the test
        test_recommendation_step_by_step(rxn_smiles, top_k)

        # Ask to continue
        continue_input = input(f"\n{Colors.BOLD}Test another reaction? (y/n): {Colors.END}").strip().lower()
        if continue_input not in ['y', 'yes', '']:
            print(f"\n{Colors.CYAN}Goodbye!{Colors.END}")
            break


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n\n{Colors.CYAN}Interrupted by user. Goodbye!{Colors.END}")
        sys.exit(0)
    except Exception as e:
        print(f"\n{Colors.RED}Unexpected error: {e}{Colors.END}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
