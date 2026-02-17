"""
Interactive test CLI for the three-tier reaction interpretation system.

Allows you to test reactions and see all three tiers of analysis side-by-side.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from reaction_agent import ReactionSMILESAnalyzer
from llmtools.clients import LLMClient
import json


class Colors:
    """ANSI color codes for terminal output."""
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    END = '\033[0m'


def print_header(text: str):
    """Print a formatted header."""
    print(f"\n{Colors.BOLD}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}  {text}{Colors.END}")
    print(f"{Colors.BOLD}{'='*80}{Colors.END}\n")


def print_tier_comparison(result: dict):
    """Print all three tiers side-by-side for comparison."""

    print(f"\n{Colors.BOLD}{'─'*80}{Colors.END}")
    print(f"{Colors.BOLD}THREE-TIER ANALYSIS COMPARISON{Colors.END}")
    print(f"{Colors.BOLD}{'─'*80}{Colors.END}\n")

    # Extract tier results
    tier1 = result.get('auto_interpretation', {}).get('interpretation', {})
    tier2 = result.get('quick_glance', {})
    tier3 = result.get('interpretation', {})

    # Tier 1
    print(f"{Colors.CYAN}🔹 TIER 1 - String Patterns (Free, <0.1s){Colors.END}")
    print(f"  Complexity: {tier1.get('reaction_complexity', 'N/A')}")
    reaction_types = tier1.get('likely_reaction_types', [])
    if reaction_types:
        print(f"  Reaction types: {', '.join(reaction_types)}")
    else:
        print(f"  Reaction types: {Colors.YELLOW}None detected{Colors.END}")
    patterns = tier1.get('patterns_detected', [])
    if patterns:
        print(f"  Patterns: {', '.join(patterns[:3])}")
        if len(patterns) > 3:
            print(f"           ... and {len(patterns)-3} more")
    print()

    # Tier 2
    print(f"{Colors.CYAN}🔹 TIER 2 - Quick LLM Glance (gpt-4o, ~1.5s, ~$0.0015){Colors.END}")
    if tier2 and tier2.get('success'):
        print(f"  Summary: {tier2.get('summary', 'N/A')}")
        t2_types = tier2.get('reaction_types', [])
        if t2_types:
            print(f"  Reaction types: {', '.join(t2_types)}")
        print(f"  Complexity: {tier2.get('complexity', 'N/A')}")
        print(f"  Confidence: {tier2.get('confidence', 0.0):.2f}")
        if 'metadata' in tier2:
            print(f"  Time: {tier2['metadata'].get('latency_ms', 0):.0f}ms")
    else:
        print(f"  {Colors.YELLOW}Not triggered{Colors.END}")
    print()

    # Tier 3
    print(f"{Colors.CYAN}🔹 TIER 3 - Deep LLM Analysis (~5-30s, ~$0.01-0.10){Colors.END}")
    print(f"  Reaction class: {tier3.get('overall_class', 'N/A')}")
    tags = tier3.get('tags', [])
    if tags:
        print(f"  Tags: {', '.join(tags)}")
    print(f"  Confidence: {tier3.get('confidence', 0.0):.2f}")
    events = tier3.get('events', [])
    if events:
        print(f"  Mechanistic events: {len(events)}")
        for event in events[:2]:
            print(f"    • {event.get('event_type', 'N/A')}")
        if len(events) > 2:
            print(f"    ... and {len(events)-2} more")
    print()

    # Comparison
    print(f"{Colors.BOLD}Comparison:{Colors.END}")

    # Check agreement
    t1_types = set(t.lower() for t in tier1.get('likely_reaction_types', []))
    t2_types = set(t.lower() for t in tier2.get('reaction_types', [])) if tier2.get('success') else set()
    t3_class = tier3.get('overall_class', '').lower()

    if t1_types:
        print(f"  Tier 1 vs Tier 2: ", end="")
        if t2_types and any(t1 in t2_item or t2_item in str(t1) for t1 in t1_types for t2_item in t2_types):
            print(f"{Colors.GREEN}✓ Agreement{Colors.END}")
        elif not t2_types:
            print(f"{Colors.YELLOW}⚠ Tier 2 not triggered{Colors.END}")
        else:
            print(f"{Colors.RED}✗ Disagreement{Colors.END}")
            print(f"    T1: {', '.join(t1_types)}")
            print(f"    T2: {', '.join(t2_types)}")

    if t1_types and t3_class:
        print(f"  Tier 1 vs Tier 3: ", end="")
        if any(t1 in t3_class or t3_class in t1 for t1 in t1_types):
            print(f"{Colors.GREEN}✓ Agreement{Colors.END}")
        else:
            print(f"{Colors.YELLOW}⚠ Different perspective{Colors.END}")

    print(f"\n{Colors.BOLD}{'─'*80}{Colors.END}\n")


def analyze_reaction(analyzer: ReactionSMILESAnalyzer, rxn_smiles: str):
    """Analyze a reaction and display results."""

    print(f"\n{Colors.BOLD}Analyzing:{Colors.END} {rxn_smiles[:80]}...")
    if len(rxn_smiles) > 80:
        print(f"           ...{rxn_smiles[80:]}")
    print()

    try:
        result = analyzer.analyze(rxn_smiles)
        print_tier_comparison(result)

        return result

    except Exception as e:
        print(f"\n{Colors.RED}❌ Analysis failed: {e}{Colors.END}\n")
        import traceback
        traceback.print_exc()
        return None


def show_help():
    """Display help message."""
    print(f"""
{Colors.BOLD}INTERACTIVE TEST CLI{Colors.END}

{Colors.BOLD}Commands:{Colors.END}
  {Colors.CYAN}<SMILES>{Colors.END}          - Analyze a reaction SMILES
  {Colors.CYAN}examples{Colors.END}           - Show example reactions to try
  {Colors.CYAN}compare <n1> <n2>{Colors.END}  - Compare two saved results
  {Colors.CYAN}save <filename>{Colors.END}    - Save last result to JSON
  {Colors.CYAN}stats{Colors.END}              - Show session statistics
  {Colors.CYAN}help{Colors.END}               - Show this help
  {Colors.CYAN}quit / exit{Colors.END}        - Exit the CLI

{Colors.BOLD}Example reactions:{Colors.END}
  {Colors.YELLOW}1. Tandem (Suzuki + Hydrolysis):{Colors.END}
     COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O

  {Colors.YELLOW}2. Simple SN2:{Colors.END}
     CCCBr.N>>CCCN

  {Colors.YELLOW}3. Esterification:{Colors.END}
     CC(=O)O.CCO>>CC(=O)OCC

  {Colors.YELLOW}4. Grignard Addition:{Colors.END}
     CC(=O)C.CCCMgBr>>CC(O)(CCC)C

{Colors.BOLD}Tips:{Colors.END}
  • Type 'examples' to see more test reactions
  • Compare how all three tiers interpret each reaction
  • Watch for agreements/disagreements between tiers
  • Save interesting results for later analysis
""")


def show_examples():
    """Show example reactions."""
    examples = [
        {
            "name": "Tandem: Suzuki Coupling + Acetal Hydrolysis",
            "smiles": "COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O",
            "expected": "Tier 1 should detect both transformations"
        },
        {
            "name": "Simple SN2 Substitution",
            "smiles": "CCCBr.N>>CCCN",
            "expected": "All tiers should agree - nucleophilic substitution"
        },
        {
            "name": "Esterification",
            "smiles": "CC(=O)O.CCO>>CC(=O)OCC",
            "expected": "Simple condensation/esterification"
        },
        {
            "name": "Grignard Addition to Ketone",
            "smiles": "CC(=O)C.CCCMgBr>>CC(O)(CCC)C",
            "expected": "Nucleophilic addition forming tertiary alcohol"
        },
        {
            "name": "Complex Suzuki Coupling",
            "smiles": "CC1(C)CC=C(B2OC(C)(C)C(C)(C)O2)CC1.CS(=O)(=O)N1CCN(C(=O)c2cnc3ccc(F)cc3c2Br)CC1>>CC1(C)CC=C(c2c(C(=O)N3CCN(S(C)(=O)=O)CC3)cnc3ccc(F)cc23)CC1",
            "expected": "Pd-catalyzed cross-coupling with complex substrate"
        },
        {
            "name": "Wittig Reaction",
            "smiles": "CC=O.C=C(C)C=P(c1ccccc1)(c1ccccc1)c1ccccc1>>CC=C(C)C",
            "expected": "Olefin formation via Wittig mechanism"
        }
    ]

    print(f"\n{Colors.BOLD}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}EXAMPLE REACTIONS{Colors.END}")
    print(f"{Colors.BOLD}{'='*80}{Colors.END}\n")

    for i, ex in enumerate(examples, 1):
        print(f"{Colors.CYAN}{i}. {ex['name']}{Colors.END}")
        print(f"   SMILES: {ex['smiles'][:60]}...")
        if len(ex['smiles']) > 60:
            print(f"           ...{ex['smiles'][60:]}")
        print(f"   {Colors.YELLOW}Expected: {ex['expected']}{Colors.END}")
        print()

    print(f"{Colors.BOLD}To test, just paste the SMILES at the prompt!{Colors.END}\n")


def show_stats(session_stats: dict):
    """Display session statistics."""
    print(f"\n{Colors.BOLD}SESSION STATISTICS{Colors.END}")
    print(f"{'─'*80}")
    print(f"Reactions analyzed: {session_stats['count']}")
    print(f"Total time: {session_stats['total_time']:.1f}s")

    if session_stats['count'] > 0:
        print(f"Average time: {session_stats['total_time']/session_stats['count']:.1f}s")

        # Tier statistics
        print(f"\nTier 2 triggers: {session_stats['tier2_triggers']}/{session_stats['count']}")

        # Agreement statistics
        agreements = session_stats['agreements']
        disagreements = session_stats['disagreements']
        if agreements + disagreements > 0:
            agree_pct = 100 * agreements / (agreements + disagreements)
            print(f"Tier 1-2 agreement: {agreements}/{agreements+disagreements} ({agree_pct:.0f}%)")

    print()


def main():
    """Run interactive test CLI."""

    # Initialize
    print(f"\n{Colors.BOLD}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}INTERACTIVE REACTION TEST CLI{Colors.END}")
    print(f"{Colors.BOLD}Three-Tier Interpretation System{Colors.END}")
    print(f"{Colors.BOLD}{'='*80}{Colors.END}\n")

    print(f"{Colors.CYAN}Initializing LLM client...{Colors.END}")
    client = LLMClient(provider="openai", model="gpt-4o-mini")
    analyzer = ReactionSMILESAnalyzer(client)
    print(f"{Colors.GREEN}✓ Ready!{Colors.END}\n")

    print(f"Type {Colors.CYAN}'help'{Colors.END} for commands or {Colors.CYAN}'examples'{Colors.END} for test reactions")
    print(f"Type {Colors.CYAN}'quit'{Colors.END} to exit\n")

    # Session tracking
    session_stats = {
        'count': 0,
        'total_time': 0.0,
        'tier2_triggers': 0,
        'agreements': 0,
        'disagreements': 0
    }
    last_result = None
    saved_results = []

    # Main loop
    while True:
        try:
            # Prompt
            user_input = input(f"{Colors.GREEN}reaction>{Colors.END} ").strip()

            if not user_input:
                continue

            # Handle commands
            cmd = user_input.lower()

            if cmd in ['quit', 'exit', 'q']:
                print(f"\n{Colors.CYAN}Thanks for testing! Goodbye!{Colors.END}\n")
                break

            elif cmd == 'help':
                show_help()

            elif cmd == 'examples':
                show_examples()

            elif cmd == 'stats':
                show_stats(session_stats)

            elif cmd.startswith('save '):
                if last_result:
                    filename = cmd[5:].strip()
                    if not filename.endswith('.json'):
                        filename += '.json'

                    import os
                    os.makedirs('results', exist_ok=True)
                    filepath = f'results/{filename}'

                    with open(filepath, 'w') as f:
                        json.dump(last_result, f, indent=2, default=str)

                    print(f"\n{Colors.GREEN}✓ Saved to {filepath}{Colors.END}\n")
                else:
                    print(f"\n{Colors.YELLOW}⚠ No result to save yet{Colors.END}\n")

            elif '>>' in user_input:
                # It's a SMILES string - analyze it
                import time
                start = time.time()

                result = analyze_reaction(analyzer, user_input)

                elapsed = time.time() - start

                if result:
                    last_result = result
                    saved_results.append(result)

                    # Update stats
                    session_stats['count'] += 1
                    session_stats['total_time'] += elapsed

                    if result.get('quick_glance', {}).get('success'):
                        session_stats['tier2_triggers'] += 1

                    # Check agreement
                    tier1 = result.get('auto_interpretation', {}).get('interpretation', {})
                    tier2 = result.get('quick_glance', {})

                    if tier1.get('likely_reaction_types') and tier2.get('success'):
                        t1_types = set(t.lower() for t in tier1['likely_reaction_types'])
                        t2_types = set(t.lower() for t in tier2.get('reaction_types', []))

                        if any(t1 in t2_item or t2_item in t1 for t1 in t1_types for t2_item in t2_types):
                            session_stats['agreements'] += 1
                        else:
                            session_stats['disagreements'] += 1

                    print(f"{Colors.CYAN}Analysis completed in {elapsed:.1f}s{Colors.END}\n")

            else:
                print(f"\n{Colors.YELLOW}Unknown command. Type 'help' for available commands.{Colors.END}\n")

        except KeyboardInterrupt:
            print(f"\n\n{Colors.CYAN}Interrupted. Type 'quit' to exit.{Colors.END}\n")

        except Exception as e:
            print(f"\n{Colors.RED}Error: {e}{Colors.END}\n")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()
