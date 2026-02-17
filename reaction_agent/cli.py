#!/usr/bin/env python
"""
Interactive CLI for Reaction SMILES Analysis Agent

A user-friendly command-line interface for analyzing chemical reactions.

Usage:
    python reaction_agent/cli.py
    python reaction_agent/cli.py --model gpt-4o
    python reaction_agent/cli.py --batch reactions.txt
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Optional, Dict, Any, List

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import ReactionSMILESAnalyzer, analyze_deterministic


# Colors for terminal output (ANSI codes)
class Colors:
    """ANSI color codes for terminal output."""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

    @staticmethod
    def disable():
        """Disable colors (for Windows compatibility or piping)."""
        Colors.HEADER = ''
        Colors.BLUE = ''
        Colors.CYAN = ''
        Colors.GREEN = ''
        Colors.YELLOW = ''
        Colors.RED = ''
        Colors.BOLD = ''
        Colors.UNDERLINE = ''
        Colors.END = ''


def print_header(text: str):
    """Print a section header."""
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'=' * 80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}  {text}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'=' * 80}{Colors.END}")


def print_success(text: str):
    """Print success message."""
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")


def print_error(text: str):
    """Print error message."""
    print(f"{Colors.RED}✗ {text}{Colors.END}")


def print_warning(text: str):
    """Print warning message."""
    print(f"{Colors.YELLOW}⚠ {text}{Colors.END}")


def print_info(text: str):
    """Print info message."""
    print(f"{Colors.BLUE}ℹ {text}{Colors.END}")


def print_result(result: Dict[str, Any], show_details: bool = True):
    """Pretty-print analysis result."""

    # Input section
    if show_details:
        print_header("INPUT")
        print(f"Raw:   {result['input']['rxn_smiles_raw']}")
        print(f"Clean: {result['input']['rxn_smiles_clean']}")
        if result['input']['spectators']:
            print(f"Spectators: {Colors.YELLOW}{', '.join(result['input']['spectators'])}{Colors.END}")
        if result['input']['parse_warnings']:
            for warning in result['input']['parse_warnings']:
                print_warning(f"Parse: {warning}")

    # Tool facts section
    if show_details:
        print_header("DETERMINISTIC ANALYSIS")
        tool_facts = result['tool_facts']

        if tool_facts.get('mapped_rxn_smiles'):
            print(f"Mapped: {tool_facts['mapped_rxn_smiles'][:80]}...")

        mapping_qc = tool_facts['mapping_qc']
        qc_ok = mapping_qc.get('ok', False)
        qc_status = f"{Colors.GREEN}✓ OK{Colors.END}" if qc_ok else f"{Colors.RED}✗ Failed{Colors.END}"
        print(f"Mapping QC: {qc_status}")

        if 'confidence' in mapping_qc:
            conf = mapping_qc['confidence']
            conf_color = Colors.GREEN if conf >= 0.7 else Colors.YELLOW if conf >= 0.5 else Colors.RED
            print(f"  Confidence: {conf_color}{conf:.2f}{Colors.END}")

        if mapping_qc.get('notes'):
            for note in mapping_qc['notes']:
                print(f"  Note: {note}")

        if tool_facts['bond_changes']:
            print(f"\nBond Changes ({len(tool_facts['bond_changes'])}):")
            for bc in tool_facts['bond_changes'][:10]:  # Show max 10
                change_color = Colors.RED if bc['change'] == 'broken' else Colors.GREEN if bc['change'] == 'formed' else Colors.YELLOW
                print(f"  {bc['id']}: {change_color}{bc['change']}{Colors.END} bond between :{bc['a1']} and :{bc['a2']} ({bc['bond']})")

            if len(tool_facts['bond_changes']) > 10:
                print(f"  ... and {len(tool_facts['bond_changes']) - 10} more")

        if tool_facts['reaction_center_atoms']:
            atoms_str = ', '.join(map(str, tool_facts['reaction_center_atoms'][:15]))
            print(f"\nReaction Center: {atoms_str}")
            if len(tool_facts['reaction_center_atoms']) > 15:
                print(f"  ... and {len(tool_facts['reaction_center_atoms']) - 15} more atoms")

    # Automatic interpretation section (deterministic pattern-based)
    if 'auto_interpretation' in result and result['auto_interpretation']:
        auto_interp = result['auto_interpretation']

        if 'error' not in auto_interp and 'interpretation' in auto_interp:
            print_header("AUTOMATIC INTERPRETATION")

            interpretation = auto_interp['interpretation']

            # Display complexity with emoji
            complexity = interpretation.get('reaction_complexity', 'unknown')
            complexity_emoji = {
                'simple': '🟢',
                'moderate': '🟡',
                'complex': '🟠',
                'tandem/multi-step': '🔴'
            }
            emoji = complexity_emoji.get(complexity, '⚪')
            print(f"{emoji} {Colors.BOLD}Complexity:{Colors.END} {Colors.CYAN}{complexity.upper()}{Colors.END}")

            # Display detected patterns
            if interpretation.get('patterns_detected'):
                print(f"\n{Colors.BOLD}Patterns Detected:{Colors.END}")
                for pattern in interpretation['patterns_detected']:
                    print(f"  • {pattern}")

            # Display likely reaction types
            if interpretation.get('likely_reaction_types'):
                print(f"\n{Colors.BOLD}Likely Reaction Type(s):{Colors.END}")
                for rxn_type in interpretation['likely_reaction_types']:
                    print(f"  • {Colors.GREEN}{rxn_type}{Colors.END}")

            # Tandem reaction flag
            if interpretation.get('tandem_reaction_suspected'):
                print(f"\n{Colors.YELLOW}⚠️  TANDEM/MULTI-STEP REACTION SUSPECTED{Colors.END}")

            # Explanation
            if interpretation.get('explanation'):
                print(f"\n{Colors.BOLD}Explanation:{Colors.END}")
                for line in interpretation['explanation'].split('\n'):
                    print(f"  {line}")

            # Recommendation
            if interpretation.get('recommendation'):
                print(f"\n{Colors.BOLD}RECOMMENDATION:{Colors.END}")
                for line in interpretation['recommendation'].split('\n'):
                    print(f"  {line}")

    # Quick glance section (Tier 2 LLM-based fast analysis)
    if 'quick_glance' in result and result['quick_glance'] and result['quick_glance'].get('success'):
        quick = result['quick_glance']

        print_header("QUICK LLM GLANCE (Tier 2)")

        # Summary
        summary = quick.get('summary', 'N/A')
        print(f"{Colors.BOLD}Summary:{Colors.END} {summary}")
        print()

        # Reaction types
        reaction_types = quick.get('reaction_types', [])
        if reaction_types:
            print(f"{Colors.BOLD}Reaction Type(s):{Colors.END}")
            for rxn_type in reaction_types:
                print(f"  • {Colors.CYAN}{rxn_type}{Colors.END}")
            print()

        # Patterns
        patterns = quick.get('patterns', [])
        if patterns:
            print(f"{Colors.BOLD}Key Patterns:{Colors.END}")
            for pattern in patterns:
                print(f"  • {pattern}")
            print()

        # All changes (if comprehensive mode)
        all_changes = quick.get('all_changes', [])
        if all_changes:
            print(f"{Colors.BOLD}All Structural Changes:{Colors.END}")
            for change in all_changes[:5]:  # Show first 5
                print(f"  • {change}")
            if len(all_changes) > 5:
                print(f"  ... and {len(all_changes)-5} more")
            print()

        # Protecting group changes (VERY IMPORTANT!)
        pg_changes = quick.get('protecting_groups', {})
        if pg_changes and (pg_changes.get('removed') or pg_changes.get('added')):
            print(f"{Colors.BOLD}{Colors.YELLOW}Protecting Group Changes:{Colors.END}")
            if pg_changes.get('removed'):
                print(f"  {Colors.RED}Removed (deprotection):{Colors.END}")
                for pg in pg_changes['removed']:
                    print(f"    • {pg}")
            if pg_changes.get('added'):
                print(f"  {Colors.GREEN}Added (protection):{Colors.END}")
                for pg in pg_changes['added']:
                    print(f"    • {pg}")
            print()

        # Side reactions
        side_rxns = quick.get('side_reactions', [])
        if side_rxns:
            print(f"{Colors.BOLD}Side/Workup Transformations:{Colors.END}")
            for rxn in side_rxns:
                print(f"  • {rxn}")
            print()

        # Complexity
        complexity = quick.get('complexity', 'unknown')
        complexity_color = {
            'simple': Colors.GREEN,
            'moderate': Colors.YELLOW,
            'complex': Colors.RED,
            'tandem': Colors.RED
        }.get(complexity, Colors.END)
        print(f"{Colors.BOLD}Complexity:{Colors.END} {complexity_color}{complexity.upper()}{Colors.END}")

        # Confidence
        llm_confidence = quick.get('confidence', 0.0)
        conf_color = Colors.GREEN if llm_confidence >= 0.7 else Colors.YELLOW if llm_confidence >= 0.5 else Colors.RED
        print(f"{Colors.BOLD}LLM Confidence:{Colors.END} {conf_color}{llm_confidence:.2f}{Colors.END}")

        # Metadata (timing)
        if 'metadata' in quick:
            meta = quick['metadata']
            latency = meta.get('latency_ms', 0)
            print(f"\n{Colors.BOLD}Analysis Time:{Colors.END} {latency:.0f}ms ({meta.get('model', 'N/A')})")

    # Interpretation section
    print_header("LLM INTERPRETATION (Tier 3)")
    interp = result.get('interpretation', {})

    if 'error' in interp:
        print_error(f"LLM Error: {interp['error']}")
        if 'raw_response' in interp:
            print(f"\nRaw response (truncated):\n{interp['raw_response'][:300]}...")
        return

    # Check for disagreement with Tier 2
    if 'quick_glance' in result and result['quick_glance'] and result['quick_glance'].get('success'):
        tier2_types = [rt.lower() for rt in result['quick_glance'].get('reaction_types', [])]
        tier3_class = interp.get('overall_class', '').lower()

        # Check if Tier 3 contradicts Tier 2 on major reaction types
        suzuki_in_t2 = any('suzuki' in rt or 'coupling' in rt for rt in tier2_types)
        substitution_in_t3 = 'substitution' in tier3_class

        if suzuki_in_t2 and substitution_in_t3:
            print(f"{Colors.YELLOW}⚠️  Warning: Classification mismatch detected!{Colors.END}")
            print(f"   Tier 2 (DeepSeek-v3.2): {', '.join(result['quick_glance']['reaction_types'])}")
            print(f"   Tier 3 (gpt-4o-mini): {tier3_class}")
            print(f"   {Colors.BOLD}→ Tier 2 is likely more accurate for this reaction.{Colors.END}\n")

    # Main classification
    rxn_class = interp.get('overall_class', 'N/A')
    print(f"{Colors.BOLD}Reaction Class:{Colors.END} {Colors.CYAN}{rxn_class}{Colors.END}")

    if interp.get('tags'):
        tags_str = ', '.join(interp['tags'])
        print(f"{Colors.BOLD}Tags:{Colors.END} {tags_str}")

    # Confidence
    confidence = interp.get('confidence', 0.0)
    conf_color = Colors.GREEN if confidence >= 0.7 else Colors.YELLOW if confidence >= 0.5 else Colors.RED
    print(f"\n{Colors.BOLD}Confidence:{Colors.END} {conf_color}{confidence:.2f}{Colors.END}")

    # Roles
    if interp.get('roles'):
        print(f"\n{Colors.BOLD}Roles:{Colors.END}")
        for role, value in interp['roles'].items():
            print(f"  {role}: {value}")

    # Events - simplified display
    if interp.get('events'):
        events_count = len(interp['events'])
        print(f"\n{Colors.BOLD}Mechanistic Events:{Colors.END} {events_count} detected")

        # Only show details if there are multiple events or if it's complex
        show_event_details = events_count > 1 or any('complex' in str(e).lower() for e in interp['events'])

        if show_event_details:
            for event in interp['events']:
                print(f"  • {event['event_id']}: {event['event_type']} ({', '.join(event.get('bond_change_refs', []))})")
        else:
            # For single simple events, just show a summary
            event = interp['events'][0]
            print(f"  • Single {event['event_type']} event")

    # Mechanism summary - more concise
    if interp.get('mechanism_summary'):
        print(f"\n{Colors.BOLD}Mechanism Summary:{Colors.END}")
        for i, step in enumerate(interp['mechanism_summary'], 1):
            print(f"  {i}. {step}")

    # Warnings
    if interp.get('warnings'):
        print(f"\n{Colors.BOLD}Warnings:{Colors.END}")
        for warning in interp['warnings']:
            print_warning(warning)

    # Metadata
    if 'metadata' in result and show_details:
        print_header("METADATA")
        meta = result['metadata']

        # Show model selection info if available
        if 'model_selected' in meta:
            print(f"Model Selected: {Colors.BOLD}{meta['model_selected']}{Colors.END}")
            if 'mode' in meta:
                print(f"Mode: {meta['mode']}")
            if 'reasoning_effort' in meta and meta['reasoning_effort']:
                print(f"Reasoning Effort: {meta['reasoning_effort']}")
        else:
            print(f"Model: {Colors.BOLD}{meta['model']}{Colors.END} ({meta['provider']})")

        print(f"Tokens: {meta['total_tokens']} (prompt: {meta['prompt_tokens']}, completion: {meta['completion_tokens']})")
        print(f"Latency: {meta['latency_ms']:.0f} ms")

        if meta['total_tokens'] > 0 and meta['latency_ms'] > 0:
            tokens_per_sec = (meta['total_tokens'] / meta['latency_ms']) * 1000
            print(f"Speed: {tokens_per_sec:.1f} tokens/sec")


def analyze_reaction_interactive(
    analyzer: ReactionSMILESAnalyzer,
    rxn_smiles: str,
    skip_mapping: bool = False,
    save_output: Optional[Path] = None,
    mode: str = "auto"
) -> Dict[str, Any]:
    """Analyze a reaction interactively with progress updates."""

    print(f"\n{Colors.BOLD}Analyzing:{Colors.END} {rxn_smiles}")
    print("-" * 80)

    # Step 1: Deterministic analysis
    print_info(f"Step 1/2: Running deterministic analysis (mode={mode})...")

    try:
        result = analyzer.analyze(rxn_smiles, skip_mapping=skip_mapping, mode=mode)
    except Exception as e:
        print_error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None

    print_success("Deterministic analysis complete")

    # Step 2: LLM interpretation (if not skipped)
    if 'interpretation' in result:
        print_info("Step 2/2: LLM interpretation complete")

    # Print results
    print_result(result)

    # Save output if requested
    if save_output:
        save_output.parent.mkdir(parents=True, exist_ok=True)
        with open(save_output, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        print_success(f"Results saved to {save_output}")

    return result


def batch_analyze(
    analyzer: ReactionSMILESAnalyzer,
    reactions: List[str],
    output_dir: Optional[Path] = None,
    mode: str = "auto"
) -> List[Dict[str, Any]]:
    """Analyze multiple reactions in batch."""

    print_header(f"BATCH ANALYSIS ({len(reactions)} reactions, mode={mode})")

    results = []
    for i, rxn in enumerate(reactions, 1):
        print(f"\n{Colors.BOLD}[{i}/{len(reactions)}]{Colors.END} {rxn[:60]}...")

        try:
            result = analyzer.analyze(rxn, mode=mode)
            results.append(result)

            # Show brief result
            interp = result.get('interpretation', {})
            rxn_class = interp.get('overall_class', 'N/A')
            confidence = interp.get('confidence', 0.0)
            print(f"  → {rxn_class} (confidence: {confidence:.2f})")

            # Save individual result
            if output_dir:
                output_file = output_dir / f"reaction_{i:03d}.json"
                with open(output_file, 'w', encoding='utf-8') as f:
                    json.dump(result, f, indent=2, ensure_ascii=False)

        except Exception as e:
            print_error(f"Failed: {e}")
            results.append({"error": str(e), "input": {"rxn_smiles_raw": rxn}})

    # Summary
    print_header("BATCH SUMMARY")
    successful = sum(1 for r in results if 'interpretation' in r and 'error' not in r.get('interpretation', {}))
    print(f"Total: {len(results)}")
    print(f"Successful: {Colors.GREEN}{successful}{Colors.END}")
    print(f"Failed: {Colors.RED}{len(results) - successful}{Colors.END}")

    if output_dir:
        summary_file = output_dir / "batch_summary.json"
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        print_success(f"Batch results saved to {output_dir}")

    return results


def deterministic_only_mode(rxn_smiles: str):
    """Run only deterministic analysis without LLM."""
    print_header("DETERMINISTIC ANALYSIS ONLY (No LLM)")
    print(f"Analyzing: {rxn_smiles}\n")

    try:
        result = analyze_deterministic(rxn_smiles)

        print_success("Analysis complete")
        print("\n" + json.dumps(result, indent=2))

        return result

    except Exception as e:
        print_error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def interactive_mode(analyzer: ReactionSMILESAnalyzer):
    """Run interactive mode with prompts."""

    print_header("INTERACTIVE MODE")
    print("\nEnter reaction SMILES to analyze (or 'quit' to exit)")
    print("Commands:")
    print("  'quit' or 'exit' - Exit the program")
    print("  'config' - Show current configuration")
    print("  'help' - Show this help message")
    print("  'batch <file>' - Analyze reactions from file")

    while True:
        try:
            # Prompt for input
            user_input = input(f"\n{Colors.BOLD}{Colors.GREEN}reaction>{Colors.END} ").strip()

            if not user_input:
                continue

            # Handle commands
            if user_input.lower() in ['quit', 'exit', 'q']:
                print("\nGoodbye!")
                break

            elif user_input.lower() == 'help':
                print("\nCommands:")
                print("  quit/exit - Exit the program")
                print("  config    - Show current configuration")
                print("  help      - Show this help message")
                print("  batch <file> - Analyze reactions from file")
                print("\nOr enter a reaction SMILES string to analyze")
                continue

            elif user_input.lower() == 'config':
                print(f"\nCurrent Configuration:")
                print(f"  Model: {analyzer.client.model}")
                print(f"  Provider: {analyzer.client.provider}")
                print(f"  Temperature: {analyzer.temperature}")
                print(f"  Max Tokens: {analyzer.max_tokens}")
                print(f"  Reasoning Effort: {analyzer.reasoning_effort or 'N/A'}")
                print(f"  Drop Spectators: {analyzer.drop_spectators}")
                continue

            elif user_input.lower().startswith('batch '):
                file_path = user_input[6:].strip()
                try:
                    with open(file_path, 'r') as f:
                        reactions = [line.strip() for line in f if line.strip() and not line.startswith('#')]

                    batch_analyze(analyzer, reactions)
                except FileNotFoundError:
                    print_error(f"File not found: {file_path}")
                except Exception as e:
                    print_error(f"Batch analysis failed: {e}")
                continue

            # Analyze the reaction
            analyze_reaction_interactive(analyzer, user_input)

        except KeyboardInterrupt:
            print("\n\nInterrupted by user. Type 'quit' to exit.")
            continue
        except EOFError:
            print("\nGoodbye!")
            break
        except Exception as e:
            print_error(f"Error: {e}")
            import traceback
            traceback.print_exc()


def main():
    """Main CLI entry point."""

    parser = argparse.ArgumentParser(
        description="Interactive CLI for Reaction SMILES Analysis Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode
  python reaction_agent/cli.py

  # Analyze single reaction (auto mode - smart switching)
  python reaction_agent/cli.py --reaction "CCBr>>CCN" --mode auto

  # Analyze with fast mode (always gpt-4o)
  python reaction_agent/cli.py --reaction "CCBr>>CCN" --mode fast

  # Analyze with deep mode (always gpt-5.2 low reasoning)
  python reaction_agent/cli.py --reaction "CCBr>>CCN" --mode deep

  # Analyze with expert mode (gpt-5.2 medium reasoning)
  python reaction_agent/cli.py --reaction "CCBr>>CCN" --mode expert

  # Batch analysis from file
  python reaction_agent/cli.py --batch reactions.txt --mode auto

  # Deterministic only (no LLM)
  python reaction_agent/cli.py --reaction "CCBr>>CCN" --no-llm

  # Always use GPT-5.2 (forces deep reasoning mode)
  python reaction_agent/cli.py --model gpt-5.2 --reaction "CCBr>>CCN"

  # Use o3-mini model
  python reaction_agent/cli.py --model o3-mini --reaction "CCBr>>CCN"

  # Save output
  python reaction_agent/cli.py --reaction "CCBr>>CCN" --output result.json
        """
    )

    # Input options
    parser.add_argument('--reaction', '-r', type=str, help='Reaction SMILES to analyze')
    parser.add_argument('--batch', '-b', type=Path, help='File with reactions (one per line)')
    parser.add_argument('--interactive', '-i', action='store_true', help='Run in interactive mode (default if no reaction provided)')

    # LLM options
    parser.add_argument('--model', '-m', type=str, default='gpt-4o-mini',
                        help='LLM model to use (default: gpt-4o-mini). Specify gpt-5.2 to force deep reasoning.')
    parser.add_argument('--provider', '-p', type=str, default='openai', help='LLM provider (default: openai)')
    parser.add_argument('--temperature', '-t', type=float, default=0.0, help='Temperature (default: 0.0)')
    parser.add_argument('--max-tokens', type=int, default=2000, help='Max tokens (default: 2000)')
    parser.add_argument('--mode', type=str, default='auto', choices=['auto', 'fast', 'deep', 'expert'],
                        help='Analysis mode: auto (smart switching), fast (always gpt-4o), deep (gpt-5.2 low), expert (gpt-5.2 medium) (default: auto)')
    parser.add_argument('--reasoning-effort', type=str, choices=['low', 'medium', 'high'],
                        help='Reasoning effort level for GPT-5.2 (overrides mode defaults)')
    parser.add_argument('--no-llm', action='store_true', help='Run deterministic analysis only (no LLM)')

    # Analysis options
    parser.add_argument('--skip-mapping', action='store_true', help='Skip atom mapping (faster)')
    parser.add_argument('--keep-spectators', action='store_true', help='Keep spectators in analysis')

    # Output options
    parser.add_argument('--output', '-o', type=Path, help='Save results to file (JSON)')
    parser.add_argument('--output-dir', type=Path, help='Directory for batch output')
    parser.add_argument('--no-color', action='store_true', help='Disable colored output')
    parser.add_argument('--quiet', '-q', action='store_true', help='Minimal output')

    args = parser.parse_args()

    # Disable colors if requested or on Windows without proper support
    if args.no_color or os.name == 'nt':
        Colors.disable()

    # Print banner
    if not args.quiet:
        print_header("Reaction SMILES Analysis Agent - Interactive CLI")
        print(f"\nModel: {Colors.BOLD}{args.model}{Colors.END} ({args.provider})")
        print(f"Mode: {Colors.BOLD}{args.mode}{Colors.END}")
        if args.reasoning_effort:
            print(f"Reasoning Effort: {Colors.BOLD}{args.reasoning_effort}{Colors.END}")

    # Check for API key unless running deterministic only
    if not args.no_llm:
        api_key_name = f"{args.provider.upper()}_API_KEY"
        if not os.getenv(api_key_name):
            print_error(f"{api_key_name} environment variable not set")
            print(f"\nPlease set it:")
            print(f"  export {api_key_name}='your-key-here'")
            print(f"\nOr run with --no-llm for deterministic analysis only")
            sys.exit(1)

    # Deterministic only mode
    if args.no_llm:
        if not args.reaction:
            print_error("--reaction required for --no-llm mode")
            sys.exit(1)

        deterministic_only_mode(args.reaction)
        return

    # Initialize LLM client and analyzer
    try:
        if not args.quiet:
            print_info("Initializing LLM client...")

        # Auto-detect if user wants to force GPT-5.2 based on model name
        is_gpt52_model = args.model.startswith("gpt-5.2") or args.model.startswith("gpt-o3") or args.model.startswith("o3-")

        # If user specified a GPT-5.2 model directly, force deep mode (unless they explicitly set a mode)
        effective_mode = args.mode
        if is_gpt52_model and args.mode == "auto":
            effective_mode = "deep"
            if not args.quiet:
                print_info(f"Detected reasoning model '{args.model}' - using deep reasoning mode")

        client = LLMClient(
            provider=args.provider,
            model=args.model,
            temperature=args.temperature,
            max_tokens=args.max_tokens
        )

        # Determine reasoning effort based on mode and explicit arg
        reasoning_effort = args.reasoning_effort
        if reasoning_effort is None:
            # Set defaults based on mode (will be overridden by mode in analyze())
            if effective_mode == "deep":
                reasoning_effort = "low"
            elif effective_mode == "expert":
                reasoning_effort = "medium"
            # auto and fast modes don't set reasoning_effort at init

        analyzer = ReactionSMILESAnalyzer(
            client=client,
            drop_spectators=not args.keep_spectators,
            temperature=args.temperature,
            max_tokens=args.max_tokens,
            reasoning_effort=reasoning_effort
        )

        if not args.quiet:
            print_success(f"Initialized: {client.model}")

    except Exception as e:
        print_error(f"Failed to initialize: {e}")
        sys.exit(1)

    # Run appropriate mode
    if args.batch:
        # Batch mode
        if not args.batch.exists():
            print_error(f"File not found: {args.batch}")
            sys.exit(1)

        with open(args.batch, 'r') as f:
            reactions = [line.strip() for line in f if line.strip() and not line.startswith('#')]

        output_dir = args.output_dir or Path('results') / args.batch.stem
        batch_analyze(analyzer, reactions, output_dir, mode=effective_mode)

    elif args.reaction:
        # Single reaction mode
        analyze_reaction_interactive(
            analyzer,
            args.reaction,
            skip_mapping=args.skip_mapping,
            save_output=args.output,
            mode=effective_mode
        )

    else:
        # Interactive mode (default)
        interactive_mode(analyzer)


if __name__ == "__main__":
    main()
