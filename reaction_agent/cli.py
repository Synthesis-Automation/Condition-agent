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

    # Interpretation section
    print_header("LLM INTERPRETATION")
    interp = result.get('interpretation', {})

    if 'error' in interp:
        print_error(f"LLM Error: {interp['error']}")
        if 'raw_response' in interp:
            print(f"\nRaw response (truncated):\n{interp['raw_response'][:300]}...")
        return

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

    # Events
    if interp.get('events'):
        print(f"\n{Colors.BOLD}Mechanistic Events ({len(interp['events'])}):{Colors.END}")
        for event in interp['events']:
            print(f"\n  {Colors.BOLD}{event['event_id']}: {event['event_type']}{Colors.END}")
            print(f"    Bond changes: {', '.join(event.get('bond_change_refs', []))}")
            print(f"    {event.get('short_rationale', 'N/A')}")
            event_conf = event.get('confidence', 0.0)
            event_color = Colors.GREEN if event_conf >= 0.7 else Colors.YELLOW
            print(f"    Confidence: {event_color}{event_conf:.2f}{Colors.END}")

    # Mechanism summary
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
    save_output: Optional[Path] = None
) -> Dict[str, Any]:
    """Analyze a reaction interactively with progress updates."""

    print(f"\n{Colors.BOLD}Analyzing:{Colors.END} {rxn_smiles}")
    print("-" * 80)

    # Step 1: Deterministic analysis
    print_info("Step 1/2: Running deterministic analysis...")

    try:
        result = analyzer.analyze(rxn_smiles, skip_mapping=skip_mapping)
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
    output_dir: Optional[Path] = None
) -> List[Dict[str, Any]]:
    """Analyze multiple reactions in batch."""

    print_header(f"BATCH ANALYSIS ({len(reactions)} reactions)")

    results = []
    for i, rxn in enumerate(reactions, 1):
        print(f"\n{Colors.BOLD}[{i}/{len(reactions)}]{Colors.END} {rxn[:60]}...")

        try:
            result = analyzer.analyze(rxn)
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

  # Analyze single reaction
  python reaction_agent/cli.py --reaction "CCBr>>CCN"

  # Batch analysis from file
  python reaction_agent/cli.py --batch reactions.txt

  # Deterministic only (no LLM)
  python reaction_agent/cli.py --reaction "CCBr>>CCN" --no-llm

  # Use different model
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
    parser.add_argument('--model', '-m', type=str, default='gpt-4o-mini', help='LLM model to use (default: gpt-4o-mini)')
    parser.add_argument('--provider', '-p', type=str, default='openai', help='LLM provider (default: openai)')
    parser.add_argument('--temperature', '-t', type=float, default=0.0, help='Temperature (default: 0.0)')
    parser.add_argument('--max-tokens', type=int, default=2000, help='Max tokens (default: 2000)')
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

        client = LLMClient(
            provider=args.provider,
            model=args.model,
            temperature=args.temperature,
            max_tokens=args.max_tokens
        )

        analyzer = ReactionSMILESAnalyzer(
            client=client,
            drop_spectators=not args.keep_spectators,
            temperature=args.temperature,
            max_tokens=args.max_tokens
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
        batch_analyze(analyzer, reactions, output_dir)

    elif args.reaction:
        # Single reaction mode
        analyze_reaction_interactive(
            analyzer,
            args.reaction,
            skip_mapping=args.skip_mapping,
            save_output=args.output
        )

    else:
        # Interactive mode (default)
        interactive_mode(analyzer)


if __name__ == "__main__":
    main()
