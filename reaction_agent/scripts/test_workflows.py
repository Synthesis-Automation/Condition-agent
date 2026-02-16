#!/usr/bin/env python
"""
Comprehensive testing of complex reactions with different workflows and parameters.

Tests reactions from llm_test_rxn.csv with:
- Multiple models (gpt-4o-mini, gpt-4o, o3-mini, o3, gpt-5.2)
- Different parameters (tokens, mapping, spectators)
- Various workflows (fast, balanced, quality, reasoning)

Finds optimal configuration for each reaction type.
"""

import os
import sys
import csv
import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any
import time

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import ReactionSMILESAnalyzer, analyze_deterministic


class Colors:
    """ANSI color codes."""
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    END = '\033[0m'


def print_header(text: str):
    """Print formatted header."""
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'=' * 80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}  {text}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'=' * 80}{Colors.END}")


def load_reactions(csv_path: Path) -> List[Dict[str, str]]:
    """Load reactions from CSV file."""
    reactions = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2 and row[1].strip():
                reactions.append({
                    'name': row[0].strip(),
                    'smiles': row[1].strip()
                })
    return reactions


def analyze_deterministic_only(rxn_smiles: str) -> Dict:
    """Run deterministic analysis only."""
    try:
        return analyze_deterministic(rxn_smiles)
    except Exception as e:
        return {"error": str(e)}


def test_configuration(
    reaction: Dict,
    model: str,
    max_tokens: int,
    skip_mapping: bool = False,
    temperature: float = 0.0,
    provider: str = "openai"
) -> Dict:
    """Test a specific configuration."""

    config_id = f"{model}_{max_tokens}t_map{not skip_mapping}"

    try:
        # Initialize
        client = LLMClient(
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=max_tokens
        )

        analyzer = ReactionSMILESAnalyzer(
            client=client,
            temperature=temperature,
            max_tokens=max_tokens
        )

        # Analyze
        start_time = time.time()
        result = analyzer.analyze(reaction['smiles'], skip_mapping=skip_mapping)
        elapsed = time.time() - start_time

        # Extract metrics
        interp = result.get('interpretation', {})
        metadata = result.get('metadata', {})
        tool_facts = result.get('tool_facts', {})

        mapping_conf = tool_facts.get('mapping_qc', {}).get('confidence', 0.0)

        return {
            'config_id': config_id,
            'model': model,
            'max_tokens': max_tokens,
            'skip_mapping': skip_mapping,
            'temperature': temperature,
            'success': True,
            'overall_class': interp.get('overall_class', 'unknown'),
            'tags': interp.get('tags', []),
            'confidence': interp.get('confidence', 0.0),
            'num_events': len(interp.get('events', [])),
            'num_steps': len(interp.get('mechanism_summary', [])),
            'warnings': interp.get('warnings', []),
            'num_warnings': len(interp.get('warnings', [])),
            'tokens': metadata.get('total_tokens', 0),
            'latency': elapsed,
            'mapping_confidence': mapping_conf,
            'full_result': result
        }

    except Exception as e:
        return {
            'config_id': config_id,
            'model': model,
            'max_tokens': max_tokens,
            'skip_mapping': skip_mapping,
            'success': False,
            'error': str(e),
            'confidence': 0.0
        }


def run_workflow_comparison(reactions: List[Dict]) -> Dict:
    """Run comprehensive workflow comparison."""

    print_header("Complex Reaction Testing - Multiple Workflows")

    # Check API key
    if not os.getenv("OPENAI_API_KEY"):
        print(f"{Colors.RED}✗ Error: OPENAI_API_KEY not set{Colors.END}")
        return {}

    # Define workflows to test
    workflows = {
        "fast_screening": {
            "description": "Fast & Cheap - High throughput screening",
            "configs": [
                {"model": "gpt-4o-mini", "max_tokens": 2000, "skip_mapping": False},
                {"model": "gpt-4o-mini", "max_tokens": 2000, "skip_mapping": True},  # Faster
            ]
        },
        "balanced": {
            "description": "Balanced - Good quality, reasonable cost",
            "configs": [
                {"model": "gpt-4o", "max_tokens": 2500, "skip_mapping": False},
                {"model": "gpt-4o", "max_tokens": 3000, "skip_mapping": False},
            ]
        },
        "high_quality": {
            "description": "High Quality - Complex reactions",
            "configs": [
                {"model": "gpt-4o", "max_tokens": 3500, "skip_mapping": False},
                {"model": "gpt-5.2", "max_tokens": 3000, "skip_mapping": False},
            ]
        },
        "reasoning": {
            "description": "Reasoning Models - Multi-stage reactions",
            "configs": [
                {"model": "o3-mini", "max_tokens": 3000, "skip_mapping": False},
                {"model": "o3-mini", "max_tokens": 4000, "skip_mapping": False},
            ]
        },
        "maximum_detail": {
            "description": "Maximum Detail - Novel/complex mechanisms",
            "configs": [
                {"model": "o3", "max_tokens": 3500, "skip_mapping": False},
                {"model": "o3", "max_tokens": 4500, "skip_mapping": False},
            ]
        }
    }

    all_results = {}

    # Test each reaction
    for rxn_idx, reaction in enumerate(reactions, 1):
        print_header(f"Reaction {rxn_idx}/{len(reactions)}: {reaction['name']}")
        print(f"SMILES: {reaction['smiles'][:80]}...")

        # First, run deterministic analysis
        print(f"\n{Colors.BLUE}→ Running deterministic analysis...{Colors.END}")
        det_result = analyze_deterministic_only(reaction['smiles'])

        if 'error' in det_result:
            print(f"{Colors.RED}✗ Deterministic analysis failed: {det_result['error']}{Colors.END}")
            continue

        # Extract deterministic metrics
        tool_facts = det_result.get('tool_facts', {})
        mapping_ok = tool_facts.get('mapping_qc', {}).get('ok', False)
        mapping_conf = tool_facts.get('mapping_qc', {}).get('confidence', 0.0)
        num_bond_changes = len(tool_facts.get('bond_changes', []))

        print(f"  Mapping: {'✓ OK' if mapping_ok else '✗ Failed'} (confidence: {mapping_conf:.3f})")
        print(f"  Bond changes: {num_bond_changes}")

        reaction_results = {
            'name': reaction['name'],
            'smiles': reaction['smiles'],
            'deterministic': det_result,
            'workflows': {}
        }

        # Test all workflows
        for workflow_name, workflow_config in workflows.items():
            print(f"\n{Colors.CYAN}▶ Workflow: {workflow_name}{Colors.END}")
            print(f"  {workflow_config['description']}")

            workflow_results = []

            for config in workflow_config['configs']:
                config_name = f"{config['model']} ({config['max_tokens']}t)"
                if config['skip_mapping']:
                    config_name += " [no-map]"

                print(f"    Testing {config_name}...", end=' ', flush=True)

                result = test_configuration(
                    reaction=reaction,
                    **config
                )

                if result['success']:
                    print(f"{Colors.GREEN}✓{Colors.END} "
                          f"conf={result['confidence']:.2f} "
                          f"events={result['num_events']} "
                          f"({result['latency']:.1f}s)")
                else:
                    print(f"{Colors.RED}✗ {result.get('error', 'Failed')}{Colors.END}")

                workflow_results.append(result)

            reaction_results['workflows'][workflow_name] = workflow_results

        all_results[reaction['name']] = reaction_results

    return all_results


def analyze_and_recommend(results: Dict) -> Dict:
    """Analyze results and provide recommendations."""

    print_header("ANALYSIS & RECOMMENDATIONS")

    recommendations = {}

    for rxn_name, rxn_results in results.items():
        print(f"\n{Colors.BOLD}Reaction: {rxn_name}{Colors.END}")
        print("-" * 80)

        # Collect all successful configs
        all_configs = []
        for workflow_name, workflow_results in rxn_results['workflows'].items():
            for result in workflow_results:
                if result['success']:
                    all_configs.append({
                        'workflow': workflow_name,
                        **result
                    })

        if not all_configs:
            print(f"{Colors.RED}✗ No successful configurations{Colors.END}")
            continue

        # Sort by confidence
        all_configs.sort(key=lambda x: x['confidence'], reverse=True)

        # Best by different criteria
        best_conf = all_configs[0]
        best_detail = max(all_configs, key=lambda x: x['num_events'])
        best_speed = min(all_configs, key=lambda x: x['latency'])
        best_value = min([c for c in all_configs if c['confidence'] > 0],
                        key=lambda x: x['tokens'] / max(c['confidence'], 0.1))

        print(f"\n{Colors.GREEN}✓ Best Confidence:{Colors.END} {best_conf['model']}")
        print(f"  Workflow: {best_conf['workflow']}")
        print(f"  Confidence: {best_conf['confidence']:.3f}")
        print(f"  Class: {best_conf['overall_class']}")
        print(f"  Tags: {', '.join(best_conf['tags'][:3])}")

        print(f"\n{Colors.YELLOW}✓ Most Detailed:{Colors.END} {best_detail['model']}")
        print(f"  Events: {best_detail['num_events']}, Steps: {best_detail['num_steps']}")

        print(f"\n{Colors.BLUE}✓ Fastest:{Colors.END} {best_speed['model']}")
        print(f"  Latency: {best_speed['latency']:.1f}s")

        print(f"\n{Colors.CYAN}✓ Best Value:{Colors.END} {best_value['model']}")
        print(f"  Conf/Token: {best_value['confidence']/best_value['tokens']*1000:.2f}")

        # Recommendation
        if best_conf['confidence'] >= 0.6:
            rec_label = "Excellent"
            rec_color = Colors.GREEN
        elif best_conf['confidence'] >= 0.4:
            rec_label = "Good"
            rec_color = Colors.YELLOW
        elif best_conf['confidence'] >= 0.2:
            rec_label = "Fair"
            rec_color = Colors.YELLOW
        else:
            rec_label = "Low"
            rec_color = Colors.RED

        print(f"\n{rec_color}→ RECOMMENDATION: {rec_label}{Colors.END}")
        print(f"  Model: {best_conf['model']}")
        print(f"  Max Tokens: {best_conf['max_tokens']}")
        print(f"  Skip Mapping: {best_conf['skip_mapping']}")
        print(f"  Workflow: {best_conf['workflow']}")

        recommendations[rxn_name] = {
            'best_config': best_conf,
            'best_detail': best_detail,
            'best_speed': best_speed,
            'best_value': best_value,
            'all_configs': all_configs
        }

    return recommendations


def generate_comparison_table(results: Dict):
    """Generate comparison table."""

    print_header("WORKFLOW COMPARISON TABLE")

    print(f"\n{Colors.BOLD}Format: Model | Conf | Events | Tokens | Time | Class{Colors.END}")
    print("-" * 100)

    for rxn_name, rxn_results in results.items():
        print(f"\n{Colors.BOLD}{rxn_name}{Colors.END}")
        print("-" * 100)

        for workflow_name, workflow_results in rxn_results['workflows'].items():
            print(f"\n  {workflow_name}:")

            for result in workflow_results:
                if not result['success']:
                    continue

                conf_color = Colors.GREEN if result['confidence'] >= 0.5 else Colors.YELLOW if result['confidence'] >= 0.3 else Colors.RED

                print(f"    {result['model']:<15} | "
                      f"{conf_color}{result['confidence']:.3f}{Colors.END} | "
                      f"{result['num_events']:2d} ev | "
                      f"{result['tokens']:4d} tok | "
                      f"{result['latency']:5.1f}s | "
                      f"{result['overall_class'][:20]}")


def save_results(results: Dict, recommendations: Dict, output_dir: Path):
    """Save detailed results."""

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save full results
    results_file = output_dir / f"workflow_comparison_{timestamp}.json"
    with open(results_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    # Save recommendations
    rec_file = output_dir / f"recommendations_{timestamp}.json"
    with open(rec_file, 'w', encoding='utf-8') as f:
        json.dump(recommendations, f, indent=2, ensure_ascii=False)

    # Generate summary markdown
    summary_file = output_dir / f"workflow_summary_{timestamp}.md"
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("# Workflow Comparison Summary\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        for rxn_name, rec in recommendations.items():
            f.write(f"## {rxn_name}\n\n")

            best = rec['best_config']
            f.write(f"**Recommended Configuration:**\n")
            f.write(f"- Model: {best['model']}\n")
            f.write(f"- Workflow: {best['workflow']}\n")
            f.write(f"- Max Tokens: {best['max_tokens']}\n")
            f.write(f"- Confidence: {best['confidence']:.3f}\n")
            f.write(f"- Classification: {best['overall_class']}\n")
            f.write(f"- Tags: {', '.join(best['tags'])}\n")
            f.write(f"- Events: {best['num_events']}\n")
            f.write(f"\n**Command:**\n")
            f.write(f"```bash\n")
            f.write(f"python reaction_agent/cli.py \\\n")
            f.write(f"  --model {best['model']} \\\n")
            f.write(f"  --max-tokens {best['max_tokens']} \\\n")
            if best['skip_mapping']:
                f.write(f"  --skip-mapping \\\n")
            f.write(f"  --reaction \"...\"\n")
            f.write(f"```\n\n")

    print(f"\n{Colors.GREEN}✓ Results saved to:{Colors.END}")
    print(f"  - {results_file}")
    print(f"  - {rec_file}")
    print(f"  - {summary_file}")


def main():
    """Main testing function."""

    # Load reactions
    csv_path = project_root / "reaction_agent" / "examples" / "llm_test_rxn.csv"

    if not csv_path.exists():
        print(f"{Colors.RED}✗ Error: {csv_path} not found{Colors.END}")
        return

    reactions = load_reactions(csv_path)

    if not reactions:
        print(f"{Colors.RED}✗ Error: No reactions found in CSV{Colors.END}")
        return

    print(f"{Colors.GREEN}✓ Loaded {len(reactions)} reactions{Colors.END}")
    for i, rxn in enumerate(reactions, 1):
        print(f"  {i}. {rxn['name']}")

    # Run comparison
    results = run_workflow_comparison(reactions)

    if not results:
        print(f"{Colors.RED}✗ Error: No results generated{Colors.END}")
        return

    # Analyze
    recommendations = analyze_and_recommend(results)

    # Generate table
    generate_comparison_table(results)

    # Save results
    output_dir = project_root / "reaction_agent" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    save_results(results, recommendations, output_dir)

    print_header("TESTING COMPLETE!")
    print(f"\n{Colors.GREEN}✓ All workflows tested successfully{Colors.END}")
    print(f"  Total reactions: {len(reactions)}")
    print(f"  Total configurations: {sum(len(r['workflows']) for r in results.values())}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n\n{Colors.YELLOW}Interrupted by user{Colors.END}")
    except Exception as e:
        print(f"\n{Colors.RED}✗ Error: {e}{Colors.END}")
        import traceback
        traceback.print_exc()
