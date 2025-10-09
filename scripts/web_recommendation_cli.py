"""
Interactive Recommendation Tester
=================================

This script prompts for a reaction SMILES and lets you choose a reaction
family (or auto-detect). It exercises three recommendation modes
(rule-based, ML, fusion), saves their JSON responses, and prints a compact
summary to the console.

Usage:
    python scripts/interactive_recommendation_cli.py

Requirements:
    - FastAPI server running (default: http://localhost:8000)
    - requests package installed
"""

from __future__ import annotations

import io
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import requests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    from scripts.recommendation_cli_utils import (
        DEFAULT_SCDB_PATH,
        FUSION_VARIANTS_DEFAULT,
        K_DEFAULT,
        LIMIT_DEFAULT,
        choose_reaction_type,
        prompt_smiles,
        save_json,
        slugify_label,
        summarize_fusion,
        summarize_ml,
        summarize_rule,
    )
except ModuleNotFoundError:
    sys.path.append(str(HERE))
    from recommendation_cli_utils import (
        DEFAULT_SCDB_PATH,
        FUSION_VARIANTS_DEFAULT,
        K_DEFAULT,
        LIMIT_DEFAULT,
        choose_reaction_type,
        prompt_smiles,
        save_json,
        slugify_label,
        summarize_fusion,
        summarize_ml,
        summarize_rule,
    )


if sys.platform == "win32":
    # Ensure UTF-8 output when running in Windows terminals.
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")


DEFAULT_BASE_URL = "http://localhost:8000"


def health_check(base_url: str) -> bool:
    """Confirm the FastAPI server is reachable."""
    try:
        response = requests.get(f"{base_url}/health", timeout=5)
        if response.status_code == 200:
            print("Server healthy.")
            return True
        print(f"Health check failed: HTTP {response.status_code}")
    except requests.RequestException as exc:
        print(f"Health check failed: {exc}")
    return False


def call_rule_based(
    base_url: str,
    reaction: str,
    db_path: Optional[str],
) -> Dict[str, Any]:
    """Execute the rule-based /match endpoint."""
    payload: Dict[str, Any] = {
        "reaction": reaction,
        "include_trace": True,
    }
    if db_path:
        payload["db"] = db_path

    try:
        response = requests.post(f"{base_url}/match", json=payload, timeout=20)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as exc:
        return {
            "error": f"Rule-based request failed: {exc}",
            "payload": payload,
        }
    except ValueError:
        return {
            "error": "Rule-based response was not valid JSON.",
            "payload": payload,
        }


def call_ml_recommendation(
    base_url: str,
    reaction: str,
    reaction_type: Optional[str],
    k: int,
    limit: int,
    rerank_strategy: str = 'rule',
    filter_unknown_reagents: bool = False,
) -> Dict[str, Any]:
    """Execute the ML recommendation endpoint."""
    payload: Dict[str, Any] = {
        "reaction": reaction,
        "reaction_type": reaction_type or None,
        "k": k,
        "limit": limit,
        "relax": {},
        "constraints": {},
        "rerank_strategy": rerank_strategy,
        "filter_unknown_reagents": filter_unknown_reagents,
    }

    try:
        response = requests.post(
            f"{base_url}/api/v1/recommend/conditions",
            json=payload,
            timeout=30,
        )
        response.raise_for_status()
        return response.json()
    except requests.RequestException as exc:
        return {
            "error": f"ML recommendation failed: {exc}",
            "payload": payload,
        }
    except ValueError:
        return {
            "error": "ML recommendation response was not valid JSON.",
            "payload": payload,
        }


def call_fusion_recommendation(
    base_url: str,
    reaction: str,
    k: int,
    max_variants: int,
) -> Dict[str, Any]:
    """
    Execute the fusion recommendation endpoint.
    
    NOTE: The /api/v1/recommend/fusion endpoint is DEPRECATED.
    Use /api/v1/recommend/conditions with rerank_strategy='rule' instead.
    """
    import warnings
    warnings.warn(
        "Fusion endpoint is deprecated. Use call_ml_recommendation() "
        "with rerank_strategy='rule' instead.",
        DeprecationWarning,
        stacklevel=2
    )
    
    payload: Dict[str, Any] = {
        "reaction": reaction,
        "k": k,
        "max_variants": max_variants,
        "relax": {},
        "constraints": {},
    }

    try:
        response = requests.post(
            f"{base_url}/api/v1/recommend/fusion",
            json=payload,
            timeout=30,
        )
        response.raise_for_status()
        return response.json()
    except requests.RequestException as exc:
        return {
            "error": f"Fusion recommendation failed: {exc}",
            "payload": payload,
        }
    except ValueError:
        return {
            "error": "Fusion recommendation response was not valid JSON.",
            "payload": payload,
        }


def main() -> None:
    """Main entry point with optional command-line arguments."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Web Recommendation Tester - Test ChemTools API via HTTP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (prompts for input):
  python scripts/web_recommendation_cli.py
  
  # Provide reaction and type via command line:
  python scripts/web_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --family Buchwald_CN
  
  # Auto-detect reaction type:
  python scripts/web_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
  
  # Custom server URL:
  python scripts/web_recommendation_cli.py --url http://myserver:8080
        """
    )
    
    parser.add_argument(
        "--url", "--base-url",
        type=str,
        default=DEFAULT_BASE_URL,
        help=f"Base URL of the FastAPI server (default: {DEFAULT_BASE_URL})"
    )
    
    parser.add_argument(
        "--rxn", "--reaction",
        type=str,
        default=None,
        help="Reaction SMILES (reactants>>products). If not provided, will prompt interactively."
    )
    
    parser.add_argument(
        "--family", "--type",
        type=str,
        default=None,
        choices=[None, "Suzuki", "Suzuki_CC", "C_N_Coupling_Cu", "Ullmann_CN", 
                 "C_N_Coupling_Pd", "Buchwald_CN", "C_N_Coupling_Ni", "Amide_formation"],
        help="Reaction family/type. If not provided, will prompt interactively or auto-detect."
    )
    
    parser.add_argument(
        "--k",
        type=int,
        default=K_DEFAULT,
        help=f"Number of precedents to retrieve for ML (default: {K_DEFAULT})"
    )
    
    parser.add_argument(
        "--limit",
        type=int,
        default=LIMIT_DEFAULT,
        help=f"Number of ML recommendations to return (default: {LIMIT_DEFAULT})"
    )
    
    parser.add_argument(
        "--fusion-variants",
        type=int,
        default=FUSION_VARIANTS_DEFAULT,
        help=f"Number of fusion recommendation variants (default: {FUSION_VARIANTS_DEFAULT})"
    )
    
    parser.add_argument(
        "--save-dir",
        type=str,
        default="results",
        help="Directory to save output JSON files (default: results)"
    )
    
    parser.add_argument(
        "--strategy",
        type=str,
        default="all",
        choices=["all", "rule", "ml", "fusion"],
        help="Which recommendation strategy to run (default: all). NOTE: 'fusion' is deprecated."
    )
    
    parser.add_argument(
        "--rerank",
        type=str,
        default="rule",
        choices=["none", "rule", "analytics"],
        help="Reranking strategy for ML recommendations: 'none' (similarity only), "
             "'rule' (boost by rule match), 'analytics' (boost by popularity). Default: rule"
    )
    
    parser.add_argument(
        "--filter-unknown",
        action="store_true",
        help="Filter out precedents with unknown base/solvent reagents (not in database)"
    )
    
    args = parser.parse_args()
    
    print("Interactive Recommendation Test")
    print("--------------------------------")

    base_url = args.url
    print(f"Base URL: {base_url}")
    if not health_check(base_url):
        print("Aborting due to failed health check.")
        return

    # Get reaction SMILES - from args or prompt
    if args.rxn:
        reaction = args.rxn.strip()
        print(f"Reaction SMILES: {reaction}")
    else:
        reaction = prompt_smiles()
    
    # Get reaction type - from args or prompt
    if args.family:
        reaction_type = args.family
        selected_label = args.family
        print(f"Reaction type: {selected_label}")
    else:
        selected_label, reaction_type = choose_reaction_type()
        print(f"Selected reaction type: {selected_label}")

    k_value = args.k
    limit_value = args.limit
    fusion_variants = args.fusion_variants

    db_path = DEFAULT_SCDB_PATH if Path(DEFAULT_SCDB_PATH).exists() else None
    if db_path is None:
        print("Warning: Default SCDB database not found; rule-based test will use server defaults.")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    type_label = reaction_type or "auto"
    label = slugify_label(type_label)

    print("\nRunning tests...\n")

    # Update output directory from args
    output_dir = Path(args.save_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a custom save function with the specified directory
    def save_to_dir(data: Dict[str, Any], filename: str) -> Path:
        """Save JSON data to the specified output directory."""
        output_path = output_dir / filename
        import json
        output_path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
        return output_path

    # Run selected strategies
    run_rule = args.strategy in ["all", "rule"]
    run_ml = args.strategy in ["all", "ml"]
    run_fusion = args.strategy in ["all", "fusion"]
    
    rule_result = None
    ml_result = None
    fusion_result = None
    
    if run_rule:
        rule_result = call_rule_based(base_url, reaction, db_path)
        rule_file = save_to_dir(rule_result, f"{timestamp}_{label}_rule.json")
    
    if run_ml:
        ml_result = call_ml_recommendation(
            base_url, 
            reaction, 
            reaction_type, 
            k_value, 
            limit_value,
            rerank_strategy=args.rerank,
            filter_unknown_reagents=args.filter_unknown
        )
        ml_file = save_to_dir(ml_result, f"{timestamp}_{label}_ml.json")
    
    if run_fusion:
        fusion_result = call_fusion_recommendation(base_url, reaction, k_value, fusion_variants)
        fusion_file = save_to_dir(fusion_result, f"{timestamp}_{label}_fusion.json")

    print("Summary\n-------")
    
    if rule_result:
        summarize_rule(rule_result)
        print()
    
    if ml_result:
        summarize_ml(ml_result)
        print()
    
    if fusion_result:
        summarize_fusion(fusion_result)

    print("\nSaved outputs:")
    if run_rule:
        print(f"  Rule JSON:   {rule_file}")
    if run_ml:
        print(f"  ML JSON:     {ml_file}")
    if run_fusion:
        print(f"  Fusion JSON: {fusion_file}")
    
    print("\nDone.")


if __name__ == "__main__":
    main()
