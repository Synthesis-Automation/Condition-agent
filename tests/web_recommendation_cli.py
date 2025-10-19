"""
Interactive Recommendation Tester
=================================

This script prompts for a reaction SMILES and lets you choose a reaction
family (or auto-detect). It exercises four recommendation modes
(rule-based, ML, fusion, protocol), saves their JSON responses, and prints
a compact summary to the console.

Usage:
    python scripts/web_recommendation_cli.py

Requirements:
    - FastAPI server running (default: http://localhost:8000)
    - requests package installed
    - Protocol index built (python -m chemtools.protocol.cli build)
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
        K_DEFAULT,
        LIMIT_DEFAULT,
        choose_reaction_type,
        choose_catalyst,
        prompt_smiles,
        save_json,
        slugify_label,
        summarize_ml,
        summarize_rule,
        summarize_protocol,
        summarize_llm_synthesis,
    )
except ModuleNotFoundError:
    sys.path.append(str(HERE))
    from recommendation_cli_utils import (
        DEFAULT_SCDB_PATH,
        K_DEFAULT,
        LIMIT_DEFAULT,
        choose_reaction_type,
        choose_catalyst,
        prompt_smiles,
        save_json,
        slugify_label,
        summarize_ml,
        summarize_rule,
        summarize_protocol,
        summarize_llm_synthesis,
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
    catalyst_preference: Optional[str] = None,
) -> Dict[str, Any]:
    """Execute the rule-based /match endpoint."""
    payload: Dict[str, Any] = {
        "reaction": reaction,
        "include_trace": True,
    }
    if db_path:
        payload["db"] = db_path
    
    # Add catalyst preference via relax parameter
    if catalyst_preference:
        payload["relax"] = {"catalyst_class": catalyst_preference}

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
    catalyst_preference: Optional[str] = None,
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
    
    # Add catalyst preference via relax parameter
    if catalyst_preference:
        payload["relax"]["catalyst_class"] = catalyst_preference

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


def call_llm_synthesis(
    base_url: str,
    reaction: str,
    family: Optional[str] = None,
    catalyst_preference: Optional[str] = None,
    timeout: int = 180,
) -> Dict[str, Any]:
    """Execute the LLM synthesis recommendation endpoint."""
    payload: Dict[str, Any] = {
        "reaction": reaction,
        "family": family,
        "catalyst_preference": catalyst_preference,
    }

    try:
        response = requests.post(
            f"{base_url}/api/v1/recommend/llm_synthesis",
            json=payload,
            timeout=timeout,
        )
        response.raise_for_status()
        return response.json()
    except requests.RequestException as exc:
        return {
            "error": f"LLM synthesis failed: {exc}",
            "payload": payload,
        }
    except ValueError:
        return {
            "error": "LLM synthesis response was not valid JSON.",
            "payload": payload,
        }


def call_protocol_recommendation(
    base_url: str,
    reaction: str,
    k: int,
    tags: Optional[list] = None,
) -> Dict[str, Any]:
    """
    Execute the protocol recommendation endpoint.
    
    Calls the local protocol module directly since there's no API endpoint yet.
    Returns standard JSON format matching other recommendation modes.
    """
    # Import here to avoid dependency if not using protocol mode
    try:
        import sys
        from pathlib import Path
        ROOT = Path(__file__).resolve().parent.parent
        if str(ROOT) not in sys.path:
            sys.path.insert(0, str(ROOT))
        
        from chemtools.protocol import ProtocolRecommender
        
        try:
            # Use local protocol recommendation
            recommender = ProtocolRecommender()
            result = recommender.recommend_with_details(
                reaction_smiles=reaction,
                k=k,
                tags=tags,
                reaction_family=None,  # Auto-detect from similarity
                use_standard_format=True
            )
            return result
        except FileNotFoundError as exc:
            return {
                "error": f"Protocol index not found. Run: python -m chemtools.protocol.cli build",
                "details": str(exc),
                "meta": {
                    "model": "Protocol-DRFP",
                    "status": "error"
                }
            }
        except Exception as exc:
            return {
                "error": f"Protocol recommendation failed: {exc}",
                "meta": {
                    "model": "Protocol-DRFP",
                    "status": "error"
                }
            }
    except ImportError:
        return {
            "error": "Protocol recommendation not available. Install with: pip install drfp",
            "meta": {
                "model": "Protocol-DRFP",
                "status": "error"
            }
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
  
  # Provide reaction, auto-detect type, specify catalyst:
  python scripts/web_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --catalyst Cu
  
  # Specify reaction type and catalyst:
  python scripts/web_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --family C_N_Coupling --catalyst Cu
  
  # Run only LLM synthesis with longer timeout:
  python scripts/web_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --strategy llm_synthesis --llm-timeout 300
  
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
        choices=[None, "Suzuki", "C_N_Coupling", "Amide_formation"],
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
        "--catalyst",
        type=str,
        choices=["None", "Pd", "Cu", "Ni", "other"],
        help="Catalyst preference: None, Pd, Cu, Ni, or other"
    )
    
    parser.add_argument(
        "--llm-timeout",
        type=int,
        default=180,
        help="Timeout in seconds for LLM API calls (default: 180)"
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
        choices=["all", "rule", "ml", "protocol", "llm_synthesis"],
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
    
    parser.add_argument(
        "--protocol-tags",
        type=str,
        default=None,
        help="Comma-separated tags to filter protocols (e.g., 'suzuki,palladium')"
    )
    
    args = parser.parse_args()
    
    print("Interactive Recommendation Test")
    print("--------------------------------")

    base_url = args.url
    print(f"Base URL: {base_url}")
    
    # Skip health check if only running protocol (it uses local module)
    if args.strategy != "protocol":
        if not health_check(base_url):
            print("Aborting due to failed health check.")
            return
    else:
        print("Protocol mode: using local module (no server required)")

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
    
    # Auto-detect reaction type if user selected "Auto-detect"
    if reaction_type is None:
        print("\nAuto-detecting reaction type...")
        from chemtools.router import detect_family_from_reaction
        family, confidence, hits = detect_family_from_reaction(reaction)
        if family:
            reaction_type = family
            print(f"✓ Detected: {family} (confidence: {confidence:.2f})")
        else:
            print("⚠ Could not auto-detect reaction type")
            reaction_type = None
    
    # Get catalyst preference - from args or prompt (MANDATORY)
    if args.catalyst:
        catalyst_label = args.catalyst
        catalyst_value = args.catalyst if args.catalyst != "None" else None
        print(f"Catalyst: {catalyst_label}")
    else:
        catalyst_label, catalyst_value = choose_catalyst()
        print(f"Selected catalyst: {catalyst_label}")

    k_value = args.k
    limit_value = args.limit
    llm_timeout = args.llm_timeout

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
    run_protocol = args.strategy in ["all", "protocol"]
    run_llm = args.strategy in ["all", "llm_synthesis"]
    
    rule_result = None
    ml_result = None
    protocol_result = None
    llm_result = None
    
    if run_rule:
        rule_result = call_rule_based(
            base_url, 
            reaction, 
            db_path,
            catalyst_preference=catalyst_value
        )
        rule_file = save_to_dir(rule_result, f"{timestamp}_{label}_rule.json")
    
    if run_ml:
        ml_result = call_ml_recommendation(
            base_url, 
            reaction, 
            reaction_type, 
            k_value, 
            limit_value,
            rerank_strategy=args.rerank,
            filter_unknown_reagents=args.filter_unknown,
            catalyst_preference=catalyst_value
        )
        ml_file = save_to_dir(ml_result, f"{timestamp}_{label}_ml.json")
    
    if run_protocol:
        # Parse protocol tags if provided
        protocol_tags = None
        if args.protocol_tags:
            protocol_tags = [tag.strip() for tag in args.protocol_tags.split(',')]
        
        protocol_result = call_protocol_recommendation(
            base_url=base_url,
            reaction=reaction,
            k=k_value,
            tags=protocol_tags
        )
        protocol_file = save_to_dir(protocol_result, f"{timestamp}_{label}_protocol.json")
    
    if run_llm:
        llm_result = call_llm_synthesis(
            base_url=base_url,
            reaction=reaction,
            family=reaction_type,
            catalyst_preference=catalyst_value,
            timeout=llm_timeout
        )
        llm_file = save_to_dir(llm_result, f"{timestamp}_{label}_llm.json")

    print("Summary\n-------")
    
    if rule_result:
        summarize_rule(rule_result)
        print()
    
    if ml_result:
        summarize_ml(ml_result)
        print()
    
    if protocol_result:
        summarize_protocol(protocol_result)
        print()
    
    if llm_result:
        summarize_llm_synthesis(llm_result)

    print("\nSaved outputs:")
    if run_rule:
        print(f"  Rule JSON:     {rule_file}")
    if run_ml:
        print(f"  ML JSON:       {ml_file}")
    if run_protocol:
        print(f"  Protocol JSON: {protocol_file}")
    if run_llm:
        print(f"  LLM JSON:      {llm_file}")
    
    print("\nDone.")


if __name__ == "__main__":
    main()
