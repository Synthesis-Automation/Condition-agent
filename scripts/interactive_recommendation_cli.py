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
) -> Dict[str, Any]:
    """Execute the ML recommendation endpoint."""
    payload: Dict[str, Any] = {
        "reaction": reaction,
        "reaction_type": reaction_type or None,
        "k": k,
        "limit": limit,
        "relax": {},
        "constraints": {},
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
    """Execute the fusion recommendation endpoint."""
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
    print("Interactive Recommendation Test")
    print("--------------------------------")

    base_url = DEFAULT_BASE_URL
    print(f"Base URL: {base_url}")
    if not health_check(base_url):
        print("Aborting due to failed health check.")
        return

    reaction = prompt_smiles()
    selected_label, reaction_type = choose_reaction_type()
    print(f"Selected reaction type: {selected_label}")

    k_value = K_DEFAULT
    limit_value = LIMIT_DEFAULT
    fusion_variants = FUSION_VARIANTS_DEFAULT

    db_path = DEFAULT_SCDB_PATH if Path(DEFAULT_SCDB_PATH).exists() else None
    if db_path is None:
        print("Warning: Default SCDB database not found; rule-based test will use server defaults.")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    type_label = reaction_type or "auto"
    label = slugify_label(type_label)

    print("\nRunning tests...\n")

    rule_result = call_rule_based(base_url, reaction, db_path)
    ml_result = call_ml_recommendation(base_url, reaction, reaction_type, k_value, limit_value)
    fusion_result = call_fusion_recommendation(base_url, reaction, k_value, fusion_variants)

    rule_file = save_json(rule_result, f"{timestamp}_{label}_rule.json")
    ml_file = save_json(ml_result, f"{timestamp}_{label}_ml.json")
    fusion_file = save_json(fusion_result, f"{timestamp}_{label}_fusion.json")

    print("Summary\n-------")
    summarize_rule(rule_result)
    print()
    summarize_ml(ml_result)
    print()
    summarize_fusion(fusion_result)

    print("\nSaved outputs:")
    print(f"  Rule JSON:   {rule_file}")
    print(f"  ML JSON:     {ml_file}")
    print(f"  Fusion JSON: {fusion_file}")
    print("\nDone.")


if __name__ == "__main__":
    main()
