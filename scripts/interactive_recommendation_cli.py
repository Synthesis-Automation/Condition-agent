"""
Interactive Recommendation Tester
=================================

This script prompts for a reaction SMILES and optional metadata, exercises
three recommendation modes (rule-based, ML, fusion), saves their JSON
responses, and prints a compact summary to the console.

Usage:
    python scripts/interactive_recommendation_cli.py

Requirements:
    - FastAPI server running (default: http://localhost:8000)
    - requests package installed
"""

from __future__ import annotations

import io
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import requests


if sys.platform == "win32":
    # Ensure UTF-8 output when running in Windows terminals.
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")


DEFAULT_BASE_URL = "http://localhost:8000"
DEFAULT_SCDB_PATH = "data/conditionDB/Suzuki_db.json"
OUTPUT_DIR = Path("results")


def prompt(text: str, default: Optional[str] = None, required: bool = False) -> str:
    """Prompt the user for input, enforcing required fields and defaults."""
    while True:
        suffix = f" [{default}]" if default else ""
        value = input(f"{text}{suffix}: ").strip()
        if value:
            return value
        if default is not None:
            return default
        if not required:
            return ""
        print("Input required. Please try again.")


def ensure_output_dir() -> None:
    """Create the results directory if it does not already exist."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


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


def format_float(value: Any, precision: int = 3, fallback: str = "N/A") -> str:
    """Format numeric values safely for console output."""
    try:
        return f"{float(value):.{precision}f}"
    except (TypeError, ValueError):
        return fallback


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


def save_json(data: Dict[str, Any], filename: str) -> Path:
    """Save JSON data to the results directory."""
    ensure_output_dir()
    output_path = OUTPUT_DIR / filename
    output_path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return output_path


def summarize_rule(data: Dict[str, Any]) -> None:
    """Print a short summary of the rule-based output."""
    if "error" in data:
        print(f"Rule-based: {data['error']}")
        return

    match_type = data.get("match_type") or "Unknown"
    entry = data.get("entry_name") or data.get("entry_id") or "N/A"
    print(f"Rule-based match type: {match_type}")
    print(f"Rule-based entry: {entry}")

    conditions = data.get("conditions")
    if isinstance(conditions, list) and conditions:
        first = conditions[0]
    else:
        first = conditions if isinstance(conditions, dict) else None

    if isinstance(first, dict):
        core = first.get("core") or first.get("catalyst")
        base = first.get("base")
        solvent = first.get("solvent")
        if core:
            print(f"  Core: {core}")
        if base:
            print(f"  Base: {base}")
        if solvent:
            print(f"  Solvent: {solvent}")


def summarize_ml(data: Dict[str, Any]) -> None:
    """Print a short summary of the ML recommendation output."""
    if "error" in data:
        print(f"ML recommendation: {data['error']}")
        return

    detection = data.get("detection", {})
    detected_type = detection.get("type", "Unknown")
    confidence = detection.get("confidence", 0.0)
    print(f"ML detected type: {detected_type} (confidence {confidence:.0%})")

    recommendations = data.get("recommendations", []) or data.get("formatted", {}).get("recommended_conditions", [])
    if not recommendations:
        print("  No ML recommendations returned.")
        return

    top = recommendations[0]
    reagents = top.get("reagents", [])
    catalyst = next((r.get("name") for r in reagents if r.get("role", "").lower() == "catalyst"), None)
    base = next((r.get("name") for r in reagents if r.get("role", "").lower() == "base"), None)
    solvent = next((r.get("name") for r in reagents if r.get("role", "").lower() == "solvent"), None)
    confidence_score = top.get("confidence")
    support = top.get("precedent_count")

    print("ML top recommendation:")
    if catalyst:
        print(f"  Catalyst: {catalyst}")
    if base:
        print(f"  Base: {base}")
    if solvent:
        print(f"  Solvent: {solvent}")
    if confidence_score is not None:
        print(f"  Confidence: {confidence_score:.0%}")
    if support is not None:
        print(f"  Support: {support} precedents")

    conditions = top.get("conditions", {})
    if isinstance(conditions, dict):
        if "temperature" in conditions:
            temp_val = conditions["temperature"]
            if isinstance(temp_val, dict):
                temp_val = temp_val.get("value")
            print(f"  Temperature: {temp_val} C")
        if "time" in conditions:
            time_val = conditions["time"]
            if isinstance(time_val, dict):
                time_val = time_val.get("value")
            print(f"  Time: {time_val} hours")


def summarize_fusion(data: Dict[str, Any]) -> None:
    """Print a short summary of the fusion recommendation output."""
    if "error" in data:
        print(f"Fusion recommendation: {data['error']}")
        return

    fusion_meta = data.get("fusion_meta", {})
    if fusion_meta and "error" not in fusion_meta:
        weights = fusion_meta.get("adaptive_weights", {})
        if weights:
            alpha = weights.get("alpha")
            beta = weights.get("beta")
            gamma = weights.get("gamma")
            delta = weights.get("delta")
            print("Fusion adaptive weights:")
            pre = format_float(alpha)
            ana = format_float(beta)
            rul = format_float(gamma)
            mlw = format_float(delta)
            print(f"  precedents={pre} analytics={ana} rules={rul} ml={mlw}")
    elif fusion_meta.get("error"):
        print(f"Fusion metadata error: {fusion_meta['error']}")

    formatted = data.get("formatted", {})
    recommendations = formatted.get("recommended_conditions") or data.get("recommended_conditions") or []
    if not recommendations:
        print("  No fusion recommendations returned.")
        return

    top = recommendations[0]
    summary = top.get("summary", {})
    core = summary.get("core", "N/A")
    base = summary.get("base")
    solvent = summary.get("solvent")
    confidence = summary.get("confidence")

    print("Fusion top recommendation:")
    print(f"  Core: {core}")
    if isinstance(base, dict):
        base = base.get("name") or base.get("id")
    if isinstance(solvent, dict):
        solvent = solvent.get("name") or solvent.get("id")
    if base:
        print(f"  Base: {base}")
    if solvent:
        print(f"  Solvent: {solvent}")
    if confidence is not None:
        print(f"  Confidence: {confidence}")


def parse_int(value: str, default: int) -> int:
    """Convert input to int, falling back to default if conversion fails."""
    try:
        return int(value)
    except ValueError:
        print(f"Invalid integer '{value}'. Using default {default}.")
        return default


def main() -> None:
    print("Interactive Recommendation Test")
    print("--------------------------------")

    base_url = prompt("Base URL", default=DEFAULT_BASE_URL, required=True)
    if not health_check(base_url):
        print("Aborting due to failed health check.")
        return

    reaction = prompt("Enter reaction SMILES", required=True)
    description = prompt("Optional description")
    reaction_type = prompt("Optional reaction type (leave blank for auto-detect)")

    k_input = prompt("Number of nearest precedents (k)", default="50")
    limit_input = prompt("Number of recommendations to return", default="5")
    fusion_variants_input = prompt("Max fusion variants", default="5")

    k_value = parse_int(k_input, 50)
    limit_value = parse_int(limit_input, 5)
    fusion_variants = parse_int(fusion_variants_input, 5)

    db_path_input = prompt(
        "SCDB path for rule-based lookup (blank to skip)",
        default=DEFAULT_SCDB_PATH if Path(DEFAULT_SCDB_PATH).exists() else "",
    )
    db_path = db_path_input if db_path_input else None
    if db_path and not Path(db_path).exists():
        print(f"Warning: SCDB file '{db_path}' not found. Rule-based test may fail.")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    label = description.replace(" ", "_") if description else "reaction"

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
