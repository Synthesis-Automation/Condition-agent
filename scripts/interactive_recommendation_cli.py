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
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import requests



if sys.platform == "win32":
    # Ensure UTF-8 output when running in Windows terminals.
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")


DEFAULT_BASE_URL = "http://localhost:8000"
DEFAULT_SCDB_PATH = "data/conditionDB/Suzuki_db.json"
OUTPUT_DIR = Path("results")
K_DEFAULT = 50
LIMIT_DEFAULT = 5
FUSION_VARIANTS_DEFAULT = 5

# Canonical reaction families supported by the backend.
REACTION_TYPE_CHOICES: Tuple[Tuple[str, Optional[str]], ...] = (
    ("Auto-detect (server decides)", None),
    ("Suzuki Coupling", "Suzuki"),
    ("Ullmann C–N (Cu)", "C_N_Coupling_Cu"),
    ("Buchwald C–N (Pd)", "C_N_Coupling_Pd"),
    ("C–N Coupling (Ni)", "C_N_Coupling_Ni"),
    ("Amide Formation", "Amide_formation"),
)


def slugify_label(text: str) -> str:
    """Create a filesystem-friendly slug from the given text."""
    slug = "".join(ch.lower() if ch.isalnum() else "_" for ch in text)
    slug = slug.strip("_")
    return slug or "reaction"


def prompt_smiles() -> str:
    """Request a reaction SMILES string from the user."""
    while True:
        value = input("Enter reaction SMILES: ").strip()
        if value:
            return value
        print("Reaction SMILES is required. Please try again.")


def choose_reaction_type() -> Tuple[str, Optional[str]]:
    """Allow the user to choose a reaction family from predefined options."""
    print("\nReaction Type Options:")
    for idx, (label, _) in enumerate(REACTION_TYPE_CHOICES, start=1):
        default_marker = " (default)" if idx == 1 else ""
        print(f"  {idx}) {label}{default_marker}")
    
    while True:
        choice = input("Select reaction type [1]: ").strip()
        if not choice:
            return REACTION_TYPE_CHOICES[0]
        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(REACTION_TYPE_CHOICES):
                return REACTION_TYPE_CHOICES[idx - 1]
        print(f"Please enter a number between 1 and {len(REACTION_TYPE_CHOICES)}.")


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


def format_condition_value(value: Any) -> Optional[str]:
    """Render condition values (temperature/time) into human-readable strings."""
    if value is None:
        return None
    if isinstance(value, dict):
        val = value.get("value")
        unit = value.get("unit")
        text = value.get("text")
        if val is not None:
            if unit:
                return f"{val} {unit}"
            return str(val)
        if text:
            return str(text)
        # Fall back to JSON-style representation for complex payloads.
        try:
            return json.dumps(value)
        except TypeError:
            return str(value)
    if isinstance(value, (int, float)):
        return str(value)
    return str(value)


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
    if not isinstance(data, dict):
        print(f"Rule-based: unexpected response type ({type(data).__name__})")
        return
    if "error" in data:
        print(f"Rule-based: {data['error']}")
        return
    
    detection = data.get("detection", {})
    family = (
        detection.get("family")
        or detection.get("detected_reaction_type")
        or "Unknown"
    )
    print(f"Rule-based family: {family}")
    
    extras = data.get("extras", {})
    match_info: Dict[str, Any] = {}
    if isinstance(extras, dict):
        match_info = extras.get("match") or {}
    match_type = match_info.get("match_type") or "Unknown"
    print(f"Rule-based match type: {match_type}")
    
    entry = match_info.get("entry_name") or match_info.get("entry_id")
    if entry:
        print(f"Rule entry: {entry}")
    
    recommendations = data.get("recommended_conditions") or []
    if not recommendations:
        print("  No rule-based recommendations returned.")
        return
    
    top = recommendations[0]
    chemicals = top.get("chemicals", [])
    
    def find_role(role: str) -> Optional[str]:
        for chem in chemicals:
            chem_role = str(chem.get("role", "")).lower()
            if chem_role == role.lower():
                return chem.get("name") or chem.get("smiles")
        return None
    
    catalyst = find_role("metal_precursor") or find_role("catalyst")
    ligand = find_role("ligand")
    base = find_role("base")
    solvent = find_role("solvent")
    
    if catalyst:
        print(f"  Catalyst: {catalyst}")
    if ligand:
        print(f"  Ligand: {ligand}")
    if base:
        print(f"  Base: {base}")
    if solvent:
        print(f"  Solvent: {solvent}")
    
    conditions = top.get("conditions", {})
    temp_text = format_condition_value(conditions.get("temperature"))
    time_text = format_condition_value(conditions.get("time"))
    if temp_text:
        print(f"  Temperature: {temp_text}")
    if time_text:
        print(f"  Time: {time_text}")


def summarize_ml(data: Dict[str, Any]) -> None:
    """Print a short summary of the ML recommendation output."""
    if not isinstance(data, dict):
        print(f"ML recommendation: unexpected response type ({type(data).__name__})")
        return
    if "error" in data:
        print(f"ML recommendation: {data['error']}")
        return

    detection = data.get("detection", {})
    detected_type = (
        detection.get("family")
        or detection.get("detected_reaction_type")
        or detection.get("type")
        or "Unknown"
    )
    confidence = detection.get("confidence")
    confidence_text = "N/A"
    if isinstance(confidence, (int, float)):
        confidence_text = f"{format_float(confidence * 100, precision=0)}%"
    print(f"ML detected family: {detected_type} (confidence {confidence_text})")

    recommendations = data.get("recommended_conditions") or data.get("recommendations") or []
    if not recommendations:
        print("  No ML recommendations returned.")
        return

    top = recommendations[0]
    chemicals = top.get("chemicals", []) or top.get("reagents", [])
    def find_role(role: str) -> Optional[str]:
        for chem in chemicals:
            if str(chem.get("role", "")).lower() == role.lower():
                return chem.get("name") or chem.get("smiles")
        return None

    catalyst = find_role("catalyst") or find_role("metal_precursor")
    base = find_role("base")
    solvent = find_role("solvent")
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
        print(f"  Confidence: {format_float(float(confidence_score) * 100, precision=0)}%")
    if support is not None:
        print(f"  Support: {support} precedents")

    conditions = top.get("conditions", {})
    temp_text = format_condition_value(conditions.get("temperature"))
    time_text = format_condition_value(conditions.get("time"))
    if temp_text:
        print(f"  Temperature: {temp_text}")
    if time_text:
        print(f"  Time: {time_text}")


def summarize_fusion(data: Dict[str, Any]) -> None:
    """Print a short summary of the fusion recommendation output."""
    if not isinstance(data, dict):
        print(f"Fusion recommendation: unexpected response type ({type(data).__name__})")
        return
    if "error" in data:
        print(f"Fusion recommendation: {data['error']}")
        return

    extras = data.get("extras", {})
    fusion_meta = {}
    if isinstance(extras, dict):
        fusion_meta = extras.get("fusion_meta") or {}
    if not fusion_meta and isinstance(data.get("fusion_meta"), dict):
        fusion_meta = data["fusion_meta"]

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

    recommendations = data.get("recommended_conditions") or data.get("recommendations") or []
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
