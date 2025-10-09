"""
Utility helpers shared by recommendation CLI tools.

These functions keep formatting, prompting, and summary behavior
consistent between HTTP-driven and local-integration scripts.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

OUTPUT_DIR = Path("results")
DEFAULT_SCDB_PATH = "data/conditionDB/Suzuki_db.json"
K_DEFAULT = 50
LIMIT_DEFAULT = 5
FUSION_VARIANTS_DEFAULT = 5

# Canonical reaction families supported by the backends.
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


def save_json(data: Dict[str, Any], filename: str) -> Path:
    """Save JSON data to the results directory."""
    ensure_output_dir()
    output_path = OUTPUT_DIR / filename
    output_path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return output_path


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
    if isinstance(value, (int, float)):
        return str(value)
    return str(value)


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

    # Use the standard 'recommended_conditions' key
    recommendations = data.get("recommended_conditions", [])
    if not recommendations:
        print("  No rule matches were returned.")
        return

    top_match = recommendations[0]
    rank = top_match.get("rank", 1)
    
    # Extract information from the match
    chemicals = top_match.get("chemicals", [])
    conditions = top_match.get("conditions", {})
    source = top_match.get("source", {})
    
    print(f"Top match (rank {rank}):")
    
    if source:
        entry_name = source.get("entry_name", "N/A")
        entry_id = source.get("entry_id", "N/A")
        match_type = source.get("match_type", "N/A")
        print(f"  Entry: {entry_name} (ID: {entry_id})")
        print(f"  Match type: {match_type}")
    
    # Show key reagents
    if chemicals:
        catalyst = next((c for c in chemicals if c.get("role") in ["catalyst", "metal_precursor"]), None)
        if catalyst:
            cat_name = catalyst.get("name") or catalyst.get("abbreviation") or "Unknown"
            print(f"  Catalyst: {cat_name}")
    
    # Show conditions
    temp = conditions.get("temperature")
    if temp:
        temp_val = temp.get("value") if isinstance(temp, dict) else temp
        print(f"  Temperature: {temp_val}°C")


def summarize_ml(data: Dict[str, Any]) -> None:
    """Print a short summary of the ML recommendation output."""
    if not isinstance(data, dict):
        print(f"ML recommendation: unexpected response type ({type(data).__name__})")
        return
    if "error" in data:
        print(f"ML recommendation: {data['error']}")
        return

    # Detection info is in the 'detection' key, not 'meta'
    detection = data.get("detection", {})
    detected_type = (
        detection.get("family")
        or detection.get("detected_reaction_type")
        or "Unknown"
    )
    confidence = detection.get("confidence")
    confidence_str = format_float(confidence, precision=2) if confidence is not None else "N/A"
    print(f"ML detected type: {detected_type} (confidence {confidence_str})")

    recommendations = data.get("recommended_conditions", [])
    if not recommendations:
        print("  No ML recommendations returned.")
        return

    top = recommendations[0]
    reagents = top.get("reagents") or []
    conditions = top.get("conditions") or {}
    confidence_score = top.get("confidence")
    support = top.get("precedent_count")

    print("Top ML recommendation:")
    print(f"  Confidence: {format_float(confidence_score, precision=2)}")
    if support is not None:
        print(f"  Precedent count: {support}")
    if reagents:
        first = reagents[0]
        name = first.get("name") or first.get("id") or "Unknown"
        print(f"  Primary reagent: {name}")

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

    # Try to extract fusion_meta from extras or top-level
    extras = data.get("extras", {})
    fusion_meta = {}
    if isinstance(extras, dict):
        fusion_meta = extras.get("fusion_meta") or {}
    if not fusion_meta and isinstance(data.get("fusion_meta"), dict):
        fusion_meta = data["fusion_meta"]

    if fusion_meta and "error" not in fusion_meta:
        weights = fusion_meta.get("adaptive_weights", {})
        if weights:
            # Weights might be under different keys: Greek letters (α,β,γ,δ) or English names
            alpha = (weights.get("alpha") or weights.get("α") or 
                    weights.get("precedents"))
            beta = (weights.get("beta") or weights.get("β") or 
                   weights.get("analytics"))
            gamma = (weights.get("gamma") or weights.get("γ") or 
                    weights.get("rules"))
            delta = (weights.get("delta") or weights.get("δ") or 
                    weights.get("ml"))
            
            print("Fusion adaptive weights:")
            pre = format_float(alpha)
            ana = format_float(beta)
            rul = format_float(gamma)
            mlw = format_float(delta)
            print(f"  precedents={pre} analytics={ana} rules={rul} ml={mlw}")
        else:
            print("Fusion weights: Not available in fusion_meta")
    elif fusion_meta.get("error"):
        print(f"Fusion metadata error: {fusion_meta['error']}")
    else:
        print("Fusion metadata: Not available")

    recommendations = data.get("recommended_conditions", [])
    if not recommendations:
        print("  No fusion recommendations returned.")
        return

    top = recommendations[0]
    
    # Try different structures for summary
    summary = top.get("summary", {})
    if summary:
        core = summary.get("core", "N/A")
        base = summary.get("base")
        solvent = summary.get("solvent")
        confidence = summary.get("confidence")
    else:
        # Fallback to extracting from chemicals
        chemicals = top.get("chemicals", [])
        core = "N/A"
        base = None
        solvent = None
        confidence = top.get("confidence")
        
        for chem in chemicals:
            role = chem.get("role", "")
            name = chem.get("name") or chem.get("abbreviation")
            if role == "base" and not base:
                base = name
            elif role == "solvent" and not solvent:
                solvent = name

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
