"""
Lightweight local rule-based matcher used by regression tests.

The original implementation lived in ``app.local_recommendation_cli``, but that
module is not present in this trimmed repo layout.  The tests only require a
deterministic matcher that:
    - Detects a reaction family
    - Loads the corresponding rule DB JSON
    - Returns a non-empty ``recommended_conditions`` list
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from chemtools import detect_reaction

_RULE_DIR = Path("data") / "rule_db_v2"
_FAMILY_TO_FILE = {
    "suzuki_miyaura": "Suzuki_db.json",
    "cn_coupling": "C_N_Coupling_Cu_db.json",
    "ullmann_cn": "C_N_Coupling_Cu_db.json",
    "buchwald_hartwig_c_n": "C_N_Coupling_Cu_db.json",
    "amide_coupling": "Amide_db.json",
}


def _load_rule_db(family: str) -> Optional[Dict[str, Any]]:
    """Load a rule DB JSON for the requested family."""
    fname = _FAMILY_TO_FILE.get(family)
    if not fname:
        return None
    path = _RULE_DIR / fname
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _build_recommendations(rule_db: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Convert base_rules into a simple recommended_conditions list."""
    recommendations: List[Dict[str, Any]] = []
    for idx, base in enumerate(rule_db.get("base_rules", []), start=1):
        entry_name = base.get("description") or base.get("id") or f"rule_{idx}"
        recommendations.append(
            {
                "score": 1.0,
                "conditions": base.get("conditions", {}),
                "source": {
                    "entry_name": entry_name,
                    "rule_id": base.get("id"),
                    "family": (rule_db.get("reaction") or {}).get("family"),
                    "db_file": rule_db.get("metadata", {}).get("name"),
                },
            }
        )
    if not recommendations and rule_db.get("default_rule"):
        default = rule_db["default_rule"]
        entry_name = default.get("description") or "default_rule"
        recommendations.append(
            {
                "score": 0.5,
                "conditions": default.get("conditions", {}),
                "source": {
                    "entry_name": entry_name,
                    "rule_id": default.get("id"),
                    "family": (rule_db.get("reaction") or {}).get("family"),
                    "db_file": rule_db.get("metadata", {}).get("name"),
                },
            }
        )
    return recommendations


def local_rule_based_match(
    reaction_smiles: str,
    catalyst_preference: Optional[str] = None,
    reaction_type: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Minimal rule-based matcher for tests.

    Args:
        reaction_smiles: reaction SMILES string
        catalyst_preference: unused placeholder
        reaction_type: optional explicit family (e.g., "cn_coupling")
    """
    detection = detect_reaction(reaction_smiles, use_ml=False) or {}
    family = reaction_type or detection.get("family")

    rule_db = _load_rule_db(family) if family else None
    if not rule_db:
        return {
            "input": {"reaction_smiles": reaction_smiles},
            "detection": {"family": family or "unknown"},
            "recommended_conditions": [],
            "error": f"No rule database found for family '{family}'",
        }

    recs = _build_recommendations(rule_db)
    return {
        "input": {"reaction_smiles": reaction_smiles},
        "detection": {"family": family, "confidence": detection.get("confidence")},
        "recommended_conditions": recs,
    }


__all__ = ["local_rule_based_match"]
