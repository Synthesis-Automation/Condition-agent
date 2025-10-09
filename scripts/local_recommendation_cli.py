"""
Local Recommendation Tester
===========================

This script mirrors `interactive_recommendation_cli.py` but calls the
ChemTools APIs directly rather than going through the FastAPI server.
It prompts for a reaction SMILES and reaction family, runs rule-based,
ML, and fusion recommendation pipelines, saves their JSON results, and
prints a compact summary to the console.

Usage:
    python scripts/local_recommendation_cli.py

Requirements:
    - ChemTools library dependencies installed (same as the FastAPI app)
    - Optional SCDB JSON database for rule-based matching
"""

from __future__ import annotations

import io
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chemtools import chem, output_formatter, recommend

try:
    from chemtools.rule_scdb_matcher.loader import SchemeConditionDBError
except Exception:  # pragma: no cover - fallback path
    class SchemeConditionDBError(RuntimeError):
        """Fallback when SCDB loader is unavailable."""
        pass

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


_ENV_SCDB_DEFAULT = (
    os.environ.get("SCDB_MATCH_DB_PATH", "cn_coupling_pd_db.json").strip()
    or "cn_coupling_pd_db.json"
)


def _resolve_rule_db(candidate: Optional[str]) -> Optional[str]:
    """Determine the best available rule DB path."""
    if candidate:
        return candidate
    if Path(DEFAULT_SCDB_PATH).exists():
        return DEFAULT_SCDB_PATH
    if Path(_ENV_SCDB_DEFAULT).exists():
        return _ENV_SCDB_DEFAULT
    return None


def local_rule_based_match(reaction: str, db_path: Optional[str]) -> Dict[str, Any]:
    """Replicate the /match endpoint using in-process ChemTools calls."""
    target_db = _resolve_rule_db(db_path) or _ENV_SCDB_DEFAULT
    start_time = time.perf_counter()

    try:
        database = chem.rules.load_database(target_db, cache=True)
        result = chem.rules.match(database, reaction)
        payload = result.to_json_dict()
        processing_time_ms = round((time.perf_counter() - start_time) * 1000, 2)
        database_label = Path(target_db).name if target_db else "SchemeConditionDB"
        return output_formatter.format_rule_match_result(
            reaction_smiles=reaction,
            match_result=payload,
            requested_type=None,
            database_name=database_label,
            processing_time_ms=processing_time_ms,
        )
    except FileNotFoundError as exc:
        return {
            "error": f"Database file not found: {target_db}",
            "exception": repr(exc),
        }
    except SchemeConditionDBError as exc:
        return {"error": f"Rule-based match failed: {exc}"}
    except Exception as exc:  # pragma: no cover - unexpected path
        return {"error": f"Unexpected error during rule-based match: {exc}"}


def _format_conditions_output(
    raw_data: Dict[str, Any],
    reaction: str,
    reaction_type: Optional[str],
    limit: int,
    elapsed_ms: float,
) -> Dict[str, Any]:
    """Apply the same formatting as /api/v1/recommend/conditions."""
    detection = raw_data.get("detection", {})
    # Try multiple keys for detected type
    detected_type = (
        detection.get("family") 
        or detection.get("detected_reaction_type")
        or detection.get("reaction_type")
        or "Unknown"
    )
    confidence = detection.get("confidence", 0.0)

    recommendations_data: List[Dict[str, Any]] = []
    for rec in raw_data.get("recommendations", [])[:limit]:
        summary = rec.get("summary", {})
        chemicals = rec.get("chemicals", [])
        conditions_info = rec.get("conditions", {})

        reagents = []
        for chemical in chemicals:
            reagents.append(
                {
                    "id": chemical.get("uid", chemical.get("cas")),
                    "role": chemical.get("role", "reagent"),
                    "name": chemical.get("name"),
                    "abbreviation": None,
                    "cas": chemical.get("cas"),
                    "smiles": chemical.get("smiles"),
                    "equivalents": None,
                }
            )

        conditions: Dict[str, Any] = {}
        if conditions_info.get("temperature") is not None:
            conditions["temperature"] = {
                "value": conditions_info["temperature"],
                "unit": "°C",
            }
        if conditions_info.get("time") is not None:
            conditions["time"] = {
                "value": conditions_info["time"],
                "unit": "hours",
            }
        if conditions_info.get("atmosphere"):
            conditions["atmosphere"] = conditions_info["atmosphere"]

        recommendations_data.append(
            {
                "rank": rec.get("rank", len(recommendations_data) + 1),
                "confidence": (
                    summary.get("confidence", 0.0) / 100.0
                    if summary.get("confidence")
                    else 0.0
                ),
                "reagents": reagents,
                "conditions": conditions,
                "precedent_count": (
                    summary.get("support", {}).get("count")
                    if isinstance(summary.get("support"), dict)
                    else summary.get("support", 0)
                ),
            }
        )

    return output_formatter.format_ml_output(
        reaction_smiles=reaction,
        requested_type=reaction_type,
        detected_type=detected_type,
        detection_confidence=confidence,
        recommendations_data=recommendations_data,
        processing_time_ms=elapsed_ms,
    )


def local_ml_recommendation(
    reaction: str,
    reaction_type: Optional[str],
    k_value: int,
    limit: int,
) -> Dict[str, Any]:
    """Replicate the /api/v1/recommend/conditions endpoint locally."""
    start_time = time.perf_counter()
    try:
        raw_data = chem.recommend.conditions(
            reaction=reaction,
            reaction_type=reaction_type,
            k=k_value,
            limit=limit,
            relax={},
            constraints={},
        )
    except Exception as exc:
        return {"error": f"Local ML recommendation failed: {exc}"}

    elapsed_ms = (time.perf_counter() - start_time) * 1000
    try:
        return _format_conditions_output(raw_data, reaction, reaction_type, limit, elapsed_ms)
    except Exception as exc:  # pragma: no cover - formatting failure
        return {"error": f"Failed to format ML recommendation: {exc}"}


def local_fusion_recommendation(
    reaction: str,
    k_value: int,
    max_variants: int,
) -> Dict[str, Any]:
    """Replicate the /api/v1/recommend/fusion endpoint locally."""
    start_time = time.perf_counter()
    try:
        result = recommend.recommend_from_reaction(
            reaction=reaction,
            k=k_value,
            use_fusion=True,
            max_variants=max_variants,
            relax={},
            constraint_rules={},
        )
    except Exception as exc:
        return {"error": f"Local fusion recommendation failed: {exc}"}

    processing_time_ms = round((time.perf_counter() - start_time) * 1000, 2)

    if "fusion_meta" not in result:
        result["fusion_meta"] = {
            "error": "Fusion metadata not available",
            "note": "Check recommend.recommend_from_reaction fusion configuration",
        }

    return output_formatter.format_fusion_output(
        reaction_smiles=reaction,
        requested_type=None,
        fusion_result=result,
        processing_time_ms=processing_time_ms,
    )


def main() -> None:
    print("Local Recommendation Test")
    print("-------------------------")

    reaction = prompt_smiles()
    selected_label, reaction_type = choose_reaction_type()
    print(f"Selected reaction type: {selected_label}")

    k_value = K_DEFAULT
    limit_value = LIMIT_DEFAULT
    fusion_variants = FUSION_VARIANTS_DEFAULT

    db_path = DEFAULT_SCDB_PATH if Path(DEFAULT_SCDB_PATH).exists() else None
    if db_path:
        print(f"Using rule database: {db_path}")
    else:
        print("Rule database: default resolver (environment or ChemTools defaults)")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    type_label = reaction_type or "auto"
    label = slugify_label(type_label)

    print("\nRunning local pipelines...\n")

    rule_result = local_rule_based_match(reaction, db_path)
    ml_result = local_ml_recommendation(reaction, reaction_type, k_value, limit_value)
    fusion_result = local_fusion_recommendation(reaction, k_value, fusion_variants)

    rule_file = save_json(rule_result, f"{timestamp}_{label}_rule_local.json")
    ml_file = save_json(ml_result, f"{timestamp}_{label}_ml_local.json")
    fusion_file = save_json(fusion_result, f"{timestamp}_{label}_fusion_local.json")

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
