"""
Simplify `chemtools/taxonomy/data/calculable_features.json`.

Goal: keep the calculable feature system maintainable by:
  - Dropping smartsrx-derived reactant-type features that are not mapped to the
    unified taxonomy reactant types (these include PFAS and heteroaryl leaving
    group "super-patterns").
  - Keeping only atomic, taxonomy-mapped reactant member detectors
    (`reactant_metadata` present).
  - Pruning derived shortcuts that reference removed tokens.

This script is intentionally conservative: it preserves all non-reactant-type
features and only removes the high-maintenance reactant-type expansion layer.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Iterable


_OPERATORS = {"AND", "OR", "NOT", ">=", "<=", ">", "<", "==", "!="}


def _extract_rule_tokens(expr: str) -> list[str]:
    """
    Extract token-like references from a boolean expression.

    Notes:
      - This is a lightweight extractor for `derive` expressions used by
        `chemtools.util.boolean_expr`.
      - It ignores numeric literals and operators and strips grouping parens.
    """
    tokens: list[str] = []
    for part in str(expr).split():
        if part in _OPERATORS:
            continue
        cleaned = part.strip().strip("()")
        if not cleaned or cleaned in {"True", "False"}:
            continue
        if re.fullmatch(r"\d+", cleaned):
            continue
        tokens.append(cleaned)
    return tokens


def _prune_reactant_type_expansion(features: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    Keep all features except reactant_types entries without reactant_metadata.

    The pruned layer corresponds to smartsrx-derived reactant type expansions,
    including PFAS and heteroaryl leaving group "super-patterns".
    """
    kept: list[dict[str, Any]] = []
    for feature in features:
        if feature.get("category") == "reactant_types" and not isinstance(
            feature.get("reactant_metadata"), dict
        ):
            continue
        kept.append(feature)
    return kept


def _prune_problematic_derived_features(features: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    Remove or simplify derived features that are invalid in the current evaluator.

    - `aryl_coupling_handle_site_count`: uses arithmetic `+` (not supported).
    - `beta_elimination_risk_high`: references tokens not present in the spec.
    - `water_sensitive_functionality_present`: drop missing `chloroformate_present`.
    """
    result: list[dict[str, Any]] = []
    for feature in features:
        token = feature.get("token")
        if token in {"aryl_coupling_handle_site_count", "beta_elimination_risk_high"}:
            continue
        if token == "water_sensitive_functionality_present":
            expr = feature.get("derive") or feature.get("derived") or ""
            expr = str(expr).replace(" OR chloroformate_present", "")
            expr = expr.replace("chloroformate_present OR ", "")
            expr = expr.strip()
            if expr:
                if "derive" in feature:
                    feature["derive"] = expr
                if "derived" in feature:
                    feature["derived"] = expr
            result.append(feature)
            continue
        result.append(feature)
    return result


def _filter_derived_shortcuts(
    derived_shortcuts: Iterable[dict[str, Any]],
    *,
    available_tokens: set[str],
) -> list[dict[str, Any]]:
    """
    Keep only derived shortcuts whose referenced tokens exist.

    This prevents dangling handle/category shortcuts (e.g. PFAS handle) after
    pruning the smartsrx-derived reactant-type layer.
    """
    kept: list[dict[str, Any]] = []
    known = set(available_tokens)

    for entry in derived_shortcuts:
        token = entry.get("token")
        expr = entry.get("derive") or ""
        if not token or not isinstance(expr, str) or not expr.strip():
            continue

        refs = _extract_rule_tokens(expr)
        if any(ref not in known for ref in refs):
            continue

        kept.append(entry)
        known.add(token)

    return kept


def simplify_spec(payload: dict[str, Any]) -> dict[str, Any]:
    """Return a simplified calculable feature spec payload."""
    features = payload.get("features", [])
    derived_shortcuts = payload.get("derived_shortcuts", [])

    if not isinstance(features, list):
        raise ValueError("Expected top-level 'features' to be a list.")
    if derived_shortcuts is None:
        derived_shortcuts = []
    if not isinstance(derived_shortcuts, list):
        raise ValueError("Expected top-level 'derived_shortcuts' to be a list.")

    features = [f for f in features if isinstance(f, dict)]
    features = _prune_reactant_type_expansion(features)
    features = _prune_problematic_derived_features(features)

    available_tokens = {f.get("token") for f in features if f.get("token")}

    pruned_shortcuts = _filter_derived_shortcuts(
        (d for d in derived_shortcuts if isinstance(d, dict)),
        available_tokens=available_tokens,
    )

    # Update metadata to reflect simplification (keep the rest intact).
    updated = dict(payload)
    updated["features"] = features
    updated["derived_shortcuts"] = pruned_shortcuts
    updated["version"] = "v6.0_simple"
    updated["description"] = (
        "Simplified, atomic SMARTS feature library for deterministic chemistry tooling. "
        "Drops smartsrx-derived reactant-type expansions (including PFAS and heteroaryl leaving-group super-patterns); "
        "keeps taxonomy-mapped reactant members plus lightweight derived logic."
    )

    schema_notes = updated.get("schema_notes")
    if isinstance(schema_notes, dict):
        notes = dict(schema_notes)
        expansion = list(notes.get("expansion_notes") or [])
        # Replace the old multi-version narrative with a single, clear statement.
        notes["expansion_notes"] = [
            "Atomic SMARTS only in the base layer; build chemistry meaning in derived logic.",
            "Removed smartsrx-derived reactant-type expansions (PFAS + heteroaryl leaving-group super-patterns).",
            "Reactant-type features are limited to unified taxonomy members (reactant_metadata present).",
        ]
        updated["schema_notes"] = notes

    return updated


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("chemtools/taxonomy/data/calculable_features.json"),
        help="Path to the source calculable_features.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output path (defaults to overwriting --input).",
    )
    args = parser.parse_args()

    source = args.input.resolve()
    target = (args.output or args.input).resolve()

    payload = json.loads(source.read_text(encoding="utf-8"))
    simplified = simplify_spec(payload)
    target.write_text(json.dumps(simplified, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote simplified spec: {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
