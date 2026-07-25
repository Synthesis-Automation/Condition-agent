"""Run weak-label condition recommendation for every reaction in a CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

from condition_recommender import recommend_conditions_from_labels


OUTPUT_FIELDNAMES = (
    "input_row",
    "reaction_type",
    "example",
    "rxn_smiles",
    "valid",
    "error",
    "warnings",
    "query_label",
    "grammar_id",
    "query_signatures",
    "candidate_count",
    "recipe_count",
    "recommendation_rank",
    "recipe_id",
    "score",
    "label_similarity",
    "signature_similarity",
    "qualifier_similarity",
    "expected_yield_pct",
    "mean_z_score",
    "support",
    "source_reaction_types",
    "source_row_numbers",
    "condition_display",
    "resolved_components",
    "temperature_c",
    "time_h",
    "concentration_m",
    "atmosphere",
    "condition_identity_uncertainty",
    "declared_absences",
    "recipe_warnings",
    "procedure_text",
    "explanation",
)

_RECIPE_BUCKETS = (
    "catalysts",
    "ligands",
    "bases",
    "acids",
    "condensation_agents",
    "oxidants",
    "reductants",
    "additives",
    "solvents",
    "other_components",
)


def _joined(values: Iterable[Any]) -> str:
    return "; ".join(str(value) for value in values)


def _resolved_components(recipe: Mapping[str, Any]) -> str:
    values = []
    for bucket in _RECIPE_BUCKETS:
        for component in recipe.get(bucket, ()):
            name = (
                component.get("canonical_name")
                or component.get("raw_identifier")
                or "unknown"
            )
            role = component.get("primary_role") or "unresolved_role"
            identity = component.get("substance_id")
            identity_suffix = f" <{identity}>" if identity else ""
            values.append(f"{name} [{role}]{identity_suffix}")
    return _joined(values)


def _base_output(
    input_row: int,
    source: Mapping[str, str],
    result: Any,
) -> Dict[str, Any]:
    return {
        "input_row": input_row,
        "reaction_type": str(source.get("reaction_type") or "").strip(),
        "example": str(source.get("example") or "").strip(),
        "rxn_smiles": str(source.get("rxn_smiles") or "").strip(),
        "valid": str(bool(result.valid)).lower(),
        "error": result.error or "",
        "warnings": _joined(result.warnings),
        "query_label": result.query_label or "",
        "grammar_id": result.grammar_id or "",
        "query_signatures": _joined(result.query_signatures),
        "candidate_count": result.candidate_count,
        "recipe_count": result.recipe_count,
    }


def recommend_rows(
    rows: Iterable[Mapping[str, str]],
    *,
    records_path: str | Path,
    top_k: int,
) -> list[Dict[str, Any]]:
    """Return long-form recommendation and diagnostic rows."""
    output = []
    for input_row, source in enumerate(rows, start=2):
        reaction_smiles = str(source.get("rxn_smiles") or "").strip()
        result = recommend_conditions_from_labels(
            reaction_smiles,
            records_path=records_path,
            top_k=top_k,
        )
        base = _base_output(input_row, source, result)
        if not result.recommendations:
            output.append(base)
            continue
        for recommendation in result.recommendations:
            recipe = recommendation.resolved_recipe
            conditions = recommendation.conditions
            output.append(
                {
                    **base,
                    "recommendation_rank": recommendation.rank,
                    "recipe_id": recommendation.recipe_id,
                    "score": recommendation.score,
                    "label_similarity": recommendation.label_similarity,
                    "signature_similarity": recommendation.signature_similarity,
                    "qualifier_similarity": recommendation.qualifier_similarity,
                    "expected_yield_pct": recommendation.expected_yield_pct,
                    "mean_z_score": recommendation.mean_z_score,
                    "support": recommendation.support,
                    "source_reaction_types": _joined(
                        recommendation.source_reaction_types
                    ),
                    "source_row_numbers": _joined(
                        recommendation.source_row_numbers
                    ),
                    "condition_display": conditions.get(
                        "condition_display", ""
                    ),
                    "resolved_components": _resolved_components(recipe),
                    "temperature_c": conditions.get("temperature_c", ""),
                    "time_h": conditions.get("time_h", ""),
                    "concentration_m": conditions.get(
                        "concentration_m", ""
                    ),
                    "atmosphere": conditions.get("atmosphere", ""),
                    "condition_identity_uncertainty": conditions.get(
                        "condition_identity_uncertainty", ""
                    ),
                    "declared_absences": _joined(
                        recipe.get("declared_absences", ())
                    ),
                    "recipe_warnings": _joined(
                        recipe.get("warnings", ())
                    ),
                    "procedure_text": conditions.get("procedure_text", ""),
                    "explanation": _joined(recommendation.explanation),
                }
            )
    return output


def recommend_csv(
    source: Path,
    destination: Path,
    *,
    records_path: Path,
    top_k: int,
) -> Dict[str, int]:
    """Read reactions, write recommendations, and return coverage counts."""
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"reaction_type", "example", "rxn_smiles"}
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"Missing reaction input columns: {missing}")
        source_rows = list(reader)

    output_rows = recommend_rows(
        source_rows,
        records_path=records_path,
        top_k=top_k,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDNAMES,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    recommended_inputs = len(
        {
            int(row["input_row"])
            for row in output_rows
            if row.get("recommendation_rank")
        }
    )
    return {
        "input_rows": len(source_rows),
        "output_rows": len(output_rows),
        "recommended_inputs": recommended_inputs,
        "unrecommended_inputs": len(source_rows) - recommended_inputs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument(
        "--records",
        type=Path,
        default=Path("datasets/reaction_label/v2.1_cleaned.csv"),
    )
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()
    if args.top_k < 1:
        parser.error("--top-k must be positive")
    stats = recommend_csv(
        args.source,
        args.destination,
        records_path=args.records,
        top_k=args.top_k,
    )
    for key in sorted(stats):
        print(f"{key}: {stats[key]}")


if __name__ == "__main__":
    main()
