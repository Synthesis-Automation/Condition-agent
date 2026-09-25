"""Evidence inspection for agent investigators; no alternative ranking or admission."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import json
from typing import Any, Sequence

from reactive_taxonomy import featurize_reaction

from .generic_indexing import GenericIndexedReaction
from .recipe_assessment import _assess_analyzed_recipe
from .signature_features import environment_tokens


@dataclass(frozen=True)
class ConditionEvidenceComparison:
    """Selected observations with structural differences and unchanged recipe evidence."""

    query_reaction_smiles: str
    query_analysis: dict[str, Any]
    precedents: tuple[dict[str, Any], ...]
    observation_count: int
    distinct_reference_count: int
    missing_reference_observation_ids: tuple[str, ...]
    recipe_support: tuple[dict[str, Any], ...]
    limitations: tuple[str, ...] = (
        "Selected observations only; no new retrieval ranking or transfer-validity claim.",
        "Distinct indexed references are not proof of independent experiments or patent families.",
        "Missing indexed fields do not establish absence from the original publication.",
    )
    schema_version: str = "condition_evidence_comparison.v1"


def _difference(left: Sequence[Any], right: Sequence[Any]) -> dict[str, Any]:
    """Compare serialized observation multisets, retaining duplicates and values."""
    a = Counter(json.dumps(value, sort_keys=True) for value in left)
    b = Counter(json.dumps(value, sort_keys=True) for value in right)
    return {key: [json.loads(value) for value in sorted(counts.elements())]
            for key, counts in (("shared", a & b), ("query_only", a - b), ("precedent_only", b - a))}


def compare_condition_evidence(
    reaction_smiles: str, precedents: Sequence[GenericIndexedReaction],
) -> ConditionEvidenceComparison:
    """Expose graph observations and canonical compatibility before agent interpretation."""
    analysis = featurize_reaction(reaction_smiles)
    query = asdict(analysis.reaction_signature) if analysis.reaction_signature else None
    rows = []
    for row in precedents:
        differences = None
        if query and row.signature:
            differences = {key: _difference(query.get(key) or (), row.signature.get(key) or ())
                           for key in ("formed_bond_types", "broken_bond_types", "order_changes", "hydrogen_changes")}
            differences["environments"] = _difference(environment_tokens(query), environment_tokens(row.signature))
        rows.append({
            "observation": asdict(row),
            "structural_comparison": differences,
            "compatibility": asdict(_assess_analyzed_recipe(analysis, row.resolved_recipe)),
            "missing_operating_fields": [key for key in ("temperature_c", "time_h", "concentration_m", "atmosphere")
                                         if row.resolved_recipe.get(key) in (None, "")],
        })
    return ConditionEvidenceComparison(
        reaction_smiles, asdict(analysis), tuple(rows),
        len({row.observation_id for row in precedents}),
        len({row.reference_id for row in precedents if row.reference_id}),
        tuple(sorted({row.observation_id for row in precedents if not row.reference_id})),
        tuple({"recipe_id": identity,
               "observation_ids": sorted({row.observation_id for row in precedents if row.recipe_id == identity}),
               "reference_ids": sorted({row.reference_id for row in precedents if row.recipe_id == identity and row.reference_id}),
               "scope": "selected_observations_only"}
              for identity in sorted({row.recipe_id for row in precedents})),
    )
