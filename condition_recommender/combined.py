"""Combine admitted recommendations without mixing source scores or procedures.

The two canonical engines remain responsible for chemistry admission. This
projection retains all source outputs, including abstentions and failures.
Only identical, registry-identified recipe payloads are coalesced; core IDs
alone are insufficient because quantities and stages may differ.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Literal, Mapping, Tuple

from .protocols import protocol_draft_for_reaction


@dataclass(frozen=True)
class RecommendationSourceResult:
    """An engine's complete output or its explicit unavailable status."""

    source: Literal["generic", "weak_label"]
    status: Literal["ok", "abstained", "unavailable", "skipped"]
    result: Mapping[str, Any]
    message: str = ""
    error_code: str | None = None


@dataclass(frozen=True)
class CombinedConditionOption:
    """One intact recipe with independently attributed evidence."""

    option_id: str
    rank: int
    evidence_kind: str
    evidence_label: str
    resolved_recipe: Mapping[str, Any]
    synthesis_protocol: Mapping[str, Any]
    evidence: Tuple[Mapping[str, Any], ...]
    cautions: Tuple[str, ...]


@dataclass(frozen=True)
class CombinedRecommendationResult:
    """Versioned shortlist and source audit for the focused application."""

    query_reaction_smiles: str
    valid: bool
    recommendations: Tuple[CombinedConditionOption, ...]
    sources: Tuple[RecommendationSourceResult, ...]
    warnings: Tuple[str, ...]
    definition_id: str = "combined_recommendation.v1"
    schema_version: str = "1.0"

    def to_dict(self) -> dict[str, Any]:
        """Serialize while retaining every engine's uncertainty and provenance."""
        return asdict(self)


def load_combined_recommendation_policy() -> dict[str, Any]:
    """Load the bounded, versioned presentation policy; no executable imports."""
    path = Path(__file__).with_name("definitions") / "combined_recommendation.v1.json"
    policy = json.loads(path.read_text(encoding="utf-8"))
    if (
        policy.get("schema_version") != "1.0"
        or policy.get("definition_id") != "combined_recommendation.v1"
        or policy.get("source_order") != ["generic", "weak_label"]
        or set(policy.get("evidence_labels", {}))
        != {"verified_signature", "structure_review", "weak_label"}
        or not all(isinstance(v, str) and v for v in policy["evidence_labels"].values())
        or not 1
        <= policy.get("default_shortlist_size", 0)
        <= policy.get("candidate_limit_per_source", 0)
        <= 50
    ):
        raise ValueError("Invalid combined recommendation policy")
    return policy


def combine_recommendation_results(
    reaction_smiles: str,
    sources: Tuple[RecommendationSourceResult, ...],
) -> CombinedRecommendationResult:
    """Preserve source order and rank; never compare scores across engines."""
    policy = load_combined_recommendation_policy()
    by_source = {source.source: source for source in sources}
    if len(by_source) != len(sources):
        raise ValueError("Duplicate recommendation source")
    options: list[CombinedConditionOption] = []
    positions: dict[str, int] = {}
    warnings: set[str] = set()
    for source_name in policy["source_order"]:
        source = by_source.get(source_name)
        if source is None:
            continue
        result = source.result
        warnings.update(result.get("warnings") or ())
        if source.status != "ok" or not result.get("valid"):
            continue
        mode = str(result.get("recommendation_mode") or "unknown")
        kind = (
            "weak_label"
            if source_name == "weak_label"
            else "verified_signature"
            if mode == "verified_signature"
            else "structure_review"
        )
        effective_query = str(
            result.get("effective_query_reaction_smiles") or reaction_smiles
        )
        for item in result.get("recommendations") or ():
            recipe = dict(item.get("resolved_recipe") or {})
            if not recipe:
                continue
            # Conservative equality includes stages, uncertainty and quantities.
            # Never deduplicate a core recipe ID or an unresolved display name.
            identity = json.dumps(
                {"recipe": recipe, "reaction_smiles": effective_query},
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            )
            if not str(recipe.get("recipe_id") or "").startswith("RCR"):
                identity += f":{source_name}:{len(options)}"
            key = hashlib.sha256(identity.encode("utf-8")).hexdigest()
            evidence = {
                "source": source_name,
                "recommendation_mode": mode,
                "query_reaction_smiles": result.get(
                    "query_reaction_smiles", reaction_smiles
                ),
                "effective_query_reaction_smiles": effective_query,
                "query_signature_id": result.get("query_signature_id"),
                "retrieval_definition_version": result.get(
                    "retrieval_definition_version"
                ),
                "source_schema_version": result.get("schema_version"),
                "source_dataset_name": result.get("source_dataset_name"),
                "warnings": list(result.get("warnings") or ()),
                "recommendation": dict(item),
            }
            cautions = tuple(dict.fromkeys(item.get("cautions") or ()))
            if key in positions:
                position = positions[key]
                previous = options[position]
                options[position] = replace(
                    previous,
                    evidence=(*previous.evidence, evidence),
                    cautions=tuple(dict.fromkeys((*previous.cautions, *cautions))),
                )
                continue
            positions[key] = len(options)
            options.append(
                CombinedConditionOption(
                    option_id=f"CCO1:{key}",
                    rank=len(options) + 1,
                    evidence_kind=kind,
                    evidence_label=policy["evidence_labels"][kind],
                    resolved_recipe=recipe,
                    synthesis_protocol=protocol_draft_for_reaction(
                        recipe, effective_query
                    ).to_dict(),
                    evidence=(evidence,),
                    cautions=cautions,
                )
            )
    return CombinedRecommendationResult(
        query_reaction_smiles=reaction_smiles,
        valid=bool(options),
        recommendations=tuple(options),
        sources=sources,
        warnings=tuple(sorted(warnings)),
    )


def build_automation_handoff(
    result: CombinedRecommendationResult,
    option_ids: Tuple[str, ...],
) -> dict[str, Any]:
    """Export selected intact recipes for a robot adapter's planning input.

    This is deliberately not a hardware command format. Relative quantities,
    reported procedures and target-specific material amounts must not silently
    become executable dispensing instructions.
    """
    options = {item.option_id: item for item in result.recommendations}
    if not option_ids or len(set(option_ids)) != len(option_ids):
        raise ValueError("Select one or more distinct condition options")
    if any(option_id not in options for option_id in option_ids):
        raise ValueError("Unknown condition option")
    experiments = []
    for option_id in option_ids:
        item = options[option_id]
        experiments.append(
            {
                "option_id": option_id,
                "protocol": dict(item.synthesis_protocol),
                "resolved_recipe": dict(item.resolved_recipe),
                "evidence_kind": item.evidence_kind,
                "evidence": list(item.evidence),
                "cautions": list(item.cautions),
            }
        )
    payload = {
        "artifact_type": "condition_automation_handoff",
        "schema_version": "1.0",
        "definition_id": result.definition_id,
        "query_reaction_smiles": result.query_reaction_smiles,
        "execution_ready": False,
        "execution_status": "requires_robot_adapter_and_review",
        "robot_target": None,
        "required_before_execution": [
            "Confirm reaction scale and all material amounts with units",
            "Confirm ordered additions, vessel, mixing, quench and workup",
            "Map materials to robot inventory and validate equipment limits",
            "Compile and validate with the target robot's command schema",
        ],
        "warnings": list(result.warnings),
        "experiments": experiments,
    }
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    return {
        "handoff_id": "CAH1:" + hashlib.sha256(canonical.encode()).hexdigest(),
        **payload,
    }
