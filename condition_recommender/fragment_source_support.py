"""Match product-only fragment requirements to reported condition components."""

from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Tuple

from condition_registry import (
    ResolvedConditionRecipe,
    resolve_substance_id,
)

from .models import (
    FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION,
    FragmentSourceSupport,
)


_RULES_PATH = (
    Path(__file__).with_name("definitions")
    / "fragment_source_capabilities.v1.json"
)
_SOURCE_FIELD_EQUIVALENTS = {
    "catalysts_json": "catalyst_cas",
    "reagents_json": "reagent_cas",
}


def _mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    if is_dataclass(value):
        return asdict(value)
    raise TypeError("fragment source requirement must be a mapping or dataclass")


@lru_cache(maxsize=1)
def load_fragment_source_capabilities() -> dict[str, Any]:
    """Load and validate curated transferable-fragment capabilities."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.1":
        raise ValueError("unsupported fragment-source capability schema")
    if str(rules.get("definition_id") or "") != (
        "fragment_source_capabilities.v1"
    ):
        raise ValueError("unexpected fragment-source capability definition")
    capabilities = tuple(dict(value) for value in rules.get("capabilities") or ())
    if not capabilities:
        raise ValueError("fragment-source capabilities must not be empty")
    seen = set()
    allowed_match_fields = {
        "element_counts",
        "maximum_atom_count",
        "rooted_fragment_smiles_any",
        "fragment_keys_any",
        "attachment_elements_any",
    }
    for capability in capabilities:
        capability_id = str(capability.get("capability_id") or "")
        if not capability_id or capability_id in seen:
            raise ValueError("fragment-source capability IDs must be unique")
        seen.add(capability_id)
        if not str(capability.get("display_name") or "").strip():
            raise ValueError(
                f"fragment-source capability {capability_id} needs a display name"
            )
        match = capability.get("match")
        if not isinstance(match, Mapping) or not match:
            raise ValueError(
                f"fragment-source capability {capability_id} needs match rules"
            )
        unknown_match = set(match) - allowed_match_fields
        if unknown_match:
            raise ValueError(
                f"unsupported fragment-source match fields: {unknown_match}"
            )
        maximum = int(match.get("maximum_atom_count") or 0)
        if maximum < 1:
            raise ValueError(
                f"fragment-source capability {capability_id} needs atom bound"
            )
        for substance_id in capability.get("substance_ids") or ():
            if resolve_substance_id(str(substance_id)).status != "resolved":
                raise ValueError(
                    f"unknown fragment-source substance: {substance_id}"
                )
        source_fields = {
            str(value) for value in capability.get("allowed_source_fields") or ()
        }
        if not source_fields or not source_fields <= {
            "catalyst_cas",
            "reagent_cas",
        }:
            raise ValueError(
                f"invalid source fields for fragment-source {capability_id}"
            )
    rules["capabilities"] = capabilities
    return rules


def _requirement_matches(
    requirement: Mapping[str, Any],
    match: Mapping[str, Any],
) -> bool:
    if int(requirement.get("atom_count") or 0) > int(
        match.get("maximum_atom_count") or 0
    ):
        return False
    configured_counts = match.get("element_counts")
    if configured_counts and {
        str(key): int(value) for key, value in configured_counts.items()
    } != {
        str(key): int(value)
        for key, value in (requirement.get("element_counts") or {}).items()
    }:
        return False
    rooted_values = {
        str(value) for value in match.get("rooted_fragment_smiles_any") or ()
    }
    if rooted_values and str(
        requirement.get("rooted_fragment_smiles") or ""
    ) not in rooted_values:
        return False
    fragment_keys = {
        str(value) for value in match.get("fragment_keys_any") or ()
    }
    if fragment_keys and str(requirement.get("fragment_key") or "") not in (
        fragment_keys
    ):
        return False
    attachment_elements = {
        str(value) for value in match.get("attachment_elements_any") or ()
    }
    if attachment_elements and str(
        requirement.get("attachment_element") or ""
    ) not in attachment_elements:
        return False
    return True


def _component_match(
    component: Any,
    capability: Mapping[str, Any],
) -> tuple[str, float] | None:
    source_fields = {
        _SOURCE_FIELD_EQUIVALENTS.get(value, value)
        for value in str(component.source_field).split("|")
    }
    if not source_fields.intersection(
        str(value) for value in capability.get("allowed_source_fields") or ()
    ):
        return None
    substance_id = str(component.substance_id or "")
    raw_identifier = str(component.raw_identifier or "")
    if substance_id and substance_id in {
        str(value) for value in capability.get("substance_ids") or ()
    }:
        return "curated_substance_capability", 0.95
    if raw_identifier in {
        str(value) for value in capability.get("raw_identifiers") or ()
    }:
        return "curated_raw_identifier_capability", 0.85
    return None


def assess_fragment_source_support(
    requirements: Iterable[Any],
    recipe: ResolvedConditionRecipe,
) -> Tuple[FragmentSourceSupport, ...]:
    """Assess every structural source requirement against one reported recipe."""
    capabilities = load_fragment_source_capabilities()["capabilities"]
    results = []
    for raw_requirement in requirements:
        requirement = _mapping(raw_requirement)
        matches: list[tuple[Any, str, str, float]] = []
        for capability in capabilities:
            if not _requirement_matches(requirement, capability["match"]):
                continue
            for component in recipe.components:
                match = _component_match(component, capability)
                if match is None:
                    continue
                evidence, confidence = match
                matches.append(
                    (
                        component,
                        str(capability["capability_id"]),
                        evidence,
                        confidence,
                    )
                )
        matches.sort(
            key=lambda value: (
                -value[3],
                str(value[0].substance_id or ""),
                str(value[0].raw_identifier or ""),
                value[1],
            )
        )
        results.append(
            FragmentSourceSupport(
                requirement_id=str(requirement.get("requirement_id") or ""),
                fragment_key=str(requirement.get("fragment_key") or ""),
                status="supported" if matches else "unsupported",
                component_substance_ids=tuple(
                    sorted(
                        {
                            str(component.substance_id)
                            for component, _, _, _ in matches
                            if component.substance_id
                        }
                    )
                ),
                component_raw_identifiers=tuple(
                    sorted(
                        {
                            str(component.raw_identifier)
                            for component, _, _, _ in matches
                            if component.raw_identifier
                        }
                    )
                ),
                capability_ids=tuple(
                    sorted({capability_id for _, capability_id, _, _ in matches})
                ),
                evidence=tuple(
                    sorted(
                        {
                            f"{evidence}:{component.raw_identifier}"
                            for component, _, evidence, _ in matches
                        }
                    )
                ),
                confidence=max(
                    (confidence for _, _, _, confidence in matches),
                    default=0.0,
                ),
                definition_version=FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION,
            )
        )
    return tuple(results)


def matching_fragment_source_capabilities(
    requirement: Any,
) -> Tuple[dict[str, Any], ...]:
    """Return deterministic curated capabilities matching one source gap."""
    requirement_mapping = _mapping(requirement)
    return tuple(
        dict(capability)
        for capability in load_fragment_source_capabilities()["capabilities"]
        if _requirement_matches(requirement_mapping, capability["match"])
    )


def fragment_source_support_is_complete(
    requirements: Iterable[Any],
    supports: Iterable[FragmentSourceSupport | Mapping[str, Any]],
) -> bool:
    """Return whether every requirement has explicit condition support."""
    requirement_ids = {
        str(_mapping(requirement).get("requirement_id") or "")
        for requirement in requirements
    }
    supported_ids = {
        str(_mapping(support).get("requirement_id") or "")
        for support in supports
        if str(_mapping(support).get("status") or "") == "supported"
    }
    return bool(requirement_ids) and requirement_ids <= supported_ids


__all__ = [
    "assess_fragment_source_support",
    "FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION",
    "fragment_source_support_is_complete",
    "load_fragment_source_capabilities",
    "matching_fragment_source_capabilities",
]
