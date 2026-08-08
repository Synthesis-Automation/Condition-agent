"""Propose and resolve condition sources for incomplete reaction queries."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from typing import Any

from condition_registry import resolve_identifier, resolve_substance_id
from reactive_taxonomy import featurize_reaction

from .fragment_source_support import matching_fragment_source_capabilities
from .models import (
    ReactionCompletionOption,
    ReactionCompletionProposal,
    ReactionCompletionRequirement,
    ReactionCompletionSelection,
)


def _proposal_id(reaction_smiles: str, requirements: tuple[Any, ...]) -> str:
    payload = {
        "query_reaction_smiles": reaction_smiles,
        "requirements": [
            {
                "requirement_id": value.requirement_id,
                "fragment_key": value.fragment_key,
                "rooted_fragment_smiles": value.rooted_fragment_smiles,
            }
            for value in requirements
        ],
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return f"RCP1:{hashlib.sha256(encoded).hexdigest()[:24]}"


def _option_id(
    requirement_id: str,
    option_kind: str,
    identity: str,
) -> str:
    encoded = f"{requirement_id}|{option_kind}|{identity}".encode("utf-8")
    return f"RCO1:{hashlib.sha256(encoded).hexdigest()[:20]}"


def build_reaction_completion_proposal(analysis: Any) -> ReactionCompletionProposal:
    """Build source choices from observed product-only fragment requirements."""
    descriptor = analysis.fallback_descriptor
    source_requirements = (
        tuple(descriptor.source_requirements) if descriptor is not None else ()
    )
    requirements = []
    for source_requirement in source_requirements:
        options = []
        capabilities = matching_fragment_source_capabilities(source_requirement)
        for capability in capabilities:
            capability_id = str(capability["capability_id"])
            display_name = str(capability["display_name"])
            options.append(
                ReactionCompletionOption(
                    option_id=_option_id(
                        source_requirement.requirement_id,
                        "compatible_source_class",
                        capability_id,
                    ),
                    option_kind="compatible_source_class",
                    display_name=f"Any compatible {display_name}",
                    capability_id=capability_id,
                )
            )
            substances = []
            for substance_id in capability.get("substance_ids") or ():
                resolved = resolve_substance_id(str(substance_id))
                if resolved.status != "resolved" or resolved.substance is None:
                    continue
                substances.append(resolved.substance)
            for substance in sorted(
                substances,
                key=lambda value: (value.canonical_name.lower(), value.substance_id),
            ):
                options.append(
                    ReactionCompletionOption(
                        option_id=_option_id(
                            source_requirement.requirement_id,
                            "registered_substance",
                            substance.substance_id,
                        ),
                        option_kind="registered_substance",
                        display_name=substance.canonical_name,
                        capability_id=capability_id,
                        substance_id=substance.substance_id,
                        canonical_name=substance.canonical_name,
                    )
                )
        options.append(
            ReactionCompletionOption(
                option_id=_option_id(
                    source_requirement.requirement_id,
                    "unresolved",
                    "unresolved",
                ),
                option_kind="unresolved",
                display_name="Leave source unresolved",
            )
        )
        requirements.append(
            ReactionCompletionRequirement(
                requirement_id=source_requirement.requirement_id,
                fragment_key=source_requirement.fragment_key,
                canonical_fragment_smiles=(
                    source_requirement.canonical_fragment_smiles
                ),
                rooted_fragment_smiles=source_requirement.rooted_fragment_smiles,
                element_counts=dict(source_requirement.element_counts),
                attachment_element=source_requirement.attachment_element,
                options=tuple(options),
            )
        )

    requirement_tuple = tuple(requirements)
    proposal_id = _proposal_id(
        str(analysis.input_reaction_smiles),
        source_requirements,
    )
    if not requirement_tuple:
        return ReactionCompletionProposal(
            query_reaction_smiles=str(analysis.input_reaction_smiles),
            proposal_id=proposal_id,
            status="not_required",
        )
    has_curated_options = all(
        any(option.option_kind != "unresolved" for option in value.options)
        for value in requirement_tuple
    )
    return ReactionCompletionProposal(
        query_reaction_smiles=str(analysis.input_reaction_smiles),
        proposal_id=proposal_id,
        status=(
            "confirmation_recommended"
            if has_curated_options
            else "no_curated_source_options"
        ),
        requirements=requirement_tuple,
        warnings=(
            "QUERY_PRODUCT_SOURCE_COMPONENT_MISSING",
            "SYSTEM_PROPOSED_SOURCES_REQUIRE_USER_CONFIRMATION",
        ),
    )


def propose_reaction_completion(
    reaction_smiles: str,
) -> ReactionCompletionProposal:
    """Featurize a reaction and propose minimal missing condition sources."""
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid:
        raise ValueError(analysis.error or "INVALID_REACTION")
    return build_reaction_completion_proposal(analysis)


def build_completion_selection(
    proposal: ReactionCompletionProposal,
    requirement_id: str,
    *,
    option_id: str | None = None,
    custom_identifier: str | None = None,
) -> ReactionCompletionSelection:
    """Resolve a confirmed proposal option or user-edited source identifier."""
    requirement = next(
        (
            value
            for value in proposal.requirements
            if value.requirement_id == requirement_id
        ),
        None,
    )
    if requirement is None:
        raise ValueError(f"unknown completion requirement: {requirement_id}")
    if custom_identifier is not None:
        raw_identifier = custom_identifier.strip()
        if not raw_identifier:
            raise ValueError("custom completion identifier must not be empty")
        resolved = resolve_identifier(raw_identifier)
        candidate_by_substance = {
            option.substance_id: option
            for option in requirement.options
            if option.substance_id
        }
        matched = (
            candidate_by_substance.get(resolved.substance.substance_id)
            if resolved.status == "resolved" and resolved.substance is not None
            else None
        )
        if matched is not None:
            return ReactionCompletionSelection(
                proposal_id=proposal.proposal_id,
                requirement_id=requirement_id,
                selection_kind="registered_substance",
                provenance="user_edited",
                display_name=matched.display_name,
                capability_id=matched.capability_id,
                substance_id=matched.substance_id,
                raw_identifier=raw_identifier,
                resolved=True,
            )
        return ReactionCompletionSelection(
            proposal_id=proposal.proposal_id,
            requirement_id=requirement_id,
            selection_kind="custom_identifier",
            provenance="user_edited",
            display_name=raw_identifier,
            raw_identifier=raw_identifier,
            resolved=False,
        )

    option = next(
        (value for value in requirement.options if value.option_id == option_id),
        None,
    )
    if option is None:
        raise ValueError(f"unknown completion option: {option_id}")
    return ReactionCompletionSelection(
        proposal_id=proposal.proposal_id,
        requirement_id=requirement_id,
        selection_kind=option.option_kind,
        provenance="user_confirmed",
        display_name=option.display_name,
        capability_id=option.capability_id,
        substance_id=option.substance_id,
        resolved=option.option_kind != "unresolved",
    )


def validate_completion_selections(
    proposal: ReactionCompletionProposal,
    selections: tuple[ReactionCompletionSelection, ...],
) -> None:
    """Reject stale, duplicate, or structurally unrelated selections."""
    requirement_ids = {value.requirement_id for value in proposal.requirements}
    selected_ids = [value.requirement_id for value in selections]
    if len(selected_ids) != len(set(selected_ids)):
        raise ValueError("duplicate completion requirement selection")
    for selection in selections:
        if selection.proposal_id != proposal.proposal_id:
            raise ValueError("completion selection belongs to a different proposal")
        if selection.requirement_id not in requirement_ids:
            raise ValueError(
                f"unknown completion requirement: {selection.requirement_id}"
            )
        requirement = next(
            value
            for value in proposal.requirements
            if value.requirement_id == selection.requirement_id
        )
        valid_pairs = {
            (option.capability_id, option.substance_id)
            for option in requirement.options
        }
        if selection.selection_kind == "compatible_source_class" and (
            selection.capability_id,
            None,
        ) not in valid_pairs:
            raise ValueError("completion capability is not valid for requirement")
        if selection.selection_kind == "registered_substance" and (
            selection.capability_id,
            selection.substance_id,
        ) not in valid_pairs:
            raise ValueError("completion substance is not valid for requirement")


def build_completed_reaction_smiles(
    reaction_smiles: str,
    selections: tuple[ReactionCompletionSelection, ...],
) -> tuple[str | None, tuple[str, ...]]:
    """Build a hypothetical reactant-completed query from exact source choices."""
    source_smiles = []
    warnings = []
    for selection in selections:
        if selection.selection_kind != "registered_substance":
            continue
        resolved = resolve_substance_id(str(selection.substance_id or ""))
        smiles = (
            str(resolved.substance.smiles or "").strip()
            if resolved.status == "resolved" and resolved.substance is not None
            else ""
        )
        if not smiles:
            warnings.append(
                "CONFIRMED_SOURCE_HAS_NO_REGISTERED_STRUCTURE:"
                f"{selection.substance_id}"
            )
            continue
        source_smiles.append(smiles)
    if not source_smiles:
        return None, tuple(warnings)

    if ">>" in reaction_smiles:
        reactants, products = reaction_smiles.split(">>", 1)
        existing_agents = None
    else:
        parts = reaction_smiles.split(">")
        if len(parts) != 3:
            raise ValueError("INVALID_REACTION_FORMAT")
        reactants, existing_agents, products = parts
    reactant_components = tuple(
        component
        for value in (reactants, *source_smiles)
        for component in value.strip().split(".")
        if component
    )
    completed_reactants = ".".join(dict.fromkeys(reactant_components))
    if existing_agents is None:
        return f"{completed_reactants}>>{products}", tuple(warnings)
    return (
        f"{completed_reactants}>{existing_agents.strip()}>{products}",
        tuple(warnings),
    )


def completion_payload(value: ReactionCompletionProposal) -> dict[str, Any]:
    """Serialize a completion proposal for recommendation results."""
    return asdict(value)


__all__ = [
    "build_completion_selection",
    "build_completed_reaction_smiles",
    "build_reaction_completion_proposal",
    "completion_payload",
    "propose_reaction_completion",
    "validate_completion_selections",
]
