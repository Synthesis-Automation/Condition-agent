"""Apply deposited molecular evidence to reported condition claims."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from condition_registry.models import ConditionComponentInput
from condition_registry.paired_identifiers import paired_identifier_evidence
from reactive_taxonomy.chemistry.source_mapping import (
    SourceMappingEvidence, assess_source_mapping, source_mapping_conflicts, source_material_graphs,
)

from .input_schema import RawReactionRecord


def structural_search_exclusions(record: Mapping[str, Any]) -> tuple[str, ...]:
    """Screen deposited graph pairs before analogue retrieval or operator learning.

    Pass the hydrated canonical source field. Historical CbyD rows lacking the
    source assessment require review. Other corpora without paired source
    representations retain their established structural admission policy.
    """
    source_value = record.get("source")
    source = source_value if isinstance(source_value, Mapping) else {}
    assessment = source.get("source_transfer_assessment")
    if assessment is None and "cbyd" not in (record.get("source_dataset"), source.get("source_dataset")):
        return ()
    mapping = assessment.get("mapping") if isinstance(assessment, Mapping) else None
    issues = mapping.get("issues") if isinstance(mapping, Mapping) else None
    if not isinstance(issues, (tuple, list)) or not all(isinstance(issue, str) for issue in issues):
        issues = ("source_mapping_unavailable_or_unreadable",)
    admitted = record.get("admission_reasons")
    flattened = admitted if isinstance(admitted, (list, tuple)) else ()
    conflicts = source_mapping_conflicts((*issues, *flattened))
    if conflicts:
        return conflicts
    return () if record.get("reaction_signature") else ("source_reaction_signature_unavailable",)


@dataclass(frozen=True)
class SourceClaimEvidence:
    """The decision for one original material reference."""

    component_key: str
    participation: str
    identifier_status: str


@dataclass(frozen=True)
class SourceTransferAssessment:
    """Traceable reasons for retaining a source observation outside transfer."""

    mapping: SourceMappingEvidence
    mapping_origin: str
    claims: tuple[SourceClaimEvidence, ...]
    chemistry_reasons: tuple[str, ...]
    condition_reasons: tuple[str, ...]
    unlisted_agent_graphs: tuple[str, ...]
    definition_version: str = "source_condition_transfer.v1"


def assess_source_claims(
    record: RawReactionRecord, evidence_quality: str,
) -> tuple[tuple[ConditionComponentInput, ...], SourceTransferAssessment | None]:
    """Exclude only unambiguous whole contributors, keeping the original inputs intact."""
    original = record.raw_fields.get("source_observation") or {}
    labels = (original.get("reaction") or {}).get("source_labels") or {}
    if "source_unmapped_reaction_smiles" not in labels:
        return record.condition_component_inputs, None
    mapping = assess_source_mapping(
        str(labels.get("source_unmapped_reaction_smiles") or ""),
        str(labels.get("source_reaction_smiles") or ""),
    )
    origin = str(labels.get("mapping_origin") or "")
    chemistry = set(mapping.issues)
    if origin and evidence_quality not in {
        "fragmented_scaffold_correspondence", "global_atom_correspondence",
        "unique_scaffold_correspondence",
    }:
        chemistry.add("source_model_mapping_requires_review")
    conditions = original.get("conditions") or {}
    review = {
        warning for warning in conditions.get("warnings") or ()
        if warning.startswith("unpaired_condition_evidence:")
        or warning in {"raw_only_condition_evidence", "unresolved_raw_condition_evidence",
                       "unresolved_source_modifiers", "unresolved_temperature_scope"}
    }
    for group in conditions.get("component_groups") or ():
        if (group.get("provenance") or {}).get("requires_review"):
            review.add(f"unpaired_condition_evidence:{group.get('group_key') or 'unknown'}")
    claims = []
    retained = []
    unlisted_agents = Counter(mapping.agent_graphs)
    contributors = set(mapping.contributor_graphs)
    available_spectators = Counter(mapping.noncontributor_graphs + mapping.agent_graphs)
    for item in record.condition_component_inputs:
        if not item.provenance.get("requires_participation_check"):
            retained.append(item)
            continue
        pair = paired_identifier_evidence(item.provenance)
        status = pair.status if pair else "unpaired"
        structures = {
            ident["value"] for ident in item.provenance.get("source_identifiers") or ()
            if ident["identifier_type"] == "smiles"
        }
        graphs = source_material_graphs(next(iter(structures))) if len(structures) == 1 else ()
        if status != "corroborated":
            review.add(f"source_material_identity_{status}")
        overlap = contributors.intersection(graphs)
        if mapping.issues or not graphs:
            participation = "unverified"
            review.add("source_material_participation_unverified")
        elif len(graphs) == 1 and overlap and not overlap.intersection(available_spectators):
            participation = "mapped_atom_contributor"
        elif overlap:
            participation = "ambiguous_contributor_overlap"
            review.add("source_material_contributor_overlap")
        elif Counter(graphs) - available_spectators:
            participation = "unverified"
            review.add("source_material_not_accounted_in_inputs")
        else:
            participation = "no_mapped_product_atoms"
            available_spectators.subtract(graphs)
            unlisted_agents -= Counter(graphs)
        claims.append(SourceClaimEvidence(
            str(item.provenance.get("component_key") or ""), participation, status,
        ))
        if participation != "mapped_atom_contributor":
            retained.append(item)
    unmatched = tuple(sorted(unlisted_agents.elements()))
    if unmatched:
        review.add("source_middle_agents_unaccounted")
    return tuple(retained), SourceTransferAssessment(
        mapping, origin, tuple(claims), tuple(sorted(chemistry)), tuple(sorted(review)), unmatched,
    )
