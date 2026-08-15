"""Evidence-faceted labels for actions observed inside synthesis routes."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Mapping, Optional

from rdkit import Chem
from reactive_taxonomy import featurize_reaction

from .chemistry import canonical_smiles, contributing_reactants, digest, split_reaction_smiles
from .generic_compiler import (
    build_generic_reaction_identity,
    compile_generic_templates,
    generic_rejection_stage,
)
from .mapping import materialize_atom_mapping
from .strategy_identity import build_strategy_id


OBSERVED_ROUTE_ACTION_LABEL_SCHEMA_VERSION = "1.0"
_DEFINITION_PATH = (
    Path(__file__).resolve().parent
    / "definitions"
    / "observed_route_action_label.v1.json"
)


@dataclass(frozen=True)
class ObservedRouteActionLabelPolicy:
    """Declarative evidence threshold for route-learning labels."""

    definition_id: str
    required_completeness_status: str
    allowed_core_statuses: tuple[str, ...]
    allowed_review_reasons: tuple[str, ...]
    required_passed_checks: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.definition_id.startswith("observed_route_action_label.v1@"):
            raise ValueError("unexpected observed route-action definition ID")
        if self.required_completeness_status != "verified":
            raise ValueError("route-action labels require verified completeness")
        if not self.allowed_core_statuses or "blocked" in self.allowed_core_statuses:
            raise ValueError("blocked reaction cores cannot produce route labels")
        if not self.required_passed_checks:
            raise ValueError("route-action policy requires structural checks")


@lru_cache(maxsize=1)
def load_observed_route_action_label_policy() -> ObservedRouteActionLabelPolicy:
    """Load and validate the versioned observed-action evidence policy."""

    value = json.loads(_DEFINITION_PATH.read_text(encoding="utf-8"))
    return ObservedRouteActionLabelPolicy(
        definition_id=str(value["definition_id"]),
        required_completeness_status=str(value["required_completeness_status"]),
        allowed_core_statuses=tuple(
            str(item) for item in value.get("allowed_core_statuses") or ()
        ),
        allowed_review_reasons=tuple(
            str(item) for item in value.get("allowed_review_reasons") or ()
        ),
        required_passed_checks=tuple(
            str(item) for item in value.get("required_passed_checks") or ()
        ),
    )


@dataclass(frozen=True)
class ObservedRouteActionLabel:
    """Independent evidence facets for one action reported in a route."""

    normalized_reaction_smiles: Optional[str]
    target_smiles: Optional[str]
    expected_precursor_smiles: Optional[str]
    disconnection_site_key: Optional[str]
    retained_operator_id: Optional[str]
    retained_operator_signature: Optional[str]
    synthon_signature: Optional[str]
    strategy_id: Optional[str]
    named_annotation: Optional[str]
    evidence_quality: str
    completeness_status: str
    core_quality_status: str
    core_review_reasons: tuple[str, ...]
    core_blocking_reasons: tuple[str, ...]
    retained_edit_count: int
    departing_edit_count: int
    departing_edit_descriptors: tuple[str, ...]
    product_site_verified: bool
    retained_edits_verified: bool
    synthon_partition_verified: bool
    exact_precursors_verified: bool
    strategy_verified: bool
    realization_verified: bool
    operator_roundtrip_verified: bool
    search_eligible: bool
    operator_admission_stage: str
    operator_admission_reason: Optional[str]
    limitations: tuple[str, ...]
    policy_definition_id: str
    schema_version: str = OBSERVED_ROUTE_ACTION_LABEL_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != OBSERVED_ROUTE_ACTION_LABEL_SCHEMA_VERSION:
            raise ValueError("unsupported observed route-action label schema")
        if self.product_site_verified and not self.disconnection_site_key:
            raise ValueError("verified product site requires SITE1 identity")
        if self.retained_edits_verified and not all(
            (self.retained_operator_id, self.retained_operator_signature)
        ):
            raise ValueError("verified retained edits require OP1 identity")
        if self.synthon_partition_verified and not self.synthon_signature:
            raise ValueError("verified synthon partition requires SYN1 identity")
        if self.exact_precursors_verified and not self.expected_precursor_smiles:
            raise ValueError("verified precursors require canonical identity")
        if self.strategy_verified and not all(
            (
                self.product_site_verified,
                self.retained_edits_verified,
                self.synthon_partition_verified,
                self.strategy_id,
            )
        ):
            raise ValueError("verified strategy requires SITE1, OP1, and SYN1")
        if self.realization_verified and not (
            self.strategy_verified and self.exact_precursors_verified
        ):
            raise ValueError("verified realization requires strategy and precursors")
        if self.operator_roundtrip_verified and not self.retained_edits_verified:
            raise ValueError("round-tripped operator requires retained edit identity")
        if self.search_eligible != (
            self.product_site_verified or self.retained_edits_verified
        ):
            raise ValueError("search eligibility contradicts learnable action facets")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible observed-action label."""

        value = asdict(self)
        for field in (
            "core_review_reasons",
            "core_blocking_reasons",
            "departing_edit_descriptors",
            "limitations",
        ):
            value[field] = list(value[field])
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ObservedRouteActionLabel":
        """Reconstruct an observed-action label from JSON data."""

        fields = dict(value)
        for field in (
            "core_review_reasons",
            "core_blocking_reasons",
            "departing_edit_descriptors",
            "limitations",
        ):
            fields[field] = tuple(str(item) for item in fields.get(field) or ())
        return cls(**fields)


def normalize_observed_reaction(reaction_smiles: str) -> Optional[str]:
    """Remove the condition middle field while preserving reaction participants."""

    if reaction_smiles.count(">>") == 1:
        reactants, products = reaction_smiles.split(">>")
    elif reaction_smiles.count(">") == 2:
        reactants, _, products = reaction_smiles.split(">")
    else:
        return None
    if not reactants or not products:
        return None
    return f"{reactants}>>{products}"


def _edit_partition(
    observation: Mapping[str, Any],
    product_maps: set[int],
) -> tuple[int, tuple[str, ...]]:
    retained_count = 0
    departing = []
    periodic_table = Chem.GetPeriodicTable()
    for raw_edit in observation.get("edits") or ():
        edit = dict(raw_edit)
        atom_1 = dict(edit.get("atom_1") or {})
        atom_2_raw = edit.get("atom_2")
        atom_2 = dict(atom_2_raw or {})
        maps = [int(atom_1.get("atom_map_number") or 0)]
        if atom_2_raw is not None:
            maps.append(int(atom_2.get("atom_map_number") or 0))
        retained = bool(maps) and all(value in product_maps for value in maps)
        if retained:
            retained_count += 1
            continue
        elements = [str(atom_1.get("element") or "?")]
        if atom_2_raw is not None:
            elements.append(str(atom_2.get("element") or "?"))
        elements.sort(
            key=lambda symbol: (
                periodic_table.GetAtomicNumber(symbol)
                if symbol != "?"
                else 10_000,
                symbol,
            )
        )
        departing.append(
            ":".join(
                (
                    str(edit.get("edit_type") or "unknown"),
                    "-".join(elements),
                    str(edit.get("old_order") or "NONE"),
                    str(edit.get("new_order") or "NONE"),
                )
            )
        )
    return retained_count, tuple(sorted(departing))


def _empty_label(
    *,
    normalized_reaction: Optional[str],
    route_product_smiles: Optional[str],
    limitation: str,
    policy: ObservedRouteActionLabelPolicy,
) -> ObservedRouteActionLabel:
    return ObservedRouteActionLabel(
        normalized_reaction_smiles=normalized_reaction,
        target_smiles=canonical_smiles(route_product_smiles or ""),
        expected_precursor_smiles=None,
        disconnection_site_key=None,
        retained_operator_id=None,
        retained_operator_signature=None,
        synthon_signature=None,
        strategy_id=None,
        named_annotation=None,
        evidence_quality="unavailable",
        completeness_status="unavailable",
        core_quality_status="unavailable",
        core_review_reasons=(),
        core_blocking_reasons=(),
        retained_edit_count=0,
        departing_edit_count=0,
        departing_edit_descriptors=(),
        product_site_verified=False,
        retained_edits_verified=False,
        synthon_partition_verified=False,
        exact_precursors_verified=False,
        strategy_verified=False,
        realization_verified=False,
        operator_roundtrip_verified=False,
        search_eligible=False,
        operator_admission_stage="not_attempted",
        operator_admission_reason=limitation,
        limitations=(limitation,),
        policy_definition_id=policy.definition_id,
    )


def build_observed_route_action_label(
    reaction_smiles: str,
    *,
    route_product_smiles: Optional[str] = None,
    reaction_id: str = "",
    reference_id: str = "",
    policy: Optional[ObservedRouteActionLabelPolicy] = None,
) -> ObservedRouteActionLabel:
    """Build task-specific supervision without relaxing operator admission."""

    resolved = policy or load_observed_route_action_label_policy()
    normalized = normalize_observed_reaction(reaction_smiles)
    if normalized is None:
        return _empty_label(
            normalized_reaction=None,
            route_product_smiles=route_product_smiles,
            limitation="invalid_reaction_format",
            policy=resolved,
        )
    materialized = materialize_atom_mapping(normalized)
    if materialized is None:
        return _empty_label(
            normalized_reaction=normalized,
            route_product_smiles=route_product_smiles,
            limitation="atom_mapping_unavailable",
            policy=resolved,
        )
    split = split_reaction_smiles(materialized.reaction_smiles)
    if split is None:
        return _empty_label(
            normalized_reaction=normalized,
            route_product_smiles=route_product_smiles,
            limitation="materialized_reaction_invalid",
            policy=resolved,
        )
    reactant_smiles, product_smiles = split
    product = Chem.MolFromSmiles(product_smiles)
    analysis = featurize_reaction(materialized.reaction_smiles)
    if product is None or not analysis.valid or analysis.reaction_core is None:
        return _empty_label(
            normalized_reaction=normalized,
            route_product_smiles=route_product_smiles,
            limitation="reaction_core_unavailable",
            policy=resolved,
        )
    analysis_value = analysis.to_dict()
    observation = dict(analysis_value.get("observation") or {})
    completeness = dict(observation.get("completeness") or {})
    quality = analysis.reaction_core.quality
    product_maps = {
        int(atom.GetAtomMapNum())
        for atom in product.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    retained_count, departing_descriptors = _edit_partition(
        observation, product_maps
    )
    identity = build_generic_reaction_identity(
        materialized.reaction_smiles,
        analysis,
    )
    expected_precursors = contributing_reactants(reactant_smiles, product_smiles)
    canonical_product = canonical_smiles(product_smiles)
    route_product = (
        canonical_smiles(route_product_smiles)
        if route_product_smiles is not None
        else canonical_product
    )
    route_product_matches = bool(
        canonical_product and route_product and canonical_product == route_product
    )
    allowed_quality = quality.status in resolved.allowed_core_statuses
    allowed_review = (
        quality.status != "review"
        or set(quality.review_reasons).issubset(resolved.allowed_review_reasons)
    )
    review_has_departing_evidence = (
        quality.status != "review" or bool(departing_descriptors)
    )
    passed_checks = set(quality.passed_checks)
    required_checks = set(resolved.required_passed_checks)
    completeness_verified = (
        str(completeness.get("status") or "")
        == resolved.required_completeness_status
    )
    structural_evidence = bool(
        route_product_matches
        and completeness_verified
        and allowed_quality
        and allowed_review
        and review_has_departing_evidence
        and required_checks.issubset(passed_checks)
        and not quality.blocking_reasons
    )
    operator_signature = identity.operator_signature if identity else ""
    operator_id = digest("OP1", operator_signature) if operator_signature else ""
    site_key = identity.disconnection_site_key if identity else ""
    synthon = identity.synthon_signature if identity else ""
    retained_verified = bool(
        structural_evidence and operator_signature and retained_count > 0
    )
    site_verified = bool(retained_verified and site_key)
    synthon_verified = bool(retained_verified and synthon)
    exact_verified = bool(structural_evidence and expected_precursors)
    strategy_verified = bool(site_verified and synthon_verified)
    strategy_id = (
        build_strategy_id(operator_id, site_key, synthon)
        if strategy_verified
        else ""
    )
    if quality.status == "pass":
        compiled = compile_generic_templates(
            {
                "reaction_id": reaction_id,
                "reference_id": reference_id,
                "reaction_smiles": normalized,
            },
            engine="reaction_core",
            levels=("L0", "L1", "L2"),
            admission_mode="data_driven",
        )
        roundtrip = bool(compiled.templates)
        admission_reason = compiled.rejection_reason
        admission_stage = (
            "accepted" if roundtrip else compiled.rejection_stage or "unknown"
        )
    else:
        roundtrip = False
        admission_reason = "materialized_core_not_verified"
        admission_stage = generic_rejection_stage(admission_reason)
    limitations = []
    if not route_product_matches:
        limitations.append("target_identity_mismatch")
    if not completeness_verified:
        limitations.append("product_completeness_not_verified")
    if not allowed_quality:
        limitations.append(f"core_quality_{quality.status}")
    if not allowed_review:
        limitations.append("unsupported_core_review_reason")
    if not review_has_departing_evidence:
        limitations.append("review_missing_departing_edit_evidence")
    if not required_checks.issubset(passed_checks):
        limitations.append("required_structural_checks_missing")
    if quality.blocking_reasons:
        limitations.append("core_has_blocking_reasons")
    if not operator_signature:
        limitations.append("retained_edit_identity_unavailable")
    if not site_key:
        limitations.append("product_site_identity_unavailable")
    if not synthon:
        limitations.append("synthon_identity_unavailable")
    if not expected_precursors:
        limitations.append("precursor_identity_unavailable")
    if not roundtrip:
        limitations.append(str(admission_reason or "operator_not_roundtripped"))
    return ObservedRouteActionLabel(
        normalized_reaction_smiles=normalized,
        target_smiles=canonical_product,
        expected_precursor_smiles=expected_precursors,
        disconnection_site_key=site_key or None,
        retained_operator_id=operator_id or None,
        retained_operator_signature=operator_signature or None,
        synthon_signature=synthon or None,
        strategy_id=strategy_id or None,
        named_annotation=(identity.named_annotation if identity else None),
        evidence_quality=str(analysis.evidence_quality),
        completeness_status=str(completeness.get("status") or "unavailable"),
        core_quality_status=str(quality.status),
        core_review_reasons=tuple(quality.review_reasons),
        core_blocking_reasons=tuple(quality.blocking_reasons),
        retained_edit_count=retained_count,
        departing_edit_count=len(departing_descriptors),
        departing_edit_descriptors=departing_descriptors,
        product_site_verified=site_verified,
        retained_edits_verified=retained_verified,
        synthon_partition_verified=synthon_verified,
        exact_precursors_verified=exact_verified,
        strategy_verified=strategy_verified,
        realization_verified=bool(strategy_verified and exact_verified),
        operator_roundtrip_verified=roundtrip,
        search_eligible=bool(site_verified or retained_verified),
        operator_admission_stage=admission_stage,
        operator_admission_reason=(None if roundtrip else str(admission_reason)),
        limitations=tuple(dict.fromkeys(limitations)),
        policy_definition_id=resolved.definition_id,
    )


__all__ = [
    "OBSERVED_ROUTE_ACTION_LABEL_SCHEMA_VERSION",
    "ObservedRouteActionLabel",
    "ObservedRouteActionLabelPolicy",
    "build_observed_route_action_label",
    "load_observed_route_action_label_policy",
    "normalize_observed_reaction",
]
