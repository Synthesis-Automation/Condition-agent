"""Deterministic admission assessment for externally proposed route steps.

External labels, citations, and claimed transformations are provenance only.
All trusted chemistry fields are reconstructed from the supplied molecular
graphs through the same mapping, signature, operator, precedent, and
compatibility paths used by the canonical retrosynthesis system.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Callable, Literal, Mapping

from forward_synthesis import ForwardOperatorLibrary, RouteStepForwardAssessment
from forward_synthesis import assess_proposed_step
from reactive_taxonomy import (
    AtomMappingProviderMetadata,
    ReactionSignature,
    featurize_reaction,
    validate_external_atom_mapping,
)

from .chemistry import canonical_smiles, digest
from .condition_ranking import RetrosynthesisConditionEvidence
from .generic_compiler import GenericReactionIdentity, analyze_generic_reaction
from .generic_models import GenericTemplateLibrary
from .mapping import materialize_atom_mapping
from .precursor_compatibility import (
    PrecursorCompatibilityResult,
    assess_precursor_compatibility,
)
from .reaction_compatibility import (
    ReactionCompatibilityResult,
    assess_candidate_reaction_compatibility,
)
from .selectivity_poc import (
    FunctionalGroupCompetitionWarning,
    detect_functional_group_competition,
)
from .step_precedents import StepPrecedentMatch, lookup_reaction_precedents


EXTERNAL_PROPOSAL_ASSESSMENT_SCHEMA_VERSION = "1.0"
EXTERNAL_PROPOSAL_ADMISSION_POLICY_PATH = (
    Path(__file__).with_name("definitions")
    / "external_proposal_admission_policy.v1.json"
)

GateStatus = Literal[
    "pass", "warn", "fail", "unresolved", "not_run", "out_of_scope"
]
ExternalProposalStatus = Literal[
    "invalid",
    "ambiguous",
    "structurally_verified",
    "operator_supported",
    "precedent_supported",
    "forward_reproduced",
    "verified_with_cautions",
    "rejected_by_compatibility",
]
OperatorMatchLevel = Literal[
    "exact_operator_signature",
    "same_edit_archetype_only",
]

_SOURCE_KINDS = frozenset(
    {
        "chemist",
        "literature",
        "web",
        "llm",
        "external_ssr",
        "imported_route",
    }
)
_EVIDENCE_TIERS = (
    "unresolved",
    "structurally_verified",
    "operator_supported",
    "precedent_supported",
    "forward_reproduced",
)


@dataclass(frozen=True)
class ExternalProposalSource:
    """Untrusted provenance for an externally proposed transformation."""

    source_id: str
    source_kind: str
    reference_id: str | None = None
    source_uri: str | None = None
    model_id: str | None = None
    retrieved_at: str | None = None
    quoted_claim: str | None = None

    def __post_init__(self) -> None:
        if not self.source_id:
            raise ValueError("external proposal source ID is required")
        if self.source_kind not in _SOURCE_KINDS:
            raise ValueError(f"unsupported external source kind: {self.source_kind}")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ExternalProposalSource":
        return cls(
            source_id=str(value["source_id"]),
            source_kind=str(value["source_kind"]),
            reference_id=_optional_string(value.get("reference_id")),
            source_uri=_optional_string(value.get("source_uri")),
            model_id=_optional_string(value.get("model_id")),
            retrieved_at=_optional_string(value.get("retrieved_at")),
            quoted_claim=_optional_string(value.get("quoted_claim")),
        )


@dataclass(frozen=True)
class ExternalRetrosynthesisProposal:
    """One concrete precursor-to-target proposal without trusted identities."""

    target_smiles: str
    precursor_smiles: str
    mapped_reaction_smiles: str | None = None
    proposed_conditions: Mapping[str, Any] | None = None
    sources: tuple[ExternalProposalSource, ...] = ()
    external_proposal_id: str | None = None
    claimed_transformation: str | None = None
    schema_version: str = "1.0"

    def to_dict(self) -> dict[str, Any]:
        return {
            **asdict(self),
            "sources": [item.to_dict() for item in self.sources],
            "proposed_conditions": (
                dict(self.proposed_conditions)
                if self.proposed_conditions is not None
                else None
            ),
        }

    @classmethod
    def from_dict(
        cls, value: Mapping[str, Any]
    ) -> "ExternalRetrosynthesisProposal":
        if str(value.get("schema_version") or "1.0") != "1.0":
            raise ValueError("unsupported external proposal schema")
        forbidden = {
            "operator_id",
            "operator_signature",
            "strategy_id",
            "template_id",
            "forward_validation_status",
            "compatibility_disposition",
            "admission_status",
            "confidence",
        }.intersection(value)
        if forbidden:
            raise ValueError(
                "external proposals cannot supply trusted fields: "
                + ", ".join(sorted(forbidden))
            )
        conditions = value.get("proposed_conditions")
        if conditions is not None and not isinstance(conditions, Mapping):
            raise ValueError("proposed conditions must be an object")
        return cls(
            target_smiles=str(value.get("target_smiles") or ""),
            precursor_smiles=str(value.get("precursor_smiles") or ""),
            mapped_reaction_smiles=_optional_string(
                value.get("mapped_reaction_smiles")
            ),
            proposed_conditions=(dict(conditions) if conditions is not None else None),
            sources=tuple(
                ExternalProposalSource.from_dict(item)
                for item in value.get("sources") or ()
            ),
            external_proposal_id=_optional_string(
                value.get("external_proposal_id")
            ),
            claimed_transformation=_optional_string(
                value.get("claimed_transformation")
            ),
        )


@dataclass(frozen=True)
class ExternalProposalAssessmentLimits:
    """Bounded work limits loaded from the admission definition."""

    maximum_operator_matches: int = 20
    maximum_precedent_matches: int = 10
    maximum_forward_operators: int = 300
    maximum_forward_products: int = 20
    maximum_route_steps: int = 40

    def __post_init__(self) -> None:
        if min(
            self.maximum_operator_matches,
            self.maximum_precedent_matches,
            self.maximum_forward_operators,
            self.maximum_forward_products,
            self.maximum_route_steps,
        ) < 1:
            raise ValueError("external proposal assessment limits must be positive")


@dataclass(frozen=True)
class ExternalProposalAdmissionPolicy:
    """Validated review-only policy for external proposal admission."""

    definition_id: str
    schema_version: str
    assessment_schema_version: str
    required_gate_order: tuple[str, ...]
    limits: ExternalProposalAssessmentLimits
    minimum_evidence_tier: str
    ranking_influence: str
    public_actionability: str


@dataclass(frozen=True)
class ExternalProposalGate:
    """One explicit deterministic evidence gate."""

    gate_id: str
    status: GateStatus
    summary: str
    evidence_ids: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    definition_id: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ExternalOperatorMatch:
    """Internally recovered operator support; never copied from the proposal."""

    match_id: str
    match_level: OperatorMatchLevel
    operator_id: str
    operator_signature: str
    template_ids: tuple[str, ...]
    observation_support: int
    independent_reference_support: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ExternalCompatibilityEvidence:
    """Existing deterministic compatibility and selectivity assessments."""

    precursor: PrecursorCompatibilityResult
    reaction: ReactionCompatibilityResult
    functional_group_competition: FunctionalGroupCompetitionWarning | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ExternalRetrosynthesisAssessment:
    """Auditable result for one externally proposed physical step."""

    assessment_id: str
    proposal_id: str
    external_proposal_id: str | None
    canonical_target_smiles: str | None
    canonical_precursor_smiles: str | None
    normalized_mapped_reaction_smiles: str | None
    status: ExternalProposalStatus
    strongest_evidence_tier: str
    admission_eligible: bool
    actionable: bool
    internal_confidence_parity: bool
    gates: tuple[ExternalProposalGate, ...]
    reaction_signature: ReactionSignature | None
    reaction_identity: GenericReactionIdentity | None
    operator_matches: tuple[ExternalOperatorMatch, ...]
    selected_operator_match_id: str | None
    selected_template_id: str | None
    precedent_matches: tuple[StepPrecedentMatch, ...]
    forward_assessment: RouteStepForwardAssessment | None
    compatibility_evidence: ExternalCompatibilityEvidence | None
    condition_evidence: RetrosynthesisConditionEvidence | None
    sources: tuple[ExternalProposalSource, ...]
    warnings: tuple[str, ...]
    schema_versions: Mapping[str, str]
    ranking_influence: str = "none_review_only"
    schema_version: str = EXTERNAL_PROPOSAL_ASSESSMENT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


ConditionEvidenceEvaluator = Callable[[str], RetrosynthesisConditionEvidence]


def _optional_string(value: Any) -> str | None:
    return str(value) if value is not None else None


def validate_external_proposal_admission_policy(
    value: Mapping[str, Any],
) -> None:
    """Validate the declarative gate order, limits, and review-only boundary."""

    if value.get("definition_id") != "external_proposal_admission_policy.v1":
        raise ValueError("unexpected external proposal admission policy ID")
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported external proposal admission policy schema")
    gates = tuple(value.get("required_gate_order") or ())
    expected = (
        "input_structure",
        "reaction_side_consistency",
        "atom_correspondence",
        "reaction_completeness",
        "reaction_signature",
        "operator_support",
        "precedent_support",
        "stereochemistry",
        "precursor_compatibility",
        "reaction_compatibility",
        "functional_group_selectivity",
        "forward_reproduction",
        "forward_competition",
        "condition_support",
    )
    if gates != expected:
        raise ValueError("external proposal gate order is incomplete or unsupported")
    raw_limits = value.get("limits")
    if not isinstance(raw_limits, Mapping):
        raise ValueError("external proposal limits must be an object")
    ExternalProposalAssessmentLimits(
        **{field: int(raw_limits[field]) for field in ExternalProposalAssessmentLimits.__dataclass_fields__}
    )
    admission = value.get("admission")
    if not isinstance(admission, Mapping):
        raise ValueError("external proposal admission policy is missing")
    if admission.get("minimum_evidence_tier") != "precedent_supported":
        raise ValueError("v1 external admission requires precedent support")
    if admission.get("ranking_influence") != "none_review_only":
        raise ValueError("external proposal assessment must remain review-only")
    if admission.get("public_actionability") != "disabled_until_release_gate":
        raise ValueError("external public actionability must remain disabled")


def load_external_proposal_admission_policy() -> ExternalProposalAdmissionPolicy:
    """Load and validate the canonical external proposal policy."""

    value = json.loads(
        EXTERNAL_PROPOSAL_ADMISSION_POLICY_PATH.read_text(encoding="utf-8")
    )
    validate_external_proposal_admission_policy(value)
    raw_limits = value["limits"]
    admission = value["admission"]
    return ExternalProposalAdmissionPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        assessment_schema_version=str(value["assessment_schema_version"]),
        required_gate_order=tuple(value["required_gate_order"]),
        limits=ExternalProposalAssessmentLimits(
            **{
                field: int(raw_limits[field])
                for field in ExternalProposalAssessmentLimits.__dataclass_fields__
            }
        ),
        minimum_evidence_tier=str(admission["minimum_evidence_tier"]),
        ranking_influence=str(admission["ranking_influence"]),
        public_actionability=str(admission["public_actionability"]),
    )


def _gate(
    gate_id: str,
    status: GateStatus,
    summary: str,
    *,
    evidence_ids: tuple[str, ...] = (),
    warnings: tuple[str, ...] = (),
    definition_id: str | None = None,
) -> ExternalProposalGate:
    return ExternalProposalGate(
        gate_id=gate_id,
        status=status,
        summary=summary,
        evidence_ids=evidence_ids,
        warnings=warnings,
        definition_id=definition_id,
    )


def _operator_archetype(operator_signature: str) -> tuple[tuple[str, str, str], ...]:
    try:
        tokens = json.loads(operator_signature)
    except (TypeError, ValueError):
        return ()
    values = []
    for token in tokens if isinstance(tokens, list) else ():
        if not isinstance(token, list) or len(token) != 4:
            return ()
        values.append((str(token[0]), str(token[2]), str(token[3])))
    return tuple(sorted(values))


def _operator_matches(
    identity: GenericReactionIdentity,
    library: GenericTemplateLibrary,
    *,
    limit: int,
) -> tuple[ExternalOperatorMatch, ...]:
    templates_by_operator: dict[str, list[Any]] = {}
    for template in library.templates:
        templates_by_operator.setdefault(template.operator_id, []).append(template)
    operator_rows: dict[str, tuple[str, int, int]] = {
        operator.operator_id: (
            operator.operator_signature,
            operator.observation_support,
            operator.independent_reference_support,
        )
        for operator in library.operators
    }
    for operator_id, templates in templates_by_operator.items():
        if operator_id in operator_rows:
            continue
        signatures = {item.operator_signature for item in templates if item.operator_signature}
        if len(signatures) != 1:
            continue
        operator_rows[operator_id] = (
            next(iter(signatures)),
            max(item.observation_support for item in templates),
            max(item.independent_reference_support for item in templates),
        )
    query_archetype = _operator_archetype(identity.operator_signature)
    matches = []
    for operator_id, (signature, observations, references) in operator_rows.items():
        if signature == identity.operator_signature:
            level: OperatorMatchLevel = "exact_operator_signature"
        elif query_archetype and _operator_archetype(signature) == query_archetype:
            level = "same_edit_archetype_only"
        else:
            continue
        template_ids = tuple(
            sorted(item.template_id for item in templates_by_operator.get(operator_id, ()))
        )
        matches.append(
            ExternalOperatorMatch(
                match_id=digest(
                    "EXTOPM1", identity.operator_signature, operator_id, level
                ),
                match_level=level,
                operator_id=operator_id,
                operator_signature=signature,
                template_ids=template_ids,
                observation_support=int(observations),
                independent_reference_support=int(references),
            )
        )
    return tuple(
        sorted(
            matches,
            key=lambda item: (
                item.match_level != "exact_operator_signature",
                -item.independent_reference_support,
                -item.observation_support,
                item.operator_id,
            ),
        )[:limit]
    )


def _precedents(
    proposal_id: str,
    target_smiles: str,
    precursor_smiles: str,
    exact_matches: tuple[ExternalOperatorMatch, ...],
    library: GenericTemplateLibrary,
    *,
    limit: int,
) -> tuple[tuple[StepPrecedentMatch, ...], str | None]:
    matches: list[StepPrecedentMatch] = []
    selected_template_id: str | None = None
    for operator_match in exact_matches:
        for template_id in operator_match.template_ids:
            lookup = lookup_reaction_precedents(
                step_id=proposal_id,
                template_id=template_id,
                operator_id=operator_match.operator_id,
                product_smiles=target_smiles,
                precursor_smiles=precursor_smiles,
                library=library,
                limit=limit,
            )
            if lookup.matches and selected_template_id is None:
                selected_template_id = template_id
            matches.extend(lookup.matches)
    unique = {item.match_id: item for item in matches}
    ranked = sorted(
        unique.values(),
        key=lambda item: (
            -item.product_similarity,
            -item.precursor_similarity,
            item.reference_id,
            item.reaction_id,
            item.template_id,
        ),
    )
    return tuple(ranked[:limit]), selected_template_id


def _finish_gates(
    policy: ExternalProposalAdmissionPolicy,
    gates: Mapping[str, ExternalProposalGate],
) -> tuple[ExternalProposalGate, ...]:
    return tuple(gates[gate_id] for gate_id in policy.required_gate_order)


def assess_external_retrosynthesis_proposal(
    proposal: ExternalRetrosynthesisProposal,
    operator_library: GenericTemplateLibrary,
    *,
    forward_library: ForwardOperatorLibrary | None = None,
    condition_evaluator: ConditionEvidenceEvaluator | None = None,
    limits: ExternalProposalAssessmentLimits | None = None,
) -> ExternalRetrosynthesisAssessment:
    """Assess one external step without trusting any external chemistry label."""

    policy = load_external_proposal_admission_policy()
    bounded = limits or policy.limits
    gates = {
        gate_id: _gate(
            gate_id,
            "not_run",
            "Gate was not run because its prerequisite evidence was unavailable.",
        )
        for gate_id in policy.required_gate_order
    }
    warnings: set[str] = {"EXTERNAL_PROPOSAL_REVIEW_ONLY"}
    target = canonical_smiles(proposal.target_smiles)
    precursors = canonical_smiles(proposal.precursor_smiles)
    if target is None or "." in target or precursors is None:
        gates["input_structure"] = _gate(
            "input_structure",
            "fail",
            "Target must be one valid molecular graph and every precursor must parse.",
            warnings=("INVALID_EXTERNAL_PROPOSAL_STRUCTURE",),
        )
        proposal_id = digest(
            "EXTPROP1", proposal.precursor_smiles, proposal.target_smiles
        )
        return _build_assessment(
            proposal=proposal,
            operator_library=operator_library,
            policy=policy,
            gates=gates,
            proposal_id=proposal_id,
            canonical_target=target,
            canonical_precursors=precursors,
            status="invalid",
            strongest_tier="unresolved",
            warnings=warnings | {"INVALID_EXTERNAL_PROPOSAL_STRUCTURE"},
        )

    proposal_id = digest("EXTPROP1", precursors, target)
    display_reaction = f"{precursors}>>{target}"
    gates["input_structure"] = _gate(
        "input_structure",
        "pass",
        "Target and precursor molecular graphs parsed and canonicalized.",
        evidence_ids=(proposal_id,),
    )

    normalized_mapped: str | None = None
    if proposal.mapped_reaction_smiles:
        mapping = validate_external_atom_mapping(
            display_reaction,
            proposal.mapped_reaction_smiles,
            provider_metadata=AtomMappingProviderMetadata(
                provider_id="external_proposal_supplied_mapping",
                provider_version="1.0",
                model_id="untrusted_external_input",
                model_sha256=None,
            ),
            mapper_confidence=None,
        )
        warnings.update(mapping.warnings)
        if not mapping.valid or not mapping.structure_preserved:
            gates["reaction_side_consistency"] = _gate(
                "reaction_side_consistency",
                "fail",
                "Supplied mapped reaction contradicts or does not validly map the displayed molecular graphs.",
                warnings=tuple(sorted(mapping.warnings + ((mapping.error,) if mapping.error else ()))),
            )
            gates["atom_correspondence"] = _gate(
                "atom_correspondence",
                "fail",
                "The supplied atom correspondence was rejected by structure-preserving validation.",
            )
            return _build_assessment(
                proposal=proposal,
                operator_library=operator_library,
                policy=policy,
                gates=gates,
                proposal_id=proposal_id,
                canonical_target=target,
                canonical_precursors=precursors,
                status="invalid",
                strongest_tier="unresolved",
                warnings=warnings | {"MAPPED_DISPLAY_STRUCTURE_CONTRADICTION"},
            )
        normalized_mapped = mapping.mapped_reaction_smiles
        gates["reaction_side_consistency"] = _gate(
            "reaction_side_consistency",
            "pass",
            "Map-stripped reaction sides exactly match the displayed proposal.",
        )
        gates["atom_correspondence"] = _gate(
            "atom_correspondence",
            "pass",
            "Externally supplied atom correspondence passed independent structural validation.",
            warnings=tuple(sorted(mapping.warnings)),
        )
    else:
        gates["reaction_side_consistency"] = _gate(
            "reaction_side_consistency",
            "pass",
            "No separate mapped reaction was supplied; the canonical display reaction is authoritative.",
        )
        materialized = materialize_atom_mapping(display_reaction)
        if materialized is None:
            gates["atom_correspondence"] = _gate(
                "atom_correspondence",
                "unresolved",
                "Deterministic atom correspondence could not be resolved without inventing a mapping.",
                warnings=("ATOM_CORRESPONDENCE_UNRESOLVED",),
            )
            return _build_assessment(
                proposal=proposal,
                operator_library=operator_library,
                policy=policy,
                gates=gates,
                proposal_id=proposal_id,
                canonical_target=target,
                canonical_precursors=precursors,
                status="ambiguous",
                strongest_tier="unresolved",
                warnings=warnings | {"ATOM_CORRESPONDENCE_UNRESOLVED"},
            )
        normalized_mapped = materialized.reaction_smiles
        gates["atom_correspondence"] = _gate(
            "atom_correspondence",
            "pass",
            f"Correspondence was resolved deterministically ({materialized.evidence}).",
        )

    analysis = featurize_reaction(normalized_mapped)
    completeness = analysis.reaction_completeness
    if completeness is None or completeness.status != "verified":
        completeness_status = completeness.status if completeness else "unavailable"
        gates["reaction_completeness"] = _gate(
            "reaction_completeness",
            "fail" if completeness_status == "incomplete" else "unresolved",
            f"Reaction completeness is {completeness_status}; complete product-atom accounting is required.",
            warnings=tuple(completeness.warnings if completeness else ()),
        )
    else:
        gates["reaction_completeness"] = _gate(
            "reaction_completeness",
            "pass",
            "Product-heavy-atom accounting is verified.",
            evidence_ids=(completeness.evidence,),
        )
    signature = analysis.reaction_signature
    if signature is None:
        gates["reaction_signature"] = _gate(
            "reaction_signature",
            "unresolved",
            "No complete generic reaction signature could be recovered.",
            warnings=tuple(sorted(analysis.warnings)),
        )
        return _build_assessment(
            proposal=proposal,
            operator_library=operator_library,
            policy=policy,
            gates=gates,
            proposal_id=proposal_id,
            canonical_target=target,
            canonical_precursors=precursors,
            normalized_mapped=normalized_mapped,
            status="ambiguous",
            strongest_tier="unresolved",
            reaction_signature=None,
            warnings=warnings | set(analysis.warnings),
        )
    gates["reaction_signature"] = _gate(
        "reaction_signature",
        "pass",
        "A complete graph-derived reaction signature was recovered.",
        evidence_ids=(signature.signature_id,),
    )
    warnings.update(analysis.warnings)
    if completeness is None or completeness.status != "verified":
        return _build_assessment(
            proposal=proposal,
            operator_library=operator_library,
            policy=policy,
            gates=gates,
            proposal_id=proposal_id,
            canonical_target=target,
            canonical_precursors=precursors,
            normalized_mapped=normalized_mapped,
            status="ambiguous",
            strongest_tier="unresolved",
            reaction_signature=signature,
            warnings=warnings | set(completeness.warnings if completeness else ()),
        )

    identity = analyze_generic_reaction(normalized_mapped)
    matches = (
        _operator_matches(identity, operator_library, limit=bounded.maximum_operator_matches)
        if identity is not None
        else ()
    )
    exact = tuple(
        item for item in matches if item.match_level == "exact_operator_signature"
    )
    if exact:
        gates["operator_support"] = _gate(
            "operator_support",
            "pass",
            f"Recovered signature exactly matches {len(exact)} admitted internal operator(s).",
            evidence_ids=tuple(item.operator_id for item in exact),
        )
    else:
        analogical = tuple(
            item for item in matches if item.match_level != "exact_operator_signature"
        )
        gates["operator_support"] = _gate(
            "operator_support",
            "unresolved",
            "No exact admitted operator signature was found; analogical matches remain non-validating evidence.",
            evidence_ids=tuple(item.operator_id for item in analogical),
        )

    precedents, selected_template_id = _precedents(
        proposal_id,
        target,
        precursors,
        exact,
        operator_library,
        limit=bounded.maximum_precedent_matches,
    )
    if precedents:
        gates["precedent_support"] = _gate(
            "precedent_support",
            "pass",
            f"Found {len(precedents)} admitted source-round-tripped precedent match(es).",
            evidence_ids=tuple(item.reaction_id for item in precedents),
        )
    elif exact:
        gates["precedent_support"] = _gate(
            "precedent_support",
            "unresolved",
            "The exact operator is admitted but no attached admitted precedent was available.",
        )
    else:
        gates["precedent_support"] = _gate(
            "precedent_support",
            "not_run",
            "Precedent retrieval requires an exact admitted operator.",
        )

    gates["stereochemistry"] = _gate(
        "stereochemistry",
        "pass",
        (
            f"Retained {len(signature.stereo_changes)} graph-derived stereochemical change(s)."
            if signature.stereo_changes
            else "No stereochemical change was asserted by the recovered signature."
        ),
    )
    precursor_compatibility = assess_precursor_compatibility(precursors)
    reaction_compatibility = assess_candidate_reaction_compatibility(
        normalized_mapped
    )
    try:
        competition = detect_functional_group_competition(normalized_mapped)
    except Exception:
        competition = None
        warnings.add("FUNCTIONAL_GROUP_SELECTIVITY_OUT_OF_SCOPE")
    compatibility = ExternalCompatibilityEvidence(
        precursor=precursor_compatibility,
        reaction=reaction_compatibility,
        functional_group_competition=competition,
    )
    gates["precursor_compatibility"] = _compatibility_gate(
        "precursor_compatibility",
        precursor_compatibility.disposition,
        "precursor",
        precursor_compatibility.policy_definition_id,
    )
    gates["reaction_compatibility"] = _compatibility_gate(
        "reaction_compatibility",
        reaction_compatibility.disposition,
        "reaction-regime",
        reaction_compatibility.policy_definition_id,
    )
    if competition is None:
        gates["functional_group_selectivity"] = _gate(
            "functional_group_selectivity",
            "pass",
            "No distinct alternative endpoint product was detected by the review-only selectivity audit.",
        )
    else:
        gates["functional_group_selectivity"] = _gate(
            "functional_group_selectivity",
            "warn",
            competition.message,
            warnings=(competition.code,),
        )
        warnings.add(competition.code)

    selected_match = exact[0] if exact else None
    forward_assessment = None
    if forward_library is None:
        gates["forward_reproduction"] = _gate(
            "forward_reproduction",
            "not_run",
            "No forward operator library was supplied.",
        )
        gates["forward_competition"] = _gate(
            "forward_competition",
            "not_run",
            "No target-hidden forward competition audit was requested.",
        )
    elif selected_match is None:
        gates["forward_reproduction"] = _gate(
            "forward_reproduction",
            "not_run",
            "Forward challenge requires an exact admitted operator match.",
        )
        gates["forward_competition"] = _gate(
            "forward_competition", "not_run", "Forward challenge was not run."
        )
    else:
        forward_assessment = assess_proposed_step(
            precursors,
            target,
            forward_library,
            recipe=proposal.proposed_conditions,
            operator_hint=selected_match.operator_id,
            top_k=bounded.maximum_forward_products,
            max_operators_to_apply=bounded.maximum_forward_operators,
        )
        reproduced = (
            forward_assessment.targeted_replay_status == "structurally_reproduced"
            or forward_assessment.intended_match in {"exact", "stereo_relaxed"}
        )
        gates["forward_reproduction"] = _gate(
            "forward_reproduction",
            "pass" if reproduced else "fail",
            (
                "The intended product was independently reproduced in the forward system."
                if reproduced
                else "The intended product was not independently reproduced in the forward system."
            ),
            warnings=tuple(forward_assessment.warnings),
        )
        competition_status: GateStatus = (
            "pass"
            if forward_assessment.disposition == "clear"
            else "warn"
            if forward_assessment.disposition in {"competitive", "unsupported", "out_of_scope"}
            else "fail"
        )
        gates["forward_competition"] = _gate(
            "forward_competition",
            competition_status,
            f"Forward competition disposition: {forward_assessment.disposition}.",
            warnings=tuple(forward_assessment.warnings),
        )

    condition_evidence = None
    if condition_evaluator is None:
        gates["condition_support"] = _gate(
            "condition_support",
            "unresolved" if proposal.proposed_conditions else "not_run",
            (
                "Proposed conditions were supplied but no canonical condition evaluator was available."
                if proposal.proposed_conditions
                else "No canonical condition evaluation was requested."
            ),
        )
    elif not exact:
        gates["condition_support"] = _gate(
            "condition_support",
            "not_run",
            "Condition evidence is evaluated only after exact operator support.",
        )
    else:
        condition_evidence = condition_evaluator(normalized_mapped)
        condition_status: GateStatus = (
            "pass"
            if condition_evidence.status
            in {"recommended_direct", "recommended_fallback"}
            else "unresolved"
        )
        gates["condition_support"] = _gate(
            "condition_support",
            condition_status,
            f"Canonical condition evidence status: {condition_evidence.status}.",
            warnings=tuple(condition_evidence.warnings),
        )

    hard_reject = (
        precursor_compatibility.disposition == "reject"
        or reaction_compatibility.disposition == "reject"
        or gates["forward_competition"].status == "fail"
    )
    strongest_tier = "structurally_verified"
    if exact:
        strongest_tier = "operator_supported"
    if precedents:
        strongest_tier = "precedent_supported"
    if (
        forward_assessment is not None
        and gates["forward_reproduction"].status == "pass"
    ):
        strongest_tier = "forward_reproduced"
    cautions = any(
        gates[gate_id].status == "warn"
        for gate_id in (
            "precursor_compatibility",
            "reaction_compatibility",
            "functional_group_selectivity",
            "forward_competition",
        )
    )
    if hard_reject:
        status: ExternalProposalStatus = "rejected_by_compatibility"
    elif cautions and _EVIDENCE_TIERS.index(strongest_tier) >= _EVIDENCE_TIERS.index("operator_supported"):
        status = "verified_with_cautions"
    else:
        status = strongest_tier  # type: ignore[assignment]
    eligible = (
        _EVIDENCE_TIERS.index(strongest_tier)
        >= _EVIDENCE_TIERS.index(policy.minimum_evidence_tier)
        and not hard_reject
    )
    parity = bool(exact and precedents and not hard_reject)
    return _build_assessment(
        proposal=proposal,
        operator_library=operator_library,
        policy=policy,
        gates=gates,
        proposal_id=proposal_id,
        canonical_target=target,
        canonical_precursors=precursors,
        normalized_mapped=normalized_mapped,
        status=status,
        strongest_tier=strongest_tier,
        admission_eligible=eligible,
        internal_confidence_parity=parity,
        reaction_signature=signature,
        reaction_identity=identity,
        operator_matches=matches,
        selected_operator_match_id=(selected_match.match_id if selected_match else None),
        selected_template_id=selected_template_id,
        precedent_matches=precedents,
        forward_assessment=forward_assessment,
        compatibility_evidence=compatibility,
        condition_evidence=condition_evidence,
        warnings=warnings,
    )


def _compatibility_gate(
    gate_id: str,
    disposition: str,
    label: str,
    definition_id: str,
) -> ExternalProposalGate:
    status: GateStatus = {
        "pass": "pass",
        "warn": "warn",
        "demote": "warn",
        "reject": "fail",
    }[disposition]
    return _gate(
        gate_id,
        status,
        f"Deterministic {label} compatibility disposition: {disposition}.",
        definition_id=definition_id,
    )


def _build_assessment(
    *,
    proposal: ExternalRetrosynthesisProposal,
    operator_library: GenericTemplateLibrary,
    policy: ExternalProposalAdmissionPolicy,
    gates: Mapping[str, ExternalProposalGate],
    proposal_id: str,
    canonical_target: str | None,
    canonical_precursors: str | None,
    status: ExternalProposalStatus,
    strongest_tier: str,
    warnings: set[str],
    normalized_mapped: str | None = None,
    admission_eligible: bool = False,
    internal_confidence_parity: bool = False,
    reaction_signature: ReactionSignature | None = None,
    reaction_identity: GenericReactionIdentity | None = None,
    operator_matches: tuple[ExternalOperatorMatch, ...] = (),
    selected_operator_match_id: str | None = None,
    selected_template_id: str | None = None,
    precedent_matches: tuple[StepPrecedentMatch, ...] = (),
    forward_assessment: RouteStepForwardAssessment | None = None,
    compatibility_evidence: ExternalCompatibilityEvidence | None = None,
    condition_evidence: RetrosynthesisConditionEvidence | None = None,
) -> ExternalRetrosynthesisAssessment:
    signature_id = reaction_signature.signature_id if reaction_signature else ""
    library_definition = str(
        operator_library.definition.get("definition_id")
        or operator_library.definition.get("algorithm_version")
        or "unspecified"
    )
    library_build = str(
        operator_library.definition.get("build_config_id")
        or operator_library.definition.get("definition_version")
        or library_definition
    )
    assessment_id = digest(
        "EXTASS1",
        proposal_id,
        signature_id,
        library_definition,
        library_build,
        operator_library.schema_version,
        policy.definition_id,
        EXTERNAL_PROPOSAL_ASSESSMENT_SCHEMA_VERSION,
    )
    all_gates = _finish_gates(policy, gates)
    for item in all_gates:
        warnings.update(item.warnings)
    return ExternalRetrosynthesisAssessment(
        assessment_id=assessment_id,
        proposal_id=proposal_id,
        external_proposal_id=proposal.external_proposal_id,
        canonical_target_smiles=canonical_target,
        canonical_precursor_smiles=canonical_precursors,
        normalized_mapped_reaction_smiles=normalized_mapped,
        status=status,
        strongest_evidence_tier=strongest_tier,
        admission_eligible=admission_eligible,
        actionable=False,
        internal_confidence_parity=internal_confidence_parity,
        gates=all_gates,
        reaction_signature=reaction_signature,
        reaction_identity=reaction_identity,
        operator_matches=operator_matches,
        selected_operator_match_id=selected_operator_match_id,
        selected_template_id=selected_template_id,
        precedent_matches=precedent_matches,
        forward_assessment=forward_assessment,
        compatibility_evidence=compatibility_evidence,
        condition_evidence=condition_evidence,
        sources=proposal.sources,
        warnings=tuple(sorted(warnings)),
        schema_versions={
            "external_proposal_assessment": EXTERNAL_PROPOSAL_ASSESSMENT_SCHEMA_VERSION,
            "external_proposal_admission_policy": policy.definition_id,
            "generic_template_library": operator_library.schema_version,
            "generic_template_definition": library_definition,
            "generic_template_build": library_build,
            "reaction_signature": (
                reaction_signature.schema_version if reaction_signature else "not_available"
            ),
        },
        ranking_influence=policy.ranking_influence,
    )


__all__ = [
    "EXTERNAL_PROPOSAL_ADMISSION_POLICY_PATH",
    "EXTERNAL_PROPOSAL_ASSESSMENT_SCHEMA_VERSION",
    "ExternalCompatibilityEvidence",
    "ExternalOperatorMatch",
    "ExternalProposalAdmissionPolicy",
    "ExternalProposalAssessmentLimits",
    "ExternalProposalGate",
    "ExternalProposalSource",
    "ExternalRetrosynthesisAssessment",
    "ExternalRetrosynthesisProposal",
    "assess_external_retrosynthesis_proposal",
    "load_external_proposal_admission_policy",
    "validate_external_proposal_admission_policy",
]
