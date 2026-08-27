"""Immutable application contracts for bounded chem-coworker assistance."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, replace
from enum import Enum
from typing import Any, Dict, Literal, Mapping, Optional, Tuple

from condition_registry import ConditionConstraint, ConditionConstraintSet
from core_retrosynthesis import RouteSearchPolicy


ASSISTANCE_SCHEMA_VERSION = "chem_coworker_assistance.v1"

AssistanceMode = Literal["conditions", "retro", "multistep"]
AssistanceStatus = Literal[
    "active",
    "completed",
    "needs_user_input",
    "abstained_insufficient_evidence",
    "blocked_by_policy",
    "budget_exhausted",
    "no_progress",
    "provider_failed",
    "domain_failed",
]
EvidenceLayer = Literal[
    "observation",
    "interpretation",
    "recommendation",
    "route",
    "user",
    "system",
]
EvidenceProvenance = Literal[
    "observed",
    "deterministic_inference",
    "user_confirmed",
    "advisory",
]
ConstraintProvenance = Literal[
    "explicit_user",
    "confirmed_user",
    "system_default",
    "model_proposed",
]
ActionStatus = Literal["completed", "rejected", "failed"]

_TERMINAL_STATUSES = frozenset(
    {
        "completed",
        "needs_user_input",
        "abstained_insufficient_evidence",
        "blocked_by_policy",
        "budget_exhausted",
        "no_progress",
        "provider_failed",
        "domain_failed",
    }
)


def _canonical_value(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if hasattr(value, "__dataclass_fields__"):
        return _canonical_value(asdict(value))
    if hasattr(value, "to_dict"):
        return _canonical_value(value.to_dict())
    if isinstance(value, Mapping):
        return {
            str(key): _canonical_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_value(item) for item in value]
    return value


def canonical_json(value: Any) -> str:
    """Serialize a contract value deterministically for IDs and replay."""

    return json.dumps(
        _canonical_value(value),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    )


def stable_assistance_id(prefix: str, value: Any) -> str:
    """Return a short stable identifier derived from normalized content."""

    digest = hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()[:20]
    return f"{prefix}:{digest}"


@dataclass(frozen=True)
class AssistanceBudget:
    """Validated upper bounds for one assistance session."""

    max_action_turns: int = 8
    max_clarification_cycles: int = 2
    max_provider_retries_per_turn: int = 1
    max_search_expansions: int = 3
    max_repair_attempts: int = 3
    max_total_tokens: int = 48_000
    max_elapsed_ms: int = 180_000

    def __post_init__(self) -> None:
        for name, value in asdict(self).items():
            if not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer")
        if self.max_action_turns < 1:
            raise ValueError("max_action_turns must be at least one")


@dataclass(frozen=True)
class AssistanceUsage:
    """Auditable resource usage accumulated by the application."""

    action_turns: int = 0
    clarification_cycles: int = 0
    search_expansions: int = 0
    repair_attempts: int = 0
    provider_attempts: int = 0
    input_tokens: int = 0
    output_tokens: int = 0
    elapsed_ms: int = 0

    def __post_init__(self) -> None:
        if any(value < 0 for value in asdict(self).values()):
            raise ValueError("assistance usage values must be non-negative")

    @property
    def total_tokens(self) -> int:
        return self.input_tokens + self.output_tokens


@dataclass(frozen=True)
class ConfirmedConstraint:
    """A typed constraint with explicit authority and ownership."""

    constraint_id: str
    owner: Literal[
        "chem_coworker",
        "condition_registry",
        "condition_recommender",
        "core_retrosynthesis",
    ]
    kind: str
    value: Any
    provenance: ConstraintProvenance
    hard: bool = False
    schema_version: str = ASSISTANCE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.constraint_id or not self.kind:
            raise ValueError("constraint_id and kind are required")
        if self.owner not in {
            "chem_coworker",
            "condition_registry",
            "condition_recommender",
            "core_retrosynthesis",
        }:
            raise ValueError("unsupported constraint owner")
        if self.provenance not in {
            "explicit_user",
            "confirmed_user",
            "system_default",
            "model_proposed",
        }:
            raise ValueError("unsupported constraint provenance")
        if self.hard and self.provenance == "model_proposed":
            raise ValueError("model-proposed constraints cannot be hard filters")


@dataclass(frozen=True)
class AssistanceRequest:
    """One user objective and the bounded authority granted to assistance."""

    objective: str
    mode: AssistanceMode
    structure_input: str
    constraints: Tuple[ConfirmedConstraint, ...] = ()
    condition_constraints: ConditionConstraintSet = field(
        default_factory=ConditionConstraintSet
    )
    route_search_policy: RouteSearchPolicy = field(default_factory=RouteSearchPolicy)
    budget: AssistanceBudget = field(default_factory=AssistanceBudget)
    provider_settings: Mapping[str, Any] = field(default_factory=dict)
    schema_version: str = ASSISTANCE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.mode not in {"conditions", "retro", "multistep"}:
            raise ValueError("unsupported assistance mode")
        if not self.objective.strip():
            raise ValueError("assistance objective is required")
        if not self.structure_input.strip():
            raise ValueError("structure_input is required")
        constraint_ids = [item.constraint_id for item in self.constraints]
        if len(constraint_ids) != len(set(constraint_ids)):
            raise ValueError("constraint IDs must be unique")
        object.__setattr__(self, "provider_settings", dict(self.provider_settings))


@dataclass(frozen=True)
class EvidenceItem:
    """One bounded model-visible fact with authority and provenance."""

    evidence_id: str
    layer: EvidenceLayer
    source_id: str
    payload_type: str
    payload: Mapping[str, Any]
    provenance: EvidenceProvenance
    schema_versions: Mapping[str, str] = field(default_factory=dict)
    uncertainty: Optional[str] = None
    schema_version: str = ASSISTANCE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.evidence_id or not self.source_id or not self.payload_type:
            raise ValueError("evidence ID, source ID, and payload type are required")
        if self.layer not in {
            "observation",
            "interpretation",
            "recommendation",
            "route",
            "user",
            "system",
        }:
            raise ValueError("unsupported evidence layer")
        if self.provenance not in {
            "observed",
            "deterministic_inference",
            "user_confirmed",
            "advisory",
        }:
            raise ValueError("unsupported evidence provenance")
        object.__setattr__(self, "payload", dict(self.payload))
        object.__setattr__(self, "schema_versions", dict(self.schema_versions))


@dataclass(frozen=True)
class AdvisoryClaim:
    """A non-authoritative statement linked to known evidence."""

    claim_type: str
    subject_id: str
    severity: Literal["info", "caution", "serious"]
    statement: str
    evidence_ids: Tuple[str, ...]
    uncertainty: str
    suggested_user_action: Optional[str] = None
    schema_version: str = ASSISTANCE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.claim_type or not self.subject_id or not self.statement.strip():
            raise ValueError("claim type, subject, and statement are required")
        if self.severity not in {"info", "caution", "serious"}:
            raise ValueError("unsupported advisory severity")
        if not self.evidence_ids:
            raise ValueError("advisory claims require at least one evidence ID")


@dataclass(frozen=True)
class ClarificationQuestion:
    """A focused question whose answer can change a deterministic rerun."""

    question_id: str
    prompt: str
    constraint_owner: str
    constraint_kind: str
    reason: str
    evidence_ids: Tuple[str, ...]
    schema_version: str = ASSISTANCE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.question_id or not self.prompt.strip() or not self.reason.strip():
            raise ValueError("clarification ID, prompt, and reason are required")


@dataclass(frozen=True)
class ActionRecord:
    """One validated controller action, without private model reasoning."""

    action_id: str
    action_name: str
    normalized_arguments: Mapping[str, Any]
    cited_evidence_ids: Tuple[str, ...]
    status: ActionStatus
    result_ref: Optional[str] = None
    added_evidence_ids: Tuple[str, ...] = ()
    decision_summary: str = ""
    error: Optional[str] = None
    provider_attempts: int = 0
    input_tokens: int = 0
    output_tokens: int = 0
    elapsed_ms: int = 0
    schema_version: str = ASSISTANCE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.action_id or not self.action_name:
            raise ValueError("action ID and action name are required")
        if self.status not in {"completed", "rejected", "failed"}:
            raise ValueError("unsupported action status")
        if min(
            self.provider_attempts,
            self.input_tokens,
            self.output_tokens,
            self.elapsed_ms,
        ) < 0:
            raise ValueError("action usage values must be non-negative")
        object.__setattr__(self, "normalized_arguments", dict(self.normalized_arguments))

    @property
    def fingerprint(self) -> str:
        return stable_assistance_id(
            "ACT",
            {"name": self.action_name, "arguments": self.normalized_arguments},
        )


@dataclass(frozen=True)
class AssistanceSessionState:
    """Application-owned replayable state for one bounded assistance loop."""

    session_id: str
    request: AssistanceRequest
    turn: int = 0
    domain_result_refs: Tuple[str, ...] = ()
    evidence: Tuple[EvidenceItem, ...] = ()
    claims: Tuple[AdvisoryClaim, ...] = ()
    unresolved_questions: Tuple[ClarificationQuestion, ...] = ()
    action_history: Tuple[ActionRecord, ...] = ()
    usage: AssistanceUsage = field(default_factory=AssistanceUsage)
    status: AssistanceStatus = "active"
    stopping_reason: Optional[str] = None
    schema_version: str = ASSISTANCE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.session_id or self.turn < 0:
            raise ValueError("valid session ID and non-negative turn are required")
        if self.status not in {"active", *_TERMINAL_STATUSES}:
            raise ValueError("unsupported assistance status")
        evidence_ids = [item.evidence_id for item in self.evidence]
        if len(evidence_ids) != len(set(evidence_ids)):
            raise ValueError("evidence IDs must be unique within a session")
        action_ids = [item.action_id for item in self.action_history]
        if len(action_ids) != len(set(action_ids)):
            raise ValueError("action IDs must be unique within a session")
        allowed_evidence = set(evidence_ids)
        referenced = {
            evidence_id
            for claim in self.claims
            for evidence_id in claim.evidence_ids
        }
        referenced.update(
            evidence_id
            for question in self.unresolved_questions
            for evidence_id in question.evidence_ids
        )
        referenced.update(
            evidence_id
            for action in self.action_history
            for evidence_id in action.cited_evidence_ids
        )
        unknown = referenced - allowed_evidence
        if unknown:
            raise ValueError(f"unknown evidence IDs: {sorted(unknown)!r}")
        if self.status in _TERMINAL_STATUSES and not self.stopping_reason:
            raise ValueError("terminal assistance states require a stopping reason")
        if self.status == "active" and self.stopping_reason:
            raise ValueError("active assistance states cannot have a stopping reason")

    @property
    def is_terminal(self) -> bool:
        return self.status in _TERMINAL_STATUSES

    @property
    def remaining_action_turns(self) -> int:
        return max(0, self.request.budget.max_action_turns - self.usage.action_turns)

    @property
    def allowed_evidence_ids(self) -> frozenset[str]:
        return frozenset(item.evidence_id for item in self.evidence)

    def to_dict(self) -> Dict[str, Any]:
        """Return a canonical JSON-compatible replay snapshot."""

        return _canonical_value(self)

    def to_json(self) -> str:
        """Return deterministic compact JSON for storage or replay comparison."""

        return canonical_json(self)


def new_session(request: AssistanceRequest) -> AssistanceSessionState:
    """Create a deterministic initial state from a normalized request."""

    return AssistanceSessionState(
        session_id=stable_assistance_id("SESSION", request),
        request=request,
    )


def add_evidence(
    state: AssistanceSessionState,
    items: Tuple[EvidenceItem, ...],
    *,
    domain_result_ref: Optional[str] = None,
) -> AssistanceSessionState:
    """Return a state with new, non-conflicting evidence appended."""

    if state.is_terminal:
        raise ValueError("cannot add evidence to a terminal session")
    existing = {item.evidence_id: item for item in state.evidence}
    for item in items:
        prior = existing.get(item.evidence_id)
        if prior is not None and prior != item:
            raise ValueError(f"evidence ID collision: {item.evidence_id}")
        existing[item.evidence_id] = item
    refs = state.domain_result_refs
    if domain_result_ref and domain_result_ref not in refs:
        refs = refs + (domain_result_ref,)
    return replace(state, evidence=tuple(existing.values()), domain_result_refs=refs)


def record_action(
    state: AssistanceSessionState,
    action: ActionRecord,
) -> AssistanceSessionState:
    """Apply one validated action record and its measured usage."""

    if state.is_terminal:
        raise ValueError("cannot record an action on a terminal session")
    if action.fingerprint in {item.fingerprint for item in state.action_history}:
        raise ValueError("identical action input has already been attempted")
    if not set(action.cited_evidence_ids).issubset(state.allowed_evidence_ids):
        raise ValueError("action cites unknown evidence")
    usage = replace(
        state.usage,
        action_turns=state.usage.action_turns + 1,
        provider_attempts=state.usage.provider_attempts + action.provider_attempts,
        input_tokens=state.usage.input_tokens + action.input_tokens,
        output_tokens=state.usage.output_tokens + action.output_tokens,
        elapsed_ms=state.usage.elapsed_ms + action.elapsed_ms,
    )
    return replace(
        state,
        turn=state.turn + 1,
        action_history=state.action_history + (action,),
        usage=usage,
    )


def finish_session(
    state: AssistanceSessionState,
    *,
    status: AssistanceStatus,
    stopping_reason: str,
    claims: Tuple[AdvisoryClaim, ...] = (),
    questions: Tuple[ClarificationQuestion, ...] = (),
) -> AssistanceSessionState:
    """Validate evidence links and return a terminal state."""

    if state.is_terminal:
        raise ValueError("cannot finish an already terminal session")
    if status not in _TERMINAL_STATUSES:
        raise ValueError("finish_session requires a terminal status")
    return replace(
        state,
        claims=claims,
        unresolved_questions=questions,
        status=status,
        stopping_reason=stopping_reason,
    )


def session_from_dict(value: Mapping[str, Any]) -> AssistanceSessionState:
    """Rebuild and fully revalidate a serialized assistance state."""

    raw_request = dict(value["request"])
    raw_budget = dict(raw_request.pop("budget", {}))
    raw_conditions = dict(raw_request.pop("condition_constraints", {}))
    raw_route_policy = dict(raw_request.pop("route_search_policy", {}))
    constraints = tuple(
        ConfirmedConstraint(**dict(item))
        for item in raw_request.pop("constraints", ())
    )
    condition_constraints = ConditionConstraintSet(
        constraints=tuple(
            ConditionConstraint(**dict(item))
            for item in raw_conditions.get("constraints", ())
        ),
        schema_version=str(
            raw_conditions.get("schema_version")
            or ConditionConstraintSet().schema_version
        ),
    )
    request = AssistanceRequest(
        **raw_request,
        constraints=constraints,
        budget=AssistanceBudget(**raw_budget),
        condition_constraints=condition_constraints,
        route_search_policy=RouteSearchPolicy(**raw_route_policy),
    )
    return AssistanceSessionState(
        session_id=str(value["session_id"]),
        request=request,
        turn=int(value.get("turn", 0)),
        domain_result_refs=tuple(value.get("domain_result_refs") or ()),
        evidence=tuple(
            EvidenceItem(**dict(item)) for item in value.get("evidence") or ()
        ),
        claims=tuple(
            AdvisoryClaim(
                **{
                    **dict(item),
                    "evidence_ids": tuple(item.get("evidence_ids") or ()),
                }
            )
            for item in value.get("claims") or ()
        ),
        unresolved_questions=tuple(
            ClarificationQuestion(
                **{
                    **dict(item),
                    "evidence_ids": tuple(item.get("evidence_ids") or ()),
                }
            )
            for item in value.get("unresolved_questions") or ()
        ),
        action_history=tuple(
            ActionRecord(
                **{
                    **dict(item),
                    "cited_evidence_ids": tuple(
                        item.get("cited_evidence_ids") or ()
                    ),
                    "added_evidence_ids": tuple(
                        item.get("added_evidence_ids") or ()
                    ),
                }
            )
            for item in value.get("action_history") or ()
        ),
        usage=AssistanceUsage(**dict(value.get("usage") or {})),
        status=value.get("status", "active"),
        stopping_reason=value.get("stopping_reason"),
        schema_version=str(value.get("schema_version") or ASSISTANCE_SCHEMA_VERSION),
    )


def session_from_json(value: str) -> AssistanceSessionState:
    """Parse a replay snapshot without trusting serialized transitions."""

    parsed = json.loads(value)
    if not isinstance(parsed, dict):
        raise ValueError("assistance session JSON must contain an object")
    return session_from_dict(parsed)
