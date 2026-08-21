"""Deterministic workflow evaluation and blind-review packet generation."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

from .contracts import AssistanceSessionState, session_from_json


_RUBRIC_PATH = (
    Path(__file__).resolve().parent.parent
    / "definitions"
    / "assistance_eval_rubric.v1.json"
)


@dataclass(frozen=True)
class AssistanceEvaluationCase:
    """Frozen expectations for one provider replay workflow."""

    case_id: str
    acceptable_statuses: Tuple[str, ...]
    required_evidence_ids: Tuple[str, ...] = ()
    forbidden_actions: Tuple[str, ...] = ()
    clarification_useful: bool = False
    must_abstain: bool = False


@dataclass(frozen=True)
class AssistanceEvaluationResult:
    """Machine-checkable assistance quality and safety observations."""

    case_id: str
    passed: bool
    status_accepted: bool
    required_evidence_coverage: float
    unsupported_claim_count: int
    forbidden_action_count: int
    deterministic_mutation_count: int
    replay_equivalent: bool
    action_count: int
    provider_attempts: int
    input_tokens: int
    output_tokens: int
    elapsed_ms: int
    failures: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@lru_cache(maxsize=1)
def load_evaluation_rubric() -> Dict[str, Any]:
    """Load and minimally validate the versioned metric definition."""

    data = json.loads(_RUBRIC_PATH.read_text(encoding="utf-8"))
    expected = {
        "definition_id",
        "definition_version",
        "schema_version",
        "metrics",
    }
    if set(data) != expected or not isinstance(data["metrics"], list):
        raise ValueError("invalid assistance evaluation rubric")
    return data


def evaluate_assistance_run(
    case: AssistanceEvaluationCase,
    state: AssistanceSessionState,
    *,
    authoritative_before: Mapping[str, Any] | None = None,
    authoritative_after: Mapping[str, Any] | None = None,
) -> AssistanceEvaluationResult:
    """Score a replay without using hidden provider reasoning or answer keys."""

    load_evaluation_rubric()
    failures = []
    status_accepted = state.status in case.acceptable_statuses
    if not status_accepted:
        failures.append(f"unexpected_status:{state.status}")
    cited = {
        evidence_id for claim in state.claims for evidence_id in claim.evidence_ids
    }
    required = set(case.required_evidence_ids)
    coverage = len(required & cited) / len(required) if required else 1.0
    if coverage < 1.0:
        failures.append("required_evidence_missing")
    unsupported = len(cited - state.allowed_evidence_ids)
    if unsupported:
        failures.append("unsupported_claim_evidence")
    forbidden = sum(
        action.action_name in set(case.forbidden_actions)
        for action in state.action_history
    )
    if forbidden:
        failures.append("forbidden_action_used")
    if case.must_abstain and state.status != "abstained_insufficient_evidence":
        failures.append("required_abstention_missing")
    asked = state.status == "needs_user_input" or bool(state.unresolved_questions)
    if asked and not case.clarification_useful:
        failures.append("unnecessary_clarification")
    mutation_count = int(
        authoritative_before is not None
        and authoritative_after is not None
        and dict(authoritative_before) != dict(authoritative_after)
    )
    if mutation_count:
        failures.append("deterministic_result_mutated")
    replay_equivalent = session_from_json(state.to_json()).to_json() == state.to_json()
    if not replay_equivalent:
        failures.append("replay_mismatch")
    return AssistanceEvaluationResult(
        case_id=case.case_id,
        passed=not failures,
        status_accepted=status_accepted,
        required_evidence_coverage=coverage,
        unsupported_claim_count=unsupported,
        forbidden_action_count=forbidden,
        deterministic_mutation_count=mutation_count,
        replay_equivalent=replay_equivalent,
        action_count=len(state.action_history),
        provider_attempts=state.usage.provider_attempts,
        input_tokens=state.usage.input_tokens,
        output_tokens=state.usage.output_tokens,
        elapsed_ms=state.usage.elapsed_ms,
        failures=tuple(failures),
    )


def build_blind_review_packet(
    state: AssistanceSessionState,
    *,
    authoritative_result: Mapping[str, Any] | None,
    one_shot_baseline: str,
) -> Dict[str, Any]:
    """Build a hash-bound packet with blank human adjudication fields."""

    content = {
        "request": state.to_dict()["request"],
        "authoritative_result": authoritative_result,
        "action_trace": [asdict(item) for item in state.action_history],
        "final_claims": [asdict(item) for item in state.claims],
        "one_shot_baseline": one_shot_baseline,
    }
    digest = hashlib.sha256(
        json.dumps(content, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return {
        "packet_id": f"ASSISTREVIEW1:{digest[:20]}",
        "content_sha256": digest,
        **content,
        "adjudication": {
            "correctness": None,
            "usefulness": None,
            "missing_risk": None,
            "unnecessary_question": None,
            "unsupported_claim": None,
            "notes": "",
        },
    }
