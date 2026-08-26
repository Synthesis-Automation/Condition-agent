"""Bounded, evidence-grounded assistance over deterministic domain services.

Importing this package does not import an LLM SDK. Provider integrations live
behind the transport boundary and are loaded only when a controller is run.
"""

from .contracts import (
    ActionRecord,
    AdvisoryClaim,
    AssistanceBudget,
    AssistanceRequest,
    AssistanceSessionState,
    AssistanceUsage,
    ClarificationQuestion,
    ConfirmedConstraint,
    EvidenceItem,
    new_session,
    session_from_dict,
    session_from_json,
)
from .evidence import (
    ConditionEvidenceProjection,
    MultistepEvidenceProjection,
    RetrosynthesisEvidenceProjection,
    project_condition_result,
    project_multistep_response,
    project_retrosynthesis_response,
)
from .capabilities import (
    ChemistryCapabilities,
    ConditionCapabilities,
    MultistepCapabilities,
    RetrosynthesisCapabilities,
)
from .controller import AssistanceController, AssistanceRunResult
from .policy import AssistancePolicy, load_assistance_policy
from .rendering import render_assistance
from .service import AssistanceApplicationService
from .evaluation import (
    AssistanceEvaluationCase,
    AssistanceEvaluationResult,
    build_blind_review_packet,
    evaluate_assistance_run,
    load_evaluation_rubric,
)
from .transport import (
    AssistanceProviderSettings,
    AssistanceTransportResult,
    OpenAICompatibleAssistanceTransport,
)

__all__ = [
    "ActionRecord",
    "AdvisoryClaim",
    "AssistanceBudget",
    "AssistanceApplicationService",
    "AssistanceController",
    "AssistanceEvaluationCase",
    "AssistanceEvaluationResult",
    "AssistancePolicy",
    "AssistanceProviderSettings",
    "AssistanceRequest",
    "AssistanceRunResult",
    "AssistanceSessionState",
    "AssistanceTransportResult",
    "AssistanceUsage",
    "ClarificationQuestion",
    "ChemistryCapabilities",
    "ConditionEvidenceProjection",
    "ConditionCapabilities",
    "MultistepCapabilities",
    "MultistepEvidenceProjection",
    "ConfirmedConstraint",
    "EvidenceItem",
    "OpenAICompatibleAssistanceTransport",
    "RetrosynthesisCapabilities",
    "RetrosynthesisEvidenceProjection",
    "load_assistance_policy",
    "load_evaluation_rubric",
    "new_session",
    "session_from_dict",
    "session_from_json",
    "project_condition_result",
    "project_multistep_response",
    "project_retrosynthesis_response",
    "render_assistance",
    "build_blind_review_packet",
    "evaluate_assistance_run",
]
