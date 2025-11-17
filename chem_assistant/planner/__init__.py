"""
Planner APIs for the auto-conditions mode.

This package exposes light-weight schemas and thin wrappers that orchestrate
existing chemtools capabilities (rule engine, DRFP precedent search, HTE
analytics) so higher-level agents can compose them deterministically.
"""

from .api import (
    ReactionInput,
    DetectedFamily,
    CandidateCondition,
    EvidenceScores,
    HteSummary,
    FusionResult,
    ProtocolOutput,
    PlannerStep,
    PlannerPlan,
    AutoConditionsResult,
    detect_family,
    fetch_rule_candidates,
    find_similar_protocols,
    fetch_hte_stats,
    score_ml_candidates,
    fuse_scores,
    build_protocol,
    plan_workflow,
    auto_conditions,
)

__all__ = [
    "ReactionInput",
    "DetectedFamily",
    "CandidateCondition",
    "EvidenceScores",
    "HteSummary",
    "FusionResult",
    "ProtocolOutput",
    "PlannerStep",
    "PlannerPlan",
    "AutoConditionsResult",
    "detect_family",
    "fetch_rule_candidates",
    "find_similar_protocols",
    "fetch_hte_stats",
    "score_ml_candidates",
    "fuse_scores",
    "build_protocol",
    "plan_workflow",
    "auto_conditions",
]
