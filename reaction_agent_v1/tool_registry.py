"""Tool registry for reaction agent foundation."""

from __future__ import annotations

from dataclasses import asdict
import time
from typing import Any, Callable, Dict, List

from poc_gpt52_reaction_v2 import analyze_reaction_general

from .contracts import ToolExecutionResult
from .coverage_advisor import CoverageAdvisor
from .validator import AgentDecisionValidator


ToolHandler = Callable[..., Dict[str, Any]]


class ToolRegistry:
    """Typed tool registry with runtime dispatch."""

    def __init__(self) -> None:
        self._handlers: Dict[str, ToolHandler] = {}

    def register(self, tool_name: str, handler: ToolHandler) -> None:
        """Register a tool handler."""
        self._handlers[tool_name] = handler

    def execute(self, tool_name: str, **kwargs: Any) -> ToolExecutionResult:
        """Execute a registered tool with kwargs."""
        handler = self._handlers.get(tool_name)
        if handler is None:
            return ToolExecutionResult(
                tool_name=tool_name,
                ok=False,
                error=f"Tool not found: {tool_name}",
            )

        start = time.perf_counter()
        try:
            payload = handler(**kwargs)
            latency_ms = (time.perf_counter() - start) * 1000
            return ToolExecutionResult(
                tool_name=tool_name,
                ok=True,
                payload=payload if isinstance(payload, dict) else {"result": payload},
                latency_ms=latency_ms,
            )
        except Exception as exc:
            latency_ms = (time.perf_counter() - start) * 1000
            return ToolExecutionResult(
                tool_name=tool_name,
                ok=False,
                error=str(exc),
                latency_ms=latency_ms,
            )


def _placeholder_payload(
    *,
    tool_name: str,
    summary: str,
    expected_inputs: List[str],
    expected_outputs: List[str],
    next_action: str,
) -> Dict[str, Any]:
    return {
        "status": "not_implemented",
        "tool": tool_name,
        "summary": summary,
        "expected_inputs": expected_inputs,
        "expected_outputs": expected_outputs,
        "next_action": next_action,
    }


def build_default_registry(*, min_confidence: float = 0.5) -> ToolRegistry:
    """Build default tool registry with deterministic chemistry tools."""
    registry = ToolRegistry()
    validator = AgentDecisionValidator(min_confidence=min_confidence)
    advisor = CoverageAdvisor()

    def reaction_diff_tool(reaction_smiles: str) -> Dict[str, Any]:
        result = analyze_reaction_general(
            reaction_smiles,
            use_llm=False,
            min_confidence=min_confidence,
        ).to_dict()
        return {
            "status": "ok",
            "analysis": result,
            "diff": {
                "principal_pair": result.get("principal_pair", {}),
                "mcs_smarts": result.get("mcs_smarts"),
                "mcs_atoms": result.get("mcs_atoms"),
                "mcs_ratio": result.get("mcs_ratio"),
                "core_formula_delta": result.get("core_formula_delta", {}),
                "side_formula_delta": result.get("side_formula_delta", {}),
                "core_reactant_formula": result.get("core_reactant_formula"),
                "core_product_formula": result.get("core_product_formula"),
                "side_reactant_formula": result.get("side_reactant_formula"),
                "side_product_formula": result.get("side_product_formula"),
            },
            "detection": {
                "detection_error": result.get("detection_error"),
                "reacted_motifs": result.get("reacted_motifs", []),
                "formed_motifs": result.get("formed_motifs", []),
                "reaction_key": result.get("reaction_key", ""),
            },
        }

    def validate_decision_tool(evidence: Dict[str, Any]) -> Dict[str, Any]:
        report = validator.validate_evidence(evidence)
        return asdict(report)

    def coverage_advice_tool(reaction_smiles: str, evidence: Dict[str, Any]) -> Dict[str, Any]:
        analysis_like = dict(evidence.get("analysis_snapshot") or {})
        if not analysis_like:
            analysis_like = {
                "decision": evidence.get("provisional_decision", {}),
                "taxonomy_candidates": evidence.get("taxonomy_candidates", []),
            }
        diff = dict(evidence.get("diff") or {})
        detection = dict(evidence.get("detection") or {})
        analysis_like.update(diff)
        analysis_like.update(detection)
        cards = advisor.suggest(reaction_smiles=reaction_smiles, analysis=analysis_like)
        return {"suggestions": cards}

    def fallback_candidate_retrieval_tool(reaction_smiles: str, evidence: Dict[str, Any]) -> Dict[str, Any]:
        return _placeholder_payload(
            tool_name="fallback_candidate_retrieval",
            summary="Fallback candidate retrieval based on bond-change and motif deltas.",
            expected_inputs=[
                "reaction_smiles",
                "evidence.detection.reaction_key",
                "evidence.diff.core_formula_delta",
                "evidence.diff.principal_pair",
            ],
            expected_outputs=[
                "candidate_reaction_types",
                "candidate_scores",
                "retrieval_evidence",
            ],
            next_action="Merge candidates into evidence.taxonomy_candidates before validate_decision.",
        )

    def confidence_calibrator_tool(evidence: Dict[str, Any], validation: Dict[str, Any]) -> Dict[str, Any]:
        return _placeholder_payload(
            tool_name="confidence_calibrator",
            summary="Calibrate confidence scores against reliability bins.",
            expected_inputs=[
                "evidence.taxonomy_candidates",
                "validation.final_decision",
                "evidence.diff.mcs_ratio",
            ],
            expected_outputs=[
                "calibrated_confidence",
                "confidence_band",
                "calibration_metadata",
            ],
            next_action="Update final decision confidence only if calibration model is available.",
        )

    def llm_rerank_constrained_tool(evidence: Dict[str, Any], validation: Dict[str, Any]) -> Dict[str, Any]:
        return _placeholder_payload(
            tool_name="llm_rerank_constrained",
            summary="Constrained LLM reranking over allowed taxonomy candidates.",
            expected_inputs=[
                "evidence.taxonomy_candidates",
                "evidence.detection.reacted_motifs",
                "evidence.detection.formed_motifs",
                "validation.final_decision",
            ],
            expected_outputs=[
                "reranked_candidates",
                "selected_reaction_type",
                "rationale",
            ],
            next_action="Apply rerank only if selected type is in allowed candidate IDs and passes validator.",
        )

    def precedent_lookup_tool(reaction_smiles: str, evidence: Dict[str, Any]) -> Dict[str, Any]:
        return _placeholder_payload(
            tool_name="precedent_lookup",
            summary="Retrieve nearest precedent reactions from internal datasets.",
            expected_inputs=[
                "reaction_smiles",
                "evidence.provisional_decision.reaction_type",
                "evidence.diff.principal_pair",
            ],
            expected_outputs=[
                "precedent_hits",
                "support_score",
                "condition_patterns",
            ],
            next_action="Attach precedent evidence to final report and confidence calibration.",
        )

    registry.register("reaction_diff", reaction_diff_tool)
    registry.register("fallback_candidate_retrieval", fallback_candidate_retrieval_tool)
    registry.register("validate_decision", validate_decision_tool)
    registry.register("confidence_calibrator", confidence_calibrator_tool)
    registry.register("llm_rerank_constrained", llm_rerank_constrained_tool)
    registry.register("precedent_lookup", precedent_lookup_tool)
    registry.register("coverage_advice", coverage_advice_tool)
    return registry
