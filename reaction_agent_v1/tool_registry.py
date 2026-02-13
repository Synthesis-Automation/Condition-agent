"""Tool registry for reaction agent foundation."""

from __future__ import annotations

from dataclasses import asdict
import time
from typing import Any, Callable, Dict, List

from poc_gpt52_reaction_v2 import analyze_reaction_general
from chemtools.taxonomy.loader import load_reaction_types_dict
from chemtools.util.rdkit_helpers import parse_smiles
from chemtools.util.smarts_cache import compile_smarts

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


def _split_reactants(reaction_smiles: str) -> List[str]:
    left = str(reaction_smiles or "").split(">>", 1)[0]
    return [token.strip() for token in left.split(".") if token.strip()]


def _count_aromatic_ring_n(smiles: str) -> int:
    mol = parse_smiles(smiles)
    if mol is None:
        return 0
    return sum(
        1
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "N" and atom.GetIsAromatic() and atom.IsInRing()
    )


def _has_any_pattern(smiles_list: List[str], smarts: str) -> bool:
    pattern = compile_smarts(smarts, validate=False)
    if not pattern:
        return False
    for smiles in smiles_list:
        mol = parse_smiles(smiles)
        if mol and mol.HasSubstructMatch(pattern):
            return True
    return False


def _fallback_rank_candidates(
    *,
    reaction_smiles: str,
    evidence: Dict[str, Any],
) -> List[Dict[str, Any]]:
    taxonomy = load_reaction_types_dict()
    principal_pair = (evidence.get("diff") or {}).get("principal_pair") or {}
    principal_reactant = str(principal_pair.get("reactant_smiles") or "")
    reacted = {str(m) for m in ((evidence.get("detection") or {}).get("reacted_motifs") or [])}
    formed = {str(m) for m in ((evidence.get("detection") or {}).get("formed_motifs") or [])}
    reaction_key = str((evidence.get("detection") or {}).get("reaction_key") or "")
    reaction_key_lower = reaction_key.lower()
    reacted_upper = {item.upper() for item in reacted}

    reactants = _split_reactants(reaction_smiles)
    has_hydrazine_like = _has_any_pattern(reactants, "[N]-[N]") or _has_any_pattern(reactants, "[N]=[N]")
    has_boron = _has_any_pattern(reactants, "[B]")
    has_grignard = _has_any_pattern(reactants, "[Mg]")
    has_amine = _has_any_pattern(reactants, "[NX3;H1,H2]")
    aromatic_ring_n = _count_aromatic_ring_n(principal_reactant)

    aryl_halide_reacted = bool({"AR-CL", "AR-BR", "AR-I"} & reacted_upper)
    c_n_event = ("lgdisp+c-n" in reaction_key_lower) or ("bond_formed: c(ar)-n" in reaction_key_lower)
    c_c_event = ("lgdisp+c-c" in reaction_key_lower) or ("bond_formed: c(ar)-c" in reaction_key_lower)

    scored: Dict[str, Dict[str, Any]] = {}

    def add_candidate(reaction_type: str, base_score: float, reason: str, signals: Dict[str, Any]) -> None:
        if reaction_type not in taxonomy:
            return
        score = max(0.0, min(0.99, base_score))
        existing = scored.get(reaction_type)
        if existing is None or score > float(existing["score"]):
            meta = taxonomy.get(reaction_type) or {}
            scored[reaction_type] = {
                "reaction_type": reaction_type,
                "score": round(score, 2),
                "taxonomy_name": str(meta.get("name") or reaction_type),
                "taxonomy_category": str(meta.get("category") or "unknown"),
                "reason": reason,
                "signals": signals,
            }

    if aryl_halide_reacted and c_n_event:
        score = 0.58
        if has_amine:
            score += 0.08
        add_candidate(
            "C_N_Coupling",
            score,
            "Aryl-halide leaving-group displacement with C-N bond formation signal.",
            {
                "aryl_halide_reacted": aryl_halide_reacted,
                "c_n_event": c_n_event,
                "has_amine": has_amine,
            },
        )

        snar_score = 0.62
        if aromatic_ring_n >= 2:
            snar_score += 0.08
        if has_hydrazine_like:
            snar_score += 0.1
        if "pyrimidine" in reaction_key_lower:
            snar_score += 0.08
        add_candidate(
            "SNAr_CN",
            snar_score,
            "Electron-poor aromatic substitution-like signature with C-N displacement.",
            {
                "aryl_halide_reacted": aryl_halide_reacted,
                "c_n_event": c_n_event,
                "principal_aromatic_ring_n": aromatic_ring_n,
                "has_hydrazine_like": has_hydrazine_like,
            },
        )

    if aryl_halide_reacted and c_c_event:
        if has_boron:
            add_candidate(
                "Suzuki_miyaura",
                0.76,
                "Aryl-halide + boron reagent with C-C formation signature.",
                {
                    "aryl_halide_reacted": aryl_halide_reacted,
                    "c_c_event": c_c_event,
                    "has_boron": has_boron,
                },
            )
        if has_grignard:
            add_candidate(
                "Kumada",
                0.74,
                "Aryl-halide + organomagnesium reagent with C-C formation signature.",
                {
                    "aryl_halide_reacted": aryl_halide_reacted,
                    "c_c_event": c_c_event,
                    "has_grignard": has_grignard,
                },
            )

    if not aryl_halide_reacted and has_boron and has_amine and c_n_event:
        add_candidate(
            "Chan_Lam_C_N_Coupling",
            0.62,
            "Boron + amine C-N coupling signature without aryl-halide leaving group.",
            {
                "has_boron": has_boron,
                "has_amine": has_amine,
                "c_n_event": c_n_event,
            },
        )

    ranked = sorted(scored.values(), key=lambda item: (-item["score"], item["reaction_type"]))
    candidates: List[Dict[str, Any]] = []
    for row in ranked:
        candidates.append(
            {
                "reaction_type": row["reaction_type"],
                "deterministic_score": row["score"],
                "detector_confidence": 0.0,
                "taxonomy_name": row["taxonomy_name"],
                "taxonomy_category": row["taxonomy_category"],
                "evidence": {
                    "fallback_rule": row["reason"],
                    "signals": row["signals"],
                    "reacted_motifs": sorted(reacted),
                    "formed_motifs": sorted(formed),
                    "reaction_key": reaction_key,
                },
            }
        )
    return candidates


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
        candidates = _fallback_rank_candidates(reaction_smiles=reaction_smiles, evidence=evidence)
        return {
            "status": "ok",
            "tool": "fallback_candidate_retrieval",
            "summary": "Recovered fallback candidates from reaction diff evidence.",
            "candidate_reaction_types": [row["reaction_type"] for row in candidates],
            "candidate_scores": {row["reaction_type"]: row["deterministic_score"] for row in candidates},
            "retrieval_evidence": {
                "reacted_motifs": list((evidence.get("detection") or {}).get("reacted_motifs") or []),
                "formed_motifs": list((evidence.get("detection") or {}).get("formed_motifs") or []),
                "reaction_key": str((evidence.get("detection") or {}).get("reaction_key") or ""),
            },
            "recovered_candidates": candidates,
        }

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
