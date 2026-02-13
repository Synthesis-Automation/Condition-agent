"""Canonical reaction evidence model for evidence-first agent workflows."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List


@dataclass
class ReactionEvidence:
    """Canonical evidence object shared across agent tools."""

    reaction_smiles: str
    diff: Dict[str, Any] = field(default_factory=dict)
    detection: Dict[str, Any] = field(default_factory=dict)
    taxonomy_candidates: List[Dict[str, Any]] = field(default_factory=list)
    provisional_decision: Dict[str, Any] = field(default_factory=dict)
    validation: Dict[str, Any] = field(default_factory=dict)
    final_decision: Dict[str, Any] = field(default_factory=dict)
    tool_artifacts: Dict[str, Any] = field(default_factory=dict)
    coverage_suggestions: List[Dict[str, Any]] = field(default_factory=list)
    analysis_snapshot: Dict[str, Any] = field(default_factory=dict)

    def has_diff(self) -> bool:
        """Return True when core transformation diff fields are present."""
        return bool(self.diff.get("principal_pair")) and "mcs_ratio" in self.diff

    def has_validation(self) -> bool:
        """Return True when validation has run."""
        return bool(self.validation)

    def decision_reaction_type(self) -> str:
        """Return best available reaction type from final/provisional decision."""
        final_rt = str((self.final_decision or {}).get("reaction_type") or "").strip()
        if final_rt:
            return final_rt
        return str((self.provisional_decision or {}).get("reaction_type") or "unknown")

    def merge_diff_payload(self, payload: Dict[str, Any]) -> None:
        """Merge payload returned by reaction_diff tool."""
        analysis = dict(payload.get("analysis") or {})
        if analysis:
            self.analysis_snapshot = analysis
            self.diff = {
                "principal_pair": analysis.get("principal_pair", {}),
                "mcs_smarts": analysis.get("mcs_smarts"),
                "mcs_atoms": analysis.get("mcs_atoms"),
                "mcs_ratio": analysis.get("mcs_ratio"),
                "core_formula_delta": analysis.get("core_formula_delta", {}),
                "side_formula_delta": analysis.get("side_formula_delta", {}),
                "core_reactant_formula": analysis.get("core_reactant_formula"),
                "core_product_formula": analysis.get("core_product_formula"),
                "side_reactant_formula": analysis.get("side_reactant_formula"),
                "side_product_formula": analysis.get("side_product_formula"),
            }
            self.detection = {
                "detection_error": analysis.get("detection_error"),
                "reacted_motifs": analysis.get("reacted_motifs", []),
                "formed_motifs": analysis.get("formed_motifs", []),
                "reaction_key": analysis.get("reaction_key", ""),
            }
            self.taxonomy_candidates = list(analysis.get("taxonomy_candidates") or [])
            self.provisional_decision = dict(analysis.get("decision") or {})
            if analysis.get("validation"):
                self.validation = dict(analysis.get("validation") or {})
        if payload.get("diff"):
            self.diff.update(dict(payload.get("diff") or {}))
        if payload.get("detection"):
            self.detection.update(dict(payload.get("detection") or {}))

    def to_analysis_view(self) -> Dict[str, Any]:
        """Build analysis-like dictionary for compatibility and reporting."""
        if self.analysis_snapshot:
            # Keep snapshot as baseline, then overlay evolving evidence.
            result = dict(self.analysis_snapshot)
        else:
            result = {"reaction_smiles": self.reaction_smiles}

        result["decision"] = dict(self.provisional_decision or result.get("decision") or {})
        result["taxonomy_candidates"] = list(self.taxonomy_candidates or result.get("taxonomy_candidates") or [])
        result["validation"] = dict(self.validation or result.get("validation") or {})
        result.setdefault("reaction_smiles", self.reaction_smiles)
        for key, value in self.diff.items():
            result[key] = value
        for key, value in self.detection.items():
            result[key] = value
        return result

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, payload: Dict[str, Any]) -> "ReactionEvidence":
        """Create evidence object from dictionary payload."""
        return cls(
            reaction_smiles=str(payload.get("reaction_smiles") or ""),
            diff=dict(payload.get("diff") or {}),
            detection=dict(payload.get("detection") or {}),
            taxonomy_candidates=list(payload.get("taxonomy_candidates") or []),
            provisional_decision=dict(payload.get("provisional_decision") or {}),
            validation=dict(payload.get("validation") or {}),
            final_decision=dict(payload.get("final_decision") or {}),
            tool_artifacts=dict(payload.get("tool_artifacts") or {}),
            coverage_suggestions=list(payload.get("coverage_suggestions") or []),
            analysis_snapshot=dict(payload.get("analysis_snapshot") or {}),
        )
