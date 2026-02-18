"""
Quality evaluation gate for the reaction SMILES pipeline.

Stage 3 of ReactionPipeline — evaluates the output of featurize_reaction()
against four chemistry-grounded criteria and decides whether the LLM fallback
(Stage 4) should be triggered.

All criteria are independently configurable via QualityConfig so the gate can
be tightened or relaxed per use case without touching pipeline logic.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Dict, FrozenSet, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Taxonomy ID loader (lazy, cached)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _load_valid_reaction_type_ids() -> FrozenSet[str]:
    """Load valid reaction type IDs from taxonomy JSON (cached after first call)."""
    try:
        from chemtools.taxonomy.loader import load_reaction_types_list
        entries = load_reaction_types_list()
        return frozenset(
            e["id"] for e in entries
            if isinstance(e, dict) and e.get("id")
        )
    except Exception:
        return frozenset()


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class QualityConfig:
    """
    Thresholds and flags that control the quality gate.

    All four criteria are applied independently; failing any one triggers
    the LLM fallback.

    Attributes:
        min_reaction_type_confidence: Minimum acceptable reaction_type_confidence
            from featurize_reaction(). Values below this indicate the taxonomy
            system is uncertain about the reaction class.
        require_reaction_key: If True, a missing or empty CRK-v1 reaction key
            triggers fallback. A missing key means product motif matching failed,
            making database lookup unreliable.
        required_reactant_motif_count: Exact number of distinct reactant motifs
            expected in reacted_motifs. For a typical 2-reactant coupling this
            must be 2 (one motif per reactant). Set to 0 to disable this check.
        allow_unclassified_reactant: If False (default), the presence of
            "Unclassified-Reactant" in reacted_motifs triggers fallback. This
            string is inserted by featurize_reaction()'s coverage guard when
            no taxonomy motif could be assigned to a reactant.
        require_known_reaction_type: If True (default), the reaction_type must
            be a valid ID from reaction_types.v4.0.json. Synthetic/fallback
            types (e.g. "Multi-Event:LGDisp+C-C") are not accepted and will
            trigger Stage 4. If the taxonomy cannot be loaded the check is
            skipped rather than failing.
        required_formed_motif_count: Exact number of formed (product) motifs
            expected. Default 1 — each reaction should produce exactly one
            product motif. Set to 0 to disable this check.
    """

    min_reaction_type_confidence: float = 0.5
    require_reaction_key: bool = True
    required_reactant_motif_count: int = 2
    allow_unclassified_reactant: bool = False
    require_known_reaction_type: bool = True
    required_formed_motif_count: int = 1


# ---------------------------------------------------------------------------
# Output dataclass
# ---------------------------------------------------------------------------

@dataclass
class QualityReport:
    """
    Result of the quality gate evaluation.

    Attributes:
        passed: True if all criteria passed (no LLM fallback needed).
        issues: Human-readable descriptions of each failing criterion.
        needs_llm_fallback: Convenience alias; equals ``not passed``.
        scores: Raw metric values for debugging / logging.
    """

    passed: bool
    issues: List[str]
    needs_llm_fallback: bool
    scores: Dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------------

class QualityEvaluator:
    """
    Evaluates featurize_reaction() output against the four quality criteria.

    Designed to be stateless — instantiate once per pipeline and call
    ``evaluate()`` for every reaction.

    Example:
        >>> config = QualityConfig(min_reaction_type_confidence=0.6)
        >>> evaluator = QualityEvaluator(config)
        >>> report = evaluator.evaluate(feat_result)
        >>> if report.needs_llm_fallback:
        ...     print("Triggering LLM fallback:", report.issues)
    """

    def __init__(self, config: Optional[QualityConfig] = None) -> None:
        self.config = config or QualityConfig()

    def evaluate(self, feat: "FeaturizationResult") -> QualityReport:  # noqa: F821
        """
        Run all quality checks against a FeaturizationResult.

        Args:
            feat: Output of ReactionPipeline.featurize()

        Returns:
            QualityReport indicating pass/fail and individual issue details
        """
        issues: List[str] = []
        cfg = self.config

        # Check 1 — featurization itself must have succeeded
        if not feat.success:
            issues.append(
                f"featurize_reaction() failed entirely: {feat.warnings}"
            )
            # No point checking further — all downstream criteria will also fail
            return _build_report(issues, feat)

        # Check 2 — reaction type confidence
        if feat.reaction_type_confidence < cfg.min_reaction_type_confidence:
            issues.append(
                f"Low reaction_type_confidence: {feat.reaction_type_confidence:.3f} "
                f"< threshold {cfg.min_reaction_type_confidence:.3f}"
                + (f" (detected: '{feat.reaction_type}')" if feat.reaction_type else "")
            )

        # Check 3 — missing or empty reaction key
        if cfg.require_reaction_key and not feat.reaction_key:
            issues.append(
                "Missing CRK-v1 reaction_key — product motif matching failed; "
                "database lookup will be unreliable"
            )

        # Check 4 — unclassified reactant
        if not cfg.allow_unclassified_reactant and feat.has_unclassified_reactant:
            issues.append(
                "'Unclassified-Reactant' found in reacted_motifs — at least one "
                "reactant has no taxonomy motif assignment"
            )

        # Check 5 — exact reactant motif count
        if cfg.required_reactant_motif_count > 0:
            actual = feat.reactant_motif_count
            expected = cfg.required_reactant_motif_count
            if actual != expected:
                issues.append(
                    f"Wrong reactant motif count: got {actual}, "
                    f"expected exactly {expected} "
                    f"(reacted_motifs={list(feat.reacted_motifs)})"
                )

        # Check 6 — reaction_type must be a known taxonomy ID, not a synthetic
        # fallback string (e.g. "Multi-Event:LGDisp+C-C")
        if cfg.require_known_reaction_type and feat.reaction_type:
            valid_ids = _load_valid_reaction_type_ids()
            if valid_ids and feat.reaction_type not in valid_ids:
                issues.append(
                    f"reaction_type '{feat.reaction_type}' is not a known taxonomy ID "
                    f"(not in reaction_types.v4.0.json). "
                    f"Synthetic/fallback types are rejected."
                )

        # Check 7 — exact formed (product) motif count
        if cfg.required_formed_motif_count > 0:
            actual_formed = len(feat.formed_motifs)
            if actual_formed != cfg.required_formed_motif_count:
                issues.append(
                    f"Wrong formed motif count: got {actual_formed}, "
                    f"expected exactly {cfg.required_formed_motif_count} "
                    f"(formed_motifs={list(feat.formed_motifs)})"
                )

        return _build_report(issues, feat)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_report(issues: List[str], feat: "FeaturizationResult") -> QualityReport:  # noqa: F821
    passed = len(issues) == 0
    valid_ids = _load_valid_reaction_type_ids()
    scores: Dict[str, Any] = {
        "reaction_type": feat.reaction_type,
        "reaction_type_confidence": feat.reaction_type_confidence,
        "reaction_type_known": (
            feat.reaction_type in valid_ids if (feat.reaction_type and valid_ids) else None
        ),
        "reaction_key_present": bool(feat.reaction_key),
        "has_unclassified_reactant": feat.has_unclassified_reactant,
        "reactant_motif_count": feat.reactant_motif_count,
        "reacted_motifs": list(feat.reacted_motifs),
        "formed_motifs": list(feat.formed_motifs),
        "formed_motif_count": len(feat.formed_motifs),
    }
    return QualityReport(
        passed=passed,
        issues=issues,
        needs_llm_fallback=not passed,
        scores=scores,
    )


__all__ = [
    "QualityConfig",
    "QualityReport",
    "QualityEvaluator",
]
