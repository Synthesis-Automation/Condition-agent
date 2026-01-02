"""
Planner-facing schemas and thin wrappers for the auto-conditions assistant.

This module keeps the deterministic spine lightweight: it normalizes inputs,
detects reaction families, pulls candidates from rules/precedents/HTE, and
returns structured payloads that higher-level agents or UIs can consume.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Tuple
import logging
from pathlib import Path

from pydantic import BaseModel, Field, field_validator

from chemtools import detect_reaction
from chemtools.taxonomy.rule_db import resolve_rule_db_v2
from chemtools.rule import RuleEngine
from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from chemtools.recommend import RecommendationResult

try:
    from chemtools.recommend import UnifiedRecommender, RecommendationResult

    UNIFIED_AVAILABLE = True
except Exception as exc:  # pragma: no cover - import-time guard
    UNIFIED_AVAILABLE = False
    UnifiedRecommender = None  # type: ignore
    RecommendationResult = Any  # type: ignore
    logger.debug("UnifiedRecommender unavailable: %s", exc)

try:
    from chemtools.HTE import HTEAnalytics

    HTE_AVAILABLE = True
except Exception as exc:  # pragma: no cover - import-time guard
    HTE_AVAILABLE = False
    HTEAnalytics = None  # type: ignore
    logger.debug("HTEAnalytics unavailable: %s", exc)

_REPO_ROOT = Path(__file__).resolve().parents[2]
_RULE_DB_SEARCH_PATHS: Tuple[Path, ...] = (
    _REPO_ROOT / "data" / "rule_db_v2",
    _REPO_ROOT.parent / "data" / "rule_db_v2",
    Path.cwd() / "data" / "rule_db_v2",
)

_FAMILY_TO_HTE_TYPE = {
    "cn_coupling": "C_N_Coupling",
    "c_n_coupling": "C_N_Coupling",
    "buchwald_cn": "C_N_Coupling",
    "buchwald_hartwig": "C_N_Coupling",
    "ullmann": "C_N_Coupling",
    "ullmann_cn": "C_N_Coupling",
    "c_n_cross_coupling": "C_N_Coupling",
    "sonogashira": "Sonogashira",
    "sonogashira_coupling": "Sonogashira",
    "suzuki": "Suzuki",
    "suzuki_miyaura": "Suzuki",
    "amide_formation": "amide_formation",
    "amide_coupling": "amide_formation",
    "amidation": "amide_formation",
    "snar": "SNAr",
    "s_nar": "SNAr",
    "snar_cn": "SNAr",
    "aromatic_nucleophilic_substitution": "SNAr",
    "c_o_coupling": "CO-Coupling",
    "co_coupling": "CO-Coupling",
    "rcm": "RCM",
}

_RULE_ENGINE_CACHE: Dict[Path, RuleEngine] = {}


class ReactionInput(BaseModel):
    """Normalized reaction payload that supports multiple input routes."""

    reaction_smiles: Optional[str] = Field(
        None, description="Full reaction SMILES (reactants>>products)."
    )
    substrates: List[str] = Field(
        default_factory=list,
        description="List of reactant SMILES strings if reaction_smiles is absent.",
    )
    products: List[str] = Field(
        default_factory=list,
        description="Optional product SMILES strings for context.",
    )
    text: Optional[str] = Field(
        None, description="Free-form description or protocol text."
    )
    context: Dict[str, Any] = Field(
        default_factory=dict, description="Arbitrary metadata (temperature, solvent)."
    )

    @field_validator("substrates", "products", mode="before")
    def _coerce_str_list(cls, value: Any) -> List[str]:
        """Ensure substrates/products are always lists of strings."""
        if value is None:
            return []
        if isinstance(value, str):
            return [value]
        return [str(item) for item in value]

    def as_reaction_smiles(self) -> Optional[str]:
        """Coerce inputs into a reaction SMILES if possible."""
        if self.reaction_smiles:
            return self.reaction_smiles
        if not self.substrates:
            return None
        reactants = ".".join(self.substrates)
        products = ".".join(self.products) if self.products else ""
        return f"{reactants}>>{products}"


class DetectedFamily(BaseModel):
    """Reaction family detection result."""

    family: Optional[str]
    confidence: Optional[float] = None
    method: Optional[str] = None
    raw: Dict[str, Any] = Field(default_factory=dict)


class CandidateCondition(BaseModel):
    """A candidate condition set from any source."""

    candidate_id: str
    components: Dict[str, Any]
    source: str
    raw_score: Optional[float] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)


class EvidenceScores(BaseModel):
    """Per-candidate evidence contributions."""

    candidate_id: str
    scores: Dict[str, float] = Field(default_factory=dict)
    notes: Dict[str, str] = Field(default_factory=dict)


class HteSummary(BaseModel):
    """Summarized HTE evidence for a reaction family."""

    reaction_type: Optional[str]
    highlights: Dict[str, Any] = Field(default_factory=dict)
    stats: Dict[str, Any] = Field(default_factory=dict)


class FusionResult(BaseModel):
    """Unified ranking across sources."""

    ranked: List[CandidateCondition]
    uncertainties: Dict[str, float] = Field(default_factory=dict)
    provenance: Dict[str, Any] = Field(default_factory=dict)


class ProtocolOutput(BaseModel):
    """Automation-friendly protocol output."""

    candidate_id: str
    additions: List[Dict[str, Any]] = Field(default_factory=list)
    notes: Optional[str] = None
    source: Optional[str] = None


class PlannerStep(BaseModel):
    """Single planner step indicating availability and intent."""

    name: str
    available: bool
    reason: Optional[str] = None


class PlannerPlan(BaseModel):
    """Ordered plan of tools the assistant intends to run."""

    steps: List[PlannerStep]


class AutoConditionsResult(BaseModel):
    """Container for the full deterministic auto-conditions run."""

    family: DetectedFamily
    plan: PlannerPlan
    rule_candidates: List[CandidateCondition] = Field(default_factory=list)
    protocol_candidates: List[CandidateCondition] = Field(default_factory=list)
    hte_summary: Optional[HteSummary] = None
    fused: FusionResult
    protocols: List[ProtocolOutput] = Field(default_factory=list)


def detect_family(reaction: ReactionInput) -> DetectedFamily:
    """Detect the reaction family using the deterministic router."""
    rxn = reaction.as_reaction_smiles()
    if not rxn:
        return DetectedFamily(family=None, method="missing_input", raw={})
    try:
        result = detect_reaction(rxn, use_ml=False)
        family = result.get("family")
        return DetectedFamily(
            family=family,
            confidence=result.get("confidence"),
            method=result.get("method"),
            raw=result,
        )
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning("Family detection failed: %s", exc)
        return DetectedFamily(family=None, method="error", raw={"error": str(exc)})


def _normalize_family_label(family: Optional[str]) -> Optional[str]:
    if not family:
        return None
    return str(family).strip().lower().replace("-", "_")


def _resolve_rule_database(family: Optional[str]) -> Optional[Path]:
    """Map a family hint to a concrete rule database path."""
    if not family:
        return None
    try:
        db_stem = resolve_rule_db_v2(family)
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug("Rule DB resolution failed: %s", exc)
        return None
    if not db_stem:
        return None
    for base in _RULE_DB_SEARCH_PATHS:
        candidate = base / f"{db_stem}.json"
        if candidate.exists():
            return candidate
    return None


def _get_rule_engine(db_path: Path) -> RuleEngine:
    """Return cached RuleEngine for the provided database path."""
    resolved = db_path.resolve()
    engine = _RULE_ENGINE_CACHE.get(resolved)
    if engine:
        return engine
    engine = RuleEngine.from_file(resolved)
    _RULE_ENGINE_CACHE[resolved] = engine
    return engine


def fetch_rule_candidates(
    reaction: ReactionInput, family: Optional[DetectedFamily] = None
) -> List[CandidateCondition]:
    """Generate rule-based candidates when a family-specific DB exists."""
    rxn = reaction.as_reaction_smiles()
    if not rxn:
        return []

    db_path = _resolve_rule_database(family.family if family else None)
    if not db_path:
        return []

    try:
        engine = _get_rule_engine(db_path)
        rec = engine.recommend(rxn)
    except Exception as exc:
        logger.warning("Rule recommendation failed: %s", exc)
        return []

    base_rule = rec.base_rule
    candidate = CandidateCondition(
        candidate_id=f"rule::{db_path.stem}",
        components=base_rule.conditions,
        source="rule",
        raw_score=base_rule.confidence,
        metadata={
            "matched_features": base_rule.matched_features,
            "database": db_path.name,
        },
    )
    return [candidate]


def _get_unified_recommender() -> Tuple[Optional[Any], Optional[str]]:
    if not UNIFIED_AVAILABLE:
        return None, "UnifiedRecommender not available"

    default_index = _REPO_ROOT / "build" / "unified_recommendation_index"
    if not default_index.exists():
        return None, f"Unified index missing at {default_index}"

    try:
        recommender = UnifiedRecommender(default_index)
        return recommender, None
    except Exception as exc:
        logger.warning("UnifiedRecommender init failed: %s", exc)
        return None, str(exc)


def _extract_unified_components(result: Any) -> Dict[str, Any]:
    if not hasattr(result, "full_data") or not result.full_data:
        return {}
    data = result.full_data
    if result.source_type == "dataset":
        return data.get("conditions", {}) or {}
    if result.source_type == "protocol":
        return data.get("conditions") or data.get("reaction_setup") or {}
    if result.source_type == "hte":
        return {
            "catalyst": data.get("Catalyst"),
            "ligand": data.get("Ligand"),
            "base": data.get("Base"),
            "solvent": data.get("Solvent"),
            "secondary_solvent": data.get("Secondary Solvent"),
            "additive": data.get("Additive"),
            "coupling_reagent": data.get("Coupling Reagent"),
        }
    return {}


def _convert_protocol_result(result: Any) -> CandidateCondition:
    components = _extract_unified_components(result)
    metadata: Dict[str, Any] = {"family": getattr(result, "family", None)}
    if hasattr(result, "full_data") and result.full_data:
        metadata["source_file"] = result.full_data.get("source_file")
    return CandidateCondition(
        candidate_id=f"{result.source_type}::{result.id}",
        components=components,
        source=result.source_type,
        raw_score=getattr(result, "similarity", None),
        metadata={
            **metadata,
            "name": getattr(result, "name", None),
            "rank": getattr(result, "rank", None),
        },
    )


def find_similar_protocols(
    reaction: ReactionInput, *, top_k: int = 5
) -> List[CandidateCondition]:
    """Retrieve unified recommendations (dataset/protocol/HTE) via DRFP similarity."""
    rxn = reaction.as_reaction_smiles()
    if not rxn:
        return []

    recommender, error = _get_unified_recommender()
    if not recommender or error:
        logger.info("Skipping DRFP similarity: %s", error)
        return []

    try:
        results: Sequence[RecommendationResult] = recommender.recommend(
            reaction_smiles=rxn, top_k=top_k, include_details=True
        )
    except Exception as exc:
        logger.warning("DRFP recommendation failed: %s", exc)
        return []

    converted: List[CandidateCondition] = []
    for result in results:
        converted.append(_convert_protocol_result(result))
    return converted


def fetch_hte_stats(
    reaction: ReactionInput, family: Optional[DetectedFamily] = None
) -> Optional[HteSummary]:
    """Summarize HTE signals for the detected family."""
    if not HTE_AVAILABLE:
        return None

    reaction_type = _FAMILY_TO_HTE_TYPE.get(
        _normalize_family_label(family.family if family else None), None
    )
    try:
        analytics = HTEAnalytics()
    except Exception as exc:
        logger.warning("HTE analytics unavailable: %s", exc)
        return None

    summary_df = analytics.get_reaction_type_summary()
    highlights: Dict[str, Any] = {}
    stats: Dict[str, Any] = {}

    if reaction_type:
        family_rows = summary_df[summary_df["Reaction_Type"] == reaction_type]
        if not family_rows.empty:
            top_row = family_rows.iloc[0].to_dict()
            highlights = {
                "reaction_type": reaction_type,
                "top_catalyst": top_row.get("Top_Catalyst"),
                "top_pair": top_row.get("Top_Reactant_Pair"),
            }
            stats = {k: v for k, v in top_row.items() if k != "Reaction_Type"}
    else:
        # Keep global context if no family match
        top_row = summary_df.iloc[0].to_dict()
        highlights = {
            "reaction_type": top_row.get("Reaction_Type"),
            "top_catalyst": top_row.get("Top_Catalyst"),
            "top_pair": top_row.get("Top_Reactant_Pair"),
        }
        stats = {k: v for k, v in top_row.items() if k != "Reaction_Type"}

    return HteSummary(reaction_type=reaction_type, highlights=highlights, stats=stats)


def score_ml_candidates(
    candidates: Sequence[CandidateCondition], reaction: ReactionInput
) -> List[EvidenceScores]:
    """
    Lightweight heuristic scorer placeholder for ML.

    This provides deterministic, fast scores so fusion can reason about ML
    signals even before a trained model is wired in.
    """
    scored: List[EvidenceScores] = []
    for cand in candidates:
        score = 0.5  # base prior
        components = cand.components or {}
        # Reward more complete condition sets
        for key in ("ligand", "base", "solvent", "catalyst"):
            if key in components and components[key]:
                score += 0.1
        score = min(score, 1.0)
        scored.append(
            EvidenceScores(
                candidate_id=cand.candidate_id,
                scores={"ml": score},
                notes={"status": "heuristic"},
            )
        )
    return scored


def fuse_scores(
    rule_candidates: Sequence[CandidateCondition],
    protocol_candidates: Sequence[CandidateCondition],
    hte_summary: Optional[HteSummary] = None,
    ml_scores: Optional[Sequence[EvidenceScores]] = None,
) -> FusionResult:
    """
    Simple fusion with optional ML weighting: rules prioritized, then protocols.
    """
    ml_map: Dict[str, float] = {}
    if ml_scores:
        for s in ml_scores:
            if "ml" in s.scores:
                ml_map[s.candidate_id] = s.scores["ml"]

    def _composite_score(cand: CandidateCondition) -> float:
        raw = cand.raw_score or 0.0
        ml = ml_map.get(cand.candidate_id, 0.0)
        source_boost = 1.5 if cand.source == "rule" else 1.0
        return source_boost * raw + ml

    combined = list(rule_candidates) + list(protocol_candidates)
    ordered = sorted(
        combined,
        key=lambda c: (_composite_score(c), 1 if c.source == "protocol" else 2),
        reverse=True,
    )

    provenance = {
        "sources": {
            "rule": len(rule_candidates),
            "protocol": len(protocol_candidates),
        },
        "composite_scores": {c.candidate_id: _composite_score(c) for c in ordered},
    }
    if hte_summary:
        provenance["hte"] = hte_summary.highlights
    return FusionResult(ranked=ordered, uncertainties={}, provenance=provenance)


def build_protocol(candidate: CandidateCondition) -> ProtocolOutput:
    """Convert a candidate into automation-friendly additions."""
    additions: List[Dict[str, Any]] = []
    notes: Optional[str] = None

    if candidate.source == "rule" and candidate.components:
        try:
            setup = rule_conditions_to_reaction_setup(
                candidate.components,
                user_substrates=None,
                scale_mmol=1.0,
                reaction_family=candidate.metadata.get("family"),
            )
            additions = setup.get("reaction_setup", []) or []
            notes = setup.get("notes")
        except Exception as exc:
            logger.warning("Rule-to-protocol formatting failed: %s", exc)

    if not additions and candidate.components:
        additions = [{"component": key, "value": value} for key, value in candidate.components.items()]

    return ProtocolOutput(
        candidate_id=candidate.candidate_id,
        additions=additions,
        notes=notes,
        source=candidate.source,
    )


def plan_workflow(reaction: ReactionInput) -> PlannerPlan:
    """Draft a minimal tool plan based on available sources."""
    family = detect_family(reaction)
    steps: List[PlannerStep] = []

    rule_available = _resolve_rule_database(family.family) is not None
    steps.append(
        PlannerStep(
            name="rule_candidates",
            available=rule_available,
            reason="rule DB found" if rule_available else "no rule DB for family",
        )
    )

    recommender, error = _get_unified_recommender()
    steps.append(
        PlannerStep(
            name="drfp_protocols",
            available=recommender is not None,
            reason="ready" if recommender else error or "unavailable",
        )
    )

    steps.append(
        PlannerStep(
            name="hte_stats",
            available=HTE_AVAILABLE,
            reason="available" if HTE_AVAILABLE else "HTE analytics missing",
        )
    )

    steps.append(
        PlannerStep(
            name="ml_score",
            available=True,
            reason="heuristic scorer available",
        )
    )

    return PlannerPlan(steps=steps)


def auto_conditions(
    reaction: ReactionInput,
    *,
    constraints: Optional[ConstraintSpec] = None,
    top_k_protocols: int = 5,
    build_protocols: bool = True,
    max_protocols: int = 3,
) -> AutoConditionsResult:
    """
    Deterministic end-to-end pipeline for auto-conditions recommendations.

    This chains detection → rules → DRFP precedents → HTE summary → fusion →
    protocol formatting. It is deliberately minimal and side-effect free.
    """
    family = detect_family(reaction)
    plan = plan_workflow(reaction)

    rule_candidates = fetch_rule_candidates(reaction, family)
    protocol_candidates = find_similar_protocols(reaction, top_k=top_k_protocols)
    hte_summary = fetch_hte_stats(reaction, family)
    
    # Apply constraints to candidates if provided
    if constraints:
        # Filter rule candidates
        rule_candidates = [
            c for c in rule_candidates 
            if _matches_constraints(c, constraints)
        ]
        # Filter protocol candidates
        protocol_candidates = [
            c for c in protocol_candidates 
            if _matches_constraints(c, constraints)
        ]

    ml_scores = score_ml_candidates(rule_candidates + protocol_candidates, reaction)

    fused = fuse_scores(rule_candidates, protocol_candidates, hte_summary, ml_scores)

    protocols: List[ProtocolOutput] = []
    if build_protocols:
        for cand in fused.ranked[:max_protocols]:
            protocols.append(build_protocol(cand))

    return AutoConditionsResult(
        family=family,
        plan=plan,
        rule_candidates=rule_candidates,
        protocol_candidates=protocol_candidates,
        hte_summary=hte_summary,
        fused=fused,
        protocols=protocols,
    )


def _matches_constraints(candidate: CandidateCondition, constraints: ConstraintSpec) -> bool:
    """Check if a candidate condition matches the provided constraints."""
    components = candidate.components
    catalyst = str(components.get("catalyst") or components.get("core") or "").upper()
    solvent = str(components.get("solvent") or "").upper()
    
    # Metal constraints
    if constraints.exclude_metals:
        if any(m.upper() in catalyst for m in constraints.exclude_metals):
            return False
    
    if constraints.allow_metals:
        if not any(m.upper() in catalyst for m in constraints.allow_metals):
            return False
            
    # Rule constraints (e.g., no_chlorinated)
    if constraints.constraint_rules:
        if constraints.constraint_rules.get("no_chlorinated"):
            chlorinated = {"DCM", "CHLOROFORM", "DCE", "CH2CL2", "CHCL3"}
            if any(c in solvent for c in chlorinated):
                return False
                
    return True
