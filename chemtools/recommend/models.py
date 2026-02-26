"""Typed request/plan/analysis models for recommendation orchestration."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple


class SourceGroup(str, Enum):
    """Canonical recommendation source groups."""

    ANY = "any"
    LITERATURE = "literature"
    PROTOCOLS = "protocols"
    EXPERIMENTS = "experiments"
    RULES = "rules"


class RunStrategy(str, Enum):
    """Execution strategy for recommendation requests."""

    SINGLE_PASS = "single_pass"
    PER_SOURCE = "per_source"
    ANALYSIS_ONLY = "analysis_only"


class OutputView(str, Enum):
    """Output layout/view semantics (separate from source selection)."""

    COMBINED = "combined"
    BY_SOURCE = "by_source"
    PRECEDENT_ONLY = "precedent_only"


class RecommendationStrategy(str, Enum):
    """Public recommendation strategy (separate from source selection)."""

    LITERATURE = "literature"
    SIMILARITY = "similarity"
    MOTIF = "motif"
    RULES = "rules"


def _coerce_enum(value: Any, enum_cls: type[Enum], default: Enum) -> Enum:
    if isinstance(value, enum_cls):
        return value
    text = str(value or "").strip().lower()
    if not text:
        return default
    aliases: Dict[str, str] = {
        "all": "any",
        "dataset": "literature",
        "datasets": "literature",
        "experiment": "experiments",
        "experiements": "experiments",
        "protocol": "literature",
        "protocols": "literature",
        "run_all": "per_source",
        "per_source_merge": "per_source",
        "merge": "per_source",
        "precedent": "precedent_only",
        "literature-based": "literature",
        "similarity-based": "similarity",
        "fingerprint": "similarity",
        "fingerprint-based": "similarity",
        "motif-based": "motif",
        "rule-based": "rules",
    }
    text = aliases.get(text, text)
    for member in enum_cls:  # type: ignore[assignment]
        if str(member.value).lower() == text:
            return member
    return default


@dataclass
class RecommendationRequest:
    """High-level recommendation request used by the facade API."""

    reaction_smiles: str
    strategy: Optional[RecommendationStrategy | str] = None
    source_group: SourceGroup | str = SourceGroup.ANY
    run_strategy: RunStrategy | str = RunStrategy.SINGLE_PASS
    output_view: OutputView | str = OutputView.COMBINED
    top_k: int = 10
    min_experiments: int = 2
    reaction_type_filter: Optional[str] = None
    catalyst_filter: Optional[str] = None
    hte_db_path: str = "data/HTE_db"
    reaction_key_only: bool = False
    use_aryl_steric_electronic_weighting: bool = False
    use_spectator_groups: bool = True
    prefer_mixfp_for_similarity: bool = True
    similarity_mixfp_weight: float = 0.75
    analysis_only: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)

    def normalized_source_group(self) -> SourceGroup:
        out = _coerce_enum(self.source_group, SourceGroup, SourceGroup.ANY)  # type: ignore[assignment]
        if out is SourceGroup.PROTOCOLS:
            return SourceGroup.LITERATURE
        return out  # type: ignore[return-value]

    def normalized_strategy(self) -> Optional[RecommendationStrategy]:
        text = str(self.strategy or "").strip()
        if not text:
            return None
        aliases = {
            "precedent": RecommendationStrategy.SIMILARITY.value,
            "precedent_only": RecommendationStrategy.SIMILARITY.value,
            "lit": RecommendationStrategy.LITERATURE.value,
            "experiments": RecommendationStrategy.MOTIF.value,
            "experiment": RecommendationStrategy.MOTIF.value,
            "experimental": RecommendationStrategy.MOTIF.value,
            "experiment-based": RecommendationStrategy.MOTIF.value,
            "rules-based": RecommendationStrategy.RULES.value,
            "similarity-based": RecommendationStrategy.SIMILARITY.value,
            "fingerprint-based": RecommendationStrategy.SIMILARITY.value,
        }
        text = aliases.get(text.lower(), text)
        out = _coerce_enum(text, RecommendationStrategy, RecommendationStrategy.MOTIF)
        return out  # type: ignore[return-value]

    def normalized_run_strategy(self) -> RunStrategy:
        if self.analysis_only:
            return RunStrategy.ANALYSIS_ONLY
        return _coerce_enum(self.run_strategy, RunStrategy, RunStrategy.SINGLE_PASS)  # type: ignore[return-value]

    def normalized_output_view(self) -> OutputView:
        return _coerce_enum(self.output_view, OutputView, OutputView.COMBINED)  # type: ignore[return-value]


@dataclass
class QueryAnalysis:
    """Featurizer-driven, dataset-independent query analysis summary."""

    reaction_smiles_input: str
    reaction_smiles_normalized: str = ""
    reactants: Tuple[str, ...] = ()
    agents: Tuple[str, ...] = ()
    products: Tuple[str, ...] = ()
    reactant_a_smiles: str = ""
    reactant_b_smiles: Optional[str] = None
    product_smiles: Optional[str] = None
    reaction_key: Optional[str] = None
    reacted_motifs: Tuple[str, ...] = ()
    formed_motifs: Tuple[str, ...] = ()
    spectator_motifs: Tuple[str, ...] = ()
    spectator_groups: Tuple[str, ...] = ()
    detected_reaction_type: Optional[str] = None
    detected_reaction_type_id: Optional[str] = None
    detected_reaction_type_name: Optional[str] = None
    detected_reaction_type_category: Optional[str] = None
    reaction_type_confidence: float = 0.0
    requested_reaction_type_filter: Optional[str] = None
    requested_reaction_type_filter_canonical: Optional[str] = None
    warnings: Tuple[str, ...] = ()
    raw_feature_summary: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SourcePlan:
    """Planner output separating source selection from execution behavior."""

    sources_to_run: Tuple[str, ...] = ()
    single_run_source_group: Optional[str] = None
    needs_hte_data: bool = False
    needs_precedent_data: bool = False
    run_strategy: RunStrategy = RunStrategy.SINGLE_PASS
    output_view: OutputView = OutputView.COMBINED
    notes: Tuple[str, ...] = ()
    recommendation_strategy: Optional[str] = None


@dataclass
class RecommendationRunResult:
    """Facade output bundling analysis, source plan, and recommendation result."""

    request: RecommendationRequest
    analysis: QueryAnalysis
    plan: SourcePlan
    recommendation: Any = None
    loaded_resources: Dict[str, Any] = field(default_factory=dict)
