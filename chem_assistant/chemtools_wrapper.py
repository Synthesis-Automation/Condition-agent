"""
LangChain tool wrappers for ChemTools functions.

This module exposes existing chemtools functionality as LangChain tools
without modifying the original chemtools codebase.

Available Tools:
    - normalize_smiles_tool: Canonicalize SMILES strings
    - normalize_reaction_tool: Canonicalize reaction SMILES
    - detect_reaction_family_tool: Detect reaction family/type
    - classify_reactant_tool: Classify reactant type (aryl halide, amine, etc.)
    - get_functional_groups_tool: Detect functional groups in a molecule
    - calculable_features_tool: Evaluate curated calculable feature library for a molecule
    - molpipeline_featurize_tool: Generate molecular features with optional MolPipeline vectors
    - analyze_bond_changes_tool: Analyze bond breaking/formation in reactions
    - reaction_similarity_tool: Compare two reactions using DRFP similarity
    - recommend_conditions_tool: Get ML-based condition recommendations
    - rule_based_conditions_tool: Deterministic rule-engine condition guidance
    - search_precedents_tool: Search for similar precedent reactions
    - protocol_recommendation_tool: Find full experimental protocols from literature
    - unified_recommender_tool: Unified DRFP-based protocol + rule search
    - rule_builder_autofill_tool: LLM-assisted drafting of rule database JSON
    - hte_recommend_tool: HTE-based condition recommendations (66K experiments, catalyst filtering)
    - hte_analytics_tool: HTE database analytics (reactant pairs, catalysts, metal usage)
    - hte_conditions_tool: Detailed conditions for specific substrate pair (top-k conditions query)
    - hte_screening_set_tool: Generate diverse condition sets for HTE screening plates (up to 24)
    - reaction_dataset_analytics_tool: Analyze reaction dataset frequency/yield statistics
    - find_reagent_tool: Look up reagent information from database
    - reagent_database_analytics_tool: Summarize reagent registry statistics
    - list_supported_cores_tool: Enumerate catalyst cores observed in precedents
    - list_all_families_tool: List all available reaction families in dataset
    - add_reagent_tool: Insert or preview reagent taxonomy entries

Usage:
    from lang_chain.chemtools_wrapper import CHEMTOOLS_TOOLS
    from langgraph.prebuilt import create_agent

    agent = create_agent(llm, CHEMTOOLS_TOOLS)
"""

from typing import Dict, Any, List, Optional, Sequence, Tuple, Union, Literal
import json
import os
import time
from datetime import date
from dataclasses import dataclass
from collections import OrderedDict
from pathlib import Path
import pandas as pd
from pydantic import BaseModel, Field
from langchain_core.tools import tool
from llmtools.clients import LLMClient, RECOMMENDED_MODELS
from llmtools.prompts import RULE_BUILDER_EXTRACTION
from chem_assistant.planner import (
    ReactionInput,
    auto_conditions,
    detect_family as planner_detect_family,
    fetch_rule_candidates,
    find_similar_protocols,
    fetch_hte_stats,
    score_ml_candidates,
    fuse_scores,
    CandidateCondition,
)

# Import chemtools functions
from chemtools.smiles import normalize, normalize_reaction
from chemtools import (
    detect_reaction,
)
from chemtools.recommend.modules import recommend_from_reaction
from chemtools.recommend.utils import canonical_family
from chemtools.precedent import knn as precedent_knn
from chemtools.dataset_analytics import (
    get_dataset_stats,
    get_common_catalysts,
    get_common_ligands,
    get_common_bases,
    get_common_solvents,
    get_common_reagents,
    get_condition_cores,
    get_all_families,
)
from chemtools.reagent import (
    find_reagent,
    classify_reactant_smiles,
    add_reagent_entry,
    ReagentAdditionError,
)
from chemtools.reagent.analytics import (
    get_database_statistics,
    get_family_statistics,
    get_missing_data_report,
)
from chemtools.util.functional_groups import (
    detect_all as detect_functional_groups,
    summarize_functional_groups,
)
from chemtools.util.reaction_center_detector import identify_changed_atoms_from_mapped_smiles
from chemtools.featurizers import reaction_pair as molecular_featurizer
from chemtools.taxonomy import reaction_catalog as _reaction_catalog
from chemtools.taxonomy.rule_db import resolve_rule_db_v2

# Rule-based recommendation engine
from chemtools.rule import RuleEngine, RuleBuilder, FeatureAnalyzer

# Global analyzer instance for tools
_FEATURE_ANALYZER = FeatureAnalyzer()

# Import bond analysis tools (NEW)
from chemtools import (
    analyze_bond_changes,
    analyze_bond_changes_hybrid,
    rxnmapper_available,
)

# Import protocol recommendation (NEW)
try:
    from chemtools.protocol import ProtocolRecommender, recommend_protocol
    PROTOCOL_AVAILABLE = True
except ImportError:
    PROTOCOL_AVAILABLE = False
    ProtocolRecommender = None
    recommend_protocol = None

# Import reaction similarity (NEW)
try:
    from chemtools.reaction_similarity import compute_drfp_similarity
    DRFP_SIMILARITY_AVAILABLE = True
except ImportError:
    DRFP_SIMILARITY_AVAILABLE = False
    compute_drfp_similarity = None

# Import unified recommender (NEW)
try:
    from chemtools.recommend import UnifiedRecommender
    UNIFIED_RECOMMENDER_AVAILABLE = True
except ImportError:
    UNIFIED_RECOMMENDER_AVAILABLE = False
    UnifiedRecommender = None

# Import HTE recommendation and analytics (NEW)
try:
    from chemtools.HTE import HTERecommender, HTEAnalytics, format_result
    HTE_AVAILABLE = True
except ImportError:
    HTE_AVAILABLE = False
    HTERecommender = None
    HTEAnalytics = None
    format_result = None

REAGENT_RESOLVER_TIMEOUT = 6.0
RULE_BUILDER_SYSTEM_PROMPT = (
    "You are an expert synthetic chemistry knowledge engineer who converts "
    "protocol notes into structured rule databases. Output JSON only."
)

from .constraint_parser import (
    ConstraintSpec,
    build_constraint_spec,
    filter_cores_by_constraints,
    format_constraints_for_prompt,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent
_RULE_DB_SEARCH_PATHS: List[Path] = [
    _REPO_ROOT / "data" / "rule_db_v2",  # Updated to v2
    _REPO_ROOT / "data",
    _REPO_ROOT,
    Path.cwd() / "data" / "rule_db_v2",  # Updated to v2
    Path.cwd() / "data",
    Path.cwd(),
]

_FAMILY_TO_RULE_DB = {
    "cn_coupling": "C_N_Coupling_Cu_db",
    "c_n_coupling": "C_N_Coupling_Pd_db",
    "buchwald_cn": "C_N_Coupling_Pd_db",
    "buchwald_hartwig": "C_N_Coupling_Pd_db",
    "buchwald_hartwig_cn": "C_N_Coupling_Pd_db",
    "buchwald_hartwig_c_n": "C_N_Coupling_Pd_db",
    "buchwald-hartwig": "C_N_Coupling_Pd_db",
    "ullmann": "C_N_Coupling_Cu_db",
    "ullmann_cn": "C_N_Coupling_Cu_db",
    "suzuki": "Suzuki_db",
    "suzuki_coupling": "Suzuki_db",
    "suzuki_miyaura": "Suzuki_db",
    "suzuki-miyaura": "Suzuki_db",
    "c_n_cross_coupling": "C_N_Coupling_Pd_db",
    "snar_cn": "SNAr_db",
    "amide_coupling": "amide_formation_db",
    "amide_formation": "amide_formation_db",
    "amide_coupling": "amide_formation_db",
    "amidation": "amide_formation_db",
    "amide": "amide_formation_db",
    "snar": "SNAr_db",
    "s_nar": "SNAr_db",
    "aromatic_nucleophilic_substitution": "SNAr_db",
    "nucleophilic_aromatic_substitution": "SNAr_db",
    "reductive_amination": "reductive_amination_db",
    # New rule databases
    "sonogashira": "sonogashira_v2",
    "sonogashira_coupling": "sonogashira_v2",
    "c_o_coupling": "C_O_coupling_db",
    "co_coupling": "C_O_coupling_db",
    "c_o": "C_O_coupling_db",
    "rcm": "RCM_db",
    "ring_closing_metathesis": "RCM_db",
    "metathesis": "RCM_db",
}

_RULE_ENGINE_CACHE: "OrderedDict[Path, RuleEngine]" = OrderedDict()
_RULE_ENGINE_CACHE_SIZE = 4

class RuleDatabaseResolutionError(FileNotFoundError):
    """Raised when a rule database cannot be located from provided identifiers."""

    def __init__(
        self,
        attempted: List[str],
        reason: Optional[str] = None,
        detection: Optional[Dict[str, Any]] = None,
    ):
        attempted_display = attempted or ["<none>"]
        message = "Could not locate rule database for identifiers: " + ", ".join(attempted_display)
        if reason:
            message += f" ({reason})"
        super().__init__(message)
        self.attempted = attempted_display
        self.reason = reason
        self.detection = detection


# ============================================================================
# Pydantic Schemas
# ============================================================================


class NormalizeSmilesInput(BaseModel):
    """Schema for SMILES normalization."""

    smiles: str = Field(..., description="SMILES string to normalize.")


class NormalizeReactionInput(BaseModel):
    """Schema for reaction SMILES normalization."""

    reaction_smiles: str = Field(
        ...,
        description="Reaction SMILES in 'reactants>>products' or "
        "'reactants>agents>products' format.",
    )


class DetectReactionFamilyInput(BaseModel):
    """Schema for reaction family detection."""

    reaction_smiles: str = Field(..., description="Reaction SMILES to analyze.")


class ClassifyReactantInput(BaseModel):
    """Schema for reactant classification."""

    smiles: str = Field(..., description="Reactant SMILES string.")


class FunctionalGroupInput(BaseModel):
    """Schema for functional group detection."""

    smiles: str = Field(..., description="SMILES string to analyze for functional groups.")


class AnalyzeBondChangesInput(BaseModel):
    """Schema for bond breaking/formation analysis."""
    
    reaction_smiles: str = Field(
        ..., 
        description="Reaction SMILES string (reactants>>products) for bond analysis."
    )
    use_hybrid: bool = Field(
        True,
        description="Use hybrid approach (Manual + RXNMapper + MCS) for best results. "
                    "If False, uses RXNMapper only."
    )


class MolPipelineFeaturizeInput(BaseModel):
    """Schema for MolPipeline-backed molecular featurization."""

    electrophile: str = Field(
        ...,
        description="Electrophile SMILES (e.g., aryl halide or activated electrophile).",
    )
    nucleophile: str = Field(
        ...,
        description="Nucleophile SMILES (amine, alcohol, etc.).",
    )
    include_molpipeline: bool = Field(
        True,
        description="Include MolPipeline fingerprint/descriptor payload when available.",
    )


class CalculableFeaturesInput(BaseModel):
    """Schema for calculable feature detection."""

    smiles: str = Field(
        ...,
        description="SMILES string of the molecule to analyze.",
    )
    feature_tokens: Optional[List[str]] = Field(
        None,
        description="Optional list of feature tokens to extract. Defaults to all features.",
    )
    only_present: bool = Field(
        False,
        description="Return only features that evaluate to True or positive counts.",
    )
    include_summary: bool = Field(
        False,
        description="Include a human-readable summary of detected features.",
    )


class RecommendConditionsInput(BaseModel):
    """Schema for condition recommendation tool."""

    reaction_smiles: str = Field(
        ..., description="Reaction SMILES string (reactants>>products)."
    )
    k: int = Field(
        25,
        ge=1,
        le=200,
        description="Maximum precedents to retrieve for DRFP similarity search.",
    )
    max_variants: int = Field(
        3,
        ge=1,
        le=10,
        description="Maximum number of recommendation variants to generate.",
    )
    rerank_strategy: str = Field(
        "rule",
        description="Condition reranking strategy: 'rule', 'analytics', or 'none'.",
    )
    constraint_text: Optional[str] = Field(
        None,
        description="Natural language constraint hints (e.g., 'Pd-free, prefer Cu').",
    )
    allow_metals: Optional[List[str]] = Field(
        None,
        description="Whitelist of allowed catalyst metals (e.g., ['Cu', 'Ni']).",
    )
    exclude_metals: Optional[List[str]] = Field(
        None,
        description="Metals that must be excluded from recommendations.",
    )
    prefer_metals: Optional[List[str]] = Field(
        None,
        description="Metals to prioritize in the recommended catalyst cores.",
    )
    search_all_families: Optional[bool] = Field(
        None,
        description="Search for precedents across all reaction families.",
    )
    constraint_rules: Optional[Dict[str, Any]] = Field(
        None,
        description="Structured constraint overrides (e.g., {'no_chlorinated': True}).",
    )


class RuleRecommendInput(BaseModel):
    """Schema for deterministic rule-engine condition recommendation."""

    reaction_smiles: str = Field(
        ...,
        description="Reaction SMILES string (reactants>>products).",
    )
    database: Optional[str] = Field(
        None,
        description="Rule database name or JSON path (e.g., 'buchwald_cn' or 'data/rule_db/buchwald_cn.json').",
    )
    family_hint: Optional[str] = Field(
        None,
        description="Optional reaction family hint to help map to a rule database (e.g., 'Buchwald_CN').",
    )
    symptoms: Optional[Union[str, List[str]]] = Field(
        None,
        description="Observed symptoms or failure modes (list or comma-separated string).",
    )
    combine_method: Literal["union", "all", "first", "separate"] = Field(
        "union",
        description="How to combine features from multiple reactants when matching rules.",
    )
    include_reasoning: bool = Field(
        True,
        description="Include matched features and modifier triggers in the output.",
    )
    auto_detect: bool = Field(
        True,
        description="Infer the rule database via reaction-family detection when database is omitted.",
    )
    include_summary: bool = Field(
        True,
        description="Include a formatted text summary of the recommendation.",
    )


class RuleBuilderAutoInput(BaseModel):
    """Schema for LLM-assisted rule database drafting."""

    family: str = Field(
        ...,
        description="Reaction family/taxonomy identifier (e.g., 'Suzuki_Miyaura').",
    )
    metadata_id: str = Field(
        ...,
        description="Unique metadata id for the rule database (lower_snake_case recommended).",
    )
    metadata_name: str = Field(
        ...,
        description="Human-readable name for the rule database.",
    )
    metadata_version: str = Field(
        ...,
        description="Version string to embed in metadata (e.g., 'v1.0-draft').",
    )
    created_date: Optional[str] = Field(
        None,
        description="Creation date (YYYY-MM-DD). Defaults to today when omitted.",
    )
    status: Optional[str] = Field(
        "draft",
        description="Metadata status tag (draft, active, deprecated, etc.).",
    )
    tags: Optional[List[str]] = Field(
        None,
        description="Optional metadata tags (e.g., ['suzuki', 'hte']).",
    )
    reference_reactions: List[str] = Field(
        ...,
        min_length=1,
        description="Representative reaction SMILES strings (reactants>>products).",
    )
    protocol_text: str = Field(
        ...,
        description="Natural-language summary of protocols, screens, or precedent trends.",
    )
    notes: Optional[str] = Field(
        None,
        description="Optional override for reaction notes (LLM output used otherwise).",
    )
    desired_focus: Optional[str] = Field(
        None,
        description="Optional focus instructions (e.g., 'stress aryl chloride cases').",
    )
    applies_if_hints: Optional[List[str]] = Field(
        None,
        description="Feature tokens to force into applies_if.all (e.g., ['sp2_halide_present']).",
    )
    modifier_hints: Optional[List[str]] = Field(
        None,
        description="Symptom triggers the model should cover (e.g., ['symptom:hydrodehalogenation_observed']).",
    )
    max_base_rules: int = Field(
        4,
        ge=1,
        le=8,
        description="Maximum number of base rules to propose.",
    )


class ListSupportedCoresInput(BaseModel):
    """Schema for listing supported catalyst cores."""

    reaction_smiles: str = Field(
        ..., description="Reaction SMILES string for precedent lookup."
    )
    k: int = Field(
        25,
        ge=1,
        le=200,
        description="Number of precedents to inspect for catalyst cores.",
    )
    search_all_families: bool = Field(
        False,
        description="Whether to search across all reaction families.",
    )


class SearchPrecedentsInput(BaseModel):
    """Schema for precedent search."""

    reaction_smiles: str = Field(
        ..., description="Reaction SMILES string used to compute DRFP embeddings."
    )
    k: int = Field(
        10,
        ge=1,
        le=200,
        description="Number of most similar precedents to return.",
    )
    family: Optional[str] = Field(
        None,
        description="Optional reaction family override if auto-detection should be skipped.",
    )


class ReactionAnalyticsInput(BaseModel):
    """Schema for reaction dataset analytics access."""

    statistic: str = Field(
        "summary",
        description=(
            "Statistic to compute: 'summary', 'top_catalysts', 'top_ligands', "
            "'top_bases', 'top_solvents', 'top_condition_cores', 'top_reagents', or 'families'."
        ),
    )
    family: Optional[str] = Field(
        None,
        description="Reaction family identifier (e.g., 'Suzuki'). Required for most statistics.",
    )
    role: Optional[str] = Field(
        None,
        description="Reagent role for 'top_reagents' analytics (e.g., 'ligand', 'base').",
    )
    top_n: int = Field(
        10,
        ge=1,
        le=100,
        description="Maximum number of ranked entries to return for top-N statistics.",
    )
    min_yield: Optional[float] = Field(
        None,
        ge=0,
        le=100,
        description="Optional minimum yield threshold (%) when ranking reagents.",
    )


class ReagentAnalyticsInput(BaseModel):
    """Schema for requesting reagent database analytics."""

    statistic: str = Field(
        "summary",
        description="Statistic to compute: 'summary', 'role_summary', or 'missing_data'.",
    )
    role: Optional[str] = Field(
        None,
        description="Reagent role identifier required when statistic is 'role_summary' (e.g., 'solvent').",
    )
    registry_dir: Optional[str] = Field(
        None,
        description="Optional override path for the reagent registry directory.",
    )
    top_families: int = Field(
        20,
        ge=1,
        le=100,
        description="Maximum number of top families to include for summary statistics.",
    )


class FindReagentInput(BaseModel):
    """Schema for reagent lookup."""

    query: str = Field(..., description="Reagent name, abbreviation, CAS, or UID.")
    reagent_type: str = Field(
        "base",
        description="Preferred reagent collection to search (base, solvent, ligand, metal, additive).",
    )


class AddReagentInput(BaseModel):
    """Schema for adding reagents to the taxonomy."""

    cas: str = Field(..., description="CAS identifier for the reagent.")
    name: Optional[str] = Field(
        None, description="Preferred reagent name (resolved automatically if omitted)."
    )
    synonyms: Optional[Any] = Field(
        None,
        description="Additional synonyms (list/tuple/set or comma-delimited string).",
    )
    role: Optional[str] = Field(
        None, description="Explicit reagent role override (e.g., base, ligand)."
    )
    family_id: Optional[str] = Field(
        None, description="Explicit taxonomy family override if required."
    )
    abbreviation: Optional[str] = Field(
        None, description="Abbreviation override (defaults to resolved name)."
    )
    smiles: Optional[str] = Field(
        None, description="Optional SMILES annotation for the reagent."
    )
    taxonomy_dir: Optional[str] = Field(
        None,
        description="Optional override for writable taxonomy directory (falls back to defaults).",
    )
    allow_default_family: bool = Field(
        True,
        description="Allow fallback to default family if inference fails.",
    )
    dry_run: bool = Field(
        False,
        description="When true, validate the payload without writing to disk.",
    )
    auto_resolve: bool = Field(
        True,
        description="Resolve missing fields via CAS resolver when available.",
    )
    resolver_timeout: float = Field(
        REAGENT_RESOLVER_TIMEOUT,
        gt=0,
        description="Timeout (seconds) for CAS resolution attempts.",
    )


class ProtocolRecommendInput(BaseModel):
    """Schema for protocol-based recommendation."""

    reaction_smiles: str = Field(
        ...,
        description="Reaction SMILES string to match against protocol database."
    )
    k: int = Field(
        5,
        gt=0,
        le=20,
        description="Number of protocol recommendations to return (1-20)."
    )
    family_filter: Optional[str] = Field(
        None,
        description="Optional reaction family filter (e.g., 'Suzuki', 'Buchwald_CN', 'Amide_formation')."
    )
    use_smarts_filter: bool = Field(
        True,
        description="Enable SMARTS-based structural pre-filtering for better matches."
    )
    min_similarity: float = Field(
        0.3,
        ge=0.0,
        le=1.0,
        description="Minimum DRFP similarity threshold (0.0-1.0). Lower values return more results."
    )


class ReactionSimilarityInput(BaseModel):
    """Schema for computing reaction similarity."""

    reaction1_smiles: str = Field(
        ...,
        description="First reaction SMILES string."
    )
    reaction2_smiles: str = Field(
        ...,
        description="Second reaction SMILES string."
    )


class ListAllFamiliesInput(BaseModel):
    """Schema for listing all reaction families in the dataset."""
    pass  # No parameters needed


class HTERecommendInput(BaseModel):
    """Schema for HTE-based condition recommendation."""
    
    reactant_a_smiles: Optional[str] = Field(
        None,
        description="SMILES string of first reactant (e.g., aryl halide, aryl boronic acid). Use either this OR reaction_smiles."
    )
    reactant_b_smiles: Optional[str] = Field(
        None,
        description="SMILES string of second reactant (e.g., amine, boronic acid). Optional for some reactions."
    )
    reaction_smiles: Optional[str] = Field(
        None,
        description="Complete reaction SMILES (reactants>>products). If provided, reactants will be auto-extracted. Use this OR reactant_a_smiles."
    )
    top_k: int = Field(
        5,
        ge=1,
        le=20,
        description="Number of condition recommendations to return (1-20)."
    )
    min_experiments: int = Field(
        2,
        ge=1,
        le=100,
        description="Minimum number of experiments required for a condition to be recommended."
    )
    reaction_type_filter: Optional[str] = Field(
        None,
        description="Optional reaction type filter (e.g., 'Suzuki', 'C-N', 'Buchwald'). Auto-detected if omitted."
    )
    catalyst_filter: Optional[str] = Field(
        None,
        description="Optional catalyst metal filter (e.g., 'Pd', 'Cu', 'Ni', 'palladium', 'copper')."
    )


class HTEConditionsInput(BaseModel):
    """Schema for querying specific substrate pair conditions from HTE database."""
    
    reactant_a_type: str = Field(
        ...,
        description="Reactant A type (e.g., 'ArI', 'ArBr', 'ArCl'). Use exact type from database."
    )
    reactant_b_type: str = Field(
        ...,
        description="Reactant B type (e.g., 'Carbamate', 'RNH2', 'arom-NH', 'Lactam')."
    )
    reaction_type: Optional[str] = Field(
        None,
        description="Reaction type filter (e.g., 'C-N', 'Suzuki'). Optional."
    )
    catalyst_filter: Optional[str] = Field(
        None,
        description="Catalyst metal filter (e.g., 'Pd', 'Cu', 'Ni'). Optional."
    )
    top_k: int = Field(
        10,
        ge=1,
        le=50,
        description="Number of top conditions to return (1-50)."
    )
    min_experiments: int = Field(
        1,
        ge=1,
        le=100,
        description="Minimum experiments per condition (default 1)."
    )
    sort_by: str = Field(
        "count",
        description="Sort by 'count' (experiments), 'success' (success rate), or 'yield' (avg yield)."
    )


class HTEAnalyticsInput(BaseModel):
    """Schema for HTE database analytics queries."""
    
    query_type: str = Field(
        ...,
        description="Type of analytics query: 'list_pairs', 'catalysts', 'reactions', 'metals', or 'similar_pairs'."
    )
    reaction_type: Optional[str] = Field(
        None,
        description="Filter by reaction type (e.g., 'Suzuki', 'C-N'). Case-insensitive substring match."
    )
    catalyst_filter: Optional[str] = Field(
        None,
        description="Filter by catalyst metal (e.g., 'Pd', 'Cu', 'palladium', 'copper')."
    )
    reactant_a_type: Optional[str] = Field(
        None,
        description="Filter by reactant A type (e.g., 'ArBr', 'ArCl')."
    )
    reactant_b_type: Optional[str] = Field(
        None,
        description="Filter by reactant B type (e.g., 'RNH2', 'ArB(OH)2')."
    )
    min_experiments: int = Field(
        5,
        ge=1,
        le=1000,
        description="Minimum number of experiments for inclusion in results. Default 5 to capture sparse data."
    )
    top_n: int = Field(
        20,
        ge=1,
        le=100,
        description="Maximum number of results to return."
    )
    sort_by: str = Field(
        "count",
        description="Sort results by 'count' (number of experiments) or 'success_rate' (percentage yield > 50%)."
    )
    similarity_criteria: Optional[str] = Field(
        None,
        description="For 'similar_pairs' query: 'reaction_type', 'catalyst', or 'both'."
    )


# ============================================================================
# Recommendation Cache Utilities
# ============================================================================


@dataclass
class _RecommendationCacheEntry:
    """Stored recommendation result with metadata about request breadth."""

    k: int
    max_variants: int
    result: Dict[str, Any]


class _RecommendationCache:
    """Simple LRU cache for expensive recommend_from_reaction calls."""

    def __init__(self, max_entries: int = 8):
        self.max_entries = max_entries
        self._store: "OrderedDict[Tuple[Any, ...], _RecommendationCacheEntry]" = OrderedDict()

    def clear(self) -> None:
        """Remove all cached entries."""
        self._store.clear()

    def __len__(self) -> int:
        return len(self._store)

    def _make_key(
        self,
        normalized_reaction: str,
        rerank_strategy: str,
        spec: ConstraintSpec,
    ) -> Tuple[Any, ...]:
        allow_metals = tuple(sorted(spec.allow_metals))
        exclude_metals = tuple(sorted(spec.exclude_metals))
        prefer_metals = tuple(sorted(spec.prefer_metals))
        rules = spec.constraint_rules or {}
        try:
            rule_items = tuple(
                sorted(
                    (key, json.dumps(value, sort_keys=True, default=str))
                    for key, value in rules.items()
                )
            )
        except TypeError:
            rule_items = tuple(sorted((key, str(value)) for key, value in rules.items()))
        return (
            normalized_reaction,
            rerank_strategy.lower() if rerank_strategy else "",
            spec.search_all_families,
            allow_metals,
            exclude_metals,
            prefer_metals,
            rule_items,
        )

    def get_or_compute(
        self,
        *,
        reaction_smiles: str,
        normalized_reaction: str,
        k: int,
        max_variants: int,
        rerank_strategy: str,
        constraint_spec: ConstraintSpec,
    ) -> Tuple[Dict[str, Any], bool]:
        """Return cached recommendation or compute and store a new entry."""
        key = self._make_key(normalized_reaction, rerank_strategy, constraint_spec)
        entry = self._store.get(key)
        if entry and entry.k >= k and entry.max_variants >= max_variants:
            self._store.move_to_end(key)
            return entry.result, True

        result = recommend_from_reaction(
            reaction=reaction_smiles,
            k=k,
            max_variants=max_variants,
            rerank_strategy=rerank_strategy,
            search_all_families=constraint_spec.search_all_families,
            constraint_rules=constraint_spec.constraint_rules or None,
        )

        self._store[key] = _RecommendationCacheEntry(
            k=k,
            max_variants=max_variants,
            result=result,
        )
        if len(self._store) > self.max_entries:
            self._store.popitem(last=False)
        return result, False


def _normalized_reaction_key(reaction_smiles: str) -> str:
    """Return a deterministic key for reaction-based cache lookups."""
    try:
        normalized = normalize_reaction(reaction_smiles)
        if isinstance(normalized, dict):
            if normalized.get("normalized"):
                return normalized["normalized"]
            if normalized.get("input"):
                return str(normalized["input"]).strip()
        if isinstance(normalized, str):
            return normalized.strip()
    except Exception:
        pass
    return reaction_smiles.strip()


_recommendation_cache = _RecommendationCache(max_entries=10)


def clear_recommendation_cache() -> None:
    """Expose cache clearing for CLI/tests."""
    _recommendation_cache.clear()


def recommendation_cache_stats() -> Dict[str, Any]:
    """Return lightweight stats describing the recommendation cache."""
    entries: List[Dict[str, Any]] = []
    for key, entry in _recommendation_cache._store.items():
        normalized_reaction = key[0]
        if len(normalized_reaction) > 80:
            normalized_reaction = normalized_reaction[:77] + "..."
        entries.append({
            "reaction": normalized_reaction,
            "search_all_families": key[2],
            "allow_metals": list(key[3]),
            "exclude_metals": list(key[4]),
            "prefer_metals": list(key[5]),
            "constraint_rules": {k: v for k, v in key[6]},
            "k": entry.k,
            "max_variants": entry.max_variants,
        })
    return {"entries": len(entries), "items": entries}


# ============================================================================
# Rule-Based Recommendation Helpers
# ============================================================================

def _normalize_family_label(value: Optional[str]) -> Optional[str]:
    """Normalize family strings for lookup/mapping."""
    if value is None:
        return None
    text = str(value).strip().lower()
    if not text:
        return None
    if text.endswith(".json"):
        text = text[:-5]
    replacements = {
        "-": "_",
        " ": "_",
        "\t": "_",
        "/": "_",
        "\\": "_",
        ":": "_",
        ";": "_",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    while "__" in text:
        text = text.replace("__", "_")
    text = text.strip("_")
    return text or None


def _map_family_to_db_name(value: Optional[str]) -> Optional[str]:
    """Map a reaction family label to the configured rule database name."""
    if not value:
        return None

    # Prefer taxonomy metadata mapping (source of truth).
    try:
        mapped = resolve_rule_db_v2(str(value))
    except Exception:
        mapped = None
    if mapped:
        return mapped

    # Legacy fallback mapping (kept for non-taxonomy aliases).
    normalized = _normalize_family_label(value)
    if not normalized:
        return None
    return _FAMILY_TO_RULE_DB.get(normalized, normalized)


def _looks_like_path_identifier(identifier: str) -> bool:
    """Heuristically determine whether the identifier is a file path."""
    identifier = identifier.strip()
    if not identifier:
        return False
    return any(sep in identifier for sep in ("/", "\\")) or identifier.startswith(".")


def _expand_with_extension(path: Path) -> List[Path]:
    """Return candidate paths including optional .json extension."""
    if path.suffix:
        return [path]
    return [path, path.with_suffix(".json")]


def _generate_identifier_variants(identifier: Optional[str]) -> List[str]:
    """Generate lookup variants for a database identifier."""
    if not identifier:
        return []
    ident = str(identifier).strip()
    if not ident or ident.lower() == "auto":
        return []
    variants: List[str] = []
    is_path = _looks_like_path_identifier(ident)
    if is_path:
        variants.append(ident)
        return variants
    if ident.lower().endswith(".json"):
        stem = Path(ident).stem
        variants.extend([ident, stem])
    else:
        variants.append(ident)
    normalized = _normalize_family_label(ident)
    if normalized and normalized not in variants:
        variants.append(normalized)
    mapped = _map_family_to_db_name(ident)
    if mapped and mapped not in variants:
        variants.append(mapped)
    if normalized:
        dashed = normalized.replace("_", "-")
        if dashed not in variants:
            variants.append(dashed)
    return variants


def _resolve_candidate_to_path(candidate: str) -> Optional[Path]:
    """Resolve a single candidate identifier to an existing JSON path."""
    candidate = candidate.strip()
    if not candidate:
        return None
    path_candidate = Path(candidate)
    targets: List[Path] = []
    if path_candidate.is_absolute() or _looks_like_path_identifier(candidate):
        targets.extend(_expand_with_extension(path_candidate if path_candidate.is_absolute() else path_candidate))
        if not path_candidate.is_absolute():
            targets.extend(
                _expand_with_extension(_REPO_ROOT / path_candidate)
            )
    else:
        names = [candidate]
        normalized = _normalize_family_label(candidate)
        if normalized and normalized not in names:
            names.append(normalized)
        mapped = _map_family_to_db_name(candidate)
        if mapped and mapped not in names:
            names.append(mapped)
        if candidate.lower().endswith(".json"):
            stem = Path(candidate).stem
            if stem and stem not in names:
                names.append(stem)
        for name in dict.fromkeys(names):
            for base in _RULE_DB_SEARCH_PATHS:
                targets.extend(_expand_with_extension(base / name))
    dedup: List[Path] = []
    seen: set[str] = set()
    for target in targets:
        key = str(target)
        if key in seen:
            continue
        seen.add(key)
        dedup.append(target)
    for target in dedup:
        try:
            if target.exists():
                return target.resolve()
        except OSError:
            continue
    return None


def _resolve_rule_database(*identifiers: Optional[str]) -> Tuple[Path, str, List[str]]:
    """Resolve the rule database path from one or more identifiers."""
    attempted: List[str] = []
    for identifier in identifiers:
        for variant in _generate_identifier_variants(identifier):
            if not variant or variant in attempted:
                continue
            attempted.append(variant)
            resolved = _resolve_candidate_to_path(variant)
            if resolved:
                return resolved, resolved.stem, attempted
    raise RuleDatabaseResolutionError(attempted)


def _select_rule_database(
    reaction_smiles: str,
    *,
    database: Optional[str],
    family_hint: Optional[str],
    auto_detect: bool,
) -> Tuple[Path, str, Optional[Dict[str, Any]], List[str]]:
    """
    Determine which rule database to use and return resolution metadata.
    
    Returns:
        (database_path, database_name, detection_payload, attempted_identifiers)
    """
    identifiers: List[str] = []
    attempted: List[str] = []
    detection_payload: Optional[Dict[str, Any]] = None
    database_clean = database.strip() if isinstance(database, str) else None

    if database_clean and database_clean.lower() != "auto":
        identifiers.append(database_clean)

    if family_hint:
        identifiers.append(str(family_hint))

    detection_error: Optional[str] = None
    if auto_detect:
        try:
            detection_payload = detect_reaction(reaction_smiles, use_ml=False)
            detected_family = detection_payload.get("family")
            if detected_family:
                identifiers.append(str(detected_family))
        except Exception as exc:
            detection_payload = {"error": str(exc)}
            detection_error = str(exc)

    if not identifiers:
        raise RuleDatabaseResolutionError(
            [],
            reason=detection_error or "no identifiers provided",
            detection=detection_payload,
        )

    try:
        db_path, db_name, attempted = _resolve_rule_database(*identifiers)
        return db_path, db_name, detection_payload, attempted
    except RuleDatabaseResolutionError as exc:
        reason = detection_error or getattr(exc, "reason", None)
        raise RuleDatabaseResolutionError(
            exc.attempted,
            reason=reason,
            detection=detection_payload or getattr(exc, "detection", None),
        ) from exc


def _get_rule_engine(db_path: Path) -> RuleEngine:
    """Return a cached RuleEngine instance for the given database path."""
    resolved = db_path.resolve()
    engine = _RULE_ENGINE_CACHE.get(resolved)
    if engine is not None:
        _RULE_ENGINE_CACHE.move_to_end(resolved)
        return engine

    engine = RuleEngine.from_file(resolved)
    _RULE_ENGINE_CACHE[resolved] = engine
    if len(_RULE_ENGINE_CACHE) > _RULE_ENGINE_CACHE_SIZE:
        _RULE_ENGINE_CACHE.popitem(last=False)
    return engine


def _normalize_symptom_list(value: Optional[Union[str, List[str]]]) -> Optional[List[str]]:
    """Coerce symptom inputs into a normalized list of strings."""
    if value is None:
        return None
    if isinstance(value, str):
        parts = [part.strip() for part in value.split(",") if part.strip()]
        return parts or None
    if isinstance(value, (list, tuple, set)):
        cleaned = [str(item).strip() for item in value if str(item).strip()]
        return cleaned or None
    coerced = str(value).strip()
    return [coerced] if coerced else None


# ============================================================================
# SMILES Normalization Tools
# ============================================================================

@tool(args_schema=NormalizeSmilesInput)
def normalize_smiles_tool(smiles: str) -> Dict[str, Any]:
    """
    Normalize (canonicalize) a SMILES string.
    
    This standardizes molecular representations for consistent comparisons.
    Returns a dictionary containing the canonical form, fragment breakdown,
    and other normalization metadata.
    
    Args:
        smiles: SMILES string to normalize (e.g., "CCO", "c1ccccc1")
    
    Returns:
        Dict[str, Any]: Structured normalization payload with success flag.
    
    Example:
        Input: "c1ccccc1"
        Output: {"success": True, "smiles_norm": "c1ccccc1", ...}
    """
    try:
        result = normalize(smiles)
        if isinstance(result, dict):
            if result.get("error"):
                return _error_response(
                    str(result.get("error")),
                    {"details": result, "smiles": smiles},
                )
            return _success_response(result)
        if result:
            return _success_response({"smiles_norm": result})
        return _error_response(f"Invalid SMILES '{smiles}'", {"smiles": smiles})
    except Exception as e:
        return _error_response(str(e), {"smiles": smiles})


@tool(args_schema=NormalizeReactionInput)
def normalize_reaction_tool(reaction_smiles: str) -> Dict[str, Any]:
    """
    Normalize a reaction SMILES string.
    
    Standardizes both reactants and products, sorting components for consistency.
    Accepts formats: "reactants>>products" or "reactants>agents>products"
    
    Args:
        reaction_smiles: Reaction SMILES (e.g., "CCBr.CCO>>CCOCC")
    
    Returns:
        Dict[str, Any]: Structured normalization payload with success flag.
    
    Example:
        Input: "c1ccccc1Br.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        Output: Normalized reaction SMILES
    """
    try:
        result = normalize_reaction(reaction_smiles)
        if not result:
            return _error_response(
                f"Invalid reaction SMILES '{reaction_smiles}'",
                {"reaction_smiles": reaction_smiles},
            )
        if isinstance(result, dict):
            return _success_response(result)
        return _success_response({"normalized": str(result)})
    except Exception as e:
        return _error_response(str(e), {"reaction_smiles": reaction_smiles})


# ============================================================================
# Reaction Analysis Tools
# ============================================================================

@tool(args_schema=DetectReactionFamilyInput)
def detect_reaction_family_tool(reaction_smiles: str) -> Dict[str, Any]:
    """
    Detect the reaction family/type from a reaction SMILES.
    
    Analyzes reactants and products to identify the reaction type
    (e.g., Suzuki, Buchwald_CN, Ullmann_CN, etc.)
    
    Args:
        reaction_smiles: Reaction SMILES string
    
    Returns:
        Dict[str, Any]: Family detection metadata with success flag.
    
    Example:
        Input: "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
        Output: {"family": "suzuki_miyaura", "confidence": 0.9}
    """
    try:
        result = detect_reaction(reaction_smiles, use_ml=False)
        detected = result.get("family")
        payload = {
            "family": detected,
            "dataset_family": canonical_family(detected),
            "confidence": result.get("confidence"),
            "method": result.get("method"),
        }
        if not payload.get("family"):
            return _error_response(
                "Could not determine reaction family",
                {"details": result, "reaction_smiles": reaction_smiles},
            )
        return _success_response(payload)
    except Exception as e:
        return _error_response(str(e), {"family": "Unknown", "reaction_smiles": reaction_smiles})


@tool(args_schema=ClassifyReactantInput)
def classify_reactant_tool(smiles: str) -> Dict[str, Any]:
    """
    Classify a reactant molecule type.
    
    Identifies reactant categories like aryl_halide, amine, boronic_acid,
    alkyne, etc. based on structural patterns.
    
    Args:
        smiles: SMILES string of the reactant
    
    Returns:
        Dict[str, Any]: Classification payload with success flag.
    
    Example:
        Input: "Brc1ccccc1"
        Output: {"category": "aryl_halide", "matches": [...]}
    """
    try:
        result = classify_reactant_smiles(smiles)
        if isinstance(result, dict):
            return _success_response(result)
        return _success_response({"classification": result})
    except Exception as e:
        return _error_response(str(e), {"category": "unknown", "smiles": smiles})


@tool(args_schema=FunctionalGroupInput)
def get_functional_groups_tool(smiles: str) -> Dict[str, Any]:
    """
    Detect functional groups in a molecule.
    
    Identifies 80+ functional groups using SMARTS patterns.
    Uses RDKit when available, falls back to text patterns.
    
    Args:
        smiles: SMILES string to analyze
    
    Returns:
        Dict[str, Any]: Functional group boolean map with success flag.
    
    Example:
        Input: "CCO"
        Output: {"alcohol": true, "primary_alcohol": true}
    """
    try:
        result = detect_functional_groups(smiles)
        if isinstance(result, dict):
            return _success_response(result)
        return _success_response({"functional_groups": result})
    except Exception as e:
        return _error_response(str(e), {"smiles": smiles})


@tool(args_schema=FunctionalGroupInput)
def summarize_functional_groups_tool(smiles: str) -> Dict[str, Any]:
    """
    Provide a human-readable summary of functional groups in a molecule.
    
    Categorizes detected groups (e.g., Oxygen-bearing, Nitrogen-bearing, etc.)
    and returns a formatted string summary.
    
    Args:
        smiles: SMILES string to analyze
    
    Returns:
        Dict[str, Any]: Formatted summary string with success flag.
    """
    try:
        summary = summarize_functional_groups(smiles)
        return _success_response({"summary": summary})
    except Exception as e:
        return _error_response(str(e), {"smiles": smiles})


class ReactionCenterInput(BaseModel):
    """Schema for reaction center identification."""
    mapped_smiles: str = Field(..., description="Atom-mapped reaction SMILES (e.g., [C:1]#[C:2].[c:3][I:4]>>[C:1]#[C:2]-[c:3])")


@tool(args_schema=ReactionCenterInput)
def identify_reaction_center_tool(mapped_smiles: str) -> Dict[str, Any]:
    """
    Identify the reaction center from an atom-mapped reaction SMILES.
    
    Detects which atoms and bonds are broken or formed during the reaction.
    This helps focus analysis on the actual transformation rather than spectator groups.
    
    Args:
        mapped_smiles: Reaction SMILES with atom mapping
    
    Returns:
        Dict[str, Any]: Lists of changed atoms, broken bonds, and formed bonds.
    """
    try:
        result = identify_changed_atoms_from_mapped_smiles(mapped_smiles)
        if "error" in result:
            return _error_response(result["error"])
        return _success_response(result)
    except Exception as e:
        return _error_response(str(e))


@tool(args_schema=MolPipelineFeaturizeInput)
def molpipeline_featurize_tool(
    electrophile: str,
    nucleophile: str,
    include_molpipeline: bool = True,
) -> Dict[str, Any]:
    """
    Featurize an electrophile/nucleophile pair with optional MolPipeline vectors.

    Returns the deterministic ChemTools substrate features (LG, bin, nucleophile class)
    and, when available and requested, MolPipeline-backed Morgan fingerprints plus
    RDKit phys-chem descriptors for each reactant. The MolPipeline payload exposes both
    list and name-value map (`physchem_map`) so properties such as MolLogP or TPSA can
    be accessed directly.

    Args:
        electrophile: Electrophile SMILES (aryl halide, activated electrophile, etc.)
        nucleophile: Nucleophile SMILES (amine, alcohol, etc.)
        include_molpipeline: Include MolPipeline fingerprint/descriptor payload.

    Returns:
        Dict[str, Any]: Feature dictionary wrapped in a success payload.
    """
    try:
        feat_result = molecular_featurizer.featurize_pair(
            electrophile,
            nucleophile,
            include_molpipeline=include_molpipeline,
        )
        features = feat_result.get("flat", {})
        return _success_response(features)
    except Exception as e:
        return _error_response(
            str(e),
            {
                "electrophile": electrophile,
                "nucleophile": nucleophile,
                "include_molpipeline": include_molpipeline,
            },
        )


@tool(args_schema=CalculableFeaturesInput)
def calculable_features_tool(
    smiles: str,
    feature_tokens: Optional[List[str]] = None,
    only_present: bool = False,
    include_summary: bool = False,
) -> Dict[str, Any]:
    """
    Detect calculable features for a molecule using the curated feature library.

    Features cover SMARTS-derived motifs, integer counts, and heuristic properties
    defined in the V2 taxonomy. When ``feature_tokens``
    is supplied the output is restricted to those tokens (missing entries are
    reported separately). Setting ``only_present`` filters the result to features
    that evaluate to ``True`` or positive integers.
    """
    try:
        all_features = _FEATURE_ANALYZER.analyze_reactant(smiles)
        filtered = dict(all_features)
        missing: List[str] = []
        if feature_tokens:
            filtered = {}
            for token in feature_tokens:
                if token in all_features:
                    filtered[token] = all_features[token]
                else:
                    missing.append(token)

        if only_present:
            filtered = {
                token: value
                for token, value in filtered.items()
                if (isinstance(value, bool) and value)
                or (isinstance(value, (int, float)) and value > 0)
            }

        present_tokens = [k for k, v in all_features.items() if v is True]

        steric_tokens = {
            "ortho_hindered": bool(all_features.get("ortho_hindered")),
            "ortho_very_hindered": bool(all_features.get("ortho_very_hindered")),
            "ortho_1": bool(all_features.get("ortho_1")),
        }
        # Lightweight textual fallback when RDKit unavailable
        if not steric_tokens["ortho_1"] and "C(C)C" in smiles:
            steric_tokens["ortho_1"] = True
        
        steric_indicators = [
            token for token, present in steric_tokens.items() if present
        ]
        steric_level = "low"
        if steric_tokens["ortho_very_hindered"]:
            steric_level = "very high"
        elif steric_tokens["ortho_hindered"]:
            steric_level = "high"
        elif steric_indicators:
            steric_level = "moderate"

        payload: Dict[str, Any] = {
            "features": filtered,
            "present_features": present_tokens,
        }
        if missing:
            payload["missing_tokens"] = missing
        if include_summary:
            payload["summary"] = ", ".join(present_tokens)
        payload["steric_hindrance"] = {
            "level": steric_level,
            "indicators": steric_indicators,
        }
        return _success_response(payload)
    except Exception as e:
        return _error_response(
            str(e),
            {
                "smiles": smiles,
                "feature_tokens": feature_tokens,
                "only_present": only_present,
            },
        )


@tool(args_schema=AnalyzeBondChangesInput)
def analyze_bond_changes_tool(reaction_smiles: str, use_hybrid: bool = True) -> Dict[str, Any]:
    """
    Analyze which bonds break and form in a chemical reaction.
    
    Uses a hybrid approach combining:
    - Manual atom mapping (if present in SMILES) - ground truth
    - RXNMapper (ML-based automatic atom mapping) - precise
    - MCS (Maximum Common Substructure) - validation fallback
    
    Detects leaving groups (like Br, I) and joining groups accurately.
    
    Args:
        reaction_smiles: Reaction SMILES (reactants>>products)
        use_hybrid: Use hybrid multi-method approach (recommended)
    
    Returns:
        Dict[str, Any]: Bond analysis with broken_bonds, formed_bonds, 
                       leaving_groups, confidence, and agreement between methods.
    
    Example:
        Input: "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
        Output: {
            "success": true,
            "method": "hybrid",
            "combined_confidence": 0.75,
            "broken_bonds": [(5, "Br (leaving group)"), (4, "B (leaving group)")],
            "formed_bonds": [(4, 5)],
            "agreement": {"rxnmapper_vs_mcs": true},
            "validation": "Both methods agree"
        }
    """
    try:
        # Check if RXNMapper is available
        if not rxnmapper_available():
            return _error_response(
                "RXNMapper not available. Install with: pip install rxnmapper",
                {"recommendation": "Install rxnmapper for accurate bond analysis"}
            )
        
        # Use hybrid approach or single method
        if use_hybrid:
            result = analyze_bond_changes_hybrid(
                reaction_smiles,
                use_rxnmapper=True,
                use_mcs=True,
                use_manual=True,
                auto_map=True
            )
        else:
            result = analyze_bond_changes(reaction_smiles, auto_map=True)
        
        if result.get('success'):
            # Format the result for better readability
            formatted = {
                'success': True,
                'method': result.get('method', 'unknown'),
                'combined_confidence': result.get('combined_confidence', result.get('mapping_confidence', 0.0)),
            }
            
            # Add recommended result if hybrid
            if use_hybrid and 'recommended_result' in result:
                rec = result['recommended_result']
                formatted['broken_bonds'] = rec.get('broken_bonds', [])
                formatted['formed_bonds'] = rec.get('formed_bonds', [])
                formatted['changed_atoms'] = list(rec.get('changed_atoms', []))
                formatted['leaving_groups'] = rec.get('leaving_groups', [])
                formatted['joining_groups'] = rec.get('joining_groups', [])
                formatted['agreement'] = result.get('agreement', {})
                formatted['validation'] = result.get('validation', 'N/A')
                
                # Add individual method results for transparency
                if result.get('manual_result'):
                    formatted['manual_detected'] = True
                if result.get('rxnmapper_result'):
                    formatted['rxnmapper_confidence'] = result['rxnmapper_result'].get('mapping_confidence')
                if result.get('mcs_result'):
                    formatted['mcs_coverage'] = result['mcs_result'].get('mcs_coverage')
            else:
                # Single method result
                formatted['broken_bonds'] = result.get('broken_bonds', [])
                formatted['formed_bonds'] = result.get('formed_bonds', [])
                formatted['changed_atoms'] = list(result.get('changed_atoms', []))
                formatted['leaving_groups'] = result.get('leaving_groups', [])
                formatted['joining_groups'] = result.get('joining_groups', [])
            
            # Add interpretation
            broken_count = len(formatted.get('broken_bonds', []))
            formed_count = len(formatted.get('formed_bonds', []))
            
            if broken_count > 0 and formed_count > 0:
                formatted['interpretation'] = f"Substitution/coupling: {broken_count} bond(s) break, {formed_count} bond(s) form"
            elif broken_count > formed_count:
                formatted['interpretation'] = "Bond-breaking dominant (fragmentation/elimination)"
            elif formed_count > broken_count:
                formatted['interpretation'] = "Bond-forming dominant (addition/condensation)"
            else:
                formatted['interpretation'] = "Balanced or rearrangement"
            
            return formatted
        else:
            return _error_response(
                result.get('error', 'Bond analysis failed'),
                {'reaction_smiles': reaction_smiles}
            )
            
    except Exception as e:
        return _error_response(str(e), {'reaction_smiles': reaction_smiles})

@tool(args_schema=RecommendConditionsInput)
def recommend_conditions_tool(
    reaction_smiles: str,
    k: int = 25,
    max_variants: int = 3,
    rerank_strategy: str = "rule",
    constraint_text: Optional[str] = None,
    allow_metals: Optional[List[str]] = None,
    exclude_metals: Optional[List[str]] = None,
    prefer_metals: Optional[List[str]] = None,
    search_all_families: Optional[bool] = None,
    constraint_rules: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Recommend reaction conditions using ML-based DRFP similarity search.
    
    Analyzes the reaction and finds the most similar precedent reactions,
    then recommends optimal conditions (catalyst, base, solvent, temperature, time).
    
    Args:
        reaction_smiles: Reaction SMILES string (reactants>>products)
        k: Number of precedents to retrieve (default: 25)
        max_variants: Maximum condition variants to generate (default: 3)
        rerank_strategy: Ranking strategy - "rule", "analytics", or "none" (default: "rule")
    
    Returns:
        Dict[str, Any]: Recommendation payload with success flag.
    
    Example:
        Input: "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        Output: {
            "family": "Buchwald_CN",
            "recommendation": {
                "core": "Pd/XPhos",
                "base": "Cs2CO3",
                "solvent": "1,4-dioxane",
                "T_C": 100,
                "time_h": 12
            },
            "alternatives": {...},
            "reasons": [...]
        }
    """
    try:
        constraint_spec = build_constraint_spec(
            text=constraint_text,
            allow_metals=allow_metals,
            exclude_metals=exclude_metals,
            prefer_metals=prefer_metals,
            search_all_families=search_all_families,
            base_constraint_rules=constraint_rules,
        )

        normalized_key = _normalized_reaction_key(reaction_smiles)
        start = time.perf_counter()
        result, cache_hit = _recommendation_cache.get_or_compute(
            reaction_smiles=reaction_smiles,
            normalized_reaction=normalized_key,
            k=k,
            max_variants=max_variants,
            rerank_strategy=rerank_strategy,
            constraint_spec=constraint_spec,
        )
        duration_ms = (time.perf_counter() - start) * 1000

        simplified = _simplify_recommendation(result)
        simplified, notes = _apply_core_constraints(simplified, constraint_spec)
        if notes:
            simplified["constraint_notes"] = notes

        prompt_hint = format_constraints_for_prompt(constraint_spec)
        if prompt_hint:
            simplified["constraint_summary"] = prompt_hint
        simplified["cache_hit"] = cache_hit
        simplified["timing_ms"] = round(duration_ms, 2)

        return _success_response(simplified)
    except Exception as e:
        return _error_response(str(e), {"family": "Unknown"})


@tool(args_schema=RuleRecommendInput)
def rule_based_conditions_tool(
    reaction_smiles: str,
    database: Optional[str] = None,
    family_hint: Optional[str] = None,
    symptoms: Optional[Union[str, List[str]]] = None,
    combine_method: Literal["union", "all", "first", "separate"] = "union",
    include_reasoning: bool = True,
    auto_detect: bool = True,
    include_summary: bool = True,
) -> Dict[str, Any]:
    """
    Recommend reaction conditions using the deterministic rule engine.
    
    Loads the upgraded rule-based database (e.g., Buchwald CN rules) and applies
    feature-driven matching to produce condition guidance along with matched
    rules and modifiers. The tool can automatically infer the appropriate
    database by detecting the reaction family or accept an explicit path/name.
    
    Args:
        reaction_smiles: Reaction SMILES string (reactants>>products)
        database: Optional rule database name or JSON path
        family_hint: Optional reaction-family hint to assist auto selection
        symptoms: Observed symptoms/failure modes (list or comma-separated)
        combine_method: How to merge features across reactants
        include_reasoning: Include matched features/modifiers in result
        auto_detect: Attempt automatic family detection when database omitted
        include_summary: Include formatted text summary of recommendation
    
    Returns:
        Dict[str, Any]: Rule-engine recommendation payload with success flag.
    """
    start = time.perf_counter()
    detection_payload: Optional[Dict[str, Any]] = None
    attempted: List[str] = []

    try:
        db_path, db_name, detection_payload, attempted = _select_rule_database(
            reaction_smiles,
            database=database,
            family_hint=family_hint,
            auto_detect=auto_detect,
        )

        engine = _get_rule_engine(db_path)
        symptom_list = _normalize_symptom_list(symptoms)
        recommendation = engine.recommend(
            reaction_smiles,
            symptoms=symptom_list or None,
            combine_method=combine_method,
            include_reasoning=include_reasoning,
        )

        payload = recommendation.to_dict()
        payload["database_name"] = db_name
        payload["database_path"] = str(db_path)
        payload["combine_method"] = combine_method
        payload["symptoms"] = symptom_list or []
        payload["attempted_identifiers"] = attempted
        payload["auto_detect"] = auto_detect
        payload["include_reasoning"] = include_reasoning
        if detection_payload is not None:
            payload["family_detection"] = detection_payload

        db_meta = getattr(getattr(engine, "database", None), "metadata", None)
        if isinstance(db_meta, dict):
            payload["database_metadata"] = {
                key: db_meta.get(key)
                for key in ("name", "reaction_type", "version", "description")
                if db_meta.get(key) is not None
            }

        if include_summary:
            payload["summary"] = recommendation.format_summary()

        payload["timing_ms"] = round((time.perf_counter() - start) * 1000, 2)
        return _success_response(payload)

    except RuleDatabaseResolutionError as exc:
        details: Dict[str, Any] = {
            "attempted": getattr(exc, "attempted", []),
            "database": database,
            "family_hint": family_hint,
            "auto_detect": auto_detect,
        }
        detection_info = getattr(exc, "detection", None)
        if detection_info is not None:
            details["family_detection"] = detection_info
        if exc.reason:
            details["reason"] = exc.reason
        return _error_response(str(exc), details)
    except Exception as exc:
        details = {
            "database": database,
            "family_hint": family_hint,
            "auto_detect": auto_detect,
        }
        if detection_payload is not None:
            details["family_detection"] = detection_payload
        details["attempted"] = attempted
        return _error_response(str(exc), details)


@tool(args_schema=RuleBuilderAutoInput)
def rule_builder_autofill_tool(
    params: RuleBuilderAutoInput,
) -> Dict[str, Any]:
    """
    Auto-draft a rule database JSON using protocol text and reference reactions.

    This tool orchestrates an LLM extraction pass plus the deterministic
    RuleBuilder to emit a validated draft. Outputs include the serialized rule
    database dictionary, validation issues, and any parsing errors.
    """
    return run_rule_builder_autofill(params)


@tool(args_schema=ListSupportedCoresInput)
def list_supported_cores_tool(
    reaction_smiles: str,
    k: int = 25,
    search_all_families: bool = False,
) -> Dict[str, Any]:
    """
    List the catalyst cores observed in precedent reactions for the query.
    
    Useful when planning constraints (e.g., deciding whether Cu or Pd options
    exist for a given substrate pairing).
    """
    try:
        constraint_spec = build_constraint_spec(
            search_all_families=search_all_families,
        )
        normalized_key = _normalized_reaction_key(reaction_smiles)
        start = time.perf_counter()
        result, cache_hit = _recommendation_cache.get_or_compute(
            reaction_smiles=reaction_smiles,
            normalized_reaction=normalized_key,
            k=k,
            max_variants=1,
            rerank_strategy="rule",
            constraint_spec=constraint_spec,
        )
        duration_ms = (time.perf_counter() - start) * 1000

        alternatives = result.get("alternatives", {}) or {}
        core_counts: Sequence[Tuple[str, int]] = alternatives.get("cores", []) or []
        cores = [
            {"core": core, "support": support}
            for core, support in core_counts
        ]

        return _success_response({
            "family": result.get("family"),
            "detected_family": result.get("detected_family"),
            "search_all_families": search_all_families,
            "core_candidates": cores,
            "precedent_count": len(result.get("precedent_pack", {}).get("precedents", [])),
            "cache_hit": cache_hit,
            "timing_ms": round(duration_ms, 2),
        })
    except Exception as e:
        return _error_response(str(e), {"core_candidates": []})


@tool(args_schema=SearchPrecedentsInput)
def search_precedents_tool(
    reaction_smiles: str,
    k: int = 10,
    family: Optional[str] = None
) -> Dict[str, Any]:
    """
    Search for similar precedent reactions using DRFP fingerprints.
    
    Finds the k most similar reactions from the precedent database
    based on reaction fingerprint similarity.
    
    Args:
        reaction_smiles: Query reaction SMILES
        k: Number of similar reactions to find (default: 10)
        family: Optional reaction family filter (e.g., "Suzuki", "Buchwald_CN")
    
    Returns:
        Dict[str, Any]: Similar precedent summary with success flag.
    
    Example:
        Input: "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
        Output: {
            "precedents": [
                {
                    "reaction": "...",
                    "similarity": 0.95,
                    "conditions": {...},
                    "yield": 85
                },
                ...
            ]
        }
    """
    try:
        detected_family_id = None
        dataset_family = None

        # Auto-detect family if not provided
        if family is None:
            detection = detect_reaction(reaction_smiles, use_ml=False)
            detected_family_id = detection.get("family", "Unknown")
            dataset_family = canonical_family(detected_family_id)
        else:
            detected_family_id = family
            dataset_family = canonical_family(family)
        family = dataset_family or detected_family_id or "Unknown"
        
        # Use relax parameter for DRFP-based search
        # Disable reagent database filtering for families with specialized coupling reagents
        # (e.g., amide formation uses HATU, COMU, EDC which may not be in general database)
        disable_filter_families = ["amide_formation", "amide_coupling", "amidation"]
        family_key = (dataset_family or detected_family_id or family or "").lower()
        filter_by_db = family_key not in disable_filter_families
        
        features: Dict[str, Any] = {}
        if detected_family_id and detected_family_id != "Unknown":
            features["reaction_type"] = detected_family_id
            resolved_id = _reaction_catalog.resolve_reaction_type(detected_family_id)
            if resolved_id is None:
                resolved_id = _reaction_catalog.resolve_reaction_type(str(detected_family_id).lower())
            if resolved_id:
                record = _reaction_catalog.get_reaction_type(resolved_id)
                if record and record.category:
                    features["reaction_category"] = record.category

        result = precedent_knn(
            family=family,
            features=features,
            k=k,
            relax={
                "use_drfp": True, 
                "reaction_smiles": reaction_smiles,
                "filter_by_reagent_database": filter_by_db
            }
        )
        
        # Simplify precedent data
        precedents = result.get("precedents", [])
        simplified_precedents = []
        for p in precedents[:k]:
            # Extract reagent information from reagents list if available
            reagents = p.get("reagents", [])
            base = p.get("base", "")
            catalyst = p.get("catalyst", "")
            solvent = p.get("solvent", "")
            coupling_reagent = ""
            additives = []
            
            # Parse reagents list for amide formation reactions
            if isinstance(reagents, list):
                for r in reagents:
                    role = r.get("role", "").upper()
                    name = r.get("name") or r.get("abbreviation") or r.get("original_name", "")
                    if role == "BASE" and not base:
                        base = name
                    elif role == "CATALYST" and not catalyst:
                        catalyst = name
                    elif role == "SOLVENT" and not solvent:
                        solvent = name
                    elif role == "COUPLING_REAGENT":
                        coupling_reagent = name
                    elif role == "ADDITIVE":
                        additives.append(name)
            
            # Also check solvents list
            if not solvent:
                solvents = p.get("solvents", [])
                if isinstance(solvents, list) and solvents:
                    solvent = solvents[0].get("name", "") if isinstance(solvents[0], dict) else str(solvents[0])
            
            # Build conditions dict
            conditions = {}
            if catalyst and catalyst != {} and str(catalyst).strip():
                conditions["catalyst"] = catalyst
            if base:
                conditions["base"] = base
            if solvent:
                conditions["solvent"] = solvent
            if coupling_reagent:
                conditions["coupling_reagent"] = coupling_reagent
            if additives:
                conditions["additives"] = ", ".join(additives)
            
            # Add temperature and time if available
            if p.get("T_C") is not None:
                conditions["temperature_C"] = p.get("T_C")
            if p.get("time_h") is not None:
                conditions["time_h"] = p.get("time_h")
            
            simplified_precedents.append({
                "similarity": round(p.get("similarity", 0), 3),
                "conditions": conditions,
                "yield": p.get("yield", ""),
                "reaction": p.get("reaction_smiles", "")[:100]  # Truncate long SMILES
            })
        
        return _success_response({
            "family": family,
            "detected_family_id": detected_family_id,
            "dataset_family": dataset_family,
            "precedent_count": len(simplified_precedents),
            "precedents": simplified_precedents
        })
    except Exception as e:
        return _error_response(str(e), {"precedents": []})


# ============================================================================
# Reaction Dataset & Reagent Database Tools
# ============================================================================

@tool(args_schema=ReactionAnalyticsInput)
def reaction_dataset_analytics_tool(
    statistic: str = "summary",
    family: Optional[str] = None,
    role: Optional[str] = None,
    top_n: int = 10,
    min_yield: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Access analytics derived from reaction precedent datasets.

    Args:
        statistic: Metric to compute ('summary', 'top_catalysts', 'top_ligands', 'top_bases',
                  'top_solvents', 'top_condition_cores', 'top_reagents', or 'families').
        family: Reaction family identifier (required for most metrics).
        role: Reagent role filter required when statistic='top_reagents'.
        top_n: Maximum number of ranked entries to return.
        min_yield: Optional yield threshold (%) when ranking reagents.

    Returns:
        Dict[str, Any]: Analytics payload with success flag.
    """
    try:
        stat_key = (statistic or "summary").strip().lower()

        def _format_ranked(items: Sequence[Tuple[str, int, Optional[float]]]) -> List[Dict[str, Any]]:
            formatted: List[Dict[str, Any]] = []
            for name, count, avg in items:
                formatted.append(
                    {
                        "name": str(name),
                        "count": int(count),
                        "avg_yield": float(avg) if avg is not None else None,
                    }
                )
            return formatted[:top_n]

        if stat_key in {"families", "list_families", "available_families"}:
            families = get_all_families()
            return _success_response({"statistic": "families", "families": families})

        if not family or not str(family).strip():
            return _error_response("family is required for reaction dataset analytics")

        family_key = str(family).strip()

        if stat_key in {"summary", "dataset_summary", "overview"}:
            stats = get_dataset_stats(family_key)
            summary = {
                "family": family_key,
                "total_reactions": stats.get("total_reactions"),
                "unique_condition_cores": stats.get("unique_condition_cores"),
                "unique_solvents": stats.get("unique_solvents"),
                "unique_bases": stats.get("unique_bases"),
                "unique_catalysts": stats.get("unique_catalysts"),
                "yield_stats": stats.get("yield_stats"),
                "temperature_stats": stats.get("temperature_stats"),
                "time_stats": stats.get("time_stats"),
            }
            return _success_response({"statistic": "summary", "summary": summary})

        if stat_key in {"top_catalysts", "catalysts"}:
            ranked = get_common_catalysts(family_key, top_n=top_n, min_yield=min_yield)
            return _success_response(
                {"statistic": "top_catalysts", "results": _format_ranked(ranked)}
            )

        if stat_key in {"top_ligands", "ligands"}:
            ranked = get_common_ligands(family_key, top_n=top_n, min_yield=min_yield)
            return _success_response(
                {"statistic": "top_ligands", "results": _format_ranked(ranked)}
            )

        if stat_key in {"top_bases", "bases"}:
            ranked = get_common_bases(family_key, top_n=top_n, min_yield=min_yield)
            return _success_response({"statistic": "top_bases", "results": _format_ranked(ranked)})

        if stat_key in {"top_solvents", "solvents"}:
            ranked = get_common_solvents(family_key, top_n=top_n, min_yield=min_yield)
            return _success_response(
                {"statistic": "top_solvents", "results": _format_ranked(ranked)}
            )

        if stat_key in {"top_condition_cores", "condition_cores", "cores"}:
            ranked = get_condition_cores(family_key, top_n=top_n, min_yield=min_yield)
            results: List[Dict[str, Any]] = []
            for core, count, avg in ranked[:top_n]:
                results.append(
                    {
                        "core": str(core),
                        "count": int(count),
                        "avg_yield": float(avg) if avg is not None else None,
                    }
                )
            return _success_response(
                {"statistic": "top_condition_cores", "results": results}
            )

        if stat_key in {"top_reagents", "reagents"}:
            if not role or not str(role).strip():
                return _error_response("role is required when statistic='top_reagents'")
            role_key = str(role).strip()
            ranked = get_common_reagents(
                family_key,
                role_key,
                top_n=top_n,
                min_yield=min_yield,
            )
            results = [
                {
                    "role": role_key,
                    "name": str(name),
                    "count": int(count),
                    "avg_yield": float(avg) if avg is not None else None,
                }
                for name, count, avg in ranked[:top_n]
            ]
            return _success_response(
                {"statistic": "top_reagents", "role": role_key, "results": results}
            )

        return _error_response(f"Unsupported statistic '{statistic}'.")
    except Exception as exc:
        return _error_response(str(exc))


@tool(args_schema=ReagentAnalyticsInput)
def reagent_database_analytics_tool(
    statistic: str = "summary",
    role: Optional[str] = None,
    registry_dir: Optional[str] = None,
    top_families: int = 20,
) -> Dict[str, Any]:
    """
    Access aggregated analytics for the reagent database.

    Args:
        statistic: Metric to compute ('summary', 'role_summary', 'missing_data').
        role: Reagent role used when requesting a role summary.
        registry_dir: Optional override for the reagent registry directory.
        top_families: Maximum number of families to include in ranked outputs.

    Returns:
        Dict[str, Any]: Analytics payload with a success flag and structured data.
    """
    try:
        stat_key = (statistic or "summary").strip().lower()
        registry_path = registry_dir.strip() if isinstance(registry_dir, str) else registry_dir
        if isinstance(registry_path, str) and not registry_path:
            registry_path = None

        if stat_key in {"summary", "database_summary", "overview"}:
            raw_stats = get_database_statistics(registry_path)

            families_by_role = {
                role_key: sorted(set(families))
                for role_key, families in raw_stats.get("families_by_role", {}).items()
            }

            family_counts_raw = raw_stats.get("family_counts", {})
            if hasattr(family_counts_raw, "most_common"):
                family_counts_list = family_counts_raw.most_common(top_families)
            elif isinstance(family_counts_raw, dict):
                family_counts_list = list(family_counts_raw.items())[:top_families]
            else:
                family_counts_list = []
            family_counts = {
                str(name): int(count)
                for name, count in family_counts_list
            }

            top_families_list = [
                {"family": str(name), "count": int(count)}
                for name, count in list(raw_stats.get("top_families", []))[:top_families]
            ]

            roles_per_reagent_entries = []
            roles_per_reagent_raw = raw_stats.get("roles_per_reagent", {})
            try:
                iterator = roles_per_reagent_raw.items()
            except AttributeError:
                iterator = []
            for num_roles, count in sorted(iterator, key=lambda item: item[0]):
                roles_per_reagent_entries.append(
                    {"roles": int(num_roles), "reagents": int(count)}
                )

            summary = {
                "registry_dir": raw_stats.get("registry_dir"),
                "total_reagents": raw_stats.get("total_reagents"),
                "by_type": dict(raw_stats.get("by_type", {})),
                "total_with_cas": raw_stats.get("total_with_cas"),
                "total_with_inchikey": raw_stats.get("total_with_inchikey"),
                "total_with_smiles": raw_stats.get("total_with_smiles"),
                "total_with_abbreviations": raw_stats.get("total_with_abbreviations"),
                "id_format_stats": dict(raw_stats.get("id_format_stats", {})),
                "multi_role_reagents": raw_stats.get("multi_role_reagents"),
                "families_by_role": families_by_role,
                "family_counts": family_counts,
                "top_families": top_families_list,
                "roles_per_reagent_breakdown": roles_per_reagent_entries,
            }

            return _success_response(
                {"statistic": "summary", "summary": summary}
            )

        if stat_key in {"role_summary", "role", "families"}:
            if not role or not role.strip():
                return _error_response("role is required when statistic='role_summary'")
            role_key = role.strip()
            role_stats = get_family_statistics(role_key, registry_path)
            families = role_stats.get("families", [])
            limited_families = families[:top_families] if top_families else families
            payload = {
                "role": role_stats.get("role", role_key),
                "total_reagents": role_stats.get("total_reagents"),
                "total_families": role_stats.get("total_families"),
                "families": [
                    {"name": str(entry.get("name", "")), "count": int(entry.get("count", 0))}
                    for entry in limited_families
                ],
            }
            return _success_response(
                {"statistic": "role_summary", "summary": payload}
            )

        if stat_key in {"missing_data", "missing", "data_gaps"}:
            report = get_missing_data_report(registry_path)
            return _success_response({"statistic": "missing_data", "report": report})

        return _error_response(f"Unsupported statistic '{statistic}'.")
    except Exception as exc:
        return _error_response(str(exc))


@tool(args_schema=FindReagentInput)
def find_reagent_tool(query: str, reagent_type: str = "base") -> Dict[str, Any]:
    """
    Look up reagent information from the reagent database.
    
    Searches by name or abbreviation to find reagent properties,
    roles, and metadata. Supports multiple reagent types.
    
    Args:
        query: Reagent name or abbreviation (e.g., "Cs2CO3", "XPhos")
        reagent_type: Type of reagent database to search 
                     (base, solvent, ligand, metal, additive, etc.)
                     Default: "base"
    
    Returns:
        Dict[str, Any]: Reagent metadata with success flag.
    
    Example:
        Input: query="Cs2CO3", reagent_type="base"
        Output: {
            "name": "Cesium carbonate",
            "cas": "534-17-8",
            "role": "base",
            "smiles": "[Cs+].[Cs+].[O-]C([O-])=O"
        }
    """
    try:
        # Try common reagent types if not found
        types_to_try = [reagent_type, "base", "solvent", "ligand", "metal", "additive"]
        
        for r_type in types_to_try:
            result = find_reagent(query, r_type)
            if result:
                # Simplify the output
                simplified = {
                    "name": result.get("name", ""),
                    "cas": result.get("cas", ""),
                    "role": result.get("role", r_type),
                    "family": result.get("family", ""),
                    "smiles": result.get("smiles", ""),
                    "source": result.get("data_source", ""),
                    "reagent_type": r_type
                }
                return _success_response(simplified)
        
        return _error_response(f"Reagent not found: {query}")
    except Exception as e:
        return _error_response(str(e), {"query": query, "reagent_type": reagent_type})


@tool(args_schema=AddReagentInput)
def add_reagent_tool(
    cas: str,
    name: Optional[str] = None,
    synonyms: Optional[Any] = None,
    role: Optional[str] = None,
    family_id: Optional[str] = None,
    abbreviation: Optional[str] = None,
    smiles: Optional[str] = None,
    taxonomy_dir: Optional[str] = None,
    allow_default_family: bool = True,
    dry_run: bool = False,
    auto_resolve: bool = True,
    resolver_timeout: float = REAGENT_RESOLVER_TIMEOUT,
) -> Dict[str, Any]:
    """
    Add a reagent entry to the taxonomy registry.

    Args:
        cas: CAS identifier for the reagent.
        name: Preferred reagent name (optional if resolver is enabled).
        synonyms: Additional synonyms (list or comma-separated string).
        role: Explicit role override.
        family_id: Explicit family override.
        abbreviation: Abbreviation override (defaults to name).
        smiles: Optional SMILES annotation.
        taxonomy_dir: Optional path to a writable taxonomy directory.
        allow_default_family: Allow default family fallback when inference fails.
        dry_run: When true, preview without saving.
        auto_resolve: Use CAS resolver to fill missing data.
        resolver_timeout: Resolver timeout in seconds.
    
    Returns:
        Dict[str, Any]: Result payload including success flag and status metadata.
    """
    try:
        normalised_synonyms = _normalise_synonyms_field(synonyms)
        taxonomy_path = taxonomy_dir.strip() if isinstance(taxonomy_dir, str) else taxonomy_dir
        if isinstance(taxonomy_path, str) and not taxonomy_path:
            taxonomy_path = None
        result = add_reagent_entry(
            cas=cas,
            taxonomy_dir=taxonomy_path,
            name=name,
            synonyms=normalised_synonyms,
            abbreviation=abbreviation,
            role=role,
            family_id=family_id,
            smiles=smiles,
            allow_default_family=allow_default_family,
            dry_run=dry_run,
            auto_resolve=auto_resolve,
            resolver_timeout=resolver_timeout,
        )
        if isinstance(result, dict) and result.get("error"):
            return _error_response(str(result["error"]), {"details": result})
        payload = result if isinstance(result, dict) else {"result": result}
        return _success_response(payload)
    except ReagentAdditionError as exc:
        return _error_response(str(exc))
    except Exception as exc:
        return _error_response(str(exc))


@tool(args_schema=ProtocolRecommendInput)
def protocol_recommendation_tool(
    reaction_smiles: str,
    k: int = 5,
    family_filter: Optional[str] = None,
    use_smarts_filter: bool = True,
    min_similarity: float = 0.3
) -> Dict[str, Any]:
    """
    Find experimental protocols for reactions similar to the query.
    
    Returns complete procedure information including:
    - Step-by-step experimental instructions
    - Reagent amounts and equivalents
    - Reaction conditions (temperature, time, atmosphere)
    - Workup and purification procedures
    - Literature references and sources
    
    This is complementary to ML-based condition recommendations - it provides
    full experimental procedures from literature rather than just reagent suggestions.
    
    Args:
        reaction_smiles: Reaction SMILES to match against protocol database
        k: Number of protocol recommendations to return (default 5, max 20)
        family_filter: Optional reaction family filter (e.g. "Suzuki", "Buchwald_CN")
        use_smarts_filter: Enable SMARTS-based structural filtering for better matches
        min_similarity: Minimum DRFP similarity threshold (0.0-1.0, default 0.3)
    
    Returns:
        Dict with ranked protocol recommendations including full procedures
    
    Example:
        >>> protocol_recommendation_tool(
        ...     reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
        ...     k=3,
        ...     family_filter='Suzuki'
        ... )
    """
    if not PROTOCOL_AVAILABLE:
        return _error_response(
            "Protocol recommendation module not available. "
            "The protocol database may not be installed or accessible.",
            {"feature": "protocol_recommendation", "available": False}
        )
    
    try:
        recommender = ProtocolRecommender()
        results = recommender.recommend(
            reaction_smiles=reaction_smiles,
            k=k,
            reaction_family=family_filter,  # Use correct parameter name
            use_smarts_filter=use_smarts_filter,
            min_similarity=min_similarity
        )
        return _success_response(results)
    except FileNotFoundError as e:
        return _error_response(
            f"Protocol database not found: {str(e)}",
            {"hint": "Run 'python -m chemtools.protocol.cli build' to create the index"}
        )
    except Exception as e:
        return _error_response(str(e), {"reaction_smiles": reaction_smiles})


@tool(args_schema=ReactionSimilarityInput)
def reaction_similarity_tool(
    reaction1_smiles: str,
    reaction2_smiles: str
) -> Dict[str, Any]:
    """
    Calculate DRFP-based similarity between two reactions.
    
    Returns Tanimoto coefficient (0.0-1.0) measuring reaction similarity:
    - 1.0 = identical reactions
    - 0.8-1.0 = very similar reactions (likely analogous)
    - 0.6-0.8 = moderately similar (same general transformation)
    - 0.4-0.6 = somewhat related
    - 0.0-0.4 = different reactions
    
    Useful for:
    - Comparing user's reaction to literature precedents
    - Assessing whether two reactions are analogous
    - Finding reaction analogs for optimization
    
    Args:
        reaction1_smiles: First reaction SMILES string
        reaction2_smiles: Second reaction SMILES string
    
    Returns:
        Dict with similarity score and interpretation
    
    Example:
        >>> reaction_similarity_tool(
        ...     'CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
        ...     'c1ccccc1Br.CCB(O)O>>CCc1ccccc1'
        ... )
    """
    if not DRFP_SIMILARITY_AVAILABLE:
        return _error_response(
            "Reaction similarity module not available. "
            "DRFP fingerprinting may not be installed.",
            {"feature": "reaction_similarity", "available": False}
        )
    
    try:
        similarity = compute_drfp_similarity(reaction1_smiles, reaction2_smiles)
        
        # Provide interpretation
        if similarity >= 0.8:
            interpretation = "Very similar reactions - likely analogous transformations"
            category = "very_similar"
        elif similarity >= 0.6:
            interpretation = "Moderately similar - same general reaction type"
            category = "moderately_similar"
        elif similarity >= 0.4:
            interpretation = "Somewhat related reactions"
            category = "somewhat_related"
        else:
            interpretation = "Different reactions"
            category = "different"
        
        return _success_response({
            "reaction1": reaction1_smiles,
            "reaction2": reaction2_smiles,
            "similarity": round(float(similarity), 4),
            "interpretation": interpretation,
            "category": category
        })
    except Exception as e:
        return _error_response(str(e), {
            "reaction1": reaction1_smiles,
            "reaction2": reaction2_smiles
        })


@tool(args_schema=ListAllFamiliesInput)
def list_all_families_tool() -> Dict[str, Any]:
    """
    List all reaction families available in the precedent dataset.
    
    Returns a comprehensive list of reaction family names that can be used
    for filtering recommendations, searching precedents, and understanding
    the scope of available data.
    
    Useful for:
    - Understanding what reaction types are supported
    - Choosing appropriate family filters for searches
    - Exploring available precedent data coverage
    
    Returns:
        Dict with list of family names and counts
    
    Example:
        >>> list_all_families_tool()
        {
            "success": True,
            "families": ["Suzuki", "Buchwald_CN", "Amide_formation", ...],
            "count": 15
        }
    """
    try:
        families = get_all_families()
        
        # Convert to list if it's a set or other iterable
        if not isinstance(families, list):
            families = sorted(list(families))
        
        return _success_response({
            "families": families,
            "count": len(families),
            "note": "Use these family names as filters in recommendation and search tools"
        })
    except Exception as e:
        return _error_response(str(e))


# ============================================================================
# Unified Recommender Tool (Protocol + Rule Search)
# ============================================================================

class UnifiedRecommenderInput(BaseModel):
    """Input schema for unified_recommender_tool."""
    reaction_smiles: str = Field(..., description="Full reaction SMILES string (reactants>>products)")
    top_k: int = Field(default=1, description="Number of recommendations to return (max 10)", ge=1, le=10)
    min_similarity: float = Field(default=0.0, description="Minimum DRFP similarity threshold (0.0-1.0)", ge=0.0, le=1.0)
    source_type: Optional[str] = Field(default=None, description="Filter by source: 'protocol' (full experimental procedures) or 'rule' (general guidelines)")
    validate_rules: bool = Field(default=True, description="Enable chemical validation: applies_if for rules, reaction_SMARTS for protocols")
    format_for_automation: bool = Field(default=False, description="Format conditions with ordered addition sequences for automated execution")
    scale_mmol: float = Field(default=1.0, description="Reaction scale in mmol (for automation format)", ge=0.1, le=1000.0)


class AutoConditionsLLMInput(BaseModel):
    """Schema for deterministic auto-conditions pipeline (LLM-callable)."""

    reaction_smiles: str = Field(..., description="Reaction SMILES (reactants>>products).")
    constraints: Optional[str] = Field(
        default=None,
        description="Natural language constraints (e.g., 'no palladium', 'prefer copper', 'avoid chlorinated solvents').",
    )
    top_k_protocols: int = Field(
        default=5,
        ge=1,
        le=10,
        description="Maximum DRFP precedent candidates to retrieve.",
    )
    max_protocols: int = Field(
        default=3,
        ge=1,
        le=5,
        description="Maximum number of protocols to format.",
    )


class PlannerDetectInput(BaseModel):
    """Schema for planner family detection."""

    reaction_smiles: str = Field(..., description="Reaction SMILES (reactants>>products).")


class PlannerRuleCandidatesInput(BaseModel):
    """Schema for planner rule candidate fetch."""

    reaction_smiles: str = Field(..., description="Reaction SMILES (reactants>>products).")


class PlannerProtocolCandidatesInput(BaseModel):
    """Schema for planner DRFP protocol candidates."""

    reaction_smiles: str = Field(..., description="Reaction SMILES (reactants>>products).")
    top_k: int = Field(5, ge=1, le=10, description="Number of DRFP candidates to retrieve.")


class PlannerHTEInput(BaseModel):
    """Schema for planner HTE summary fetch."""

    reaction_smiles: str = Field(..., description="Reaction SMILES (reactants>>products).")


class PlannerScoreInput(BaseModel):
    """Schema for ML scoring of arbitrary candidates."""

    reaction_smiles: str = Field(..., description="Reaction SMILES (reactants>>products).")
    candidates: List[Dict[str, Any]] = Field(
        ...,
        description="List of candidate condition dicts (candidate_id, source, components, raw_score optional).",
    )


class PlannerFuseInput(BaseModel):
    """Schema for fuse step on auto-generated candidates."""

    reaction_smiles: str = Field(..., description="Reaction SMILES (reactants>>products).")
    top_k_protocols: int = Field(5, ge=1, le=10, description="DRFP candidates to fetch.")


@tool(args_schema=UnifiedRecommenderInput)
def unified_recommender_tool(
    reaction_smiles: str,
    top_k: int = 1,
    min_similarity: float = 0.0,
    source_type: Optional[str] = None,
    validate_rules: bool = True,
    format_for_automation: bool = False,
    scale_mmol: float = 1.0
) -> Dict[str, Any]:
    """
    Find similar reactions and condition recommendations using DRFP similarity with chemical validation.
    
    This tool provides a unified search across both:
    - **Protocols**: Full experimental procedures from literature
    - **Rules**: General reaction guidelines and best practices
    
    Uses DRFP (Differentiable Reaction Fingerprints) for reaction similarity,
    with optional post-similarity validation for chemical appropriateness:
    - **Protocols**: Validated using reaction_SMARTS patterns (exact transformation matching)
    - **Rules**: Validated using applies_if criteria (functional group detection)
    
    **When to use this tool:**
    - Need recommendations for an unfamiliar reaction
    - Want to find similar precedents from literature
    - Looking for general reaction guidelines
    - Comparing different condition options
    
    **Source types:**
    - `protocol`: Returns only full experimental procedures (e.g., from Org. Synth.)
    - `rule`: Returns only general reaction guidelines
    - None (default): Returns both protocols and rules, ranked by similarity
    
    **Validation (default: enabled):**
    - Filters out chemically inappropriate recommendations
    - Protocols: Checks if reaction_SMARTS patterns match the query transformation
    - Rules: Checks if detected features meet applies_if criteria
    - Recommended to keep enabled for chemical accuracy
    
    Args:
        reaction_smiles: Full reaction SMILES (reactants>>products)
        top_k: Number of results to return (1-10, default 1)
        min_similarity: Minimum similarity threshold (0.0-1.0, default 0.0)
        source_type: Filter by 'protocol', 'rule', or None for both
        validate_rules: Enable chemical validation (default True, recommended)
    
    Returns:
        Dict with ranked recommendations including similarity scores, metadata,
        and source attribution
    
    Example:
        >>> unified_recommender_tool(
        ...     reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        ...     top_k=3,
        ...     min_similarity=0.3,
        ...     validate_rules=True
        ... )
        {
            "success": True,
            "recommendations": [
                {
                    "rank": 1,
                    "name": "Buchwald–Hartwig C–N Coupling",
                    "source_type": "rule",
                    "family": "Buchwald–Hartwig_C–N",
                    "similarity": 0.82,
                    "tags": ["C-N-coupling", "amination", "cross-coupling"],
                    "validated": true
                },
                ...
            ],
            "query": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            "count": 3,
            "validation_enabled": true
        }
    """
    if not UNIFIED_RECOMMENDER_AVAILABLE:
        return _error_response(
            "Unified recommender not available. Install DRFP: pip install drfp"
        )
    
    try:
        # Initialize recommender (uses default index)
        recommender = UnifiedRecommender()
        
        # Convert source_type to list format
        source_types = None
        if source_type:
            source_types = [source_type.lower()]
        
        # Get recommendations with validation
        results = recommender.recommend(
            reaction_smiles=reaction_smiles,
            top_k=top_k,
            min_similarity=min_similarity,
            source_types=source_types,
            validate_rules=validate_rules,  # Enable chemical validation
            format_for_automation=format_for_automation,
            scale_mmol=scale_mmol
        )
        
        # Format results for LLM
        formatted_results = []
        for result in results:
            result_dict = {
                "rank": result.rank,
                "id": result.id,
                "name": result.name,
                "source_type": result.source_type,
                "family": result.family,
                "similarity": round(result.similarity, 3),
                "tags": result.tags,
                "version": result.version,
                "validated": validate_rules  # Indicate if validation was applied
            }
            
            # Add automation format if requested
            if format_for_automation and result.full_data:
                result_dict["reaction_setup"] = result.full_data.get("reaction_setup", [])
                result_dict["metadata"] = result.full_data.get("metadata", {})
            
            formatted_results.append(result_dict)
        
        return _success_response({
            "recommendations": formatted_results,
            "query": reaction_smiles,
            "count": len(formatted_results),
            "validation_enabled": validate_rules,
            "automation_format": format_for_automation,
            "filters": {
                "top_k": top_k,
                "min_similarity": min_similarity,
                "source_type": source_type or "all",
                "scale_mmol": scale_mmol if format_for_automation else None
            }
        })
    
    except Exception as e:
        return _error_response(f"Unified recommendation failed: {str(e)}")


# ============================================================================
# Auto-Conditions Deterministic Pipeline (LLM-callable)
# ============================================================================

@tool(args_schema=AutoConditionsLLMInput)
def auto_conditions_llm_tool(
    reaction_smiles: str,
    constraints: Optional[str] = None,
    top_k_protocols: int = 5,
    max_protocols: int = 3,
) -> Dict[str, Any]:
    """
    Deterministic auto-conditions pipeline for LLM agents.

    Chains family detection → rule candidates → DRFP precedents → HTE summary →
    heuristic/ML scoring → fusion → automation-ready protocol formatting.
    
    Supports natural language constraints like 'no palladium' or 'prefer copper'.
    """
    try:
        spec = build_constraint_spec(text=constraints) if constraints else None
        rxn_input = ReactionInput(reaction_smiles=reaction_smiles)
        result = auto_conditions(
            rxn_input,
            constraints=spec,
            top_k_protocols=top_k_protocols,
            max_protocols=max_protocols,
            build_protocols=True,
        )
        payload = result.model_dump()
        return _success_response(payload)
    except Exception as e:
        return _error_response(str(e))


# ============================================================================
# Planner Building Blocks (LLM-callable)
# ============================================================================


@tool(args_schema=PlannerDetectInput)
def planner_detect_family_tool(reaction_smiles: str) -> Dict[str, Any]:
    """Detect reaction family using the deterministic planner router."""
    result = planner_detect_family(ReactionInput(reaction_smiles=reaction_smiles))
    return _success_response(result.model_dump())


@tool(args_schema=PlannerRuleCandidatesInput)
def planner_rule_candidates_tool(reaction_smiles: str) -> Dict[str, Any]:
    """Fetch rule-based candidates for a reaction (auto-detected family)."""
    rxn = ReactionInput(reaction_smiles=reaction_smiles)
    family = planner_detect_family(rxn)
    candidates = fetch_rule_candidates(rxn, family)
    return _success_response(
        {
            "family": family.model_dump(),
            "candidates": [c.model_dump() for c in candidates],
        }
    )


@tool(args_schema=PlannerProtocolCandidatesInput)
def planner_protocol_candidates_tool(reaction_smiles: str, top_k: int = 5) -> Dict[str, Any]:
    """Fetch DRFP precedent/protocol candidates for a reaction."""
    rxn = ReactionInput(reaction_smiles=reaction_smiles)
    candidates = find_similar_protocols(rxn, top_k=top_k)
    return _success_response({"candidates": [c.model_dump() for c in candidates]})


@tool(args_schema=PlannerHTEInput)
def planner_hte_summary_tool(reaction_smiles: str) -> Dict[str, Any]:
    """Fetch HTE summary for the detected family."""
    rxn = ReactionInput(reaction_smiles=reaction_smiles)
    family = planner_detect_family(rxn)
    summary = fetch_hte_stats(rxn, family)
    return _success_response(
        {
            "family": family.model_dump(),
            "hte_summary": summary.model_dump() if summary else None,
        }
    )


@tool(args_schema=PlannerScoreInput)
def planner_score_candidates_tool(reaction_smiles: str, candidates: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Score arbitrary candidates with the planner heuristic/ML stub."""
    rxn = ReactionInput(reaction_smiles=reaction_smiles)
    candidate_objs: List[CandidateCondition] = []
    for c in candidates:
        candidate_objs.append(
            CandidateCondition(
                candidate_id=c.get("candidate_id", f"cand_{len(candidate_objs)}"),
                components=c.get("components", {}),
                source=c.get("source", "unknown"),
                raw_score=c.get("raw_score"),
                metadata=c.get("metadata", {}),
            )
        )
    scores = score_ml_candidates(candidate_objs, rxn)
    return _success_response({"scores": [s.model_dump() for s in scores]})


@tool(args_schema=PlannerFuseInput)
def planner_fuse_tool(reaction_smiles: str, top_k_protocols: int = 5) -> Dict[str, Any]:
    """
    Generate candidates (rules + DRFP), score, and fuse into a ranked list.

    Returns fused ranking and provenance without formatting protocol steps.
    """
    rxn = ReactionInput(reaction_smiles=reaction_smiles)
    family = planner_detect_family(rxn)
    rule_candidates = fetch_rule_candidates(rxn, family)
    protocol_candidates = find_similar_protocols(rxn, top_k=top_k_protocols)
    hte_summary = fetch_hte_stats(rxn, family)
    ml_scores = score_ml_candidates(rule_candidates + protocol_candidates, rxn)
    fused = fuse_scores(rule_candidates, protocol_candidates, hte_summary, ml_scores)
    return _success_response(
        {
            "family": family.model_dump(),
            "fused": [c.model_dump() for c in fused.ranked],
            "provenance": fused.provenance,
        }
    )


# ============================================================================
# HTE Recommendation and Analytics Tools
# ============================================================================

@tool(args_schema=HTERecommendInput)
def hte_recommend_tool(
    reactant_a_smiles: Optional[str] = None,
    reactant_b_smiles: Optional[str] = None,
    reaction_smiles: Optional[str] = None,
    top_k: int = 5,
    min_experiments: int = 1,
    reaction_type_filter: Optional[str] = None,
    catalyst_filter: Optional[str] = None
) -> Dict[str, Any]:
    """
    Recommend reaction conditions based on HTE (High-Throughput Experimentation) data.
    
    This tool provides condition recommendations based on **66,308 experimental results**
    across **41 reaction types**. Can work with either individual reactants OR complete
    reaction SMILES (reactants will be auto-extracted).
    
    **Key Features:**
    - Fast: <100ms query time
    - Works with reaction SMILES OR individual reactants
    - Z-score based ranking (primary metric)
    - Statistical confidence with success rates & sample sizes
    - Covers Suzuki, C-N coupling, Buchwald-Hartwig, amide formation, and 37+ more
    
    **Use when:**
    - Planning a new reaction and need starting conditions
    - Comparing catalyst options (Pd vs Cu vs Ni)
    - Need conditions for common reaction types
    - Want data-backed recommendations with success statistics
    
    Args:
        reactant_a_smiles: First reactant SMILES (e.g., aryl bromide). Use this OR reaction_smiles.
        reactant_b_smiles: Second reactant SMILES (e.g., amine). Optional.
        reaction_smiles: Complete reaction SMILES (reactants>>products). Reactants auto-extracted.
        top_k: Number of recommendations to return (1-20, default 5)
        min_experiments: Minimum experiments required per condition (default 2)
        reaction_type_filter: Optional reaction filter ("Suzuki", "C-N", "C_N_Coupling", "Buchwald", etc.)
        catalyst_filter: Optional catalyst metal filter ("Pd", "Cu", "Ni", "copper", "palladium", etc.)
    
    Returns:
        Dict with 'success' boolean, 'matching_experiments' count, 'recommendations' list,
        and reactant type information. When successful with data, 'matching_experiments' > 0
        and 'recommendations' contains ranked conditions with Z-scores, yields, and reagents.
        
        CRITICAL: If 'success' is True and 'matching_experiments' > 0, DATA WAS FOUND!
        The 'recommendations' list contains actual experimental conditions to use.
    
    Examples:
        >>> # With individual reactants
        >>> hte_recommend_tool(
        ...     reactant_a_smiles="Brc1ccccc1",
        ...     reactant_b_smiles="CCN",
        ...     top_k=3
        ... )
        
        >>> # With reaction SMILES (auto-extracts reactants)
        >>> hte_recommend_tool(
        ...     reaction_smiles="Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1",
        ...     reaction_type_filter="C_N_Coupling",
        ...     catalyst_filter="Cu",
        ...     top_k=10,
        ...     min_experiments=1
        ... )
    """
    if not HTE_AVAILABLE:
        return _error_response(
            "HTE tools not available. Install with: pip install chemtools[hte]",
            {"recommendation": "Install HTE dependencies"}
        )
    
    try:
        # Parse reaction SMILES if provided
        if reaction_smiles:
            from chemtools.featurizers.analysis.smiles import _split_reaction_smiles
            
            parts = _split_reaction_smiles(reaction_smiles.strip())
            reactants_str = parts[0] if len(parts) > 0 else ""
            
            if not reactants_str:
                return _error_response(
                    "Could not extract reactants from reaction SMILES. Check format (should be reactants>>products)."
                )
            
            # Split reactants by dot notation
            reactant_list = [r.strip() for r in reactants_str.split(".") if r.strip()]
            
            if len(reactant_list) >= 1:
                reactant_a_smiles = reactant_list[0]
            if len(reactant_list) >= 2:
                reactant_b_smiles = reactant_list[1]
            
            # If more than 2 reactants, log but use first two
            if len(reactant_list) > 2:
                print(f"Note: Reaction has {len(reactant_list)} reactants, using first two for HTE lookup.")
        
        # Validate we have at least one reactant
        if not reactant_a_smiles:
            return _error_response(
                "Must provide either reactant_a_smiles or reaction_smiles."
            )
        
        recommender = HTERecommender()
        
        result = recommender.recommend(
            reactant_a_smiles=reactant_a_smiles,
            reactant_b_smiles=reactant_b_smiles,
            top_k=top_k,
            min_experiments=min_experiments,
            reaction_type_filter=reaction_type_filter,
            catalyst_filter=catalyst_filter
        )
        
        # Convert to serializable dict
        recommendations = []
        for i, rec in enumerate(result.recommendations, 1):
            rec_dict = {
                "rank": i,
                "avg_z_score": round(rec.avg_z_score, 2),
                "catalyst": rec.catalyst,
                "ligand": rec.ligand,
                "base": rec.base,
                "solvent": rec.solvent,
                "success_rate": round(rec.success_rate, 1),
                "avg_yield": round(rec.avg_yield, 1),
                "median_yield": round(rec.median_yield, 1),
                "num_experiments": rec.num_experiments,
                "confidence_score": round(rec.confidence_score, 1)
            }
            
            # Add optional components
            if rec.secondary_solvent:
                rec_dict["secondary_solvent"] = rec.secondary_solvent
            if rec.additive:
                rec_dict["additive"] = rec.additive
            if rec.coupling_reagent:
                rec_dict["coupling_reagent"] = rec.coupling_reagent
            
            recommendations.append(rec_dict)
        
        # Check if we have any recommendations
        if result.total_matching_experiments == 0 or len(recommendations) == 0:
            # Provide helpful error with detected types
            error_msg = (
                f"No HTE data found for reactant combination: {result.reactant_a_type} + {result.reactant_b_type}"
            )
            
            # Add filter information if applicable
            filters_applied = []
            if reaction_type_filter:
                filters_applied.append(f"reaction type '{reaction_type_filter}'")
            if catalyst_filter:
                filters_applied.append(f"catalyst '{catalyst_filter}'")
            
            if filters_applied:
                error_msg += f" with {' and '.join(filters_applied)}"
            
            # Suggest checking available data
            suggestion = (
                f"The reactants were successfully classified as {result.reactant_a_type} and {result.reactant_b_type}, "
                f"but this specific combination has no matching experiments in the HTE database. "
                f"Try: (1) removing filters to see if data exists without restrictions, "
                f"(2) using the hte_analytics_tool to explore available reactant pairs for this reaction type, "
                f"or (3) checking if the reactant types are common in the database."
            )
            
            return _error_response(
                error_msg,
                {
                    "reactant_a_type": result.reactant_a_type,
                    "reactant_b_type": result.reactant_b_type,
                    "reactant_a_smiles": reactant_a_smiles,
                    "reactant_b_smiles": reactant_b_smiles,
                    "predicted_reaction_type": result.predicted_reaction_type,
                    "matching_experiments": 0,
                    "suggestion": suggestion,
                    "filters_applied": {
                        "reaction_type": reaction_type_filter,
                        "catalyst": catalyst_filter
                    }
                }
            )
        
        return _success_response({
            "reactant_a_type": result.reactant_a_type,
            "reactant_b_type": result.reactant_b_type,
            "reactant_a_smiles": reactant_a_smiles,
            "reactant_b_smiles": reactant_b_smiles,
            "predicted_reaction_type": result.predicted_reaction_type,
            "reaction_confidence": round(result.reaction_type_confidence * 100, 1),
            "matching_experiments": result.total_matching_experiments,
            "recommendations": recommendations
        })
    
    except Exception as e:
        return _error_response(f"HTE recommendation failed: {str(e)}")


@tool(args_schema=HTEAnalyticsInput)
def hte_analytics_tool(
    query_type: str,
    reaction_type: Optional[str] = None,
    catalyst_filter: Optional[str] = None,
    reactant_a_type: Optional[str] = None,
    reactant_b_type: Optional[str] = None,
    min_experiments: int = 5,
    top_n: int = 20,
    sort_by: str = "count",
    similarity_criteria: Optional[str] = None
) -> Dict[str, Any]:
    """
    Analyze the HTE database to explore reactant pairs, catalysts, and reaction patterns.
    
    This tool provides analytics on **66,308 experimental results** to help understand:
    - Which reactant pairs work for specific reactions/catalysts
    - Which catalysts are most commonly used
    - Success rates and experimental coverage
    - Metal usage patterns (Pd vs Cu vs Ni)
    
    **Query Types:**
    - `list_pairs`: List reactant type pairs (e.g., ArBr + RNH2) with statistics
    - `catalysts`: Analyze catalyst usage and performance
    - `reactions`: Summarize all reaction types in database
    - `metals`: Analyze metal distribution (Pd, Cu, Ni, etc.)
    - `similar_pairs`: Find similar reactant combinations
    
    **Use when:**
    - Exploring what substrates work with specific catalysts
    - Comparing Pd vs Cu catalyst usage
    - Understanding reaction scope and coverage
    - Finding alternative substrate combinations
    
    Args:
        query_type: Type of analysis ('list_pairs', 'catalysts', 'reactions', 'metals', 'similar_pairs')
        reaction_type: Filter by reaction type. Supported formats:
            - 'Suzuki', 'Buchwald', 'Sonogashira' (common name)
            - 'C-N', 'C-N coupling' → C_N_Coupling
            - 'C-O', 'C-O coupling' → CO-Coupling (Ullmann)
            - 'C-C', 'C-C coupling' → CC-Coupling
            - 'C-S', 'C-S coupling' → CS-Coupling
        catalyst_filter: Filter by catalyst metal (e.g., 'Pd', 'Cu', 'Ni')
        reactant_a_type: Filter by reactant A type (e.g., 'ArBr')
        reactant_b_type: Filter by reactant B type (e.g., 'RNH2')
        min_experiments: Minimum experiments for inclusion (default 5)
        top_n: Maximum results to return (default 20)
        sort_by: Sort by 'count' or 'success_rate' (default 'count')
        similarity_criteria: For similar_pairs: 'reaction_type', 'catalyst', or 'both'
    
    Returns:
        Dict with query results including statistics, rankings, and metadata
    
    Example:
        >>> hte_analytics_tool(
        ...     query_type="list_pairs",
        ...     reaction_type="Suzuki",
        ...     catalyst_filter="Pd",
        ...     min_experiments=50,
        ...     top_n=10
        ... )
        {
            "success": True,
            "query_type": "list_pairs",
            "filters": {"reaction_type": "Suzuki", "catalyst": "Pd"},
            "results": [
                {
                    "reactant_a": "ArCl",
                    "reactant_b": "ArB(OR)2",
                    "reaction": "Suzuki",
                    "experiments": 2528,
                    "avg_yield": 30.4,
                    "success_rate": 25.9,
                    "top_catalyst": "dtbpfPdCl2"
                },
                ...
            ],
            "total_results": 16
        }
    """
    if not HTE_AVAILABLE:
        return _error_response(
            "HTE analytics not available. Install with: pip install chemtools[hte]",
            {"recommendation": "Install HTE dependencies"}
        )
    
    try:
        analytics = HTEAnalytics()
        
        # Normalize reaction type name to handle common variations
        # The HTE database uses underscores (e.g., "C_N_Coupling", "CO-Coupling")
        # Users might input "C-N" or "C-N coupling", so normalize to match
        if reaction_type:
            # Common mappings for user-friendly names
            reaction_type_map = {
                'c-n': 'C_N',
                'cn': 'C_N',
                'c-n coupling': 'C_N_Coupling',
                'cn coupling': 'C_N_Coupling',
                'buchwald': 'Buchwald',
                'buchwald-hartwig': 'Buchwald',
                'c-o': 'CO',
                'co': 'CO',
                'c-o coupling': 'CO-Coupling',
                'co coupling': 'CO-Coupling',
                'ullmann': 'CO-Coupling',
                'c-c': 'CC',
                'cc': 'CC',
                'c-c coupling': 'CC-Coupling',
                'cc coupling': 'CC-Coupling',
                'c-s': 'CS',
                'cs': 'CS',
                'c-s coupling': 'CS-Coupling',
                'cs coupling': 'CS-Coupling',
            }
            
            # Try exact match first (case-insensitive)
            reaction_lower = reaction_type.lower().strip()
            if reaction_lower in reaction_type_map:
                reaction_type = reaction_type_map[reaction_lower]
            # Otherwise use as-is and rely on str.contains() matching
        
        if query_type == "list_pairs":
            df = analytics.list_reactant_pairs(
                reaction_type=reaction_type,
                catalyst_filter=catalyst_filter,
                min_experiments=min_experiments,
                sort_by=sort_by
            )
            
            results = []
            for _, row in df.head(top_n).iterrows():
                results.append({
                    "reactant_a": row["Reactant_A_Type"],
                    "reactant_b": row["Reactant_B_Type"],
                    "reaction": row["Reaction_Type"],
                    "experiments": int(row["Num_Experiments"]),
                    "avg_yield": round(row["Avg_Yield"], 1),
                    "median_yield": round(row["Median_Yield"], 1),
                    "success_rate": round(row["Success_Rate"], 1),
                    "top_catalyst": row["Top_Catalyst"]
                })
            
            return _success_response({
                "query_type": "list_pairs",
                "filters": {
                    "reaction_type": reaction_type,
                    "catalyst": catalyst_filter,
                    "min_experiments": min_experiments
                },
                "results": results,
                "total_results": len(df)
            })
        
        elif query_type == "catalysts":
            df = analytics.get_catalyst_stats(
                reaction_type=reaction_type,
                reactant_a_type=reactant_a_type,
                reactant_b_type=reactant_b_type
            )
            
            results = []
            for _, row in df.head(top_n).iterrows():
                results.append({
                    "catalyst": row["Catalyst"],
                    "metal": row["Metal"],
                    "experiments": int(row["Num_Experiments"]),
                    "avg_yield": round(row["Avg_Yield"], 1),
                    "success_rate": round(row["Success_Rate"], 1),
                    "reaction_types": row["Reaction_Types"]
                })
            
            return _success_response({
                "query_type": "catalysts",
                "filters": {
                    "reaction_type": reaction_type,
                    "reactant_a": reactant_a_type,
                    "reactant_b": reactant_b_type
                },
                "results": results,
                "total_results": len(df)
            })
        
        elif query_type == "reactions":
            df = analytics.get_reaction_type_summary()
            
            results = []
            for _, row in df.head(top_n).iterrows():
                results.append({
                    "reaction_type": row["Reaction_Type"],
                    "experiments": int(row["Num_Experiments"]),
                    "reactant_pairs": int(row["Num_Reactant_Pairs"]),
                    "catalysts": int(row["Num_Catalysts"]),
                    "avg_yield": round(row["Avg_Yield"], 1),
                    "success_rate": round(row["Success_Rate"], 1),
                    "top_catalyst": row["Top_Catalyst"],
                    "top_pair": row["Top_Reactant_Pair"]
                })
            
            return _success_response({
                "query_type": "reactions",
                "results": results,
                "total_reactions": len(df),
                "database_total": 66308
            })
        
        elif query_type == "metals":
            result = analytics.analyze_metal_usage()
            
            metals = []
            for _, row in result["metal_distribution"].iterrows():
                metal_data = {
                    "metal": row["Metal"],
                    "experiments": int(row["Num_Experiments"]),
                    "percentage": round(row["Percentage"], 1)
                }
                
                # Add top reactions for this metal
                if row["Metal"] in result["by_reaction_type"]:
                    reactions = result["by_reaction_type"][row["Metal"]]
                    top_3 = sorted(reactions.items(), key=lambda x: x[1], reverse=True)[:3]
                    metal_data["top_reactions"] = [
                        {"reaction": rxn, "count": count}
                        for rxn, count in top_3
                    ]
                
                metals.append(metal_data)
            
            return _success_response({
                "query_type": "metals",
                "total_experiments": result["total_experiments"],
                "metals": metals
            })
        
        elif query_type == "similar_pairs":
            if not reactant_a_type or not reactant_b_type:
                return _error_response(
                    "similar_pairs query requires both reactant_a_type and reactant_b_type"
                )
            
            df = analytics.find_similar_pairs(
                reactant_a_type=reactant_a_type,
                reactant_b_type=reactant_b_type,
                similarity_criteria=similarity_criteria or "reaction_type"
            )
            
            results = []
            for _, row in df.head(top_n).iterrows():
                results.append({
                    "reactant_a": row["Reactant_A_Type"],
                    "reactant_b": row["Reactant_B_Type"],
                    "reaction": row["Reaction_Type"],
                    "experiments": int(row["Num_Experiments"]),
                    "avg_yield": round(row["Avg_Yield"], 1),
                    "success_rate": round(row["Success_Rate"], 1),
                    "top_catalyst": row["Top_Catalyst"]
                })
            
            return _success_response({
                "query_type": "similar_pairs",
                "query_pair": f"{reactant_a_type} + {reactant_b_type}",
                "similarity_criteria": similarity_criteria or "reaction_type",
                "results": results,
                "total_results": len(df)
            })
        
        else:
            return _error_response(
                f"Unknown query_type '{query_type}'. "
                "Valid types: list_pairs, catalysts, reactions, metals, similar_pairs"
            )
    
    except Exception as e:
        return _error_response(f"HTE analytics failed: {str(e)}")


@tool(args_schema=HTEConditionsInput)
def hte_conditions_tool(
    reactant_a_type: str,
    reactant_b_type: str,
    reaction_type: Optional[str] = None,
    catalyst_filter: Optional[str] = None,
    top_k: int = 10,
    min_experiments: int = 1,
    sort_by: str = "count"
) -> Dict[str, Any]:
    """
    Get detailed experimental conditions for a specific substrate pair from HTE database.
    
    This tool retrieves all tested conditions (catalyst, ligand, base, solvent combinations)
    for a specific pair of reactant types. Use this when you need to see what specific
    conditions were tested for a particular substrate combination.
    
    **Use when:**
    - User asks for "top conditions" for specific substrates
    - Need detailed breakdown of tested catalyst/ligand/base/solvent combinations
    - Comparing different conditions for the same substrate pair
    - Finding best performing conditions for specific reactants
    
    **Example queries this tool answers:**
    - "What are the top 10 conditions for ArI + Carbamate with copper catalyst?"
    - "Show me all conditions tested for ArBr + RNH2 in C-N coupling"
    - "What ligands work best with palladium for ArCl + Boronic acid?"
    
    Args:
        reactant_a_type: Reactant A type (e.g., 'ArI', 'ArBr', 'ArCl')
        reactant_b_type: Reactant B type (e.g., 'Carbamate', 'RNH2', 'arom-NH', 'Lactam')
        reaction_type: Optional reaction type filter (e.g., 'C-N', 'Suzuki')
        catalyst_filter: Optional catalyst metal filter (e.g., 'Pd', 'Cu', 'Ni')
        top_k: Number of top conditions to return (default 10)
        min_experiments: Minimum experiments per condition (default 1)
        sort_by: Sort by 'count', 'success', or 'yield' (default 'count')
    
    Returns:
        Dict with top conditions including catalyst, ligand, base, solvent, statistics
    
    Example:
        >>> hte_conditions_tool(
        ...     reactant_a_type="ArI",
        ...     reactant_b_type="Carbamate",
        ...     catalyst_filter="Cu",
        ...     top_k=10
        ... )
        {
            "success": True,
            "reactant_pair": "ArI + Carbamate",
            "filters": {"catalyst": "Cu"},
            "total_experiments": 736,
            "total_conditions": 45,
            "conditions": [
                {
                    "rank": 1,
                    "catalyst": "CuI",
                    "ligand": "L1",
                    "base": "Cs2CO3",
                    "solvent": "Dioxane",
                    "experiments": 24,
                    "avg_yield": 68.5,
                    "success_rate": 75.0
                },
                ...
            ]
        }
    """
    if not HTE_AVAILABLE:
        return _error_response(
            "HTE tools not available. Install with: pip install chemtools[hte]",
            {"recommendation": "Install HTE dependencies"}
        )
    
    try:
        from chemtools.HTE import HTERecommender
        
        recommender = HTERecommender()
        df = recommender.df.copy()
        
        # Apply filters
        # Filter by reactant types
        df = df[
            (df['Reactant_A_Type'] == reactant_a_type) &
            (df['Reactant_B_Type'] == reactant_b_type)
        ]
        
        if len(df) == 0:
            return _error_response(
                f"No experiments found for {reactant_a_type} + {reactant_b_type}",
                {
                    "suggestion": "Check reactant type names. Use hte_analytics_tool with query_type='list_pairs' to see available pairs."
                }
            )
        
        # Normalize reaction type if provided
        if reaction_type:
            reaction_type_map = {
                'c-n': 'C_N', 'cn': 'C_N',
                'c-n coupling': 'C_N_Coupling', 'cn coupling': 'C_N_Coupling',
                'c-o': 'CO', 'co': 'CO',
                'c-o coupling': 'CO-Coupling', 'co coupling': 'CO-Coupling',
                'ullmann': 'CO-Coupling',
                'buchwald': 'Buchwald', 'buchwald-hartwig': 'Buchwald',
            }
            reaction_lower = reaction_type.lower().strip()
            if reaction_lower in reaction_type_map:
                reaction_type = reaction_type_map[reaction_lower]
            
            df = df[df['Reaction_Type_Standardized'].str.contains(reaction_type, case=False, na=False)]
        
        # Filter by catalyst if provided
        if catalyst_filter:
            metal_map = {
                'palladium': 'Pd', 'copper': 'Cu', 'nickel': 'Ni',
                'iridium': 'Ir', 'rhodium': 'Rh', 'platinum': 'Pt'
            }
            filter_lower = catalyst_filter.lower()
            search_term = metal_map.get(filter_lower, catalyst_filter)
            df = df[df['Catalyst'].str.contains(search_term, case=False, na=False)]
        
        if len(df) == 0:
            return _error_response(
                f"No experiments found with applied filters",
                {
                    "reactant_pair": f"{reactant_a_type} + {reactant_b_type}",
                    "filters": {
                        "reaction_type": reaction_type,
                        "catalyst": catalyst_filter
                    },
                    "suggestion": "Try removing filters or use different catalyst type."
                }
            )
        
        # Group by condition components
        condition_cols = ['Catalyst', 'Ligand', 'Base', 'Solvent']
        
        # Calculate statistics for each condition
        grouped = df.groupby(condition_cols).agg({
            'AREA_TOTAL_REDUCED': ['count', 'mean', 'median', 'std']
        }).reset_index()
        
        grouped.columns = condition_cols + ['Count', 'Avg_Yield', 'Median_Yield', 'Std_Yield']
        
        # Calculate success rate (yield > 50)
        success_rates = df.groupby(condition_cols)['AREA_TOTAL_REDUCED'].apply(
            lambda x: (x > 50).sum() / len(x) * 100
        ).reset_index(name='Success_Rate')
        
        grouped = grouped.merge(success_rates, on=condition_cols)
        
        # Filter by min_experiments
        grouped = grouped[grouped['Count'] >= min_experiments]
        
        # Sort based on user preference
        if sort_by == "success":
            grouped = grouped.sort_values(['Success_Rate', 'Count'], ascending=[False, False])
        elif sort_by == "yield":
            grouped = grouped.sort_values(['Avg_Yield', 'Count'], ascending=[False, False])
        else:  # count (default)
            grouped = grouped.sort_values('Count', ascending=False)
        
        # Format results
        conditions = []
        for i, (_, row) in enumerate(grouped.head(top_k).iterrows(), 1):
            condition = {
                "rank": i,
                "catalyst": row['Catalyst'],
                "ligand": row['Ligand'] if pd.notna(row['Ligand']) else None,
                "base": row['Base'] if pd.notna(row['Base']) else None,
                "solvent": row['Solvent'] if pd.notna(row['Solvent']) else None,
                "experiments": int(row['Count']),
                "avg_yield": round(row['Avg_Yield'], 1),
                "median_yield": round(row['Median_Yield'], 1),
                "success_rate": round(row['Success_Rate'], 1)
            }
            
            # Add optional components from original data
            condition_df = df[
                (df['Catalyst'] == row['Catalyst']) &
                (df['Ligand'] == row['Ligand']) &
                (df['Base'] == row['Base']) &
                (df['Solvent'] == row['Solvent'])
            ]
            
            # Check for secondary solvent, additive, coupling reagent
            if 'Secondary Solvent' in condition_df.columns:
                sec_solv = condition_df['Secondary Solvent'].dropna().mode()
                if len(sec_solv) > 0:
                    condition['secondary_solvent'] = sec_solv.iloc[0]
            
            if 'Additive' in condition_df.columns:
                additive = condition_df['Additive'].dropna().mode()
                if len(additive) > 0:
                    condition['additive'] = additive.iloc[0]
            
            conditions.append(condition)
        
        return _success_response({
            "reactant_pair": f"{reactant_a_type} + {reactant_b_type}",
            "filters": {
                "reaction_type": reaction_type,
                "catalyst": catalyst_filter
            },
            "total_experiments": len(df),
            "total_conditions": len(grouped),
            "conditions": conditions,
            "sort_by": sort_by
        })
    
    except Exception as e:
        return _error_response(f"Failed to query HTE conditions: {str(e)}")


@tool
def hte_screening_set_tool(
    reactant_a_smiles: Optional[str] = None,
    reactant_b_smiles: Optional[str] = None,
    reaction_smiles: Optional[str] = None,
    num_conditions: int = 24,
    min_experiments: int = 1,
    reaction_type_filter: Optional[str] = None,
    catalyst_filter: Optional[str] = None,
    diversity_strategy: str = "balanced"
) -> Dict[str, Any]:
    """
    Generate a diverse set of conditions for HTE screening plates (up to 24 conditions).
    
    **PRIMARY HTE USE CASE**: Generate a group of diverse conditions to test in parallel
    on a screening plate (e.g., 4x6 = 24 wells, 8x12 = 96 wells).
    
    This tool is specifically designed for:
    - Planning HTE screening experiments
    - Generating diverse condition sets for parallel testing
    - Balancing top performers with exploratory conditions
    - Optimizing plate layouts with reagent diversity
    
    **Key Features:**
    - Generates up to 24 conditions (standard 4x6 plate) or custom amounts
    - Three diversity strategies: balanced, top_performers, diverse
    - Ensures reagent variation across catalyst, ligand, base, solvent
    - Based on 66,308 experiments with z-score ranking
    
    **Diversity Strategies:**
    - "balanced" (DEFAULT): Mix of top ~8 performers + 16 diverse alternatives
    - "top_performers": Focus on highest z-score conditions (best for optimization)
    - "diverse": Maximize reagent diversity (best for broad exploration)
    
    **Use when:**
    - Setting up HTE screening plates
    - Need 12-24 conditions for parallel testing
    - Want to explore reagent space systematically
    - Designing hit-finding experiments
    
    Args:
        reactant_a_smiles: SMILES of first reactant
        reactant_b_smiles: SMILES of second reactant (optional)
        reaction_smiles: Complete reaction SMILES (alternative to separate reactants)
        num_conditions: Number of conditions to generate (default 24 for 4x6 plate)
        min_experiments: Minimum experiments for a condition to be included (default 1)
        reaction_type_filter: Optional filter for reaction type (e.g., 'C_N_Coupling', 'Suzuki')
        catalyst_filter: Optional filter by metal (e.g., 'Pd', 'Cu', 'Ni')
        diversity_strategy: Selection strategy ('balanced', 'top_performers', 'diverse')
    
    Returns:
        Dict with diverse condition set including z-scores, yields, and statistics
    
    Example:
        >>> hte_screening_set_tool(
        ...     reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        ...     num_conditions=24,
        ...     reaction_type_filter="C_N_Coupling",
        ...     catalyst_filter="Cu",
        ...     diversity_strategy="balanced"
        ... )
        {
            "success": True,
            "reactant_a_type": "ArBr",
            "reactant_b_type": "ArNH2",
            "matching_experiments": 112,
            "num_conditions": 24,
            "diversity_strategy": "balanced",
            "screening_conditions": [
                {
                    "rank": 1,
                    "plate_position": "A1",
                    "catalyst": "CuI",
                    "ligand": "PPBO",
                    "base": "NaOtBu",
                    "solvent": "tAmOH",
                    "avg_z_score": 2.61,
                    "avg_yield": 49.8,
                    "confidence_score": 63.8
                },
                ...
            ]
        }
    """
    if not HTE_AVAILABLE:
        return _error_response(
            "HTE tools not available. Install with: pip install chemtools[hte]",
            {"recommendation": "Install HTE dependencies"}
        )
    
    try:
        # Parse reaction SMILES if provided
        if reaction_smiles:
            from chemtools.featurizers.analysis.smiles import _split_reaction_smiles
            
            parts = _split_reaction_smiles(reaction_smiles.strip())
            reactants_str = parts[0] if len(parts) > 0 else ""
            
            if not reactants_str:
                return _error_response(
                    "Could not extract reactants from reaction SMILES. Check format (should be reactants>>products)."
                )
            
            # Split reactants by dot notation
            reactant_list = [r.strip() for r in reactants_str.split(".") if r.strip()]
            
            if len(reactant_list) >= 1:
                reactant_a_smiles = reactant_list[0]
            if len(reactant_list) >= 2:
                reactant_b_smiles = reactant_list[1]
            
            if len(reactant_list) > 2:
                print(f"Note: Reaction has {len(reactant_list)} reactants, using first two for HTE lookup.")
        
        # Validate we have at least one reactant
        if not reactant_a_smiles:
            return _error_response(
                "Must provide either reactant_a_smiles or reaction_smiles."
            )
        
        recommender = HTERecommender()
        
        result = recommender.generate_screening_set(
            reactant_a_smiles=reactant_a_smiles,
            reactant_b_smiles=reactant_b_smiles,
            num_conditions=num_conditions,
            min_experiments=min_experiments,
            reaction_type_filter=reaction_type_filter,
            catalyst_filter=catalyst_filter,
            diversity_strategy=diversity_strategy
        )
        
        # Check if we have any recommendations
        if result.total_matching_experiments == 0 or len(result.recommendations) == 0:
            error_msg = (
                f"No HTE data found for reactant combination: {result.reactant_a_type} + {result.reactant_b_type}"
            )
            
            filters_applied = []
            if reaction_type_filter:
                filters_applied.append(f"reaction type '{reaction_type_filter}'")
            if catalyst_filter:
                filters_applied.append(f"catalyst '{catalyst_filter}'")
            
            if filters_applied:
                error_msg += f" with {' and '.join(filters_applied)}"
            
            suggestion = (
                f"The reactants were successfully classified as {result.reactant_a_type} and {result.reactant_b_type}, "
                f"but this specific combination has no matching experiments in the HTE database. "
                f"Try: (1) removing filters, (2) using hte_analytics_tool to explore available pairs, "
                f"or (3) checking if the reactant types are common in the database."
            )
            
            return _error_response(
                error_msg,
                {
                    "reactant_a_type": result.reactant_a_type,
                    "reactant_b_type": result.reactant_b_type,
                    "matching_experiments": 0,
                    "suggestion": suggestion
                }
            )
        
        # Convert to screening format with plate positions
        # Generate plate positions (A1-D6 for 24-well, can extend to A1-H12 for 96-well)
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        cols = list(range(1, 13))  # 1-12
        
        plate_positions = []
        for row in rows:
            for col in cols:
                plate_positions.append(f"{row}{col}")
                if len(plate_positions) >= num_conditions:
                    break
            if len(plate_positions) >= num_conditions:
                break
        
        screening_conditions = []
        for i, rec in enumerate(result.recommendations, 1):
            condition = {
                "rank": i,
                "plate_position": plate_positions[i-1] if i <= len(plate_positions) else None,
                "catalyst": rec.catalyst,
                "ligand": rec.ligand,
                "base": rec.base,
                "solvent": rec.solvent,
                "avg_z_score": round(rec.avg_z_score, 2),
                "avg_yield": round(rec.avg_yield, 1),
                "median_yield": round(rec.median_yield, 1),
                "success_rate": round(rec.success_rate, 1),
                "num_experiments": rec.num_experiments,
                "confidence_score": round(rec.confidence_score, 1)
            }
            
            # Add optional components
            if rec.secondary_solvent:
                condition["secondary_solvent"] = rec.secondary_solvent
            if rec.additive:
                condition["additive"] = rec.additive
            if rec.coupling_reagent:
                condition["coupling_reagent"] = rec.coupling_reagent
            
            screening_conditions.append(condition)
        
        return _success_response({
            "reactant_a_type": result.reactant_a_type,
            "reactant_b_type": result.reactant_b_type,
            "reactant_a_smiles": reactant_a_smiles,
            "reactant_b_smiles": reactant_b_smiles,
            "predicted_reaction_type": result.predicted_reaction_type,
            "reaction_confidence": round(result.reaction_type_confidence * 100, 1),
            "matching_experiments": result.total_matching_experiments,
            "num_conditions": len(screening_conditions),
            "diversity_strategy": diversity_strategy,
            "screening_conditions": screening_conditions,
            "plate_format": f"{len([p for p in plate_positions if p[0] <= 'D'])//6}x6" if num_conditions <= 24 else "8x12"
        })
    
    except Exception as e:
        return _error_response(f"HTE screening set generation failed: {str(e)}")


# ============================================================================
# Tool Collection
# ============================================================================

CHEMTOOLS_TOOLS = [
    # SMILES tools
    normalize_smiles_tool,
    normalize_reaction_tool,
    
    # Analysis tools
    detect_reaction_family_tool,
    classify_reactant_tool,
    molpipeline_featurize_tool,
    get_functional_groups_tool,
    summarize_functional_groups_tool,  # NEW: Human-readable functional group summary
    calculable_features_tool,
    analyze_bond_changes_tool,  # Bond breaking/formation analysis
    identify_reaction_center_tool,  # NEW: Reaction center identification from atom mapping
    reaction_similarity_tool,  # NEW: DRFP-based reaction comparison
    
    # Recommendation tools
    recommend_conditions_tool,
    rule_based_conditions_tool,
    search_precedents_tool,
    protocol_recommendation_tool,  # NEW: Protocol-based recommendations with full procedures
    unified_recommender_tool,  # NEW: Unified DRFP-based protocol + rule search
    rule_builder_autofill_tool,  # NEW: LLM-assisted rule database drafting
    hte_recommend_tool,  # NEW: HTE-based condition recommendations with catalyst filtering
    hte_analytics_tool,  # NEW: HTE database analytics (pairs, catalysts, metals)
    hte_conditions_tool,  # NEW: Detailed conditions for specific substrate pairs
    hte_screening_set_tool,  # NEW: Generate diverse condition sets for HTE screening plates (up to 24)
    
    # Database tools
    reaction_dataset_analytics_tool,
    find_reagent_tool,
    reagent_database_analytics_tool,
    list_supported_cores_tool,
    list_all_families_tool,  # NEW: List all available reaction families
    add_reagent_tool,
    auto_conditions_llm_tool,
    planner_detect_family_tool,
    planner_rule_candidates_tool,
    planner_protocol_candidates_tool,
    planner_hte_summary_tool,
    planner_score_candidates_tool,
    planner_fuse_tool,
]


# ============================================================================
# Helper Functions
# ============================================================================

def _success_response(data: Any) -> Dict[str, Any]:
    """Return a standardized success payload merged with tool output."""
    payload: Dict[str, Any] = {"success": True}
    if isinstance(data, dict):
        payload.update(data)
    else:
        payload["data"] = data
    return payload


def _error_response(message: str, extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Return a standardized error payload with optional contextual fields."""
    payload: Dict[str, Any] = {"success": False, "error": message}
    if extra:
        if isinstance(extra, dict):
            payload.update(extra)
        else:
            payload["details"] = extra
    return payload


def _schema_properties(tool_obj: Any) -> Dict[str, Any]:
    """Return JSON-schema properties for a tool's input model."""
    args_schema = getattr(tool_obj, "args_schema", None)
    if not args_schema:
        return {}
    try:
        return args_schema.model_json_schema().get("properties", {})
    except AttributeError:
        # Pydantic <2 compatibility
        return args_schema.schema().get("properties", {})


def _schema_required(tool_obj: Any) -> List[str]:
    """Return required parameter names for a tool."""
    args_schema = getattr(tool_obj, "args_schema", None)
    if not args_schema:
        return []
    try:
        return list(args_schema.model_json_schema().get("required", []))
    except AttributeError:
        return list(args_schema.schema().get("required", []))


def get_tool_descriptions() -> List[Dict[str, Any]]:
    """Get structured descriptions (name, docstring, parameters) for tools."""
    descriptions: List[Dict[str, Any]] = []
    for tool_obj in CHEMTOOLS_TOOLS:
        properties = _schema_properties(tool_obj)
        required = set(_schema_required(tool_obj))
        params: List[Dict[str, Any]] = []
        for param_name, metadata in properties.items():
            params.append({
                "name": param_name,
                "required": param_name in required,
                "description": metadata.get("description", ""),
                "type": metadata.get("type"),
                "default": metadata.get("default", None),
            })
        descriptions.append({
            "name": tool_obj.name,
            "description": tool_obj.description or "",
            "parameters": params,
        })
    return descriptions


def print_tool_summary():
    """Print a summary of available tools."""
    print("=" * 70)
    print("ChemTools LangChain Tools")
    print("=" * 70)
    descriptions = get_tool_descriptions()
    for i, entry in enumerate(descriptions, 1):
        print(f"\n{i}. {entry['name']}")
        desc_line = (entry["description"] or "").split("\n")[0].strip()
        print(f"   {desc_line if desc_line else 'No description'}")
        params = entry.get("parameters") or []
        if not params:
            print("   Parameters: none")
        else:
            print("   Parameters:")
            for param in params:
                requirement = "required" if param.get("required") else "optional"
                param_type = param.get("type") or "any"
                print(f"     - {param['name']} ({param_type}, {requirement})")
                if param.get("description"):
                    print(f"       {param['description']}")
    print("\n" + "=" * 70)


if __name__ == "__main__":
    # Print tool summary when run directly
    print_tool_summary()


# ============================================================================
# Internal helpers (kept at end to avoid cluttering the main tool definitions)
# ============================================================================

def _normalise_synonyms_field(value: Any) -> Optional[List[str]]:
    """Normalize synonyms input from the tool into a list of strings."""
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        try:
            parsed = json.loads(stripped)
        except json.JSONDecodeError:
            parsed = None
        if isinstance(parsed, (list, tuple, set)):
            return [str(item).strip() for item in parsed if str(item).strip()]
        return [part.strip() for part in stripped.split(",") if part.strip()]
    if isinstance(value, (list, tuple, set)):
        return [str(item).strip() for item in value if str(item).strip()]
    coerced = str(value).strip()
    return [coerced] if coerced else None


def _simplify_recommendation(raw: Dict[str, Any]) -> Dict[str, Any]:
    """Extract the subset of recommendation fields most useful to the LLM."""
    reasons_src = raw.get("reasons")
    reasons_list: List[str] = []
    if isinstance(reasons_src, dict):
        reasons_list = [str(v) for v in reasons_src.values()]
    elif isinstance(reasons_src, (list, tuple)):
        reasons_list = [str(v) for v in reasons_src]
    elif reasons_src:
        reasons_list = [str(reasons_src)]

    # Get base recommendation
    recommendation = dict(raw.get("recommendation", {}) or {})
    
    # Resolve base and solvent names from UIDs
    base_uid = recommendation.get("base_uid")
    solvent_uid = recommendation.get("solvent_uid")
    
    if base_uid and not recommendation.get("base"):
        try:
            base_info = find_reagent(str(base_uid), "base")
            if isinstance(base_info, dict) and not base_info.get("error"):
                recommendation["base"] = base_info.get("name") or base_info.get("abbreviation")
        except Exception:
            pass
    
    if solvent_uid and not recommendation.get("solvent"):
        try:
            solvent_info = find_reagent(str(solvent_uid), "solvent")
            if isinstance(solvent_info, dict) and not solvent_info.get("error"):
                recommendation["solvent"] = solvent_info.get("name") or solvent_info.get("abbreviation")
        except Exception:
            pass
    
    # For amide formation reactions, extract coupling reagents and additives from precedents
    family = str(raw.get("family", "")).lower()
    is_amide = any(term in family for term in ["amide", "amidation"])
    
    if is_amide and not recommendation.get("coupling_reagent"):
        # Extract reagent information from top precedents
        precedent_pack = raw.get("precedent_pack", {})
        precedents = precedent_pack.get("precedents", [])
        
        if precedents:
            # Analyze top precedents to find most common coupling reagents and additives
            coupling_reagents = {}
            additives = {}
            
            for prec in precedents[:10]:  # Check top 10 precedents
                reagents = prec.get("reagents", [])
                if isinstance(reagents, list):
                    for r in reagents:
                        role = r.get("role", "").upper()
                        name = r.get("name") or r.get("abbreviation") or r.get("original_name", "")
                        
                        if role == "COUPLING_REAGENT" and name:
                            coupling_reagents[name] = coupling_reagents.get(name, 0) + 1
                        elif role == "ADDITIVE" and name:
                            additives[name] = additives.get(name, 0) + 1
            
            # Add most common coupling reagent to recommendation
            if coupling_reagents:
                top_coupling = max(coupling_reagents.items(), key=lambda x: x[1])[0]
                recommendation["coupling_reagent"] = top_coupling
            
            # Add most common additives to recommendation
            if additives:
                top_additives = sorted(additives.items(), key=lambda x: x[1], reverse=True)
                recommendation["additives"] = ", ".join([name for name, _ in top_additives[:2]])

    return {
        "family": raw.get("family", "Unknown"),
        "detected_family": raw.get("detected_family"),
        "search_all_families": raw.get("search_all_families", False),
        "recommendation": recommendation,
        "alternatives": raw.get("alternatives", {}) or {},
        "reasons": reasons_list[:5],
        "precedent_count": len(raw.get("precedent_pack", {}).get("precedents", [])),
        "detection": raw.get("detection", {}),
        "constraint_filters": raw.get("constraint_filters"),
    }


def _apply_core_constraints(
    simplified: Dict[str, Any],
    spec: ConstraintSpec,
) -> Tuple[Dict[str, Any], Optional[str]]:
    """
    Apply catalyst-core constraints using the structured specification.
    
    Returns the updated simplified dict and a note describing any adjustments.
    """
    if not spec.allow_metals and not spec.exclude_metals and not spec.prefer_metals:
        return simplified, None

    recommendation = dict(simplified.get("recommendation") or {})
    if not recommendation:
        return simplified, None

    base_core = recommendation.get("core")

    alt_list = simplified.get("alternatives", {}) or {}
    core_counts: Sequence[Tuple[str, int]] = alt_list.get("cores", []) or []
    ordered_cores: List[str] = []

    if base_core:
        ordered_cores.append(base_core)
    ordered_cores.extend(
        core for core, _ in core_counts
        if core and core not in ordered_cores
    )

    updated_order = filter_cores_by_constraints(
        ordered_cores,
        allow_metals=spec.allow_metals,
        exclude_metals=spec.exclude_metals,
        prefer_metals=spec.prefer_metals,
    )

    if not updated_order:
        return simplified, None

    new_core = updated_order[0]
    note_bits: List[str] = []

    if new_core != base_core:
        recommendation["core"] = new_core
        note_bits.append(f"core updated to {new_core}")

    simplified["recommendation"] = recommendation
    simplified["alternatives"] = dict(alt_list)
    simplified["alternatives"]["cores"] = [
        (core, next((count for c, count in core_counts if c == core), None))
        for core in updated_order
    ]

    if spec.allow_metals:
        note_bits.append(f"allowed metals: {', '.join(sorted(spec.allow_metals))}")
    if spec.exclude_metals:
        note_bits.append(f"excluded metals: {', '.join(sorted(spec.exclude_metals))}")
    if spec.prefer_metals:
        note_bits.append(f"preferred metals: {', '.join(sorted(spec.prefer_metals))}")

    note = "; ".join(note_bits) if note_bits else None
    return simplified, note


def run_rule_builder_autofill(params: RuleBuilderAutoInput) -> Dict[str, Any]:
    """Generate a rule database draft from protocol text programmatically."""
    builder = RuleBuilder.new(params.family)
    tags = params.tags if params.tags else [params.family.lower()]
    created_date = params.created_date or date.today().isoformat()
    builder.set_metadata(
        id=params.metadata_id,
        name=params.metadata_name,
        version=params.metadata_version,
        created_date=created_date,
        status=params.status or "draft",
        tags=tags,
    )
    builder.add_reference_reactions(params.reference_reactions)
    if params.notes:
        builder.set_notes(params.notes.strip())

    focus_bits: List[str] = []
    if params.desired_focus:
        focus_bits.append(f"Focus: {params.desired_focus}")
    if params.applies_if_hints:
        focus_bits.append(
            "Required applies_if features: " + ", ".join(params.applies_if_hints)
        )
    if params.modifier_hints:
        focus_bits.append(
            "Modifier triggers to include: " + ", ".join(params.modifier_hints)
        )
    focus_block = "\n".join(focus_bits) if focus_bits else "None specified"
    reference_block = "\n".join(f"- {rxn}" for rxn in params.reference_reactions) or "- None provided"

    prompt = RULE_BUILDER_EXTRACTION.format(
        family=params.family,
        reference_block=reference_block,
        protocol_text=params.protocol_text.strip(),
        focus=focus_block,
        max_base_rules=params.max_base_rules,
    )

    try:
        llm_response = _invoke_rule_builder_llm(prompt)
    except Exception as exc:  # pragma: no cover - network failures handled at runtime
        return _error_response(f"LLM call failed: {exc}")

    payload_text = _strip_json_fences(llm_response)
    try:
        structured = json.loads(payload_text)
    except json.JSONDecodeError as exc:
        return _error_response(
            "LLM response could not be parsed as JSON.",
            {"raw_response": llm_response, "parse_error": str(exc)},
        )

    if not isinstance(structured, dict):
        return _error_response(
            "LLM response did not produce a JSON object.",
            {"raw_response": llm_response},
        )

    _apply_rule_builder_payload(builder, structured, params)
    if params.applies_if_hints:
        _ensure_applies_if_hints(builder, params.applies_if_hints)

    issues = builder.validate(strict=False)
    result = {
        "rule_database": builder.to_dict(),
        "issues": _issues_to_dicts(issues),
    }
    result["message"] = "Draft created" if not issues else "Draft created with warnings"
    return _success_response(result)


def _apply_rule_builder_payload(
    builder: RuleBuilder,
    payload: Dict[str, Any],
    params: RuleBuilderAutoInput,
) -> None:
    """Map LLM output into the deterministic builder."""
    scope = payload.get("scope")
    if isinstance(scope, dict):
        builder.set_scope(
            scope_type=scope.get("scope_type"),
            compatible=scope.get("compatible_functional_groups"),
            incompatible=scope.get("incompatible_functional_groups"),
        )

    notes = payload.get("notes")
    if notes and not params.notes:
        builder.set_notes(str(notes))

    applies = payload.get("applies_if")
    if isinstance(applies, dict):
        builder.set_applies_if(raw=applies)

    default_rule = payload.get("default_rule") or {}
    conditions = default_rule.get("conditions") or {}
    if conditions:
        builder.set_default_rule(
            rule_id=default_rule.get("id"),
            description=default_rule.get("description"),
            conditions=conditions,
        )

    base_rules = payload.get("base_rules") or []
    for idx, rule in enumerate(base_rules, 1):
        if not isinstance(rule, dict):
            continue
        rule_id = rule.get("id") or f"auto_rule_{idx}"
        builder.upsert_base_rule(
            rule_id,
            name=rule.get("name") or rule_id,
            description=rule.get("description") or "",
            reactant_features=rule.get("reactant_features") or {},
            conditions=rule.get("conditions") or {},
            priority=rule.get("priority"),
        )

    modifiers = payload.get("modifiers") or []
    for idx, modifier in enumerate(modifiers, 1):
        if not isinstance(modifier, dict):
            continue
        mod_id = modifier.get("id") or f"auto_modifier_{idx}"
        when = modifier.get("when") or []
        suggestion = modifier.get("suggest") or modifier.get("suggestion")
        if not when or not suggestion:
            continue
        builder.upsert_modifier(
            mod_id,
            when=when,
            suggestion=suggestion,
            rationale=modifier.get("rationale"),
        )

    # Fallback: ensure at least one base rule exists
    if not builder.data.get("base_rules"):
        default_conditions = (builder.data.get("default_rule") or {}).get("conditions") or {}
        if default_conditions:
            builder.upsert_base_rule(
                "auto_base_rule",
                name="Auto Base Rule",
                description="Generated from default rule conditions.",
                reactant_features=builder.data.get("applies_if") or {},
                conditions=default_conditions,
            )


def _ensure_applies_if_hints(builder: RuleBuilder, hints: List[str]) -> None:
    applies = builder.data.get("applies_if") or {}
    existing_all = applies.get("all") or []
    merged = list(dict.fromkeys(existing_all + hints))
    applies["all"] = merged
    builder.set_applies_if(raw=applies)


def _issues_to_dicts(issues: List[Any]) -> List[Dict[str, Any]]:
    return [
        {
            "field": getattr(issue, "field", ""),
            "message": getattr(issue, "message", ""),
            "severity": getattr(issue, "severity", "error"),
        }
        for issue in issues
    ]


def _strip_json_fences(text: str) -> str:
    stripped = text.strip()
    if stripped.startswith("```"):
        lines = stripped.splitlines()
        if lines and lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].startswith("```"):
            lines = lines[:-1]
        stripped = "\n".join(lines).strip()
    return stripped


def _invoke_rule_builder_llm(prompt: str) -> str:
    client = _get_rule_builder_llm_client()
    response = client.chat(
        prompt=prompt,
        system=RULE_BUILDER_SYSTEM_PROMPT,
        temperature=0.1,
        max_tokens=2200,
    )
    return response.content.strip()


def _get_rule_builder_llm_client() -> LLMClient:
    provider = os.getenv("LLM_PROVIDER", "openai")
    model = os.getenv("LLM_MODEL")
    if not model:
        model = (
            RECOMMENDED_MODELS.get(provider, {}).get("reasoning")
            or RECOMMENDED_MODELS.get(provider, {}).get("balanced")
            or "gpt-4o"
        )
    return LLMClient(provider=provider, model=model, temperature=0.1, max_tokens=2000)
