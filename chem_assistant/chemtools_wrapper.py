"""
LangChain tool wrappers for ChemTools functions.

This module exposes existing chemtools functionality as LangChain tools
without modifying the original chemtools codebase.

Available Tools:
    - normalize_smiles_tool: Canonicalize SMILES strings
    - normalize_reaction_tool: Canonicalize reaction SMILES
    - detect_reaction_family_tool: Detect reaction family/type
    - recommend_conditions_tool: Get ML-based condition recommendations
    - rule_based_conditions_tool: Deterministic rule-engine condition guidance (NEW)
    - enhanced_cross_family_recommend_tool: Cross-family precedent search with mechanism filters
    - search_precedents_tool: Search for similar precedent reactions
    - reaction_dataset_analytics_tool: Analyze reaction dataset frequency/yield statistics
    - find_reagent_tool: Look up reagent information from database
    - reagent_database_analytics_tool: Summarize reagent registry statistics
    - list_supported_cores_tool: Enumerate catalyst cores observed in precedents
    - add_reagent_tool: Insert or preview reagent taxonomy entries
    - classify_reactant_tool: Classify reactant type (aryl halide, amine, etc.)
    - get_functional_groups_tool: Detect functional groups in a molecule
    - calculable_features_tool: Evaluate curated calculable feature library for a molecule
    - molpipeline_featurize_tool: Generate molecular features with optional MolPipeline vectors
    - analyze_bond_changes_tool: Analyze bond breaking/formation in reactions (NEW)

Usage:
    from lang_chain.chemtools_wrapper import CHEMTOOLS_TOOLS
    from langgraph.prebuilt import create_react_agent
    
    agent = create_react_agent(llm, CHEMTOOLS_TOOLS)
"""

from typing import Dict, Any, List, Optional, Sequence, Tuple, Union, Literal
import json
import time
from dataclasses import dataclass
from collections import OrderedDict
from pathlib import Path
from pydantic import BaseModel, Field
from langchain_core.tools import tool

# Import chemtools functions
from chemtools.smiles import normalize, normalize_reaction
from chemtools import detect_reaction
from chemtools.output_formatter import ensure_standard_output
from chemtools.recommend.modules import recommend_from_reaction
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
from chemtools.util.functional_groups import detect_all as detect_functional_groups
from chemtools.featurizers import molecular as molecular_featurizer
from chemtools.featurizers import calculable as calculable_features

# Rule-based recommendation engine
from chemtools.rule import RuleEngine

# Import bond analysis tools (NEW)
from chemtools import (
    analyze_bond_changes,
    analyze_bond_changes_hybrid,
    rxnmapper_available,
)

REAGENT_RESOLVER_TIMEOUT = 6.0

from .constraint_parser import (
    ConstraintSpec,
    build_constraint_spec,
    filter_cores_by_constraints,
    format_constraints_for_prompt,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent
_RULE_DB_SEARCH_PATHS: List[Path] = [
    _REPO_ROOT / "data" / "rule_db",
    _REPO_ROOT / "data",
    _REPO_ROOT,
    Path.cwd() / "data" / "rule_db",
    Path.cwd() / "data",
    Path.cwd(),
]

_FAMILY_TO_RULE_DB = {
    "cn_coupling": "buchwald_cn",
    "c_n_coupling": "buchwald_cn",
    "buchwald_cn": "buchwald_cn",
    "buchwald_hartwig": "buchwald_cn",
    "buchwald_hartwig_cn": "buchwald_cn",
    "buchwald_hartwig_c_n": "buchwald_cn",
    "buchwald-hartwig": "buchwald_cn",
    "suzuki": "suzuki",
    "suzuki_coupling": "suzuki",
    "suzuki_miyaura": "suzuki",
    "suzuki-miyaura": "suzuki",
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


class EnhancedCrossFamilyInput(BaseModel):
    """Schema for enhanced cross-family recommendation tool."""

    reaction_smiles: str = Field(
        ...,
        description="Reaction SMILES string for cross-family comparison (reactants>>products).",
    )
    k: int = Field(
        50,
        ge=1,
        le=250,
        description="Number of cross-family precedents to retrieve.",
    )
    reaction_type_threshold: float = Field(
        0.15,
        ge=0.0,
        le=1.0,
        description="Minimum fraction of precedents a reaction type must represent.",
    )
    mechanism_similarity_threshold: float = Field(
        0.4,
        ge=0.0,
        le=1.0,
        description="Minimum allowed mechanism similarity for precedents.",
    )
    mechanism_weight: float = Field(
        0.3,
        ge=0.0,
        le=1.0,
        description="Weighting for mechanism similarity in cross-family scoring.",
    )
    filter_unknown_reagents: bool = Field(
        True,
        description="Exclude precedents missing base/solvent annotations.",
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
        payload = {
            "family": result.get("family"),
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
        features = molecular_featurizer.featurize(
            electrophile,
            nucleophile,
            include_molpipeline=include_molpipeline,
        )
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
    defined in ``chemtools/featurizers/calculable_features.json``. When ``feature_tokens``
    is supplied the output is restricted to those tokens (missing entries are
    reported separately). Setting ``only_present`` filters the result to features
    that evaluate to ``True`` or positive integers.
    """
    try:
        all_features = calculable_features.detect_all_features(smiles)
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
                or (isinstance(value, int) and value > 0)
            }

        present_tokens = calculable_features.get_present_features(smiles)

        steric_tokens = {
            "ortho_substitution_present": bool(all_features.get("ortho_substitution_present")),
            "tert_butyl_present": bool(all_features.get("tert_butyl_present")),
            "isopropyl_present": bool(all_features.get("isopropyl_present")),
        }
        # Lightweight textual fallback when RDKit unavailable
        if not steric_tokens["isopropyl_present"] and "C(C)C" in smiles:
            steric_tokens["isopropyl_present"] = True
        if not steric_tokens["tert_butyl_present"] and "C(C)(C)C" in smiles:
            steric_tokens["tert_butyl_present"] = True
        steric_indicators = [
            token for token, present in steric_tokens.items() if present
        ]
        steric_level = "low"
        if steric_tokens["ortho_substitution_present"] or len(steric_indicators) >= 2:
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
            payload["summary"] = calculable_features.feature_summary(smiles)
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


# ============================================================================
# Recommendation Tools
# ============================================================================

def _build_cross_family_relax_config(
    reaction_type_threshold: float,
    mechanism_similarity_threshold: float,
    mechanism_weight: float,
) -> Dict[str, float]:
    """Compose the relax configuration for enhanced cross-family search."""
    return {
        "reaction_type_threshold": reaction_type_threshold,
        "mechanism_similarity_threshold": mechanism_similarity_threshold,
        "mechanism_weight": mechanism_weight,
        "use_drfp": True,
    }


def _extract_cross_family_metrics(
    recommended_conditions: Sequence[Dict[str, Any]]
) -> Dict[str, Any]:
    """Summarize family diversity and mechanism compatibility across precedents."""
    if not recommended_conditions:
        return {}

    family_counts: Dict[str, int] = {}
    mechanism_scores: List[float] = []
    compatibility_stats = {
        "compatible": 0,
        "moderate": 0,
        "incompatible": 0,
        "well_represented": 0,
        "underrepresented": 0,
    }

    for rec in recommended_conditions:
        cond_section = rec.get("conditions") if isinstance(rec, dict) else {}
        if not isinstance(cond_section, dict):
            cond_section = {}
        cross_meta = cond_section.get("cross_family_metadata")
        if not isinstance(cross_meta, dict):
            cross_meta = {}

        family = (
            cond_section.get("reaction_family")
            or cond_section.get("rxn_type")
            or cross_meta.get("precedent_family")
            or "Unknown"
        )
        family_counts[family] = family_counts.get(family, 0) + 1

        mech_sim = cross_meta.get("mechanism_similarity")
        if isinstance(mech_sim, (int, float)):
            mechanism_scores.append(float(mech_sim))

        mech_status = cross_meta.get("mechanism_status")
        if mech_status in compatibility_stats:
            compatibility_stats[mech_status] += 1

        type_status = cross_meta.get("reaction_type_status")
        if type_status in compatibility_stats:
            compatibility_stats[type_status] += 1

    total = len(recommended_conditions)
    family_diversity = len(family_counts) / total if total else 0.0
    metrics: Dict[str, Any] = {
        "total_recommendations": total,
        "family_diversity_score": round(family_diversity, 3),
        "unique_families": len(family_counts),
        "family_distribution": dict(sorted(family_counts.items(), key=lambda kv: kv[1], reverse=True)),
        "compatibility_breakdown": {
            "mechanism_compatibility": {
                "compatible": compatibility_stats["compatible"],
                "moderate": compatibility_stats["moderate"],
                "incompatible": compatibility_stats["incompatible"],
            },
            "reaction_type_representation": {
                "well_represented": compatibility_stats["well_represented"],
                "underrepresented": compatibility_stats["underrepresented"],
            },
        },
    }

    if mechanism_scores:
        metrics.update({
            "avg_mechanism_similarity": round(sum(mechanism_scores) / len(mechanism_scores), 3),
            "min_mechanism_similarity": round(min(mechanism_scores), 3),
            "max_mechanism_similarity": round(max(mechanism_scores), 3),
        })

    return metrics


def _format_chemical_label(entry: Any) -> Optional[str]:
    """Return the most human-readable identifier for a chemical entry."""
    if isinstance(entry, dict):
        return (
            entry.get("name")
            or entry.get("abbreviation")
            or entry.get("cas")
            or entry.get("smiles")
        )
    if entry is None:
        return None
    return str(entry)


def _format_condition_value(value: Any) -> Optional[str]:
    """Normalize condition values (temperature/time) for summaries."""
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        compact = [str(v) for v in value if v is not None]
        return "-".join(compact) if compact else None
    return str(value)


def _summarize_top_cross_family_condition(
    recommended_conditions: Sequence[Dict[str, Any]]
) -> Dict[str, Any]:
    """Capture key details from the lead cross-family recommendation."""
    if not recommended_conditions:
        return {}

    first = recommended_conditions[0] or {}
    summary = first.get("summary", {}) if isinstance(first, dict) else {}
    conditions = first.get("conditions", {}) if isinstance(first, dict) else {}
    if not isinstance(summary, dict):
        summary = {}
    if not isinstance(conditions, dict):
        conditions = {}
    cross_meta = conditions.get("cross_family_metadata")
    if not isinstance(cross_meta, dict):
        cross_meta = {}

    support = summary.get("support")
    support_count: Optional[int] = None
    if isinstance(support, dict):
        support_count = support.get("count") or support.get("reference_population")
    elif isinstance(support, (int, float)):
        support_count = int(support)
    precedents = summary.get("precedents")
    if support_count is None and isinstance(precedents, list):
        support_count = len(precedents)

    top_payload = {
        "rank": summary.get("rank") or first.get("rank"),
        "core": summary.get("core"),
        "base": _format_chemical_label(summary.get("base")),
        "ligand": _format_chemical_label(summary.get("ligand")),
        "solvent": _format_chemical_label(summary.get("solvent")),
        "temperature": _format_condition_value(conditions.get("temperature")),
        "time": _format_condition_value(conditions.get("time")),
        "confidence": summary.get("confidence"),
        "supporting_precedents": support_count,
        "reaction_family": conditions.get("reaction_family"),
        "precedent_family": cross_meta.get("precedent_family") or cross_meta.get("detected_family"),
        "mechanism_similarity": cross_meta.get("mechanism_similarity") or conditions.get("mechanism_similarity"),
        "mechanism_status": cross_meta.get("mechanism_status") or conditions.get("mechanism_status"),
        "reaction_type_status": cross_meta.get("reaction_type_status") or conditions.get("reaction_type_status"),
    }

    return {key: value for key, value in top_payload.items() if value is not None}


def _generate_cross_family_insights(
    recommended_conditions: Sequence[Dict[str, Any]],
    metrics: Dict[str, Any],
) -> List[str]:
    """Produce concise insights for LLM consumption."""
    if not recommended_conditions:
        return ["No cross-family precedents met the specified thresholds."]

    insights: List[str] = []
    total = metrics.get("total_recommendations") or len(recommended_conditions)
    family_distribution = metrics.get("family_distribution") or {}

    if family_distribution:
        sorted_families = list(family_distribution.items())
        top_family, top_count = sorted_families[0]
        insights.append(f"Top precedent family: {top_family} ({top_count}/{total}).")
        if len(sorted_families) > 1:
            others = ", ".join(
                f"{fam} ({count})" for fam, count in sorted_families[1:3]
            )
            if others:
                insights.append(f"Secondary families represented: {others}.")

    mech_avg = metrics.get("avg_mechanism_similarity")
    if mech_avg is not None:
        mech_min = metrics.get("min_mechanism_similarity", mech_avg)
        mech_max = metrics.get("max_mechanism_similarity", mech_avg)
        insights.append(
            f"Mechanism similarity averages {mech_avg:.2f} "
            f"(range {mech_min:.2f}-{mech_max:.2f})."
        )

    compat = metrics.get("compatibility_breakdown", {})
    mech_compat = compat.get("mechanism_compatibility", {}) if isinstance(compat, dict) else {}
    if mech_compat:
        incompatible = mech_compat.get("incompatible", 0)
        moderate = mech_compat.get("moderate", 0)
        if incompatible:
            insights.append(f"{incompatible} precedent(s) flagged as mechanism incompatible were filtered out.")
        if moderate:
            insights.append(f"{moderate} precedent(s) show moderate mechanism similarity and may need scrutiny.")

    repr_stats = compat.get("reaction_type_representation", {}) if isinstance(compat, dict) else {}
    underrepresented = repr_stats.get("underrepresented", 0)
    if underrepresented:
        insights.append(
            f"{underrepresented} reaction family/families are underrepresented and may require manual validation."
        )

    lead_summary = _summarize_top_cross_family_condition(recommended_conditions)
    if lead_summary:
        components: List[str] = []
        core_val = lead_summary.get("core")
        if core_val:
            components.append(str(core_val))
        core_text_lower = str(core_val).lower() if core_val else ""

        base_val = lead_summary.get("base")
        if base_val:
            base_text = f"Base: {base_val}"
            if "base" not in core_text_lower and not any(str(base_val) in comp for comp in components):
                components.append(base_text)

        ligand_val = lead_summary.get("ligand")
        if ligand_val:
            ligand_text = f"Ligand: {ligand_val}"
            if not any(str(ligand_val) in comp for comp in components):
                components.append(ligand_text)

        solvent_val = lead_summary.get("solvent")
        if solvent_val:
            solvent_text = f"Solvent: {solvent_val}"
            if not any(str(solvent_val) in comp for comp in components):
                components.append(solvent_text)

        if lead_summary.get("temperature"):
            components.append(f"T: {lead_summary['temperature']}")
        if lead_summary.get("time"):
            components.append(f"t: {lead_summary['time']}")
        if components:
            insights.append("Lead cross-family condition -> " + "; ".join(components))

        mech_sim = lead_summary.get("mechanism_similarity")
        if isinstance(mech_sim, (int, float)):
            status = lead_summary.get("mechanism_status", "unknown")
            insights.append(f"Lead precedent mechanism similarity {float(mech_sim):.2f} ({status}).")

    return insights


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


@tool(args_schema=EnhancedCrossFamilyInput)
def enhanced_cross_family_recommend_tool(
    reaction_smiles: str,
    k: int = 50,
    reaction_type_threshold: float = 0.15,
    mechanism_similarity_threshold: float = 0.4,
    mechanism_weight: float = 0.3,
    filter_unknown_reagents: bool = True,
) -> Dict[str, Any]:
    """
    Run the enhanced cross-family recommendation workflow with mechanism-aware filtering.

    Provides DRFP-backed precedent suggestions drawn from all families, while enforcing
    reaction-type representation and mechanism-similarity thresholds. The output mirrors
    the CLI but is tailored for LLM consumption (metrics + concise insights).
    """
    try:
        relax_config = _build_cross_family_relax_config(
            reaction_type_threshold=reaction_type_threshold,
            mechanism_similarity_threshold=mechanism_similarity_threshold,
            mechanism_weight=mechanism_weight,
        )

        raw_result = recommend_from_reaction(
            reaction=reaction_smiles,
            k=k,
            relax=relax_config,
            rerank_strategy="none",
            filter_unknown_reagents=filter_unknown_reagents,
            search_all_families=True,
        )

        formatted_source = raw_result.get("formatted")
        if not isinstance(formatted_source, dict):
            raise ValueError("Recommendation engine did not return formatted output.")

        formatted = ensure_standard_output(
            formatted_source,
            default_model="Enhanced-Cross-Family-Search",
            fallback_reaction_smiles=reaction_smiles,
            extras={"raw_recommendation": raw_result},
        )

        recommended_conditions = formatted.get("recommended_conditions", [])
        metrics = _extract_cross_family_metrics(recommended_conditions)
        if metrics:
            formatted.setdefault("meta", {})["cross_family_metrics"] = metrics

        insights = _generate_cross_family_insights(recommended_conditions, metrics)
        top_summary = _summarize_top_cross_family_condition(recommended_conditions)

        payload = {
            "input_reaction": reaction_smiles,
            "parameters": {
                "k": k,
                "reaction_type_threshold": reaction_type_threshold,
                "mechanism_similarity_threshold": mechanism_similarity_threshold,
                "mechanism_weight": mechanism_weight,
                "filter_unknown_reagents": filter_unknown_reagents,
            },
            "recommended_conditions": recommended_conditions,
            "recommendations": recommended_conditions,
            "top_recommendation": top_summary,
            "cross_family_metrics": metrics,
            "insights": insights,
            "meta": formatted.get("meta", {}),
            "formatted": formatted,
        }
        return _success_response(payload)
    except Exception as exc:
        return _error_response(str(exc), {"reaction_smiles": reaction_smiles})


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
        # Auto-detect family if not provided
        if family is None:
            detection = detect_reaction(reaction_smiles, use_ml=False)
            family = detection.get("family", "Unknown")
        
        # Use relax parameter for DRFP-based search
        result = precedent_knn(
            family=family,
            features={},  # DRFP uses reaction fingerprint, not substrate features
            k=k,
            relax={"use_drfp": True, "reaction_smiles": reaction_smiles}
        )
        
        # Simplify precedent data
        precedents = result.get("precedents", [])
        simplified_precedents = []
        for p in precedents[:k]:
            simplified_precedents.append({
                "similarity": round(p.get("similarity", 0), 3),
                "conditions": {
                    "catalyst": p.get("catalyst", ""),
                    "base": p.get("base", ""),
                    "solvent": p.get("solvent", ""),
                    "temperature_C": p.get("T_C", ""),
                    "time_h": p.get("time_h", "")
                },
                "yield": p.get("yield", ""),
                "reaction": p.get("rxn_smiles", "")[:100]  # Truncate long SMILES
            })
        
        return _success_response({
            "family": family,
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
    calculable_features_tool,
    analyze_bond_changes_tool,  # NEW: Bond breaking/formation analysis
    
    # Recommendation tools
    recommend_conditions_tool,
    rule_based_conditions_tool,
    enhanced_cross_family_recommend_tool,
    search_precedents_tool,
    
    # Database tools
    reaction_dataset_analytics_tool,
    find_reagent_tool,
    reagent_database_analytics_tool,
    list_supported_cores_tool,
    add_reagent_tool,
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

    return {
        "family": raw.get("family", "Unknown"),
        "detected_family": raw.get("detected_family"),
        "search_all_families": raw.get("search_all_families", False),
        "recommendation": raw.get("recommendation", {}) or {},
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
