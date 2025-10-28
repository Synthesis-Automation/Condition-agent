"""
LangChain tool wrappers for ChemTools functions.

This module exposes existing chemtools functionality as LangChain tools
without modifying the original chemtools codebase.

Available Tools:
    - normalize_smiles_tool: Canonicalize SMILES strings
    - normalize_reaction_tool: Canonicalize reaction SMILES
    - detect_reaction_family_tool: Detect reaction family/type
    - recommend_conditions_tool: Get ML-based condition recommendations
    - search_precedents_tool: Search for similar precedent reactions
    - find_reagent_tool: Look up reagent information from database
    - classify_reactant_tool: Classify reactant type (aryl halide, amine, etc.)
    - get_functional_groups_tool: Detect functional groups in a molecule

Usage:
    from lang_chain.chemtools_wrapper import CHEMTOOLS_TOOLS
    from langgraph.prebuilt import create_react_agent
    
    agent = create_react_agent(llm, CHEMTOOLS_TOOLS)
"""

from typing import Dict, Any, List, Optional, Sequence, Tuple
import json
import time
from dataclasses import dataclass
from collections import OrderedDict
from pydantic import BaseModel, Field
from langchain_core.tools import tool

# Import chemtools functions
from chemtools.smiles import normalize, normalize_reaction
from chemtools.router import detect_family_from_reaction
from chemtools.recommend.modules import recommend_from_reaction
from chemtools.precedent import knn as precedent_knn
from chemtools.reagent import (
    find_reagent,
    classify_reactant_smiles,
    add_reagent_entry,
    ReagentAdditionError,
)
from chemtools.util.functional_groups import detect_all as detect_functional_groups

REAGENT_RESOLVER_TIMEOUT = 6.0

from .constraint_parser import (
    ConstraintSpec,
    build_constraint_spec,
    filter_cores_by_constraints,
    format_constraints_for_prompt,
)


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
        Output: {"family": "Suzuki", "confidence": "high"}
    """
    try:
        result = detect_family_from_reaction(reaction_smiles)
        payload = dict(result or {})
        if not payload.get("family"):
            return _error_response(
                "Could not determine reaction family",
                {"details": payload, "reaction_smiles": reaction_smiles},
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


# ============================================================================
# Recommendation Tools
# ============================================================================

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
            family_result = detect_family_from_reaction(reaction_smiles)
            family = family_result.get("family", "Unknown")
        
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
# Reagent Database Tools
# ============================================================================

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
    get_functional_groups_tool,
    
    # Recommendation tools
    recommend_conditions_tool,
    search_precedents_tool,
    
    # Database tools
    find_reagent_tool,
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
