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
    - reaction_dataset_analytics_tool: Analyze reaction dataset frequency/yield statistics
    - find_reagent_tool: Look up reagent information from database
    - reagent_database_analytics_tool: Summarize reagent registry statistics
    - list_supported_cores_tool: Enumerate catalyst cores observed in precedents
    - add_reagent_tool: Insert or preview reagent taxonomy entries
    - classify_reactant_tool: Classify reactant type (aryl halide, amine, etc.)
    - get_functional_groups_tool: Detect functional groups in a molecule
    - molpipeline_featurize_tool: Generate molecular features with optional MolPipeline vectors
    - analyze_bond_changes_tool: Analyze bond breaking/formation in reactions (NEW)

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
from chemtools import detect_reaction
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
    analyze_bond_changes_tool,  # NEW: Bond breaking/formation analysis
    
    # Recommendation tools
    recommend_conditions_tool,
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
