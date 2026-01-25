"""
LangChain tool wrappers for ChemTools featurization, analysis, and reagent registry access.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

from langchain_core.tools import tool
from pydantic import BaseModel, Field

from chemtools.featurizers import analysis as analysis_tools
from chemtools.featurizers import calculable as calculable_tools
from chemtools.featurizers import molecule as motif_molecule
from chemtools.featurizers import reaction as motif_reaction
from chemtools.featurizers import reaction_detection as reaction_detection_tools
from chemtools.featurizers import reaction_pair as reaction_pair_tools
from chemtools.featurizers import unified as unified_tools
from chemtools.featurizers.analysis import reaction_context as reaction_context_tools


# ============================================================================
# Helper utilities
# ============================================================================

def _to_jsonable(value: Any) -> Any:
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, dict):
        return {k: _to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        return _to_jsonable(to_dict())
    to_list = getattr(value, "tolist", None)
    if callable(to_list):
        return _to_jsonable(to_list())
    return value


def _success_response(data: Any) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"success": True}
    if isinstance(data, dict):
        payload.update(_to_jsonable(data))
    else:
        payload["data"] = _to_jsonable(data)
    return payload


def _error_response(message: str, extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"success": False, "error": message}
    if extra:
        payload.update(_to_jsonable(extra))
    return payload


def _extract_smiles_from_normalized(items: List[Dict[str, Any]]) -> List[str]:
    smiles_list: List[str] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        for key in ("smiles_norm", "largest_smiles", "input"):
            value = item.get(key)
            if value:
                smiles_list.append(str(value).strip())
                break
    return [smiles for smiles in smiles_list if smiles]


def _summarize_reagent_entry(entry: Dict[str, Any], include_properties: bool = False) -> Dict[str, Any]:
    summary = {
        "name": entry.get("name"),
        "cas": entry.get("cas"),
        "abbreviation": entry.get("abbreviation"),
        "smiles": entry.get("smiles"),
        "family_id": entry.get("family_id"),
        "tag": entry.get("tag"),
        "roles": sorted((entry.get("roles") or {}).keys()),
    }
    if include_properties:
        summary["properties"] = entry.get("properties") or {}
    return summary


def _canonical_reagent_roles() -> List[str]:
    try:
        from chemtools import reagent as reagent_tools

        roles = list(reagent_tools.ROLE_FILES.keys())
        role_priority = getattr(reagent_tools, "ROLE_PRIORITY", {})
        roles.sort(key=lambda role: (role_priority.get(role, 99), role))
        return roles
    except Exception:
        return []


def _is_reagent_exact_match(query_norm: str, entry: Dict[str, Any]) -> bool:
    if not query_norm:
        return False
    from chemtools import reagent as reagent_tools

    name_norm = reagent_tools.normalize_name(entry.get("name", ""))
    if name_norm and name_norm == query_norm:
        return True
    for abbr in entry.get("abbreviation", []) or []:
        if not abbr:
            continue
        abbr_norm = reagent_tools.normalize_name(str(abbr))
        if abbr_norm == query_norm:
            return True
    cas = entry.get("cas")
    if cas:
        cas_norm = reagent_tools.normalize_name(str(cas))
        if cas_norm == query_norm:
            return True
    return False


def _is_reagent_partial_match(query_norm: str, entry: Dict[str, Any]) -> bool:
    if not query_norm:
        return False
    from chemtools import reagent as reagent_tools

    name_norm = reagent_tools.normalize_name(entry.get("name", ""))
    if name_norm and (query_norm in name_norm or name_norm in query_norm):
        return True
    for abbr in entry.get("abbreviation", []) or []:
        if not abbr:
            continue
        abbr_norm = reagent_tools.normalize_name(str(abbr))
        if abbr_norm and (query_norm in abbr_norm or abbr_norm in query_norm):
            return True
    return False


def _build_reagent_info(entry: Dict[str, Any], role: str, match_kind: str) -> Dict[str, Any]:
    abbr = None
    abbreviations = entry.get("abbreviation") or []
    if isinstance(abbreviations, list) and abbreviations:
        abbr = abbreviations[0]
    elif isinstance(abbreviations, str):
        abbr = abbreviations
    return {
        "name": entry.get("name"),
        "type": role,
        "cas": entry.get("cas"),
        "abbreviation": abbr,
        "smiles": entry.get("smiles"),
        "inchi_key": entry.get("inchi_key"),
        "aliases": entry.get("aliases", []),
        "roles": entry.get("roles", {}),
        "properties": entry.get("properties", {}),
        "found": True,
        "match_kind": match_kind,
    }


def _resolve_hte_inputs(
    reaction_smiles: Optional[str],
    reactant_a_smiles: Optional[str],
    reactant_b_smiles: Optional[str],
    product_smiles: Optional[str],
) -> tuple[str, Optional[str], Optional[str]]:
    if reaction_smiles:
        try:
            from chemtools.smiles import normalize_reaction
        except Exception:
            normalize_reaction = None
        if normalize_reaction:
            normalized = normalize_reaction(reaction_smiles)
            reactants = _extract_smiles_from_normalized(normalized.get("reactants", []))
            products = _extract_smiles_from_normalized(normalized.get("products", []))
            if not reactant_a_smiles and reactants:
                reactant_a_smiles = reactants[0]
            if reactant_b_smiles is None and len(reactants) > 1:
                reactant_b_smiles = reactants[1]
            if product_smiles is None and products:
                product_smiles = products[0]

    reactant_a_smiles = (reactant_a_smiles or "").strip()
    reactant_b_smiles = (reactant_b_smiles or "").strip() or None
    product_smiles = (product_smiles or "").strip() or None
    return reactant_a_smiles, reactant_b_smiles, product_smiles


# ============================================================================
# Input schemas
# ============================================================================

class SmilesInput(BaseModel):
    """Schema for molecule SMILES."""

    smiles: str = Field(..., description="SMILES string.")


class SmilesBatchInput(BaseModel):
    """Schema for a batch of SMILES strings."""

    smiles_list: List[str] = Field(..., description="List of SMILES strings.")


class ReactionSmilesInput(BaseModel):
    """Schema for reaction SMILES."""

    reaction_smiles: str = Field(
        ...,
        description="Reaction SMILES in 'reactants>>products' or 'reactants>agents>products' format.",
    )


class ReactionSmilesWithMaxHitsInput(BaseModel):
    """Schema for reaction SMILES with motif hit caps."""

    reaction_smiles: str = Field(..., description="Reaction SMILES string to analyze.")
    max_hits_per_compound: Optional[int] = Field(
        None,
        description="Optional cap for motif hits per compound.",
    )


class ReactionListDetectionInput(BaseModel):
    """Schema for reaction detection from reactant/product lists."""

    reactant_smiles: List[str] = Field(
        ..., description="Reactant SMILES list for detection."
    )
    product_smiles: Optional[List[str]] = Field(
        None, description="Optional product SMILES list."
    )
    max_hits_per_compound: Optional[int] = Field(
        None,
        description="Optional cap for motif hits per compound.",
    )


class ReactionMotifIdsInput(BaseModel):
    """Schema for motif-id detection from reactant lists."""

    reactant_smiles: List[str] = Field(
        ..., description="Reactant SMILES list for motif detection."
    )
    max_hits_per_compound: Optional[int] = Field(
        None,
        description="Optional cap for motif hits per compound.",
    )


class AnalyzeReactionInput(BaseModel):
    """Schema for analysis.analyze_reaction."""

    reaction_smiles: str = Field(..., description="Reaction SMILES string.")
    use_rxn_insight: bool = Field(
        True, description="Reserved for future reaction insight overrides."
    )


class ReactionContextInput(BaseModel):
    """Schema for context-aware reactant classification."""

    reaction_smiles: str = Field(..., description="Reaction SMILES string.")
    reaction_type: Optional[str] = Field(
        None,
        description="Optional reaction type override (taxonomy id).",
    )
    auto_detect: bool = Field(
        True,
        description="Auto-detect reaction type when reaction_type is not provided.",
    )


class LabelInput(BaseModel):
    """Schema for label normalization helpers."""

    label: str = Field(..., description="Label to normalize.")


class FamilyInput(BaseModel):
    """Schema for reaction family normalization."""

    family: str = Field(..., description="Reaction family or alias.")


class UnifiedFeaturizeMoleculeInput(BaseModel):
    """Schema for unified molecule featurization."""

    smiles: str = Field(..., description="Molecule SMILES string.")
    registry_paths: Optional[Dict[str, str]] = Field(
        None,
        description="Optional registry overrides: {'groups': path, 'compounds': path}.",
    )
    options: Optional[Dict[str, Any]] = Field(
        None,
        description="Optional feature toggles (rdkit props, etc.).",
    )


class UnifiedFeaturizeReactionInput(BaseModel):
    """Schema for unified reaction featurization."""

    reaction_smiles: str = Field(..., description="Reaction SMILES string.")
    registry_paths: Optional[Dict[str, str]] = Field(
        None,
        description="Optional registry overrides: {'groups': path, 'compounds': path}.",
    )
    options: Optional[Dict[str, Any]] = Field(
        None,
        description="Optional feature toggles (roles, agent roles, rdkit props, etc.).",
    )


class MotifFeaturizeMoleculeInput(BaseModel):
    """Schema for motif-based molecule featurization."""

    smiles: str = Field(..., description="Molecule SMILES string.")
    registry_paths: Optional[Dict[str, str]] = Field(
        None,
        description="Optional registry overrides: {'groups': path, 'compounds': path}.",
    )
    options: Optional[Dict[str, Any]] = Field(
        None, description="Optional feature toggles for motif analysis."
    )


class MotifFeaturizeReactionInput(BaseModel):
    """Schema for motif-based reaction featurization."""

    reaction_smiles: str = Field(..., description="Reaction SMILES string.")
    registry_paths: Optional[Dict[str, str]] = Field(
        None,
        description="Optional registry overrides: {'groups': path, 'compounds': path}.",
    )
    options: Optional[Dict[str, Any]] = Field(
        None, description="Optional feature toggles for motif analysis."
    )


class ReactionPairInput(BaseModel):
    """Schema for electrophile/nucleophile pair featurization."""

    electrophile: str = Field(..., description="Electrophile SMILES string.")
    nucleophile: str = Field(..., description="Nucleophile SMILES string.")
    include_calculable: Optional[bool] = Field(
        None, description="Include calculable features when available."
    )
    include_structural: Optional[bool] = Field(
        None, description="Include motif-based structural payloads."
    )


class ReactionPairFlatInput(BaseModel):
    """Schema for flat electrophile/nucleophile features."""

    electrophile: str = Field(..., description="Electrophile SMILES string.")
    nucleophile: str = Field(..., description="Nucleophile SMILES string.")
    include_calculable: Optional[bool] = Field(
        None, description="Include calculable features when available."
    )


class HTERecommendInput(BaseModel):
    """Schema for HTE recommender queries."""

    reaction_smiles: Optional[str] = Field(
        None,
        description="Optional reaction SMILES (A.B>>P) to parse reactants from.",
    )
    reactant_a_smiles: Optional[str] = Field(
        None, description="Reactant A SMILES string."
    )
    reactant_b_smiles: Optional[str] = Field(
        None, description="Optional reactant B SMILES string."
    )
    product_smiles: Optional[str] = Field(
        None, description="Optional product SMILES string."
    )
    top_k: int = Field(10, ge=1, le=200, description="Number of conditions to return.")
    min_experiments: int = Field(
        1, ge=1, le=200, description="Minimum experiments per condition."
    )
    reaction_type_filter: Optional[str] = Field(
        None, description="Optional reaction type/category filter."
    )
    catalyst_filter: Optional[str] = Field(
        None, description="Optional catalyst metal/name filter."
    )
    source_group: Optional[str] = Field(
        None,
        description="Optional source group filter (literature, experiments, rules under data/HTE_db/rules).",
    )
    use_aryl_steric_electronic_weighting: bool = Field(
        False,
        description="Apply aryl steric/electronic similarity weighting.",
    )
    hte_db_path: Optional[str] = Field(
        None, description="Path to HTE database folder or file."
    )
    use_spectator_groups: bool = Field(
        True, description="Apply spectator group weighting when available."
    )


class HTEDatasetSummaryInput(BaseModel):
    """Schema for HTE dataset summary queries."""

    reaction_type_filter: Optional[str] = Field(
        None, description="Reaction type/category filter (e.g., suzuki_miyaura)."
    )
    reactant_type_filters: Optional[List[str]] = Field(
        None,
        description="Reactant motif filters (e.g., ['Ar-Cl']).",
    )
    match_all_reactants: bool = Field(
        False, description="Require all reactant_type_filters to be present."
    )
    source_group: Optional[str] = Field(
        "literature",
        description="Source group filter (literature, experiments, rules under data/HTE_db/rules).",
    )
    top_k: int = Field(10, ge=1, le=200, description="Number of conditions to return.")
    min_experiments: int = Field(
        2, ge=1, le=200, description="Minimum experiments per condition."
    )
    hte_db_path: Optional[str] = Field(
        None, description="Path to HTE database folder or file."
    )


class HTEStatsInput(BaseModel):
    """Schema for HTE database statistics."""

    hte_db_path: Optional[str] = Field(
        None, description="Path to HTE database folder or file."
    )


class ReagentLookupInput(BaseModel):
    """Schema for reagent lookup by name, abbreviation, or CAS."""

    query: str = Field(..., description="Reagent name, abbreviation, or CAS number.")
    role: Optional[str] = Field(
        None,
        description="Optional reagent role filter (ligand, base, solvent, etc.).",
    )
    include_all: bool = Field(
        False,
        description="If True, return all matching roles (up to max_results).",
    )
    max_results: int = Field(
        5, ge=1, le=50, description="Maximum matches to return when include_all is True."
    )


class ReagentRolesInput(BaseModel):
    """Schema for listing available reagent roles."""

    include_counts: bool = Field(
        True, description="Include per-role reagent counts."
    )


class ReagentListByRoleInput(BaseModel):
    """Schema for listing reagents by role."""

    role: str = Field(..., description="Reagent role to list (ligand, base, solvent, etc.).")
    limit: int = Field(25, ge=1, le=200, description="Maximum number of reagents to return.")
    include_properties: bool = Field(
        False,
        description="Include extra CSV properties (formula, bp, mp, etc.) when available.",
    )


class ReagentFamilySearchInput(BaseModel):
    """Schema for listing reagents by family within a role."""

    role: str = Field(..., description="Reagent role to filter by.")
    family_id: str = Field(..., description="Family identifier to filter by.")
    limit: int = Field(25, ge=1, le=200, description="Maximum number of reagents to return.")
    include_properties: bool = Field(
        False,
        description="Include extra CSV properties (formula, bp, mp, etc.) when available.",
    )


class RagSearchInput(BaseModel):
    """Schema for knowledge base RAG search."""

    query: str = Field(..., description="Natural language query to search.")
    top_k: int = Field(3, ge=1, le=10, description="Number of results to return.")
    root: Optional[str] = Field(
        None, description="Optional knowledge base root path override."
    )
    include_text: bool = Field(True, description="Include chunk text in results.")
    max_chars: int = Field(
        1200, ge=200, le=5000, description="Max characters per chunk text."
    )
    chunk_size: int = Field(
        300, ge=50, le=800, description="Approximate words per chunk."
    )
    chunk_overlap: int = Field(
        40, ge=0, le=200, description="Approximate overlap between chunks."
    )


class KBConditionSearchInput(BaseModel):
    """Schema for knowledge base condition summaries."""

    query: str = Field(..., description="Natural language query to match summary tables.")
    top_k: int = Field(5, ge=1, le=20, description="Number of condition rows to return.")
    root: Optional[str] = Field(
        None, description="Optional knowledge base root path override."
    )


# ============================================================================
# Analysis tools
# ============================================================================

@tool(args_schema=SmilesInput)
def analysis_normalize_smiles(smiles: str) -> Dict[str, Any]:
    """Normalize a SMILES string using the analysis layer."""
    try:
        return _success_response(analysis_tools.normalize(smiles))
    except Exception as exc:
        return _error_response("Failed to normalize SMILES.", {"details": str(exc)})


@tool(args_schema=ReactionSmilesInput)
def analysis_normalize_reaction(reaction_smiles: str) -> Dict[str, Any]:
    """Normalize a reaction SMILES string into reactants/agents/products."""
    try:
        return _success_response(analysis_tools.normalize_reaction(reaction_smiles))
    except Exception as exc:
        return _error_response("Failed to normalize reaction SMILES.", {"details": str(exc)})


@tool(args_schema=AnalyzeReactionInput)
def analysis_analyze_reaction(reaction_smiles: str, use_rxn_insight: bool = True) -> Dict[str, Any]:
    """Analyze a reaction: normalization, reactant taxonomy, and family detection."""
    try:
        return _success_response(
            analysis_tools.analyze_reaction(reaction_smiles, use_rxn_insight=use_rxn_insight)
        )
    except Exception as exc:
        return _error_response("Failed to analyze reaction.", {"details": str(exc)})


@tool(args_schema=SmilesInput)
def analysis_classify_reactant_smiles(smiles: str) -> Dict[str, Any]:
    """Return the best reactant match for a SMILES input."""
    try:
        match = analysis_tools.classify_reactant_smiles(smiles)
        return _success_response(_to_jsonable(match))
    except Exception as exc:
        return _error_response("Failed to classify reactant SMILES.", {"details": str(exc)})


@tool(args_schema=SmilesInput)
def analysis_classify_reactant_category(smiles: str) -> Dict[str, Any]:
    """Return the reactant category id for a SMILES input."""
    try:
        category = analysis_tools.classify_reactant_category(smiles)
        return _success_response({"category": category})
    except Exception as exc:
        return _error_response("Failed to classify reactant category.", {"details": str(exc)})


@tool(args_schema=SmilesInput)
def analysis_classify_reactant_group(smiles: str) -> Dict[str, Any]:
    """Return the reactant group label for a SMILES input."""
    try:
        group = analysis_tools.classify_reactant_group(smiles)
        return _success_response({"group": group})
    except Exception as exc:
        return _error_response("Failed to classify reactant group.", {"details": str(exc)})


@tool(args_schema=SmilesBatchInput)
def analysis_classify_reactant_batch(smiles_list: List[str]) -> Dict[str, Any]:
    """Classify a batch of reactant SMILES strings."""
    try:
        matches = analysis_tools.classify_reactant_batch(smiles_list)
        return _success_response([_to_jsonable(m) for m in matches])
    except Exception as exc:
        return _error_response("Failed to classify reactant batch.", {"details": str(exc)})


@tool(args_schema=SmilesInput)
def analysis_get_reactant_category_matches(smiles: str) -> Dict[str, Any]:
    """Return all reactant categories matched by the SMILES input."""
    try:
        categories = analysis_tools.get_reactant_category_matches(smiles)
        return _success_response({"categories": categories})
    except Exception as exc:
        return _error_response("Failed to read reactant category matches.", {"details": str(exc)})


@tool(args_schema=SmilesInput)
def analysis_get_all_reactant_matches(smiles: str) -> Dict[str, Any]:
    """Return all reactant SMARTS matches for the SMILES input."""
    try:
        matches = analysis_tools.get_all_reactant_matches(smiles)
        return _success_response([_to_jsonable(m) for m in matches])
    except Exception as exc:
        return _error_response("Failed to read reactant matches.", {"details": str(exc)})


@tool(args_schema=LabelInput)
def analysis_normalize_reactant_identifier(label: str) -> Dict[str, Any]:
    """Normalize a reactant identifier to its canonical category id."""
    try:
        normalized = analysis_tools.normalize_reactant_identifier(label)
        return _success_response({"reactant_id": normalized})
    except Exception as exc:
        return _error_response("Failed to normalize reactant identifier.", {"details": str(exc)})


@tool(args_schema=LabelInput)
def analysis_normalize_reaction_type(label: str) -> Dict[str, Any]:
    """Normalize a reaction type label to canonical taxonomy id."""
    try:
        normalized = analysis_tools.normalize_reaction_type(label)
        return _success_response({"reaction_type": normalized})
    except Exception as exc:
        return _error_response("Failed to normalize reaction type label.", {"details": str(exc)})


@tool(args_schema=FamilyInput)
def analysis_resolve_reaction_family(family: str) -> Dict[str, Any]:
    """Resolve reaction family aliases to canonical ids."""
    try:
        resolved = analysis_tools.resolve_reaction_family(family)
        return _success_response({"reaction_family": resolved})
    except Exception as exc:
        return _error_response("Failed to resolve reaction family.", {"details": str(exc)})


@tool(args_schema=FamilyInput)
def analysis_canonical_family_label(family: str) -> Dict[str, Any]:
    """Resolve reaction family labels to canonical taxonomy identifiers."""
    try:
        resolved = analysis_tools.canonical_family_label(family)
        return _success_response({"reaction_family": resolved})
    except Exception as exc:
        return _error_response("Failed to resolve canonical family label.", {"details": str(exc)})


@tool(args_schema=ReactionContextInput)
def analysis_classify_reactants_with_context(
    reaction_smiles: str,
    reaction_type: Optional[str] = None,
    auto_detect: bool = True,
) -> Dict[str, Any]:
    """Classify reactants with reaction-type context."""
    try:
        result = reaction_context_tools.classify_reactants_with_context(
            reaction_smiles,
            reaction_type=reaction_type,
            auto_detect=auto_detect,
        )
        payload = _to_jsonable(result)
        payload["summary"] = reaction_context_tools.get_reactant_summary(result)
        return _success_response(payload)
    except Exception as exc:
        return _error_response("Failed to classify reactants with context.", {"details": str(exc)})


@tool(args_schema=ReactionContextInput)
def analysis_reactant_summary(
    reaction_smiles: str,
    reaction_type: Optional[str] = None,
    auto_detect: bool = True,
) -> Dict[str, Any]:
    """Return a summary of context-aware reactant classification."""
    try:
        result = reaction_context_tools.classify_reactants_with_context(
            reaction_smiles,
            reaction_type=reaction_type,
            auto_detect=auto_detect,
        )
        summary = reaction_context_tools.get_reactant_summary(result)
        return _success_response(summary)
    except Exception as exc:
        return _error_response("Failed to summarize reactant roles.", {"details": str(exc)})


# ============================================================================
# Reaction detection tools
# ============================================================================

@tool(args_schema=ReactionSmilesWithMaxHitsInput)
def detection_detect_reaction_types(
    reaction_smiles: str,
    max_hits_per_compound: Optional[int] = None,
) -> Dict[str, Any]:
    """Detect reaction types from a reaction SMILES string."""
    try:
        result = reaction_detection_tools.detect_reaction_types(
            reaction_smiles,
            max_hits_per_compound=max_hits_per_compound,
        )
        return _success_response(result.to_dict())
    except Exception as exc:
        return _error_response("Failed to detect reaction types.", {"details": str(exc)})


@tool(args_schema=ReactionListDetectionInput)
def detection_detect_reaction_types_from_smiles(
    reactant_smiles: List[str],
    product_smiles: Optional[List[str]] = None,
    max_hits_per_compound: Optional[int] = None,
) -> Dict[str, Any]:
    """Detect reaction types from reactant/product SMILES lists."""
    try:
        result = reaction_detection_tools.detect_reaction_types_from_smiles(
            reactant_smiles,
            product_smiles=product_smiles,
            max_hits_per_compound=max_hits_per_compound,
        )
        return _success_response(result.to_dict())
    except Exception as exc:
        return _error_response("Failed to detect reaction types from lists.", {"details": str(exc)})


@tool(args_schema=ReactionMotifIdsInput)
def detection_detect_motif_ids_from_smiles(
    reactant_smiles: List[str],
    max_hits_per_compound: Optional[int] = None,
) -> Dict[str, Any]:
    """Return detected motif IDs for a list of reactant SMILES strings."""
    try:
        motifs = reaction_detection_tools.detect_motif_ids_from_smiles(
            reactant_smiles,
            max_hits_per_compound=max_hits_per_compound,
        )
        return _success_response({"motif_ids": sorted(motifs)})
    except Exception as exc:
        return _error_response("Failed to detect motif ids.", {"details": str(exc)})


# ============================================================================
# Featurizer tools
# ============================================================================

@tool(args_schema=UnifiedFeaturizeMoleculeInput)
def unified_featurize_molecule(
    smiles: str,
    registry_paths: Optional[Dict[str, str]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Return the unified molecule feature bundle."""
    try:
        return _success_response(
            unified_tools.featurize_molecule(
                smiles,
                registry_paths=registry_paths,
                options=options,
            )
        )
    except Exception as exc:
        return _error_response("Unified molecule featurization failed.", {"details": str(exc)})


@tool(args_schema=UnifiedFeaturizeReactionInput)
def unified_featurize_reaction(
    reaction_smiles: str,
    registry_paths: Optional[Dict[str, str]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Return the unified reaction feature bundle."""
    try:
        return _success_response(
            unified_tools.featurize_reaction(
                reaction_smiles,
                registry_paths=registry_paths,
                options=options,
            )
        )
    except Exception as exc:
        return _error_response("Unified reaction featurization failed.", {"details": str(exc)})


@tool(args_schema=MotifFeaturizeMoleculeInput)
def motif_featurize_molecule(
    smiles: str,
    registry_paths: Optional[Dict[str, str]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Return motif-based molecule features (sterics/electronics)."""
    try:
        return _success_response(
            motif_molecule.featurize_molecule(
                smiles,
                registry_paths=registry_paths,
                options=options,
            )
        )
    except Exception as exc:
        return _error_response("Motif molecule featurization failed.", {"details": str(exc)})


@tool(args_schema=MotifFeaturizeReactionInput)
def motif_featurize_reaction(
    reaction_smiles: str,
    registry_paths: Optional[Dict[str, str]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Return motif-based reaction features."""
    try:
        return _success_response(
            motif_reaction.featurize_reaction(
                reaction_smiles,
                registry_paths=registry_paths,
                options=options,
            )
        )
    except Exception as exc:
        return _error_response("Motif reaction featurization failed.", {"details": str(exc)})


@tool(args_schema=ReactionPairInput)
def reaction_pair_featurize_pair(
    electrophile: str,
    nucleophile: str,
    include_calculable: Optional[bool] = None,
    include_structural: Optional[bool] = None,
) -> Dict[str, Any]:
    """Return structured electrophile/nucleophile pair features."""
    try:
        return _success_response(
            reaction_pair_tools.featurize_pair(
                electrophile,
                nucleophile,
                include_calculable=include_calculable,
                include_structural=include_structural,
            )
        )
    except Exception as exc:
        return _error_response("Reaction-pair featurization failed.", {"details": str(exc)})


@tool(args_schema=ReactionPairFlatInput)
def reaction_pair_featurize_flat(
    electrophile: str,
    nucleophile: str,
    include_calculable: Optional[bool] = None,
) -> Dict[str, Any]:
    """Return flat electrophile/nucleophile feature dictionary."""
    try:
        return _success_response(
            reaction_pair_tools.featurize_flat(
                electrophile,
                nucleophile,
                include_calculable=include_calculable,
            )
        )
    except Exception as exc:
        return _error_response("Reaction-pair flat featurization failed.", {"details": str(exc)})


@tool(args_schema=SmilesInput)
def calculable_detect_all_features(smiles: str) -> Dict[str, Any]:
    """Detect all calculable features for a SMILES string."""
    try:
        return _success_response(calculable_tools.detect_all_features(smiles))
    except Exception as exc:
        return _error_response("Failed to detect calculable features.", {"details": str(exc)})


@tool(args_schema=SmilesInput)
def calculable_get_reactant_type_features(smiles: str) -> Dict[str, Any]:
    """Return calculable reactant-type features for a SMILES string."""
    try:
        return _success_response(calculable_tools.get_reactant_type_features(smiles))
    except Exception as exc:
        return _error_response("Failed to read reactant-type features.", {"details": str(exc)})


@tool(args_schema=SmilesInput)
def calculable_classify_reactant_smiles(smiles: str) -> Dict[str, Any]:
    """Return the best reactant match using calculable features."""
    try:
        return _success_response(calculable_tools.classify_reactant_smiles(smiles))
    except Exception as exc:
        return _error_response("Failed to classify reactant via calculable features.", {"details": str(exc)})


# ============================================================================
# HTE recommendation tools
# ============================================================================

@tool(args_schema=HTERecommendInput)
def hte_recommend_conditions(
    reaction_smiles: Optional[str] = None,
    reactant_a_smiles: Optional[str] = None,
    reactant_b_smiles: Optional[str] = None,
    product_smiles: Optional[str] = None,
    top_k: int = 10,
    min_experiments: int = 1,
    reaction_type_filter: Optional[str] = None,
    catalyst_filter: Optional[str] = None,
    source_group: Optional[str] = None,
    use_aryl_steric_electronic_weighting: bool = False,
    hte_db_path: Optional[str] = None,
    use_spectator_groups: bool = True,
) -> Dict[str, Any]:
    """Recommend conditions from the HTE database."""
    try:
        from chemtools.recommend import HTERecommender

        reactant_a, reactant_b, product = _resolve_hte_inputs(
            reaction_smiles,
            reactant_a_smiles,
            reactant_b_smiles,
            product_smiles,
        )
        if not reactant_a:
            return _error_response(
                "HTE recommendation requires at least one reactant SMILES.",
                {"details": "Provide reaction_smiles or reactant_a_smiles."},
            )

        recommender = HTERecommender(hte_db_path or "data/HTE_db")
        result = recommender.recommend(
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            product_smiles=product,
            top_k=top_k,
            min_experiments=min_experiments,
            reaction_type_filter=reaction_type_filter or None,
            catalyst_filter=catalyst_filter or None,
            source_group=source_group or None,
            use_aryl_steric_electronic_weighting=use_aryl_steric_electronic_weighting,
            use_spectator_groups=use_spectator_groups,
        )
        payload = _to_jsonable(result)
        payload["diagnostics"] = {
            "reactant_a_type": result.reactant_a_type,
            "reactant_b_type": result.reactant_b_type,
            "reactant_a_category": result.reactant_a_category,
            "reactant_b_category": result.reactant_b_category,
            "matched_motifs": result.matched_motifs,
            "predicted_reaction_type": result.predicted_reaction_type,
            "reaction_type_confidence": result.reaction_type_confidence,
            "total_matching_experiments": result.total_matching_experiments,
            "database_coverage": result.database_coverage,
            "is_fallback_match": result.is_fallback_match,
        }
        payload["hte_stats"] = recommender.get_statistics()
        payload["input"] = {
            "reaction_smiles": reaction_smiles,
            "reactant_a_smiles": reactant_a,
            "reactant_b_smiles": reactant_b,
            "product_smiles": product,
            "source_group": source_group,
            "use_aryl_steric_electronic_weighting": use_aryl_steric_electronic_weighting,
        }
        return _success_response(payload)
    except Exception as exc:
        return _error_response("HTE recommendation failed.", {"details": str(exc)})


@tool(args_schema=HTEDatasetSummaryInput)
def hte_dataset_summary(
    reaction_type_filter: Optional[str] = None,
    reactant_type_filters: Optional[List[str]] = None,
    match_all_reactants: bool = False,
    source_group: Optional[str] = "literature",
    top_k: int = 10,
    min_experiments: int = 2,
    hte_db_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Summarize top conditions from a filtered HTE dataset slice."""
    try:
        from chemtools.recommend import HTERecommender

        recommender = HTERecommender(hte_db_path or "data/HTE_db")
        payload = recommender.summarize_conditions(
            reaction_type_filter=reaction_type_filter or None,
            reactant_type_filters=reactant_type_filters,
            match_all_reactants=match_all_reactants,
            source_group=source_group or None,
            top_k=top_k,
            min_experiments=min_experiments,
        )
        return _success_response(_to_jsonable(payload))
    except Exception as exc:
        return _error_response("HTE dataset summary failed.", {"details": str(exc)})


@tool(args_schema=HTEStatsInput)
def hte_database_stats(hte_db_path: Optional[str] = None) -> Dict[str, Any]:
    """Return summary statistics for the HTE database."""
    try:
        from chemtools.recommend import HTERecommender

        recommender = HTERecommender(hte_db_path or "data/HTE_db")
        stats = recommender.get_statistics()
        return _success_response({"hte_stats": stats})
    except Exception as exc:
        return _error_response("Failed to read HTE database stats.", {"details": str(exc)})


# ============================================================================
# Reagent registry tools
# ============================================================================

@tool(args_schema=ReagentLookupInput)
def reagent_lookup(
    query: str,
    role: Optional[str] = None,
    include_all: bool = False,
    max_results: int = 5,
) -> Dict[str, Any]:
    """Lookup a reagent by name, abbreviation, or CAS."""
    try:
        from chemtools import reagent as reagent_tools

        query = (query or "").strip()
        if not query:
            return _error_response("Reagent lookup requires a non-empty query.")

        query_norm = reagent_tools.normalize_name(query)
        roles = _canonical_reagent_roles() or reagent_tools.get_all_reagent_types()

        def search_role(role_name: str, allow_partial: bool) -> Optional[Dict[str, Any]]:
            db = reagent_tools.load_reagent_database(role_name)
            for entry in db:
                if not isinstance(entry, dict):
                    continue
                if _is_reagent_exact_match(query_norm, entry):
                    return _build_reagent_info(entry, role_name, "exact")
            if not allow_partial:
                return None
            for entry in db:
                if not isinstance(entry, dict):
                    continue
                if _is_reagent_partial_match(query_norm, entry):
                    return _build_reagent_info(entry, role_name, "partial")
            return None

        if role:
            match = search_role(role, allow_partial=True)
            matches = [match] if match else []
            return _success_response({
                "query": query,
                "role": role,
                "matches": matches,
                "searched_roles": [role],
            })

        matches: List[Dict[str, Any]] = []
        for reagent_type in roles:
            match = search_role(reagent_type, allow_partial=False)
            if match:
                matches.append(match)
                if not include_all:
                    break
                if len(matches) >= max_results:
                    break

        if not matches:
            for reagent_type in roles:
                match = search_role(reagent_type, allow_partial=True)
                if match:
                    matches.append(match)
                    if not include_all:
                        break
                    if len(matches) >= max_results:
                        break

        return _success_response({
            "query": query,
            "matches": matches,
            "searched_roles": roles,
        })
    except Exception as exc:
        return _error_response("Reagent lookup failed.", {"details": str(exc)})


@tool(args_schema=ReagentRolesInput)
def reagent_list_roles(include_counts: bool = True) -> Dict[str, Any]:
    """List available reagent roles in the registry."""
    try:
        from chemtools import reagent as reagent_tools

        roles = _canonical_reagent_roles() or reagent_tools.get_all_reagent_types()
        payload = {"roles": roles, "total_roles": len(roles)}
        if include_counts:
            payload["role_counts"] = {
                role: reagent_tools.count_reagents_by_type(role) for role in roles
            }
        return _success_response(payload)
    except Exception as exc:
        return _error_response("Failed to list reagent roles.", {"details": str(exc)})


@tool(args_schema=ReagentListByRoleInput)
def reagent_list_by_role(
    role: str,
    limit: int = 25,
    include_properties: bool = False,
) -> Dict[str, Any]:
    """List reagents for a specific role."""
    try:
        from chemtools import reagent as reagent_tools

        role = (role or "").strip()
        if not role:
            return _error_response("Role is required to list reagents.")

        entries = reagent_tools.get_all_reagents_by_type(role)
        entries = sorted(entries, key=lambda item: (item.get("name") or "").lower())
        total = len(entries)
        items = [
            _summarize_reagent_entry(entry, include_properties=include_properties)
            for entry in entries[:limit]
        ]
        return _success_response({
            "role": role,
            "total": total,
            "returned": len(items),
            "items": items,
        })
    except Exception as exc:
        return _error_response("Failed to list reagents by role.", {"details": str(exc)})


@tool(args_schema=ReagentFamilySearchInput)
def reagent_list_by_family(
    role: str,
    family_id: str,
    limit: int = 25,
    include_properties: bool = False,
) -> Dict[str, Any]:
    """List reagents for a role/family combination."""
    try:
        from chemtools import reagent as reagent_tools

        role = (role or "").strip()
        family_id = (family_id or "").strip()
        if not role or not family_id:
            return _error_response("Role and family_id are required.")

        entries = reagent_tools.find_reagents_by_family(role, family_id)
        entries = sorted(entries, key=lambda item: (item.get("name") or "").lower())
        total = len(entries)
        items = [
            _summarize_reagent_entry(entry, include_properties=include_properties)
            for entry in entries[:limit]
        ]
        return _success_response({
            "role": role,
            "family_id": family_id,
            "total": total,
            "returned": len(items),
            "items": items,
        })
    except Exception as exc:
        return _error_response("Failed to list reagents by family.", {"details": str(exc)})


# ============================================================================
# RAG tools
# ============================================================================

@tool(args_schema=KBConditionSearchInput)
def kb_recommend_conditions(
    query: str,
    top_k: int = 5,
    root: Optional[str] = None,
) -> Dict[str, Any]:
    """Return condition summary rows extracted from the knowledge base."""
    try:
        from chem_assistant.kb_conditions import search_condition_summaries

        payload = search_condition_summaries(
            query,
            top_k=top_k,
            root=root,
        )
        return _success_response(payload)
    except Exception as exc:
        return _error_response(
            "Knowledge base condition lookup failed.", {"details": str(exc)}
        )


@tool(args_schema=RagSearchInput)
def rag_search(
    query: str,
    top_k: int = 3,
    root: Optional[str] = None,
    include_text: bool = True,
    max_chars: int = 1200,
    chunk_size: int = 300,
    chunk_overlap: int = 40,
) -> Dict[str, Any]:
    """Search the curated knowledge base for relevant snippets."""
    try:
        from chem_assistant.rag import search_knowledge_base

        payload = search_knowledge_base(
            query,
            top_k=top_k,
            root=root,
            include_text=include_text,
            max_chars=max_chars,
            chunk_size=chunk_size,
            chunk_overlap=chunk_overlap,
        )
        return _success_response(payload)
    except Exception as exc:
        return _error_response("RAG search failed.", {"details": str(exc)})


# ============================================================================
# Tool registry and helpers
# ============================================================================

CHEMTOOLS_TOOLS = [
    # Analysis
    analysis_normalize_smiles,
    analysis_normalize_reaction,
    analysis_analyze_reaction,
    analysis_classify_reactant_smiles,
    analysis_classify_reactant_category,
    analysis_classify_reactant_group,
    analysis_classify_reactant_batch,
    analysis_get_reactant_category_matches,
    analysis_get_all_reactant_matches,
    analysis_normalize_reactant_identifier,
    analysis_normalize_reaction_type,
    analysis_resolve_reaction_family,
    analysis_canonical_family_label,
    analysis_classify_reactants_with_context,
    analysis_reactant_summary,
    # Reaction detection
    detection_detect_reaction_types,
    detection_detect_reaction_types_from_smiles,
    detection_detect_motif_ids_from_smiles,
    # Featurizers
    unified_featurize_molecule,
    unified_featurize_reaction,
    motif_featurize_molecule,
    motif_featurize_reaction,
    reaction_pair_featurize_pair,
    reaction_pair_featurize_flat,
    calculable_detect_all_features,
    calculable_get_reactant_type_features,
    calculable_classify_reactant_smiles,
    # HTE recommendations
    hte_recommend_conditions,
    hte_dataset_summary,
    hte_database_stats,
    # Reagent registry
    reagent_lookup,
    reagent_list_roles,
    reagent_list_by_role,
    reagent_list_by_family,
    # RAG
    rag_search,
    kb_recommend_conditions,
]


def _schema_properties(tool_obj: Any) -> Dict[str, Any]:
    args_schema = getattr(tool_obj, "args_schema", None)
    if not args_schema:
        return {}
    try:
        return args_schema.model_json_schema().get("properties", {})
    except AttributeError:
        return args_schema.schema().get("properties", {})


def _schema_required(tool_obj: Any) -> List[str]:
    args_schema = getattr(tool_obj, "args_schema", None)
    if not args_schema:
        return []
    try:
        return list(args_schema.model_json_schema().get("required", []))
    except AttributeError:
        return list(args_schema.schema().get("required", []))


def get_tool_descriptions() -> List[Dict[str, Any]]:
    """Return structured descriptions (name, docstring, parameters) for tools."""
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


def print_tool_summary() -> None:
    """Print a summary of available tools."""
    print("=" * 70)
    print("ChemTools Featurization/Analysis Tools")
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
    print_tool_summary()
