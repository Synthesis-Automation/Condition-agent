"""
LangChain tool wrappers for ChemTools featurization and analysis.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from langchain_core.tools import tool
from pydantic import BaseModel, Field

from chemtools.featurizers import analysis as analysis_tools
from chemtools.featurizers import calculable as calculable_tools
from chemtools.featurizers import molpipeline as molpipeline_tools
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
        description="Optional feature toggles (rdkit props, functional groups, etc.).",
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
    include_molpipeline: Optional[bool] = Field(
        None, description="Include MolPipeline descriptors when available."
    )
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
    include_molpipeline: Optional[bool] = Field(
        None, description="Include MolPipeline descriptors when available."
    )
    include_calculable: Optional[bool] = Field(
        None, description="Include calculable features when available."
    )


class MolpipelineMorganInput(BaseModel):
    """Schema for Morgan fingerprints via MolPipeline."""

    smiles: Union[str, List[str]] = Field(
        ..., description="SMILES string or list of SMILES strings."
    )
    n_bits: int = Field(2048, ge=128, le=4096, description="Fingerprint size.")
    radius: int = Field(2, ge=1, le=6, description="Morgan radius.")
    n_jobs: int = Field(1, ge=1, le=32, description="Parallel worker count.")
    return_sparse: bool = Field(
        False,
        description="Return sparse matrices (default False to keep JSON-friendly output).",
    )


class MolpipelinePhyschemInput(BaseModel):
    """Schema for MolPipeline physchem descriptors."""

    smiles: Union[str, List[str]] = Field(
        ..., description="SMILES string or list of SMILES strings."
    )
    descriptor_list: Optional[List[str]] = Field(
        None,
        description="Optional list of descriptor names (default MolPipeline set).",
    )
    n_jobs: int = Field(1, ge=1, le=32, description="Parallel worker count.")


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
    include_molpipeline: Optional[bool] = None,
    include_calculable: Optional[bool] = None,
    include_structural: Optional[bool] = None,
) -> Dict[str, Any]:
    """Return structured electrophile/nucleophile pair features."""
    try:
        return _success_response(
            reaction_pair_tools.featurize_pair(
                electrophile,
                nucleophile,
                include_molpipeline=include_molpipeline,
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
    include_molpipeline: Optional[bool] = None,
    include_calculable: Optional[bool] = None,
) -> Dict[str, Any]:
    """Return flat electrophile/nucleophile feature dictionary."""
    try:
        return _success_response(
            reaction_pair_tools.featurize_flat(
                electrophile,
                nucleophile,
                include_molpipeline=include_molpipeline,
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


@tool(args_schema=MolpipelineMorganInput)
def molpipeline_morgan_fingerprint(
    smiles: Union[str, List[str]],
    n_bits: int = 2048,
    radius: int = 2,
    n_jobs: int = 1,
    return_sparse: bool = False,
) -> Dict[str, Any]:
    """Return Morgan fingerprints via MolPipeline."""
    try:
        fingerprint = molpipeline_tools.morgan_fingerprint(
            smiles,
            n_bits=n_bits,
            radius=radius,
            n_jobs=n_jobs,
            return_sparse=return_sparse,
        )
        return _success_response(_to_jsonable(fingerprint))
    except Exception as exc:
        return _error_response("MolPipeline Morgan fingerprint failed.", {"details": str(exc)})


@tool(args_schema=MolpipelinePhyschemInput)
def molpipeline_physchem_features(
    smiles: Union[str, List[str]],
    descriptor_list: Optional[List[str]] = None,
    n_jobs: int = 1,
) -> Dict[str, Any]:
    """Return physchem descriptors via MolPipeline."""
    try:
        features = molpipeline_tools.physchem_features(
            smiles,
            descriptor_list=descriptor_list,
            n_jobs=n_jobs,
        )
        return _success_response(_to_jsonable(features))
    except Exception as exc:
        return _error_response("MolPipeline physchem features failed.", {"details": str(exc)})


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
    molpipeline_morgan_fingerprint,
    molpipeline_physchem_features,
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
