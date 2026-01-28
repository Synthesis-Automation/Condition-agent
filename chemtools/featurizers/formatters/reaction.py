"""
Reaction-level feature formatting and bundling.

Handles reaction type detection, reactant/product processing, and reaction key generation.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from chemtools.util import rdkit_helpers
from chemtools.smiles import normalize_reaction

from ..detection import detect_reaction_types
from ..analysis.reaction_context import classify_reactants_with_context, get_reactant_summary
from ..analysis.feasibility import analyze_snar_feasibility

from .molecule import build_molecule_bundle, to_bool
from .aggregation import aggregate_reaction_features, infer_intramolecular
from .utils import extract_motif_ids


def format_reaction_type_summary(detection: Any) -> Dict[str, Any]:
    """Extract reaction type information with alternatives."""
    matches = detection.matches if detection else []
    if not matches:
        return {"reaction_type": "Unknown", "confidence": 0.0, "slot_evidence": {}}
    
    best = matches[0]
    result = {
        "reaction_type": best.reaction_type,
        "name": best.name,
        "category": best.category,
        "confidence": best.confidence,
        "slot_evidence": {slot: list(values) for slot, values in best.slot_evidence.items()},
    }
    
    # Add alternatives if available (top 3 total)
    if len(matches) > 1:
        alts = []
        for alt in matches[1:3]:
            alts.append({
                "reaction_type": alt.reaction_type,
                "name": alt.name,
                "confidence": alt.confidence
            })
        result["alternatives"] = alts
        
    return result


def format_reaction_key(
    reacted: Iterable[str],
    formed: Iterable[str],
    spectators: Iterable[str]
) -> str:
    """
    Format motif sets into standardized Reaction_Key string.
    
    Format: Reacted -> Formed || Spectators
    Example: Ar-Br|R_acidic-H -> Ar-Alkyl || Ar-COR|RCH2-COR|R_acidic-H
    
    Args:
        reacted: Motifs consumed in the reaction
        formed: Motifs created in the product
        spectators: Motifs present but unchanged
    
    Returns:
        Formatted Reaction_Key string
    """
    reacted_list = sorted(reacted) if reacted else []
    formed_list = sorted(formed) if formed else []
    spectators_list = sorted(spectators) if spectators else []
    
    reacted_str = "|".join(reacted_list) if reacted_list else "[]"
    formed_str = "|".join(formed_list) if formed_list else "[]"
    spectators_str = "|".join(spectators_list) if spectators_list else "[]"
    
    return f"{reacted_str} -> {formed_str} || {spectators_str}"


def select_primary_reacted_motifs(
    reactant_entries: Iterable[Any],
    reacted_set: set[str],
) -> List[str]:
    """Select primary reacted motif from each reactant."""
    primary: List[str] = []
    for entry in reactant_entries or []:
        motifs = []
        if isinstance(entry, dict):
            motifs = extract_motif_ids(entry.get("motifs", []))
        else:
            motifs = extract_motif_ids(entry)
        reacted_here = [m for m in motifs if m in reacted_set]
        if reacted_here:
            primary.append(reacted_here[0])
    return primary


def select_primary_formed_motif(
    product_motifs: Iterable[Any],
    formed_set: set[str],
) -> Optional[str]:
    """Select the first formed motif from products."""
    for motif in extract_motif_ids(product_motifs):
        if motif in formed_set:
            return motif
    return None


def classify_agent_roles(agents: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    """Classify reagents/solvents by role using reagent taxonomy."""
    from functools import lru_cache
    from chemtools.reagent.reagent_v2 import ReagentTaxonomyV2
    
    @lru_cache(maxsize=1)
    def load_taxonomy() -> Optional[ReagentTaxonomyV2]:
        try:
            return ReagentTaxonomyV2.from_path()
        except Exception:
            return None
    
    def get_agent_smiles(agent: Dict[str, Any]) -> str:
        """Extract SMILES from agent dict."""
        for key in ("smiles", "smiles_norm", "largest_smiles"):
            value = agent.get(key)
            if value:
                return str(value)
        return ""
    
    taxonomy = load_taxonomy()
    if not taxonomy:
        return {"agent_count": 0, "role_counts": {}, "family_counts": {}, "role_flags": {}, "flags": {}}
    
    role_flags = [
        "metal_catalyst", "organo_catalyst", "enzyme", "ligand", "base",
        "acid", "solvent", "additive", "oxidant", "reductant",
        "condensation_agent", "other_reagent"
    ]
    
    entries: List[Dict[str, Any]] = []
    role_counts: Dict[str, int] = {}
    family_counts: Dict[str, int] = {}
    flags: Dict[str, bool] = {role: False for role in role_flags}
    
    for agent in agents or []:
        smiles = get_agent_smiles(agent)
        if not smiles:
            continue
        
        reagent = taxonomy.lookup_reagent(smiles)
        if not reagent:
            continue
        
        role = reagent.get("role") or "other_reagent"
        family = reagent.get("family") or "Unknown"
        
        role_counts[role] = role_counts.get(role, 0) + 1
        family_counts[family] = family_counts.get(family, 0) + 1
        
        if role in flags:
            flags[role] = True
        
        entries.append({
            "smiles": smiles,
            "role": role,
            "family": family,
            "name": reagent.get("name"),
        })
    
    return {
        "entries": entries,
        "agent_count": len(entries),
        "role_counts": role_counts,
        "family_counts": family_counts,
        "role_flags": {k: v for k, v in flags.items() if v},
        "flags": flags,
    }


def featurize_reaction(
    reaction_smiles: str,
    *,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Return a unified reaction feature bundle.
    
    Args:
        reaction_smiles: Reaction SMILES with >> separator
        registry_paths: Custom taxonomy paths
        options: Featurization options
        
    Returns:
        Complete reaction bundle with detection, roles, aggregates
    """
    options = options or {}
    include_roles = to_bool(options.get("include_roles"), default=True)
    include_agent_roles = to_bool(options.get("include_agent_roles"), default=True)
    confirm_suzuki_products = to_bool(options.get("confirm_suzuki_products"), default=False)
    confirm_coupling_products = to_bool(options.get("confirm_coupling_products"), default=False)

    # Detect reaction type
    detection = detect_reaction_types(
        reaction_smiles,
        max_hits_per_compound=options.get("max_hits_per_compound"),
        confirm_suzuki_products=confirm_suzuki_products,
        confirm_coupling_products=confirm_coupling_products,
    )
    detection_payload = detection.to_dict()
    reaction_type = format_reaction_type_summary(detection)

    # Normalize reaction SMILES
    normalized = normalize_reaction(reaction_smiles)
    reactant_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in (normalized.get("reactants") or [])
    ]
    reactant_smiles = [s for s in reactant_smiles if s]

    # Classify agents/reagents
    agent_roles = None
    if include_agent_roles:
        agent_roles = classify_agent_roles(normalized.get("agents") or [])

    # Featurize reactants
    reactant_bundles = [
        build_molecule_bundle(smiles, registry_paths=registry_paths, options=options)
        for smiles in reactant_smiles
    ]

    # Featurize products
    product_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in (normalized.get("products") or [])
    ]
    product_smiles = [s for s in product_smiles if s]
    product_bundles: List[Dict[str, Any]] = []
    product_motif_ids: List[str] = []
    for smiles in product_smiles:
        try:
            bundle = build_molecule_bundle(smiles, registry_paths=registry_paths, options=options)
        except Exception:
            continue
        product_bundles.append(bundle)
        product_motif_ids.extend(extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", [])))

    # Aggregate features
    aggregates = aggregate_reaction_features(
        reactant_bundles, product_motif_ids=product_motif_ids
    )

    # Classify reactant roles
    roles_summary = None
    if include_roles:
        try:
            roles = classify_reactants_with_context(reaction_smiles)
            roles_summary = get_reactant_summary(roles)
        except Exception:
            roles_summary = None

    intramolecular = infer_intramolecular(reactant_smiles, product_smiles, roles_summary)

    # Generate Reaction_Key
    reaction_key = None
    if product_bundles and reactant_bundles:
        reactant_primary: List[str] = []
        for bundle in reactant_bundles:
            reactant_primary.extend(extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", [])))

        product_primary: List[str] = []
        for bundle in product_bundles:
            product_primary.extend(extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", [])))

        reactant_counts = Counter(reactant_primary)
        product_counts = Counter(product_primary)

        reacted: set[str] = set()
        formed: set[str] = set()
        spectators: set[str] = set()
        all_motifs = set(reactant_counts.keys()) | set(product_counts.keys())
        for motif in all_motifs:
            rc = reactant_counts.get(motif, 0)
            pc = product_counts.get(motif, 0)
            if pc > rc:
                formed.add(motif)
                if rc > 0:
                    spectators.add(motif)
            elif pc < rc:
                reacted.add(motif)
                if pc > 0:
                    spectators.add(motif)
            else:
                if rc > 0:
                    spectators.add(motif)

        primary_reacted = select_primary_reacted_motifs(reactant_bundles, reacted)
        primary_formed = select_primary_formed_motif(product_primary, formed)

        reaction_key = format_reaction_key(
            primary_reacted if primary_reacted else reacted,
            [primary_formed] if primary_formed else formed,
            spectators,
        )

    reaction = {
        "reaction_smiles": reaction_smiles,
        "normalized": normalized,
        "detection": detection_payload,
        "reaction_type": reaction_type,
        "reactants": reactant_bundles,
        "products": product_bundles,
        "aggregates": aggregates,
        "reaction_key": reaction_key,
        "roles": roles_summary,
        "agent_roles": agent_roles,
        "intramolecular": intramolecular,
    }

    # Add feasibility analysis for specific reaction types
    rt_id = reaction_type.get("reaction_type")
    if rt_id == "snar_cn" or rt_id == "c_n_cross_coupling":
        reaction["feasibility"] = analyze_snar_feasibility(reaction)

    meta = {
        "rdkit_available": rdkit_helpers.rdkit_available(),
        "errors": [detection_payload.get("error")] if detection_payload.get("error") else [],
    }
    return {
        "schema_version": "v1",
        "kind": "reaction",
        "reaction": reaction,
        "meta": meta,
    }
