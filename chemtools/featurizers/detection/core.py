"""
Main API for reaction type detection.

Combines motif-based matching, bond change analysis, and coupling confirmation.
"""

from __future__ import annotations

from typing import Iterable, List, Optional

from chemtools.util.rdkit_helpers import rdkit_available
from ...taxonomy.reaction_catalog import load_reaction_catalog

from .models import ReactionMatch, ReactionDetectionResult
from .bond_changes import detect_reaction_types_by_bond_changes, clear_bond_change_cache
from .matchers import _detect_motif_profile, detect_motif_ids_from_smiles, match_reaction_definition
from .coupling import confirm_coupling_product_by_attachment, get_coupling_confirmation_specs


# Re-export for backward compatibility
__all__ = [
    "detect_reaction_types",
    "detect_reaction_types_from_smiles",
    "detect_motif_ids_from_smiles",
    "clear_bond_change_cache",
    "ReactionMatch",
    "ReactionDetectionResult",
]


def detect_reaction_types(
    reaction_smiles: str,
    *,
    max_hits_per_compound: Optional[int] = None,
    use_bond_changes: bool = False,
    bond_change_threshold: float = 0.4,
    confirm_coupling_products: bool = False,
) -> ReactionDetectionResult:
    """
    Detect reaction types from a reaction SMILES string.
    
    Args:
        reaction_smiles: Reaction SMILES with >> separator
        max_hits_per_compound: Limit motif matches per compound
        use_bond_changes: Try bond change analysis first
        bond_change_threshold: Minimum similarity for bond changes
        confirm_coupling_products: Verify coupling products by attachment (default False)
            DEPRECATED: Reacted motifs analysis is now sufficient for validation
        
    Returns:
        ReactionDetectionResult with ranked matches
        
    Note:
        Validation now relies on reacted/formed motifs analysis instead of
        SMARTS-based coupling confirmation, which had issues with implicit hydrogens.
    """
    if not reaction_smiles:
        return ReactionDetectionResult(matches=[], error="empty_reaction")
    
    if use_bond_changes:
        bond_change_result = detect_reaction_types_by_bond_changes(
            reaction_smiles,
            min_similarity=bond_change_threshold,
        )
        if bond_change_result.matches:
            return bond_change_result
    
    # Lazy import to avoid circular dependency
    from chemtools.smiles import normalize_reaction
    
    normalized = normalize_reaction(reaction_smiles)
    reactants = normalized.get("reactants") or []
    reactant_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in reactants
    ]
    reactant_smiles = [s for s in reactant_smiles if s]

    # Only consider products if the reaction SMILES explicitly contains them (via > or >>)
    product_smiles = None
    if ">" in reaction_smiles:
        product_smiles = [
            item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
            for item in (normalized.get("products") or [])
        ]
        product_smiles = [s for s in product_smiles if s]

    if not reactant_smiles:
        return ReactionDetectionResult(matches=[], error="no_reactants")
    
    return detect_reaction_types_from_smiles(
        reactant_smiles,
        product_smiles=product_smiles,
        max_hits_per_compound=max_hits_per_compound,
        confirm_coupling_products=confirm_coupling_products,
    )


def detect_reaction_types_from_smiles(
    reactant_smiles: Iterable[str],
    *,
    product_smiles: Optional[Iterable[str]] = None,
    max_hits_per_compound: Optional[int] = None,
    confirm_coupling_products: bool = False,
) -> ReactionDetectionResult:
    """
    Detect reaction types from a list of reactant SMILES strings.
    
    Args:
        reactant_smiles: List of reactant SMILES
        product_smiles: Optional list of product SMILES
        max_hits_per_compound: Limit motif matches per compound
        confirm_coupling_products: Verify coupling products by attachment (default False)
            DEPRECATED: Reacted motifs analysis is now sufficient for validation
        
    Returns:
        ReactionDetectionResult with ranked matches
        
    Note:
        Validation now relies on reacted/formed motifs analysis instead of
        SMARTS-based coupling confirmation, which had issues with implicit hydrogens.
    """
    if not rdkit_available():
        return ReactionDetectionResult(matches=[], error="rdkit_unavailable")

    detected_profile = _detect_motif_profile(
        reactant_smiles, max_hits_per_compound=max_hits_per_compound
    )
    if not detected_profile:
        return ReactionDetectionResult(matches=[])

    product_profile = (
        _detect_motif_profile(product_smiles, max_hits_per_compound=max_hits_per_compound)
        if product_smiles is not None
        else None
    )

    definitions, _ = load_reaction_catalog()
    matches: List[ReactionMatch] = []

    for definition in definitions.values():
        match = match_reaction_definition(definition, detected_profile, product_profile)
        if match is not None:
            matches.append(match)

    # Validate coupling products by analyzing bond formation patterns
    if confirm_coupling_products and product_smiles is not None:
        coupling_specs = get_coupling_confirmation_specs()
        confirmed: List[ReactionMatch] = []
        for match in matches:
            # Skip if reaction type doesn't support coupling confirmation
            if match.reaction_type not in coupling_specs:
                confirmed.append(match)
                continue
            
            # Validate coupling product using general bond formation analysis
            ok, reason = confirm_coupling_product_by_attachment(
                reactant_smiles,
                product_smiles,
                match.reaction_type,
            )
            if not ok:
                continue
            
            # Add confirmation metadata to slot evidence
            slot_evidence = {**match.slot_evidence}
            if reason:
                slot_evidence["product_confirmation"] = [reason]
            confirmed.append(
                ReactionMatch(
                    reaction_type=match.reaction_type,
                    name=match.name,
                    category=match.category,
                    slot_evidence=slot_evidence,
                    matched_slots=match.matched_slots,
                    required_slots=match.required_slots,
                )
            )
        matches = confirmed

    # Priority map to break ties between overlapping definitions (e.g., Suzuki vs C-H Arylation)
    # Higher values come first.
    _PRIORITIES = {
        "Suzuki_miyaura": 100,
        "C_N_Coupling": 100,
        "C_O_Coupling": 100,
        "C_S_Coupling": 100,
        "Amide_formation": 100,
        "Azide_coupling": 100,  # Specific transformation - prioritize over generic C-H arylation
        "C_H_arylation": 10,
        "Arylation_acidic_C_H": 10,
    }

    matches.sort(
        key=lambda m: (
            -_PRIORITIES.get(m.reaction_type, 50),  # Priority first!
            -m.matched_slots,
            -sum(len(v) for v in m.slot_evidence.values()),
            m.reaction_type,
        )
    )
    return ReactionDetectionResult(matches=matches)
