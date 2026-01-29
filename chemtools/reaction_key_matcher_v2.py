"""
Reaction type detection using reaction keys against existing reaction_types.v4.0.json.

This matcher uses the original taxonomy definitions but matches via reaction key
(reacted/formed motifs) instead of complex slot-filling logic.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


@dataclass
class ReactionKeyMatch:
    """Result of reaction type matching from reaction key."""
    reaction_type: str
    confidence: float
    matched_reacted: Dict[str, List[str]] = field(default_factory=dict)  # slot -> matched motifs
    matched_formed: Dict[str, List[str]] = field(default_factory=dict)
    priority: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "reaction_type": self.reaction_type,
            "confidence": self.confidence,
            "matched_reacted": self.matched_reacted,
            "matched_formed": self.matched_formed,
            "priority": self.priority,
        }


@lru_cache(maxsize=1)
def _load_reaction_types() -> Dict[str, Any]:
    """Load the original reaction_types.v4.0.json taxonomy."""
    path = Path(__file__).resolve().parent / "taxonomy" / "data" / "reaction_types.v4.0.json"
    if not path.exists():
        return {"reaction_types": []}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


@lru_cache(maxsize=1)
def _load_compound_logic() -> Dict[str, Any]:
    """Load compound logic for motif set expansion."""
    path = Path(__file__).resolve().parent / "taxonomy" / "data" / "compound_logic.json"
    if not path.exists():
        return {"motif_sets": {}}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def expand_pattern(pattern: Any, motif_sets: Dict[str, Any]) -> Set[str]:
    """
    Expand a pattern to concrete motif IDs.
    
    Handles:
    - "@set_name" → expand from motif_sets
    - "motif_id" → literal motif
    - ["a", "b", ...] → expand each and union
    """
    if pattern is None:
        return set()
    
    if isinstance(pattern, str):
        if pattern.startswith("@"):
            set_name = pattern[1:]
            set_data = motif_sets.get(set_name, {})
            members = set_data.get("members", [])
            return set(members) if members else set()
        else:
            return {pattern}
    
    if isinstance(pattern, list):
        result: Set[str] = set()
        for item in pattern:
            result.update(expand_pattern(item, motif_sets))
        return result
    
    return set()


def match_reaction_type_v2(
    reacted_motifs: Set[str],
    formed_motifs: Set[str],
    reaction_def: Dict[str, Any],
    motif_sets: Dict[str, Any],
) -> Optional[ReactionKeyMatch]:
    """
    Check if reaction key matches a reaction type definition from v4.0 taxonomy.
    
    Matching logic:
    - For each slot in "reactants", check if any reacted motif is in the expanded set
    - For each slot in "products", check if any formed motif is in the expanded set
    - All slots must have at least one match
    
    Returns match result or None if no match.
    """
    reaction_id = reaction_def.get("id", "")
    reactants_def = reaction_def.get("reactants", {})
    products_def = reaction_def.get("products", {})
    
    matched_reacted: Dict[str, List[str]] = {}
    matched_formed: Dict[str, List[str]] = {}
    
    # Match each reactant slot
    for slot_name, slot_pattern in reactants_def.items():
        slot_set = expand_pattern(slot_pattern, motif_sets)
        if not slot_set:
            continue  # Skip empty patterns
        
        matches = reacted_motifs & slot_set
        if not matches:
            return None  # Required slot not matched
        matched_reacted[slot_name] = sorted(matches)
    
    # Match at least one product slot
    product_matched = False
    for slot_name, slot_pattern in products_def.items():
        slot_set = expand_pattern(slot_pattern, motif_sets)
        if not slot_set:
            continue
        
        matches = formed_motifs & slot_set
        if matches:
            matched_formed[slot_name] = sorted(matches)
            product_matched = True
    
    # Require at least one product match if products are defined
    if products_def and not product_matched:
        return None
    
    # Calculate confidence based on match quality
    total_slots = len(reactants_def) + len(products_def)
    matched_slots = len(matched_reacted) + len(matched_formed)
    
    if total_slots > 0:
        confidence = matched_slots / total_slots
    else:
        confidence = 0.5
    
    # Boost for matching all reactant slots (more specific match)
    if len(matched_reacted) == len(reactants_def) and len(reactants_def) > 1:
        confidence = min(1.0, confidence + 0.15)
    
    return ReactionKeyMatch(
        reaction_type=reaction_id,
        confidence=round(confidence, 2),
        matched_reacted=matched_reacted,
        matched_formed=matched_formed,
        priority=100,  # All from same taxonomy
    )


def detect_from_reaction_key_v2(
    reacted_motifs: List[str],
    formed_motifs: List[str],
    spectator_motifs: Optional[List[str]] = None,
) -> Tuple[Optional[ReactionKeyMatch], List[ReactionKeyMatch]]:
    """
    Detect reaction type from reaction key motifs using v4.0 taxonomy.
    
    Args:
        reacted_motifs: Motifs consumed in the reaction
        formed_motifs: Motifs formed in the product
        spectator_motifs: Unchanged motifs (optional, for future use)
        
    Returns:
        Tuple of (top_match, all_matches) sorted by confidence descending.
    """
    taxonomy = _load_reaction_types()
    logic = _load_compound_logic()
    motif_sets = logic.get("motif_sets", {})
    
    reacted_set = set(reacted_motifs)
    formed_set = set(formed_motifs)
    
    matches: List[ReactionKeyMatch] = []
    
    for reaction_def in taxonomy.get("reaction_types", []):
        match = match_reaction_type_v2(reacted_set, formed_set, reaction_def, motif_sets)
        if match:
            matches.append(match)
    
    # Sort by confidence descending, then by number of matched slots
    matches.sort(key=lambda m: (
        -m.confidence,
        -(len(m.matched_reacted) + len(m.matched_formed))
    ))
    
    top_match = matches[0] if matches else None
    return top_match, matches


def detect_reaction_type_from_smiles_v2(
    reaction_smiles: str,
) -> Tuple[Optional[ReactionKeyMatch], List[ReactionKeyMatch]]:
    """
    Detect reaction type from a reaction SMILES string using v4.0 taxonomy.
    
    Uses the featurizer to extract reacted/formed motifs, then matches.
    """
    from .featurizers.formatters.reaction import featurize_reaction
    
    try:
        result = featurize_reaction(reaction_smiles)
        aggregates = result.get("aggregates", {})
        reacted = aggregates.get("reacted_motifs", [])
        formed = aggregates.get("formed_motifs", [])
        spectators = aggregates.get("spectator_motifs", [])
        
        return detect_from_reaction_key_v2(reacted, formed, spectators)
    except Exception as e:
        return None, []


__all__ = [
    "ReactionKeyMatch",
    "detect_from_reaction_key_v2",
    "detect_reaction_type_from_smiles_v2",
]
