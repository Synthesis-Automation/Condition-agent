"""
Reaction type detection based purely on reaction keys (reacted/formed motifs).

This is a simplified detection system that maps reaction keys directly to
reaction types without complex slot-matching logic.
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
    matched_reacted: List[str] = field(default_factory=list)
    matched_formed: List[str] = field(default_factory=list)
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
def _load_reaction_key_taxonomy() -> Dict[str, Any]:
    """Load the reaction key-based taxonomy."""
    path = Path(__file__).resolve().parent / "taxonomy" / "data" / "reaction_key_taxonomy.json"
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


def expand_motif_set(pattern: str, motif_sets: Dict[str, Any]) -> Set[str]:
    """
    Expand a pattern to concrete motif IDs.
    
    If pattern starts with @, it references a motif set.
    Otherwise it's a literal motif ID.
    """
    if not pattern.startswith("@"):
        return {pattern}
    
    set_name = pattern[1:]  # Remove @
    set_data = motif_sets.get(set_name, {})
    members = set_data.get("members", [])
    return set(members) if members else set()


def expand_patterns(patterns: List[str], motif_sets: Dict[str, Any]) -> Set[str]:
    """Expand a list of patterns to all concrete motif IDs."""
    result: Set[str] = set()
    for pattern in patterns:
        result.update(expand_motif_set(pattern, motif_sets))
    return result


def match_reaction_type(
    reacted_motifs: Set[str],
    formed_motifs: Set[str],
    reaction_def: Dict[str, Any],
    motif_sets: Dict[str, Any],
) -> Optional[ReactionKeyMatch]:
    """
    Check if reaction key matches a reaction type definition.
    
    Matching logic:
    - reacted.any_of: At least one reacted motif must be in this set
    - reacted.all_of: At least one reacted motif from each set must be present
    - formed.any_of: At least one formed motif must be in this set
    - exclude_reacted: None of these motifs should be in reacted set
    
    Returns match result or None if no match.
    """
    reaction_id = reaction_def.get("id", "")
    priority = reaction_def.get("priority", 50)
    
    reacted_def = reaction_def.get("reacted", {})
    formed_def = reaction_def.get("formed", {})
    
    matched_reacted: List[str] = []
    matched_formed: List[str] = []
    
    # Check exclude_reacted first (disqualifiers)
    exclude_patterns = reaction_def.get("exclude_reacted", [])
    if exclude_patterns:
        excluded_set = expand_patterns(exclude_patterns, motif_sets)
        if reacted_motifs & excluded_set:
            return None
    
    # Check reacted.any_of - need at least one match
    any_of_reacted = reacted_def.get("any_of", [])
    if any_of_reacted:
        any_of_set = expand_patterns(any_of_reacted, motif_sets)
        matches = reacted_motifs & any_of_set
        if not matches:
            return None
        matched_reacted.extend(sorted(matches))
    
    # Check reacted.all_of - need at least one from each set
    all_of_reacted = reacted_def.get("all_of", [])
    for pattern in all_of_reacted:
        pattern_set = expand_motif_set(pattern, motif_sets)
        matches = reacted_motifs & pattern_set
        if not matches:
            return None
        # Add only one representative from each required set
        matched_reacted.append(sorted(matches)[0])
    
    # Check formed.any_of - need at least one match
    any_of_formed = formed_def.get("any_of", [])
    if any_of_formed:
        any_of_set = expand_patterns(any_of_formed, motif_sets)
        matches = formed_motifs & any_of_set
        if not matches:
            return None
        matched_formed.extend(sorted(matches))
    
    # Calculate confidence based on match quality
    total_required = len(any_of_reacted) + len(all_of_reacted) + len(any_of_formed)
    total_matched = len(set(matched_reacted)) + len(set(matched_formed))
    
    if total_required > 0:
        confidence = min(1.0, total_matched / max(1, total_required))
    else:
        confidence = 0.5  # Weak match if no requirements
    
    # Boost confidence for more specific matches
    if len(all_of_reacted) > 0:
        confidence = min(1.0, confidence + 0.1)  # Bonus for stricter matching
    
    return ReactionKeyMatch(
        reaction_type=reaction_id,
        confidence=round(confidence, 2),
        matched_reacted=list(set(matched_reacted)),
        matched_formed=list(set(matched_formed)),
        priority=priority,
    )


def detect_from_reaction_key(
    reacted_motifs: List[str],
    formed_motifs: List[str],
    spectator_motifs: Optional[List[str]] = None,
) -> Tuple[Optional[ReactionKeyMatch], List[ReactionKeyMatch]]:
    """
    Detect reaction type from reaction key motifs.
    
    Args:
        reacted_motifs: Motifs consumed in the reaction
        formed_motifs: Motifs formed in the product
        spectator_motifs: Unchanged motifs (optional, for context)
        
    Returns:
        Tuple of (top_match, all_matches) where matches are sorted by
        (priority descending, confidence descending).
    """
    taxonomy = _load_reaction_key_taxonomy()
    logic = _load_compound_logic()
    motif_sets = logic.get("motif_sets", {})
    
    reacted_set = set(reacted_motifs)
    formed_set = set(formed_motifs)
    
    matches: List[ReactionKeyMatch] = []
    
    for reaction_def in taxonomy.get("reaction_types", []):
        match = match_reaction_type(reacted_set, formed_set, reaction_def, motif_sets)
        if match:
            matches.append(match)
    
    # Sort by priority (desc), then confidence (desc)
    matches.sort(key=lambda m: (-m.priority, -m.confidence))
    
    top_match = matches[0] if matches else None
    return top_match, matches


def detect_reaction_type_from_key(reaction_key: str) -> Tuple[Optional[ReactionKeyMatch], List[ReactionKeyMatch]]:
    """
    Parse a reaction key string and detect reaction type.
    
    Reaction key format: "reacted1|reacted2 -> formed1|formed2 || spectator1|spectator2"
    """
    # Parse the reaction key
    parts = reaction_key.split(" || ")
    main_part = parts[0] if parts else ""
    spectator_part = parts[1] if len(parts) > 1 else ""
    
    arrow_parts = main_part.split(" -> ")
    reacted_str = arrow_parts[0] if arrow_parts else ""
    formed_str = arrow_parts[1] if len(arrow_parts) > 1 else ""
    
    # Parse motif lists
    def parse_motifs(s: str) -> List[str]:
        if not s or s == "[]":
            return []
        return [m.strip() for m in s.split("|") if m.strip()]
    
    reacted = parse_motifs(reacted_str)
    formed = parse_motifs(formed_str)
    spectators = parse_motifs(spectator_part)
    
    return detect_from_reaction_key(reacted, formed, spectators)


# Convenience function for direct SMILES-based detection
def detect_reaction_type_from_smiles(
    reaction_smiles: str,
) -> Tuple[Optional[ReactionKeyMatch], List[ReactionKeyMatch]]:
    """
    Detect reaction type from a reaction SMILES string.
    
    Uses the featurizer to extract reacted/formed motifs, then matches.
    """
    from .featurizers.formatters.reaction import featurize_reaction
    
    try:
        result = featurize_reaction(reaction_smiles)
        aggregates = result.get("aggregates", {})
        reacted = aggregates.get("reacted_motifs", [])
        formed = aggregates.get("formed_motifs", [])
        spectators = aggregates.get("spectator_motifs", [])
        
        return detect_from_reaction_key(reacted, formed, spectators)
    except Exception:
        return None, []


__all__ = [
    "ReactionKeyMatch",
    "detect_from_reaction_key",
    "detect_reaction_type_from_key",
    "detect_reaction_type_from_smiles",
]
