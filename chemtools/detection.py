"""
Reaction type detection based on reaction key analysis.

This is a clean, simple detection system that:
1. Extracts reacted/formed motifs from reaction SMILES
2. Matches against reaction type definitions
3. Returns ranked matches

No backward compatibility baggage - just clean, focused logic.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


# =============================================================================
# Data Models
# =============================================================================

@dataclass
class ReactionMatch:
    """A matched reaction type with evidence."""
    reaction_type: str
    confidence: float
    electrophile: List[str] = field(default_factory=list)
    nucleophile: List[str] = field(default_factory=list)
    product: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "reaction_type": self.reaction_type,
            "confidence": self.confidence,
            "slots": {
                "electrophile": self.electrophile,
                "nucleophile": self.nucleophile,
                "product": self.product,
            },
        }


@dataclass
class DetectionResult:
    """Detection result with all matches."""
    matches: List[ReactionMatch]
    reacted_motifs: List[str] = field(default_factory=list)
    formed_motifs: List[str] = field(default_factory=list)
    reaction_key: str = ""
    error: Optional[str] = None
    
    @property
    def top_match(self) -> Optional[ReactionMatch]:
        return self.matches[0] if self.matches else None
    
    def to_dict(self) -> Dict[str, Any]:
        result = {
            "matches": [m.to_dict() for m in self.matches],
            "reacted_motifs": self.reacted_motifs,
            "formed_motifs": self.formed_motifs,
            "reaction_key": self.reaction_key,
        }
        if self.error:
            result["error"] = self.error
        return result


# =============================================================================
# Data Loading
# =============================================================================

@lru_cache(maxsize=1)
def _load_reaction_types() -> List[Dict[str, Any]]:
    """Load reaction type definitions."""
    path = Path(__file__).resolve().parent / "taxonomy" / "data" / "reaction_types.v4.0.json"
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    return data.get("reaction_types", [])


@lru_cache(maxsize=1)
def _load_motif_sets() -> Dict[str, Set[str]]:
    """Load and expand motif sets from compound_logic.json."""
    path = Path(__file__).resolve().parent / "taxonomy" / "data" / "compound_logic.json"
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    
    motif_sets = {}
    for set_name, set_data in data.get("motif_sets", {}).items():
        members = set_data.get("members", [])
        motif_sets[set_name] = set(members)
    return motif_sets


def _expand_pattern(pattern: Any, motif_sets: Dict[str, Set[str]]) -> Set[str]:
    """Expand pattern to concrete motif IDs."""
    if pattern is None:
        return set()
    
    if isinstance(pattern, str):
        if pattern.startswith("@"):
            return motif_sets.get(pattern[1:], set())
        return {pattern}
    
    if isinstance(pattern, list):
        result: Set[str] = set()
        for item in pattern:
            result.update(_expand_pattern(item, motif_sets))
        return result
    
    return set()


# =============================================================================
# Reaction Key Extraction
# =============================================================================

def extract_reaction_key(reaction_smiles: str) -> Tuple[List[str], List[str], List[str], str]:
    """
    Extract reacted/formed/spectator motifs from reaction SMILES.
    
    Returns:
        Tuple of (reacted, formed, spectator, reaction_key_string). The key is CRK-v1.
    """
    from .featurizers.formatters.reaction import featurize_reaction
    
    # Skip role classification to avoid circular dependency
    result = featurize_reaction(reaction_smiles, options={"include_roles": False})
    aggregates = result.get("aggregates", {})
    
    reacted = aggregates.get("reacted_motifs", [])
    formed = aggregates.get("formed_motifs", [])
    spectator = aggregates.get("spectator_motifs", [])
    reaction_key = result.get("reaction_key", "")
    
    return reacted, formed, spectator, reaction_key


# =============================================================================
# Matching Logic
# =============================================================================

def _match_reaction_type(
    reacted: Set[str],
    formed: Set[str],
    reaction_def: Dict[str, Any],
    motif_sets: Dict[str, Set[str]],
) -> Optional[ReactionMatch]:
    """
    Match reaction key against a reaction type definition.
    
    Logic:
    - Each slot in "reactants" must have at least one match in reacted motifs
    - At least one slot in "products" must have a match in formed motifs
    """
    reaction_id = reaction_def.get("id", "")
    reactants_def = reaction_def.get("reactants", {})
    products_def = reaction_def.get("products", {})

    # Exclude reactions if any disqualifying motifs are present
    exclude_patterns = reaction_def.get("exclude_reacted") or []
    exclude_set = _expand_pattern(exclude_patterns, motif_sets) if exclude_patterns else set()
    if exclude_set and reacted & exclude_set:
        return None
    
    electrophile_matches = []
    nucleophile_matches = []
    product_matches = []
    
    # Match reactant slots
    for slot_name, slot_pattern in reactants_def.items():
        slot_set = _expand_pattern(slot_pattern, motif_sets)
        if not slot_set:
            continue
        
        matches = reacted & slot_set
        if not matches:
            return None  # Required slot not matched
        
        if slot_name == "electrophile":
            electrophile_matches = sorted(matches)
        elif slot_name == "nucleophile":
            nucleophile_matches = sorted(matches)
        elif slot_name == "substrate":
            electrophile_matches = sorted(matches)  # Substrate acts like electrophile
    
    # Match product slots
    for slot_name, slot_pattern in products_def.items():
        slot_set = _expand_pattern(slot_pattern, motif_sets)
        if not slot_set:
            continue
        
        matches = formed & slot_set
        if matches:
            product_matches = sorted(matches)
            break  # Only need one product match
    
    # Require product match if products defined
    if products_def and not product_matches:
        return None
    
    # Calculate confidence
    total_slots = len(reactants_def) + (1 if products_def else 0)
    matched_slots = (1 if electrophile_matches else 0) + \
                   (1 if nucleophile_matches else 0) + \
                   (1 if product_matches else 0)
    
    confidence = matched_slots / max(total_slots, 1)
    
    # Boost for matching both reactant slots
    if electrophile_matches and nucleophile_matches:
        confidence = min(1.0, confidence + 0.15)
    
    return ReactionMatch(
        reaction_type=reaction_id,
        confidence=round(confidence, 2),
        electrophile=electrophile_matches,
        nucleophile=nucleophile_matches,
        product=product_matches,
    )


# =============================================================================
# Public API
# =============================================================================

def detect_reaction_type(reaction_smiles: str) -> DetectionResult:
    """
    Detect reaction type from reaction SMILES.
    
    This is the main entry point. It:
    1. Extracts reacted/formed motifs (reaction key)
    2. Matches against reaction type definitions
    3. Returns ranked matches
    
    Args:
        reaction_smiles: Reaction SMILES with >> separator
        
    Returns:
        DetectionResult with ranked matches and reaction key info
    """
    if not reaction_smiles or ">>" not in reaction_smiles:
        return DetectionResult(matches=[], error="invalid_reaction_smiles")
    
    try:
        reacted, formed, spectator, reaction_key = extract_reaction_key(reaction_smiles)
    except Exception as e:
        return DetectionResult(matches=[], error=str(e))
    
    if not reacted and not formed:
        return DetectionResult(
            matches=[],
            reacted_motifs=reacted,
            formed_motifs=formed,
            reaction_key=reaction_key,
            error="no_motif_changes",
        )
    
    reaction_types = _load_reaction_types()
    motif_sets = _load_motif_sets()
    
    reacted_set = set(reacted)
    formed_set = set(formed)
    
    matches: List[ReactionMatch] = []
    for reaction_def in reaction_types:
        match = _match_reaction_type(reacted_set, formed_set, reaction_def, motif_sets)
        if match:
            matches.append(match)
    
    # Sort by confidence, then by specificity (number of matched slots)
    matches.sort(key=lambda m: (
        -m.confidence,
        -(len(m.electrophile) + len(m.nucleophile) + len(m.product)),
    ))
    
    return DetectionResult(
        matches=matches,
        reacted_motifs=reacted,
        formed_motifs=formed,
        reaction_key=reaction_key,
    )


def detect_reaction_types(reaction_smiles: str) -> DetectionResult:
    """Alias for detect_reaction_type (for compatibility)."""
    return detect_reaction_type(reaction_smiles)


__all__ = [
    "ReactionMatch",
    "DetectionResult",
    "detect_reaction_type",
    "detect_reaction_types",
    "extract_reaction_key",
]
