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

from .synthon import classify_reactant_synthons
from .taxonomy.reaction_catalog import motif_tokens_compatible


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
    slot_sources: Dict[str, str] = field(default_factory=dict)
    synthon_slot_evidence: Dict[str, List[str]] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "reaction_type": self.reaction_type,
            "confidence": self.confidence,
            "slots": {
                "electrophile": self.electrophile,
                "nucleophile": self.nucleophile,
                "product": self.product,
            },
        }
        if self.slot_sources:
            payload["slot_sources"] = dict(self.slot_sources)
        if self.synthon_slot_evidence:
            payload["synthon_slot_evidence"] = dict(self.synthon_slot_evidence)
        return payload


@dataclass
class DetectionResult:
    """Detection result with all matches."""
    matches: List[ReactionMatch]
    reacted_motifs: List[str] = field(default_factory=list)
    formed_motifs: List[str] = field(default_factory=list)
    reaction_key: str = ""
    synthon_evidence: Dict[str, Any] = field(default_factory=dict)
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
        if self.synthon_evidence:
            result["synthon_evidence"] = self.synthon_evidence
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
    from .featurizers.formatters.reaction import featurize_reaction, get_crk_options
    
    # Skip role classification to avoid circular dependency
    result = featurize_reaction(reaction_smiles, options=get_crk_options())
    aggregates = result.get("aggregates", {})
    
    reacted = aggregates.get("reacted_motifs", [])
    formed = aggregates.get("formed_motifs", [])
    spectator = aggregates.get("spectator_motifs", [])
    reaction_key = result.get("reaction_key", "")
    
    return reacted, formed, spectator, reaction_key


def _extract_reactant_smiles(reaction_smiles: str) -> List[str]:
    if not reaction_smiles or ">>" not in reaction_smiles:
        return []
    try:
        from .smiles import normalize_reaction
    except Exception:
        normalize_reaction = None  # type: ignore[assignment]
    if normalize_reaction is None:
        left = reaction_smiles.split(">>", 1)[0]
        return [tok for tok in left.split(".") if tok]
    try:
        normalized = normalize_reaction(reaction_smiles)
    except Exception:
        normalized = None
    reactants: List[str] = []
    if isinstance(normalized, dict):
        for entry in normalized.get("reactants", []) or []:
            if not isinstance(entry, dict):
                continue
            smiles = (
                entry.get("smiles_norm")
                or entry.get("largest_smiles")
                or entry.get("input")
            )
            if smiles:
                reactants.append(str(smiles))
    if reactants:
        return reactants
    left = reaction_smiles.split(">>", 1)[0]
    return [tok for tok in left.split(".") if tok]


def _build_synthon_context(reaction_smiles: str) -> Dict[str, Any]:
    reactants = _extract_reactant_smiles(reaction_smiles)
    role_to_motifs: Dict[str, Set[str]] = {"electrophile": set(), "nucleophile": set()}
    role_to_indices: Dict[str, Set[int]] = {"electrophile": set(), "nucleophile": set()}
    reactant_assignments: List[Dict[str, Any]] = []
    all_motifs: Set[str] = set()

    for idx, smiles in enumerate(reactants):
        assignments = classify_reactant_synthons(smiles)
        formatted_assignments: List[Dict[str, Any]] = []
        for hit in assignments:
            matched = sorted({str(m) for m in hit.matched_motifs if str(m)})
            if not matched:
                continue
            role = str(hit.role or "").strip().lower()
            if role in role_to_motifs:
                role_to_motifs[role].update(matched)
                role_to_indices[role].add(idx)
            all_motifs.update(matched)
            formatted_assignments.append(
                {
                    "synthon_id": hit.synthon_id,
                    "role": role,
                    "matched_motifs": matched,
                    "priority": int(hit.priority),
                }
            )
        reactant_assignments.append(
            {
                "index": idx,
                "smiles": smiles,
                "assignments": formatted_assignments,
            }
        )

    role_motifs = {
        role: sorted(values) for role, values in role_to_motifs.items() if values
    }
    role_indices = {
        role: sorted(values) for role, values in role_to_indices.items() if values
    }

    distinct_roles = False
    electrophiles = role_to_indices.get("electrophile", set())
    nucleophiles = role_to_indices.get("nucleophile", set())
    if electrophiles and nucleophiles:
        distinct_roles = any(e != n for e in electrophiles for n in nucleophiles)

    return {
        "reactants": reactants,
        "reactant_assignments": reactant_assignments,
        "role_motifs": role_motifs,
        "role_indices": role_indices,
        "all_motifs": sorted(all_motifs),
        "distinct_roles": distinct_roles,
    }


def _slot_roles(slot_name: str) -> Tuple[str, ...]:
    text = str(slot_name or "").strip().lower()
    if not text:
        return tuple()
    if text == "electrophile" or text.startswith("electrophile_"):
        return ("electrophile",)
    if text == "nucleophile" or text.startswith("nucleophile_"):
        return ("nucleophile",)
    if text == "substrate":
        return ("electrophile",)
    if "electrophile" in text:
        return ("electrophile",)
    if "nucleophile" in text or "partner" in text:
        return ("nucleophile",)
    return tuple()


def _slot_synthon_matches(
    slot_name: str,
    slot_set: Set[str],
    synthon_context: Dict[str, Any],
) -> Set[str]:
    if not slot_set or not synthon_context:
        return set()
    role_motifs = synthon_context.get("role_motifs", {}) or {}
    role_matches: Set[str] = set()
    for role in _slot_roles(slot_name):
        role_matches.update(set(role_motifs.get(role, []) or []))
    if role_matches:
        return _compatible_observed_matches(role_matches, slot_set)
    all_motifs = set(synthon_context.get("all_motifs", []) or [])
    return _compatible_observed_matches(all_motifs, slot_set)


def _compatible_observed_matches(observed: Set[str], allowed: Set[str]) -> Set[str]:
    """Return observed motif tokens that match an allowed token under scope rules."""
    if not observed or not allowed:
        return set()
    matches: Set[str] = set()
    for observed_token in observed:
        if any(motif_tokens_compatible(observed_token, allowed_token) for allowed_token in allowed):
            matches.add(observed_token)
    return matches


# =============================================================================
# Matching Logic
# =============================================================================

def _match_reaction_type(
    reacted: Set[str],
    formed: Set[str],
    reaction_def: Dict[str, Any],
    motif_sets: Dict[str, Set[str]],
    synthon_context: Optional[Dict[str, Any]] = None,
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
    slot_sources: Dict[str, str] = {}
    synthon_slot_evidence: Dict[str, List[str]] = {}
    reactant_slot_total = 0
    
    # Match reactant slots
    for slot_name, slot_pattern in reactants_def.items():
        slot_set = _expand_pattern(slot_pattern, motif_sets)
        if not slot_set:
            continue

        reactant_slot_total += 1
        motif_matches = _compatible_observed_matches(reacted, slot_set)
        synthon_matches = _slot_synthon_matches(
            slot_name,
            slot_set,
            synthon_context or {},
        )
        if motif_matches:
            matches = set(motif_matches)
            slot_sources[slot_name] = "motif"
        elif synthon_matches:
            matches = set(synthon_matches)
            slot_sources[slot_name] = "synthon"
            synthon_slot_evidence[slot_name] = sorted(synthon_matches)
        else:
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
        
        matches = _compatible_observed_matches(formed, slot_set)
        if matches:
            product_matches = sorted(matches)
            break  # Only need one product match
    
    # Require product match if products defined
    if products_def and not product_matches:
        return None
    
    # Calculate confidence
    total_slots = reactant_slot_total + (1 if products_def else 0)
    matched_slots = reactant_slot_total + (1 if product_matches else 0)
    
    confidence = matched_slots / max(total_slots, 1)
    
    # Boost for matching both reactant slots
    if electrophile_matches and nucleophile_matches:
        confidence = min(1.0, confidence + 0.15)
    if synthon_context and synthon_context.get("distinct_roles"):
        confidence = min(1.0, confidence + 0.05)
    
    return ReactionMatch(
        reaction_type=reaction_id,
        confidence=round(confidence, 2),
        electrophile=electrophile_matches,
        nucleophile=nucleophile_matches,
        product=product_matches,
        slot_sources=slot_sources,
        synthon_slot_evidence=synthon_slot_evidence,
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
    synthon_context = _build_synthon_context(reaction_smiles)
    
    reacted_set = set(reacted)
    formed_set = set(formed)
    
    matches: List[ReactionMatch] = []
    for reaction_def in reaction_types:
        match = _match_reaction_type(
            reacted_set,
            formed_set,
            reaction_def,
            motif_sets,
            synthon_context=synthon_context,
        )
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
        synthon_evidence={
            "reactants": synthon_context.get("reactant_assignments", []),
            "role_motifs": synthon_context.get("role_motifs", {}),
            "distinct_roles": bool(synthon_context.get("distinct_roles")),
        },
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
