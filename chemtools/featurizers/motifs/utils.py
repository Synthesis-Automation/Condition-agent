"""SMARTS manipulation and validation utilities for motif detection."""

from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Mapping

from .models import _MAP_RE


def _format_compound_smarts(
    *,
    template_format: str,
    a_smarts: str,
    b_smarts: str,
) -> str:
    """Format compound SMARTS from template using group patterns."""
    return template_format.format(A=a_smarts, B=b_smarts)


def _has_atom_map(smarts: str) -> bool:
    """Check if SMARTS contains any atom map numbers."""
    return bool(_MAP_RE.search(smarts))


def _has_map(smarts: str, *, map_num: int) -> bool:
    """Check if SMARTS contains a specific atom map number."""
    if not smarts:
        return False
    return bool(re.search(rf":{map_num}(?=\])", smarts))


def _inject_map_on_first_atom(smarts: str, *, map_num: int) -> str:
    """Add atom map number to the first atom in a SMARTS pattern."""
    if not smarts:
        return smarts
    if smarts.startswith("["):
        end = smarts.find("]")
        if end == -1:
            return smarts
        return f"{smarts[:end]}:{map_num}{smarts[end:]}"

    match = re.match(r"([A-Za-z][a-z]?)", smarts)
    if not match:
        return smarts
    atom = match.group(1)
    return f"[{atom}:{map_num}]{smarts[len(atom):]}"


def _extract_compound_smarts(entry: Mapping[str, Any]) -> List[str]:
    """Extract SMARTS patterns from compound entry.
    
    Handles both direct "smarts" field and "smarts_any" list field.
    """
    direct_smarts = entry.get("smarts")
    smarts_any = entry.get("smarts_any")
    smarts_list: List[str] = []
    if isinstance(direct_smarts, str):
        smarts_list = [direct_smarts]
    elif isinstance(direct_smarts, list):
        smarts_list = [s for s in direct_smarts if isinstance(s, str)]
    elif isinstance(smarts_any, list):
        smarts_list = [s for s in smarts_any if isinstance(s, str)]
    return [s.strip() for s in smarts_list if isinstance(s, str) and s.strip()]


def _validate_group_maps(groups: Mapping[str, Mapping[str, Any]]) -> None:
    """Validate that group SMARTS have required atom map numbers.
    
    Scaffolds (kind="scaffold") must have :1 map.
    Substituents (kind="substituent") must have :2 map.
    
    Raises:
        ValueError: If any groups are missing required maps
    """
    errors: List[str] = []
    for group_id, group in groups.items():
        if not isinstance(group, Mapping):
            continue
        kind = str(group.get("kind") or "")
        smarts = str(group.get("smarts") or "")
        if not kind or not smarts:
            continue
        expected_map = 1 if kind == "scaffold" else 2
        if not _has_map(smarts, map_num=expected_map):
            errors.append(f"{group_id} (kind={kind}, expected :{expected_map})")
    if errors:
        joined = ", ".join(sorted(errors))
        raise ValueError(f"Group SMARTS missing required attach maps: {joined}")


def _validate_compound_templates(
    compounds: Iterable[Dict[str, Any]],
    templates: Mapping[str, str],
) -> None:
    """Validate that all compound templates are defined.
    
    Raises:
        ValueError: If any compound references undefined templates
    """
    missing: set[str] = set()
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        if _extract_compound_smarts(entry):
            continue
        template_id = entry.get("template")
        if template_id and template_id not in templates:
            missing.add(str(template_id))
    if missing:
        joined = ", ".join(sorted(missing))
        raise ValueError(f"Compound templates missing from registry: {joined}")


def calculate_smarts_complexity(query_mol: Any) -> int:
    """Calculate a structural complexity score for a SMARTS query molecule.
    
    Higher scores indicate more specific patterns. Scoring factors:
    - Base score per atom/bond
    - Explicit element specification
    - Constraints (hybridization, ring membership, H-count)
    - Stereochemistry
    - Multiple bonds
    - Negations (increase specificity)
    - OR operators (decrease specificity)
    
    Args:
        query_mol: RDKit query molecule from compiled SMARTS
        
    Returns:
        Complexity score (higher = more specific)
    """
    if not query_mol:
        return 0
    score = 0
    # Score atoms
    for atom in query_mol.GetAtoms():
        score += 10  # Base atom
        symbol = atom.GetSymbol()
        # Non-wildcard/Non-H atoms are more specific
        if symbol != "*" and symbol != "H":
            score += 10  # Heavy atom bonus
        if symbol != "*":
            score += 5  # Explicit element

        # Parse SMARTS constraints for specificity
        smarts = atom.GetSmarts()
        score += len(smarts)
        if "X" in smarts:
            score += 5  # Hybridization constraint
        if "R" in smarts or "r" in smarts:
            score += 5  # Ring constraint
        if "H" in smarts or "h" in smarts:
            score += 3  # Explicit H-count
        if "!" in smarts:
            score += 2  # Negations
        score += smarts.count(";") * 2  # AND connections
        score -= smarts.count(",") * 2  # OR connections (genericness)

    # Score bonds
    for bond in query_mol.GetBonds():
        score += 5  # Base bond
        bond_smarts = bond.GetSmarts()
        score += len(bond_smarts)
        if bond.GetBondTypeAsDouble() > 1:
            score += 3  # Multiple bonds are more specific
        if bond_smarts and "@" in bond_smarts:
            score += 5  # Stereochemical constraints

    return score


__all__ = [
    "_format_compound_smarts",
    "_has_atom_map",
    "_has_map",
    "_inject_map_on_first_atom",
    "_extract_compound_smarts",
    "_validate_group_maps",
    "_validate_compound_templates",
    "calculate_smarts_complexity",
]
