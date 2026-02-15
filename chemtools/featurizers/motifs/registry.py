"""Registry loading and compilation for motif detection system."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

from chemtools.taxonomy.substituent_composer import load_organic_groups_with_compositions
from chemtools.util.smarts_cache import compile_smarts

from .models import CompoundPattern, _DEFAULT_TEMPLATES
from .utils import (
    _extract_compound_smarts,
    _format_compound_smarts,
    _has_atom_map,
    _inject_map_on_first_atom,
    _validate_compound_templates,
    _validate_group_maps,
    calculate_smarts_complexity,
)


def _registry_cache_key(
    registry_paths: Mapping[str, str | Path],
) -> Tuple[Tuple[str, str], ...]:
    """Create a stable cache key from registry paths."""
    key_parts: List[Tuple[str, str]] = []
    for name, path_value in registry_paths.items():
        resolved = str(Path(path_value).resolve())
        key_parts.append((str(name), resolved))
    key_parts.sort(key=lambda item: item[0])
    return tuple(key_parts)


@lru_cache(maxsize=32)
def _build_compound_registry_cached(
    cache_key: Tuple[Tuple[str, str], ...],
) -> Dict[str, Any]:
    path_map = {name: Path(path_str) for name, path_str in cache_key}
    return _build_compound_registry_uncached(path_map)


def build_compound_registry(registry_paths: Mapping[str, str | Path]) -> Dict[str, Any]:
    """Load and compile motif registry from taxonomy files.

    Results are cached by registry path set to avoid repeated JSON loading and
    SMARTS compilation across per-molecule/per-reaction calls.
    """
    return _build_compound_registry_cached(_registry_cache_key(registry_paths))


def _build_compound_registry_uncached(
    registry_paths: Mapping[str, str | Path],
) -> Dict[str, Any]:
    """Load and compile motif registry from taxonomy files.
    
    This is the main entry point for building the motif detection system.
    Loads groups, templates, and compounds, then compiles SMARTS patterns
    and calculates priorities and complexity scores.
    
    Args:
        registry_paths: Paths to registry files:
            - "groups": organic_groups.v1.3.json (required)
            - "compounds": organic_compounds.v1.3.json (required)
            - "templates": template patterns (optional, uses defaults if missing)
            - "scaffold_motifs": scaffold_motifs.v1.3.json (optional)
            - "logic": group_logic.json (optional, for priority overrides)
            
    Returns:
        Dictionary with compiled registry:
        - groups: Raw group definitions
        - templates: Template formats
        - compiled_compounds: List of CompoundPattern (sorted by priority/complexity)
        - compound_map: compound_id -> CompoundPattern
        - combination_map: (group_a, group_b) -> compound_id
        - priorities: group_id -> priority score
        - compiled_groups: Compiled group queries
    """
    groups_path = Path(registry_paths["groups"])
    templates_path = (
        Path(registry_paths["templates"]) if "templates" in registry_paths else None
    )
    compounds_path = Path(registry_paths["compounds"])
    logic_path = Path(registry_paths["logic"]) if "logic" in registry_paths else None

    groups = _load_groups(groups_path)
    templates = _load_templates(templates_path)
    compounds = _load_compounds(compounds_path)

    if "scaffold_motifs" in registry_paths:
        scaffold_path = Path(registry_paths["scaffold_motifs"])
        compounds.extend(_load_compounds(scaffold_path))

    # Load priorities from groups and logic sets
    priorities = {g_id: g.get("priority", 1) for g_id, g in groups.items()}
    if logic_path and logic_path.exists():
        logic_data = _load_json(logic_path)
        group_sets = logic_data.get("group_sets", {})
        for set_id, set_data in group_sets.items():
            if "priority" in set_data:
                priorities[set_id] = set_data["priority"]

    _validate_group_maps(groups)
    _validate_compound_templates(compounds, templates)

    compiled: List[CompoundPattern] = []
    compound_map: Dict[str, CompoundPattern] = {}
    combination_map: Dict[Tuple[str, str], str] = {}

    for entry in compounds:
        compound_id = entry.get("id")
        group_a = str(entry.get("A") or "")
        group_b = str(entry.get("B") or "")

        # Auto-generate ID as A+B when missing (v1.3 simplified format)
        if not compound_id:
            if group_a and group_b:
                # Avoid double dash if group_b starts with dash
                if group_b.startswith("-"):
                    compound_id = f"{group_a}{group_b}"
                else:
                    compound_id = f"{group_a}-{group_b}"
            else:
                continue
        if group_a and group_b:
            combination_map[(group_a, group_b)] = compound_id

        # Calculate priority: Priority(A) + Priority(B), unless overridden
        priority = entry.get("priority")
        if priority is None:
            priority = priorities.get(group_a, 1) + priorities.get(group_b, 1)
        reactivity_weight = float(entry.get("reactivity_weight") or 0.0)

        smarts_list = _extract_compound_smarts(entry)
        if smarts_list:
            for smarts in smarts_list:
                if not _has_atom_map(smarts):
                    smarts = _inject_map_on_first_atom(smarts, map_num=1)
                query = compile_smarts(smarts, validate=False)
                if query is None:
                    continue
                pattern = CompoundPattern(
                    compound_id=compound_id,
                    smarts=smarts,
                    query=query,
                    group_a=group_a,
                    group_b=group_b,
                    b_tags=list(entry.get("tags") or []),
                    priority=priority,
                    complexity=calculate_smarts_complexity(query),
                    reactivity_weight=reactivity_weight,
                )
                compiled.append(pattern)
                compound_map[compound_id] = pattern
            continue

        template_id = entry.get("template", "")
        template_format = templates.get(template_id)
        if not template_format:
            continue

        if not group_a or not group_b:
            continue
        group_a_record = groups.get(group_a)
        group_b_record = groups.get(group_b)

        if not group_a_record or not group_b_record:
            continue

        a_smarts = group_a_record.get("smarts", "")
        b_smarts = group_b_record.get("smarts", "")
        if not a_smarts or not b_smarts:
            continue

        compound_smarts = _format_compound_smarts(
            template_format=template_format,
            a_smarts=a_smarts,
            b_smarts=b_smarts,
        )
        query = compile_smarts(compound_smarts, validate=False)
        if query is None:
            continue

        pattern = CompoundPattern(
            compound_id=compound_id,
            smarts=compound_smarts,
            query=query,
            group_a=group_a,
            group_b=group_b,
            b_tags=list(group_b_record.get("tags") or []),
            priority=priority,
            complexity=calculate_smarts_complexity(query),
            reactivity_weight=reactivity_weight,
        )
        compiled.append(pattern)
        compound_map[compound_id] = pattern

    # Re-sort compiled patterns by priority and complexity
    # Higher priority first. If priority tied, more complex (narrower) SMARTS first.
    def sort_key(p: CompoundPattern) -> tuple[int, int, str]:
        return (-p.priority, -p.complexity, p.compound_id)

    compiled.sort(key=sort_key)

    compiled_groups: Dict[str, Dict[str, Any]] = {}
    for group_id, group in groups.items():
        smarts = group.get("smarts")
        if smarts:
            query = compile_smarts(smarts, validate=False)
            if query:
                compiled_groups[group_id] = {
                    "id": group_id,
                    "kind": group.get("kind"),
                    "query": query,
                    "priority": priorities.get(group_id, 1),
                    "complexity": calculate_smarts_complexity(query),
                    "tags": list(group.get("tags") or []),
                }

    return {
        "groups": groups,
        "templates": templates,
        "compiled_compounds": compiled,
        "compound_map": compound_map,
        "combination_map": combination_map,
        "priorities": priorities,
        "compiled_groups": compiled_groups,
    }


def clear_compound_registry_cache() -> None:
    """Clear in-memory registry cache (useful for tests or taxonomy reloads)."""
    _build_compound_registry_cached.cache_clear()


def _load_groups(path: Path) -> Dict[str, Dict[str, Any]]:
    """Load group definitions from JSON file."""
    payload = load_organic_groups_with_compositions(path)
    groups = payload.get("groups", [])
    return {
        group["id"]: group
        for group in groups
        if isinstance(group, dict) and group.get("id")
    }


def _load_templates(path: Optional[Path]) -> Dict[str, str]:
    """Load template formats or return defaults."""
    if path is None or not path.exists():
        return dict(_DEFAULT_TEMPLATES)

    payload = _load_json(path)
    templates = payload.get("templates", {})
    formatted: Dict[str, str] = dict(_DEFAULT_TEMPLATES)
    for template_id, entry in templates.items():
        if isinstance(entry, dict) and "format" in entry:
            formatted[template_id] = str(entry["format"])
    return formatted


def _load_compounds(path: Path) -> List[Dict[str, Any]]:
    """Load compound definitions from JSON file.
    
    Auto-generates compound IDs from A+B if not specified.
    """
    payload = _load_json(path)
    compounds = list(payload.get("compounds") or [])
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        if "anchors" not in entry:
            entry["anchors"] = {"scaffold": "A", "substituent": "B"}
        # Auto-generate compound ID from A+B
        if "id" not in entry:
            group_a = str(entry.get("A") or "")
            group_b = str(entry.get("B") or "")
            if group_a and group_b:
                # Add hyphen separator if B doesn't already start with "-"
                # Substituents have "-" prefix (e.g., "-Cl"), scaffolds don't (e.g., "Alkenyl")
                separator = "" if group_b.startswith("-") else "-"
                entry["id"] = f"{group_a}{separator}{group_b}"
    return compounds


def _load_json(path: Path) -> Dict[str, Any]:
    """Load and parse JSON file."""
    if not path.exists():
        raise FileNotFoundError(f"Motif registry file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _default_registry_paths() -> Dict[str, Path]:
    """Return default paths to taxonomy files."""
    base = Path(__file__).resolve().parents[2] / "taxonomy" / "data"
    return {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
        "scaffold_motifs": base / "scaffold_motifs.v1.3.json",
        "logic": base / "group_logic.json",
    }


__all__ = [
    "build_compound_registry",
    "clear_compound_registry_cache",
    "_load_groups",
    "_load_templates",
    "_load_compounds",
    "_load_json",
    "_default_registry_paths",
]
