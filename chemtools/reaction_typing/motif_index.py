"""
Motif query index utilities for reaction typing.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from chemtools.featurizers.motif_registry import build_compound_detect_registry
from chemtools.util.rdkit_helpers import rdkit_available


def build_query_index(
    *,
    groups_path: str | Path | None = None,
    compounds_path: str | Path | None = None,
    templates_path: str | Path | None = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
) -> Dict[str, List[Any]]:
    """
    Build a detect-SMARTS query index from groups/compounds.
    """
    if not rdkit_available():
        return {}

    if registry_paths is None:
        registry_paths = _default_registry_paths()
        if groups_path is not None:
            registry_paths = dict(registry_paths)
            registry_paths["groups"] = groups_path
        if compounds_path is not None:
            registry_paths = dict(registry_paths)
            registry_paths["compounds"] = compounds_path
        if templates_path is not None:
            registry_paths = dict(registry_paths)
            registry_paths["templates"] = templates_path
    else:
        registry_paths = dict(registry_paths)
        if "templates" not in registry_paths:
            registry_paths["templates"] = _default_registry_paths()["templates"]

    registry = build_compound_detect_registry(registry_paths)
    return dict(registry.get("compiled_compounds") or {})


def type_molecule(mol: Any, query_index: Mapping[str, Iterable[Any]]) -> List[str]:
    """
    Return motif IDs matched by the molecule.
    """
    hits: List[str] = []
    for motif_id, queries in query_index.items():
        if any(mol.HasSubstructMatch(query) for query in queries):
            hits.append(motif_id)
    return hits


def _default_registry_paths() -> Dict[str, Path]:
    base = Path(__file__).resolve().parents[1] / "taxonomy" / "v2_data"
    return {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
        "templates": base / "smarts_templates.v1.json",
    }
