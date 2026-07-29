"""Interpretable, site-relative molecular environment descriptors."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List

from .descriptors.profiles import build_site_reactivity_profile
from .models import FunctionalGroup, ReactiveSite, SiteEnvironment

_RULES_PATH = Path(__file__).with_name("definitions") / "descriptor_rules.v1.json"


@lru_cache(maxsize=1)
def _environment_rules() -> Dict[str, Any]:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        return dict((json.load(handle).get("site_environment") or {}))


def _distance(mol: Any, start: int, targets: Iterable[int]) -> int | None:
    from rdkit import Chem
    values = []
    for target in targets:
        if int(start) == int(target):
            values.append(0)
            continue
        try:
            values.append(len(Chem.GetShortestPath(mol, int(start), int(target))) - 1)
        except Exception:
            continue
    return min(values) if values else None


def build_site_environment(mol: Any, site: ReactiveSite, groups: Iterable[FunctionalGroup]) -> SiteEnvironment:
    """Build raw local descriptors without assuming a reaction mechanism."""
    roles = site.details.get("atom_roles") or {}
    preferred = roles.get("anchor") or roles.get("electrophile") or roles.get("center") or roles.get("heteroatom")
    center = int((preferred or site.atom_indices)[0])
    group_values = tuple(groups)
    nearby: List[Dict[str, Any]] = []
    rules = _environment_rules()
    radius = int(rules.get("local_group_radius", 3))
    for group in group_values:
        distance = _distance(mol, center, group.atom_indices)
        if distance is None:
            continue
        if distance <= radius:
            nearby.append({"group_id": group.group_id, "label": group.chemist_label, "distance": distance, "tags": list(group.tags)})
    atom = mol.GetAtomWithIdx(center)
    nearby_groups = tuple(
        sorted(nearby, key=lambda item: (item["distance"], item["group_id"]))
    )
    reactivity_profile = build_site_reactivity_profile(
        mol,
        site,
        group_values,
        nearby_groups,
        center_atom_index=center,
    )
    return SiteEnvironment(
        site_id=site.site_id,
        center_atom_index=center,
        first_shell=tuple(sorted(n.GetSymbol() for n in atom.GetNeighbors())),
        nearby_groups=nearby_groups,
        reactivity_profile=reactivity_profile,
    )


__all__ = ["build_site_environment"]
