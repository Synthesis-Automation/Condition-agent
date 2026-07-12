"""Interpretable, site-relative molecular environment descriptors."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

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


def _steric(mol: Any, center: int) -> Dict[str, Any]:
    atom = mol.GetAtomWithIdx(center)
    heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
    local = {center}
    frontier = {center}
    for _ in range(int(_environment_rules().get("steric_radius", 2))):
        frontier = {
            n.GetIdx() for idx in frontier for n in mol.GetAtomWithIdx(idx).GetNeighbors()
            if n.GetAtomicNum() > 1
        } - local
        local.update(frontier)
    if atom.GetIsAromatic():
        ring_neighbors = [n for n in heavy_neighbors if n.GetIsAromatic()]
        ortho_substituted = sum(
            1 for neighbor in ring_neighbors
            if any(x.GetAtomicNum() > 1 and not x.GetIsAromatic() for x in neighbor.GetNeighbors())
        )
        steric_class = "ortho_hindered" if ortho_substituted >= 2 else ("ortho_substituted" if ortho_substituted else "open")
        return {"class": steric_class, "ortho_substituent_count": ortho_substituted, "local_heavy_atoms_r2": len(local), "method": "graph_local_v1"}
    carbon_neighbors = sum(n.GetAtomicNum() == 6 for n in heavy_neighbors)
    classes = {0: "methyl", 1: "primary", 2: "secondary"}
    return {"class": classes.get(carbon_neighbors, "tertiary"), "carbon_neighbor_count": carbon_neighbors, "local_heavy_atoms_r2": len(local), "method": "graph_local_v1"}


def build_site_environment(mol: Any, site: ReactiveSite, groups: Iterable[FunctionalGroup]) -> SiteEnvironment:
    """Build raw local descriptors without assuming a reaction mechanism."""
    roles = site.details.get("atom_roles") or {}
    preferred = roles.get("anchor") or roles.get("electrophile") or roles.get("center") or roles.get("heteroatom")
    center = int((preferred or site.atom_indices)[0])
    nearby: List[Dict[str, Any]] = []
    electronic_sum = 0.0
    rules = _environment_rules()
    radius = int(rules.get("local_group_radius", 3))
    tag_weights = dict(rules.get("electronic_tag_weights") or {})
    for group in groups:
        distance = _distance(mol, center, group.atom_indices)
        if distance is None:
            continue
        if distance <= radius:
            nearby.append({"group_id": group.group_id, "label": group.chemist_label, "distance": distance, "tags": list(group.tags)})
        if 0 < distance <= radius:
            attenuation = 1.0 / distance
            electronic_sum += attenuation * sum(float(tag_weights.get(tag, 0.0)) for tag in group.tags)
    threshold = float(rules.get("electronic_threshold", 0.35))
    electronic_class = "electron_poor" if electronic_sum > threshold else ("electron_rich" if electronic_sum < -threshold else "neutral")
    atom = mol.GetAtomWithIdx(center)
    return SiteEnvironment(
        site_id=site.site_id,
        center_atom_index=center,
        first_shell=tuple(sorted(n.GetSymbol() for n in atom.GetNeighbors())),
        nearby_groups=tuple(sorted(nearby, key=lambda item: (item["distance"], item["group_id"]))),
        steric=_steric(mol, center),
        electronic={"class": electronic_class, "qualitative_sum": round(electronic_sum, 3), "method": "tag_distance_v1"},
    )


__all__ = ["build_site_environment"]
