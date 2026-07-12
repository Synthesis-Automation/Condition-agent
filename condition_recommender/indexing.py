"""Loading and indexing of verified converted reaction records."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .condition_normalization import split_identifiers
from .models import ConditionIdentity


@dataclass(frozen=True)
class IndexedReaction:
    reaction_id: str
    reaction_smiles: str
    yield_pct: float
    conditions: ConditionIdentity
    recipe_id: str
    electrophile: Dict[str, Any]
    transfer_partner: Dict[str, Any]
    spectator_group_ids: Tuple[str, ...]
    family_flags: Tuple[str, ...]
    product_connection_label: str


def load_verified_index(path: str | Path) -> Tuple[IndexedReaction, ...]:
    """Load only converter-verified rows from a flat pilot table."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    indexed: List[IndexedReaction] = []
    for row in rows:
        if row.get("admission_tier") != "verified":
            continue
        environment = json.loads(row.get("family_environment_json") or "{}")
        partners = {item.get("role"): item for item in environment.get("partners") or []}
        if not {"electrophile", "transfer_partner"} <= set(partners):
            continue
        indexed.append(IndexedReaction(
            reaction_id=str(row["reaction_id"]),
            reaction_smiles=str(row["reaction_smiles"]),
            yield_pct=float(row["yield_pct"]),
            conditions=ConditionIdentity(
                catalyst_cas=split_identifiers(row.get("catalyst_cas")),
                reagent_cas=split_identifiers(row.get("reagent_cas")),
                solvent_cas=split_identifiers(row.get("solvent_cas")),
            ),
            recipe_id=str(row["condition_recipe_id"]),
            electrophile=dict(partners["electrophile"]),
            transfer_partner=dict(partners["transfer_partner"]),
            spectator_group_ids=split_identifiers(row.get("spectator_group_ids")),
            family_flags=tuple(environment.get("flags") or ()),
            product_connection_label=str(row.get("product_connection_label") or ""),
        ))
    return tuple(indexed)


__all__ = ["IndexedReaction", "load_verified_index"]
