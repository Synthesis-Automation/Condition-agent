"""Read the package-owned substance and role/family definitions."""

from __future__ import annotations

import csv
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterator, Tuple

from .models import RoleAssignment, Substance
from .normalization import normalize_cas

DEFINITIONS_DIR = Path(__file__).with_name("definitions")
SUBSTANCES_PATH = DEFINITIONS_DIR / "substances.v1.csv"
PENDING_PATH = DEFINITIONS_DIR / "pending_substances.csv"
TAXONOMY_PATH = DEFINITIONS_DIR / "roles_families.v1.json"


@lru_cache(maxsize=1)
def load_taxonomy() -> Dict[str, Any]:
    with TAXONOMY_PATH.open("r", encoding="utf-8") as handle:
        return dict(json.load(handle))


def iter_substance_rows() -> Iterator[Dict[str, str]]:
    with SUBSTANCES_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle)


def _role_assignments(row: Dict[str, str]) -> Tuple[RoleAssignment, ...]:
    roles = []
    for index in (1, 2):
        role = str(row.get(f"role_{index}") or "").strip()
        if role:
            roles.append(RoleAssignment(
                role_id=role,
                family_id=str(row.get(f"family_{index}") or "").strip() or None,
                tag=str(row.get(f"tag_{index}") or "").strip() or None,
            ))
    return tuple(roles)


def row_to_substance(row: Dict[str, str], row_number: int) -> Substance:
    name = str(row.get("name") or "").strip()
    cas_raw = str(row.get("cas") or "").strip()
    cas = normalize_cas(cas_raw)
    abbreviation = str(row.get("abbreviation") or "").strip().strip('"')
    identity = f"cas:{cas}" if cas else f"registry-row:{row_number}"
    core = {"name", "abbreviation", "cas", "smiles", "role_1", "family_1", "tag_1", "role_2", "family_2", "tag_2"}
    properties = {key: str(value).strip() for key, value in row.items() if key and key not in core and value and str(value).strip()}
    return Substance(
        substance_id=identity, canonical_name=name, cas=cas,
        smiles=str(row.get("smiles") or "").strip() or None,
        aliases=(abbreviation,) if abbreviation else (), roles=_role_assignments(row),
        properties=properties,
    )


__all__ = ["DEFINITIONS_DIR", "PENDING_PATH", "SUBSTANCES_PATH", "TAXONOMY_PATH", "iter_substance_rows", "load_taxonomy", "row_to_substance"]
