"""Read the package-owned substance and role/family definitions."""

from __future__ import annotations

import csv
import json
from dataclasses import replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Mapping, Tuple

from .models import RoleAssignment, Substance, SubstanceIdentifier
from .normalization import normalize_cas

DEFINITIONS_DIR = Path(__file__).with_name("definitions")
SUBSTANCES_PATH = DEFINITIONS_DIR / "substances.v1.csv"
ADDITIONS_PATH = DEFINITIONS_DIR / "substance_additions.v1.csv"
IDENTIFIERS_PATH = DEFINITIONS_DIR / "substance_identifiers.v1.csv"
PENDING_PATH = DEFINITIONS_DIR / "pending_substances.csv"
TAXONOMY_PATH = DEFINITIONS_DIR / "roles_families.v1.json"


@lru_cache(maxsize=1)
def load_taxonomy() -> Dict[str, Any]:
    with TAXONOMY_PATH.open("r", encoding="utf-8") as handle:
        return dict(json.load(handle))


def iter_substance_rows(
    path: str | Path = SUBSTANCES_PATH,
) -> Iterator[Dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle)


def iter_identifier_rows(
    path: str | Path = IDENTIFIERS_PATH,
) -> Iterator[Dict[str, str]]:
    """Yield supplemental one-to-many identifier definition rows."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle)


def iter_addition_rows(
    path: str | Path = ADDITIONS_PATH,
) -> Iterator[Dict[str, str]]:
    """Yield explicitly curated compound additions."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
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


def _legacy_identifiers(
    *,
    substance_id: str,
    canonical_name: str,
    cas: str | None,
    abbreviation: str,
    source_definition: str,
) -> Tuple[SubstanceIdentifier, ...]:
    identifiers = [
        SubstanceIdentifier(
            identifier_id=f"legacy:{substance_id}:canonical_name",
            substance_id=substance_id,
            identifier_type="canonical_name",
            value=canonical_name,
            is_preferred=True,
            source=f"{source_definition}:name",
            normalization_profile="chemical_name_v1",
        )
    ]
    if cas:
        identifiers.append(
            SubstanceIdentifier(
                identifier_id=f"legacy:{substance_id}:cas",
                substance_id=substance_id,
                identifier_type="cas",
                value=cas,
                is_preferred=True,
                source=f"{source_definition}:cas",
                normalization_profile="cas_v1",
            )
        )
    if abbreviation:
        identifiers.append(
            SubstanceIdentifier(
                identifier_id=f"legacy:{substance_id}:abbreviation",
                substance_id=substance_id,
                identifier_type="abbreviation",
                value=abbreviation,
                source=f"{source_definition}:abbreviation",
                normalization_profile="abbreviation_v1",
            )
        )
    return tuple(identifiers)


def row_to_identifier(row: Mapping[str, str]) -> SubstanceIdentifier:
    """Convert one declarative supplemental identifier row."""
    confidence_raw = str(row.get("confidence") or "1.0").strip()
    return SubstanceIdentifier(
        identifier_id=str(row.get("identifier_id") or "").strip(),
        substance_id=str(row.get("substance_id") or "").strip(),
        identifier_type=str(row.get("identifier_type") or "").strip(),
        value=str(row.get("value") or "").strip(),
        language=str(row.get("language") or "").strip() or None,
        is_preferred=str(row.get("is_preferred") or "").strip().casefold()
        in {"1", "true", "yes"},
        source=str(row.get("source") or "").strip() or None,
        confidence=float(confidence_raw),
        status=str(row.get("status") or "active").strip(),
        normalization_profile=(
            str(row.get("normalization_profile") or "").strip() or None
        ),
        allow_ambiguous=str(row.get("allow_ambiguous") or "")
        .strip()
        .casefold()
        in {"1", "true", "yes"},
    )


def row_to_substance(
    row: Dict[str, str],
    row_number: int,
    *,
    supplemental_identifiers: Iterable[SubstanceIdentifier] = (),
    source_definition: str = "substances.v1.csv",
) -> Substance:
    name = str(row.get("name") or "").strip()
    cas_raw = str(row.get("cas") or "").strip()
    cas = normalize_cas(cas_raw)
    abbreviation = str(row.get("abbreviation") or "").strip().strip('"')
    explicit_identity = str(row.get("substance_id") or "").strip()
    identity = explicit_identity or (
        f"cas:{cas}" if cas else f"registry-row:{row_number}"
    )
    core = {"name", "abbreviation", "cas", "smiles", "role_1", "family_1", "tag_1", "role_2", "family_2", "tag_2"}
    properties = {key: str(value).strip() for key, value in row.items() if key and key not in core and value and str(value).strip()}
    return Substance(
        substance_id=identity,
        canonical_name=name,
        cas=cas,
        smiles=str(row.get("smiles") or "").strip() or None,
        identifiers=(
            *_legacy_identifiers(
                substance_id=identity,
                canonical_name=name,
                cas=cas,
                abbreviation=abbreviation,
                source_definition=source_definition,
            ),
            *tuple(supplemental_identifiers),
        ),
        roles=_role_assignments(row),
        properties=properties,
    )


def load_substances(
    *,
    substances_path: str | Path = SUBSTANCES_PATH,
    additions_path: str | Path = ADDITIONS_PATH,
    identifiers_path: str | Path = IDENTIFIERS_PATH,
) -> Tuple[Substance, ...]:
    """Load substances and merge arbitrary supplemental identifiers."""
    legacy_substances = tuple(
        row_to_substance(row, row_number)
        for row_number, row in enumerate(
            iter_substance_rows(substances_path), start=2
        )
    )
    added_substances = tuple(
        row_to_substance(
            row,
            row_number,
            source_definition=Path(additions_path).name,
        )
        for row_number, row in enumerate(
            iter_addition_rows(additions_path), start=2
        )
        if str(row.get("status") or "active").strip() == "active"
    )
    substances = (*legacy_substances, *added_substances)
    by_id = {substance.substance_id: substance for substance in substances}
    additions: Dict[str, list[SubstanceIdentifier]] = {}
    for row in iter_identifier_rows(identifiers_path):
        identifier = row_to_identifier(row)
        if identifier.substance_id not in by_id:
            raise ValueError(
                "Supplemental identifier references unknown substance: "
                f"{identifier.identifier_id}:{identifier.substance_id}"
            )
        additions.setdefault(identifier.substance_id, []).append(identifier)
    return tuple(
        replace(
            substance,
            identifiers=(
                *substance.identifiers,
                *additions.get(substance.substance_id, ()),
            ),
        )
        for substance in substances
    )


__all__ = [
    "ADDITIONS_PATH",
    "DEFINITIONS_DIR",
    "IDENTIFIERS_PATH",
    "PENDING_PATH",
    "SUBSTANCES_PATH",
    "TAXONOMY_PATH",
    "iter_addition_rows",
    "iter_identifier_rows",
    "iter_substance_rows",
    "load_substances",
    "load_taxonomy",
    "row_to_identifier",
    "row_to_substance",
]
