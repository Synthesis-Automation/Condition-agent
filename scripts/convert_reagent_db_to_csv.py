"""
Convert reagent registry JSON files into a single curated CSV.

Usage:
  python scripts/convert_reagent_db_to_csv.py
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


BASE_DIR = Path("data/reagent_db")
OUTPUT_CSV = BASE_DIR / "reagents.csv"
AUDIT_JSON = BASE_DIR / "reagents_audit.json"
EXCLUDE_FILES = {"not_determined_reagents.json", "reagents_audit.json"}
EXCLUDE_PREFIXES = ("taxonomy_",)

CSV_FIELDS = [
    "name",
    "abbreviation",
    "cas",
    "smile",
    "role",
    "family_id",
    "tag",
]

_SLUG_RE = re.compile(r"[^A-Za-z0-9]+")


def _first_value(value: Any) -> str:
    if isinstance(value, (list, tuple, set)):
        for item in value:
            if item or item == 0:
                return str(item)
        return ""
    if value is None:
        return ""
    return str(value)


def _slugify(text: Optional[str]) -> str:
    slug = _SLUG_RE.sub("-", (text or "").strip()).strip("-")
    return slug.lower() or "unknown"


def _pick_role(entry: Dict[str, Any], fallback_role: str) -> str:
    roles = entry.get("roles")
    if isinstance(roles, dict) and roles:
        if fallback_role in roles:
            return fallback_role
        if len(roles) == 1:
            return next(iter(roles))
        return sorted(roles.keys())[0]
    role = entry.get("role")
    return str(role) if role else fallback_role


def _family_id(roles: Dict[str, Any], role: str) -> str:
    payload = roles.get(role) if isinstance(roles, dict) else None
    if not isinstance(payload, dict):
        return ""
    families = payload.get("families")
    if isinstance(families, (list, tuple)) and families:
        return str(families[0])
    if isinstance(families, str):
        return families
    return ""


def _build_tag(role: str, family_id: str, payload: Dict[str, Any]) -> str:
    role_label = role.replace("_", " ")
    family_lower = family_id.lower()

    def pick(value: Any) -> str:
        return _first_value(value)

    if role == "acid":
        acidity = pick(payload.get("acidity"))
        if "lewis" in family_lower:
            return f"{acidity} Lewis acid".strip() if acidity else "Lewis acid"
        if acidity:
            return f"{acidity} acid"
        return "acid"

    if role == "base":
        basicity = pick(payload.get("basicity") or payload.get("strength_band"))
        sterics = pick(payload.get("sterics"))
        nucleophilicity = pick(payload.get("nucleophilicity"))
        base_tag = f"{basicity} base" if basicity else "base"
        qualifiers: List[str] = []
        if sterics in {"hindered", "bulky"}:
            qualifiers.append(sterics)
        if nucleophilicity == "low":
            qualifiers.append("non-nucleophilic")
        if qualifiers:
            return f"{base_tag}; {', '.join(qualifiers)}"
        return base_tag

    if role == "solvent":
        polarity = pick(payload.get("polarity_band") or payload.get("polarity"))
        proticity = pick(payload.get("proticity"))
        tokens = [tok for tok in (polarity, proticity) if tok]
        return f"{' '.join(tokens)} solvent" if tokens else "solvent"

    if role in {"oxidant", "reductant", "condensation_agent"}:
        strength = pick(payload.get("strength_band") or payload.get("strength"))
        return f"{strength} {role_label}".strip() if strength else role_label

    if role == "ligand":
        denticity = pick(payload.get("denticity"))
        donors = pick(payload.get("donors"))
        tokens = [tok for tok in (denticity, donors) if tok]
        return f"{' '.join(tokens)} ligand" if tokens else "ligand"

    if role == "metal_catalyst":
        metal = pick(payload.get("metal"))
        oxidation = pick(payload.get("oxidation_states"))
        if metal and oxidation:
            return f"{metal}({oxidation}) catalyst"
        if metal:
            return f"{metal} catalyst"
        return "metal catalyst"

    if role in {"additive", "other_reagent", "organo_catalyst", "enzyme"}:
        return role_label

    return role_label


def _ensure_unique_cas(raw: str, used: set[str]) -> str:
    if raw not in used:
        used.add(raw)
        return raw
    suffix = 2
    while True:
        candidate = f"{raw}-{suffix}"
        if candidate not in used:
            used.add(candidate)
            return candidate
        suffix += 1


def _artificial_cas(entry: Dict[str, Any], role: str) -> str:
    entry_id = _first_value(entry.get("id")).strip()
    if entry_id:
        return f"${entry_id}"
    name = _first_value(entry.get("name"))
    return f"${_slugify(name)}-{role}"


def _row_from_entry(
    entry: Dict[str, Any],
    fallback_role: str,
    used_artificial_cas: set[str],
    audit: Dict[str, Any],
) -> Dict[str, str]:
    roles = entry.get("roles") or {}
    role = _pick_role(entry, fallback_role)
    payload = roles.get(role, {}) if isinstance(roles, dict) else {}

    entry_id = _first_value(entry.get("id"))
    name = _first_value(entry.get("name"))
    smiles = _first_value(entry.get("smiles"))
    family_id = _family_id(roles, role)

    abbr = _first_value(entry.get("abbreviation"))
    if not abbr:
        abbr = _first_value(entry.get("abbr"))

    cas = _first_value(entry.get("cas"))
    if not cas:
        audit["missing_cas"] += 1
        cas = _artificial_cas(entry, role)
        cas = _ensure_unique_cas(cas, used_artificial_cas)
        audit["generated_cas"] += 1

    if not entry_id:
        audit["missing_id"] += 1
        entry_id = f"cas-{cas}" if cas else f"name-{_slugify(name)}"

    if not family_id:
        audit["missing_family_id"] += 1

    tag = _build_tag(role, family_id, payload)
    if (not tag or tag == role.replace("_", " ")) and family_id:
        tag = family_id.replace("_", " ")

    return {
        "name": name,
        "abbreviation": abbr,
        "cas": cas,
        "smile": smiles,
        "role": role,
        "family_id": family_id,
        "tag": tag,
    }


def _load_entries(path: Path) -> Iterable[Dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError(f"{path} must contain a JSON list.")
    return payload


def main() -> None:
    base = BASE_DIR
    files = [p for p in base.glob("*.json") if p.is_file()]
    files = [
        p for p in files
        if p.name not in EXCLUDE_FILES and not p.name.startswith(EXCLUDE_PREFIXES)
    ]

    used_artificial_cas: set[str] = set()
    rows: List[Dict[str, str]] = []

    audit: Dict[str, Any] = {
        "source_files": [p.name for p in sorted(files)],
        "excluded_files": sorted(EXCLUDE_FILES),
        "total_entries": 0,
        "written_entries": 0,
        "missing_cas": 0,
        "generated_cas": 0,
        "missing_id": 0,
        "missing_family_id": 0,
        "role_counts": Counter(),
        "duplicate_cas": {},
        "samples": defaultdict(list),
    }

    for path in sorted(files):
        fallback_role = path.stem
        entries = list(_load_entries(path))
        audit["total_entries"] += len(entries)
        for entry in entries:
            row = _row_from_entry(entry, fallback_role, used_artificial_cas, audit)
            audit["role_counts"][row["role"]] += 1
            rows.append(row)
            if not row["family_id"] and len(audit["samples"]["missing_family_id"]) < 5:
                audit["samples"]["missing_family_id"].append(row["name"])

    rows.sort(key=lambda item: (item["role"], item["name"].lower(), item["cas"]))

    cas_counts = Counter(row["cas"] for row in rows if row["cas"])
    duplicates = {cas: count for cas, count in cas_counts.items() if count > 1}
    audit["duplicate_cas"] = dict(sorted(duplicates.items(), key=lambda item: item[1], reverse=True))
    audit["written_entries"] = len(rows)
    audit["role_counts"] = dict(audit["role_counts"])
    audit["samples"] = dict(audit["samples"])

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    AUDIT_JSON.write_text(
        json.dumps(audit, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
