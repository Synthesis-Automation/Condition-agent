"""Dataset-driven condition identity and role coverage auditing."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

from .api import resolve_substance

_SPLIT_RE = re.compile(r"[,;|]")
_GENERIC_ROLES = {"organic_compound", "inorganic_compound", "other_reagent"}
_EXPECTED_ROLES = {
    "catalyst_cas": {"metal_catalyst", "ligand", "organo_catalyst", "enzyme", "additive"},
    "solvent_cas": {"solvent"},
    # Reagent columns are deliberately broad; reaction context resolves these later.
    "reagent_cas": {"base", "acid", "condensation_agent", "oxidant", "reductant", "additive", "other_reagent", "organic_compound", "inorganic_compound", "ligand", "metal_catalyst", "solvent"},
}

DEFAULT_DATASETS = {
    "suzuki_miyaura": "examples/reaction_dataset_300/suzuki_miyaura.csv",
    "c_n_coupling": "examples/reaction_dataset_300/C_N_Coupling.csv",
    "c_o_coupling": "examples/reaction_dataset_300/C_O_Coupling.csv",
    "c_s_coupling": "examples/reaction_dataset_300/C_S_Coupling.csv",
    "amide_formation": "examples/reaction_dataset_300/Amide_formation.csv",
}


@dataclass
class _ObservedIdentity:
    identifier: str
    source_field: str
    occurrences: int = 0
    datasets: set[str] = field(default_factory=set)
    dataset_occurrences: Counter[str] = field(default_factory=Counter)
    reaction_ids: set[str] = field(default_factory=set)


def _identifiers(value: Any) -> Tuple[str, ...]:
    return tuple(sorted({item.strip() for item in _SPLIT_RE.split(str(value or "")) if item.strip()}))


def _classify(identifier: str, source_field: str) -> Dict[str, Any]:
    result = resolve_substance(cas=identifier)
    base: Dict[str, Any] = {
        "identifier": identifier,
        "source_field": source_field,
        "resolution_status": result.status,
        "match_kind": result.match_kind or "",
        "substance_id": "",
        "canonical_name": "",
        "role_ids": "",
        "family_ids": "",
        "coverage_category": "",
        "coverage_flags": "",
    }
    if result.status == "invalid_identifier":
        base["coverage_category"] = "invalid_identifier"
        return base
    if result.status != "resolved" or result.substance is None:
        base["coverage_category"] = "missing_substance"
        return base
    substance = result.substance
    roles = {role.role_id for role in substance.roles if role.role_id}
    families = {role.family_id for role in substance.roles if role.family_id}
    base.update({
        "substance_id": substance.substance_id,
        "canonical_name": substance.canonical_name,
        "role_ids": "|".join(sorted(roles)),
        "family_ids": "|".join(sorted(families)),
    })
    specific_roles = roles - _GENERIC_ROLES
    expected = _EXPECTED_ROLES[source_field]
    flags = []
    if not specific_roles:
        category = "generic_only"
    elif not roles.intersection(expected):
        category = "role_conflict"
    elif len(roles) > 1:
        category = "multiple_roles"
    else:
        category = "resolved"
    if len(roles) > 1:
        flags.append("MULTIPLE_POSSIBLE_ROLES")
    if roles.intersection(_GENERIC_ROLES) and specific_roles:
        flags.append("GENERIC_AND_SPECIFIC_ROLES")
    base["coverage_category"] = category
    base["coverage_flags"] = "|".join(flags)
    return base


def analyze_coverage(
    datasets: Mapping[str, str | Path],
    output_dir: str | Path,
) -> Dict[str, Any]:
    """Resolve condition identifiers and write mutually exclusive review queues."""
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    observed: Dict[Tuple[str, str], _ObservedIdentity] = {}
    source_rows = 0
    for family_id, raw_path in datasets.items():
        path = Path(raw_path)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        source_rows += len(rows)
        for row in rows:
            reaction_id = str(row.get("reaction_id") or "")
            for source_field in _EXPECTED_ROLES:
                for identifier in _identifiers(row.get(source_field)):
                    key = (identifier, source_field)
                    entry = observed.setdefault(key, _ObservedIdentity(identifier, source_field))
                    entry.occurrences += 1; entry.datasets.add(family_id); entry.dataset_occurrences[family_id] += 1
                    if reaction_id: entry.reaction_ids.add(reaction_id)
    records: List[Dict[str, Any]] = []
    for entry in observed.values():
        record = _classify(entry.identifier, entry.source_field)
        record.update({
            "occurrence_count": entry.occurrences,
            "dataset_count": len(entry.datasets),
            "datasets": "|".join(sorted(entry.datasets)),
            "dataset_occurrences_json": json.dumps(dict(sorted(entry.dataset_occurrences.items())), sort_keys=True),
            "reaction_count": len(entry.reaction_ids),
            "reaction_ids": "|".join(sorted(entry.reaction_ids)),
        })
        records.append(record)
    records.sort(key=lambda row: (row["coverage_category"], -int(row["occurrence_count"]), row["identifier"], row["source_field"]))
    fieldnames = list(records[0]) if records else []
    files = {
        "resolved": "resolved.csv",
        "multiple_roles": "ambiguous_roles.csv",
        "generic_only": "generic_only.csv",
        "role_conflict": "role_conflicts.csv",
        "missing_substance": "missing_substances.csv",
        "invalid_identifier": "invalid_identifiers.csv",
    }
    by_category = {category: [row for row in records if row["coverage_category"] == category] for category in files}
    for category, filename in files.items():
        with (output_dir / filename).open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames); writer.writeheader(); writer.writerows(by_category[category])
    occurrence_counts = Counter()
    dataset_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    source_field_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    for row in records:
        occurrence_counts[row["coverage_category"]] += int(row["occurrence_count"])
        source_field_counts[row["source_field"]][row["coverage_category"]] += int(row["occurrence_count"])
        for dataset, count in json.loads(row["dataset_occurrences_json"]).items():
            dataset_counts[dataset][row["coverage_category"]] += int(count)
    total_occurrences = sum(occurrence_counts.values())
    identity_resolved = total_occurrences - occurrence_counts["missing_substance"] - occurrence_counts["invalid_identifier"]
    report = {
        "schema_version": "1.0",
        "datasets": {family: str(path) for family, path in datasets.items()},
        "source_rows": source_rows,
        "unique_identifier_source_pairs": len(records),
        "unique_identifiers": len({row["identifier"] for row in records}),
        "category_counts": {category: len(by_category[category]) for category in files},
        "occurrence_counts": dict(sorted(occurrence_counts.items())),
        "by_dataset_occurrence_counts": {
            dataset: dict(sorted(counts.items())) for dataset, counts in sorted(dataset_counts.items())
        },
        "by_source_field_occurrence_counts": {
            field: dict(sorted(counts.items())) for field, counts in sorted(source_field_counts.items())
        },
        "identity_resolution_rate": round(identity_resolved / total_occurrences, 6) if total_occurrences else 0.0,
        "recommendation_ready_role_rate": round(occurrence_counts["resolved"] / total_occurrences, 6) if total_occurrences else 0.0,
    }
    (output_dir / "coverage_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


__all__ = ["DEFAULT_DATASETS", "analyze_coverage"]
