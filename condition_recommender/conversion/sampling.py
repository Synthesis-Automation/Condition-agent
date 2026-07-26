"""Deterministic, reference-safe sampling for large reaction corpora."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple

from .input_schema import discover_csv_datasets, iter_csv_records
from .references import normalize_reference

SAMPLER_SCHEMA_VERSION = "1.0"
SAMPLER_VERSION = "generic_reference_safe_sampling.v1"
_PRIMARY_PARTITIONS = ("development", "validation", "untouched_test")
_METADATA_FIELDS = (
    "source_dataset",
    "source_path",
    "source_row_number",
    "sample_partition",
    "sample_is_smoke",
    "sampling_group_id",
    "reference_id",
    "reaction_text_id",
    "sampling_strata",
    "sampling_priority",
)


@dataclass
class _Candidate:
    candidate_id: str
    source_dataset: str
    source_path: str
    source_row_number: int
    reaction_id: str
    reference_id: str
    reaction_text_id: str
    raw_recipe_id: str
    stage_class: str
    condition_class: str
    yield_class: str
    procedure_class: str
    priority: str
    group_id: str = ""
    partition: str = ""
    reference_class: str = ""

    @property
    def strata(self) -> Tuple[str, ...]:
        return (
            self.reference_class,
            self.stage_class,
            self.condition_class,
            self.yield_class,
            self.procedure_class,
        )


def _sha256_text(prefix: str, value: str) -> str:
    return prefix + hashlib.sha256(value.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _reaction_text_id(reaction_smiles: str) -> str:
    normalized = "".join(str(reaction_smiles or "").split())
    return _sha256_text("RXNTEXT1:", normalized) if normalized else ""


def _raw_recipe_id(
    catalyst_cas: Iterable[str],
    reagent_cas: Iterable[str],
    solvent_cas: Iterable[str],
) -> str:
    payload = json.dumps(
        {
            "catalyst_cas": tuple(sorted(catalyst_cas)),
            "reagent_cas": tuple(sorted(reagent_cas)),
            "solvent_cas": tuple(sorted(solvent_cas)),
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return _sha256_text("SAMPLERCOND1:", payload)


def _stage_class(stages: str) -> str:
    value = str(stages or "").strip()
    if not value:
        return "stage:unknown"
    try:
        count = int(value)
    except ValueError:
        return "stage:invalid"
    return "stage:multi" if count > 1 else "stage:single"


def _condition_class(
    catalyst_cas: Tuple[str, ...],
    reagent_cas: Tuple[str, ...],
    solvent_cas: Tuple[str, ...],
) -> str:
    present = sum(bool(values) for values in (catalyst_cas, reagent_cas, solvent_cas))
    if present == 3:
        return "conditions:all_source_fields"
    if present:
        return "conditions:partial"
    return "conditions:missing"


def _yield_class(value: Optional[float]) -> str:
    if value is None:
        return "yield:missing"
    if not 0.0 <= value <= 100.0:
        return "yield:invalid"
    if value < 20.0:
        return "yield:0_20"
    if value < 50.0:
        return "yield:20_50"
    if value < 80.0:
        return "yield:50_80"
    return "yield:80_100"


def _priority(seed: int, *tokens: object) -> str:
    identity = "\0".join((str(seed), *(str(token) for token in tokens)))
    return hashlib.sha256(identity.encode("utf-8")).hexdigest()


def _scan_candidates(
    paths: Tuple[Path, ...],
    *,
    seed: int,
) -> tuple[
    list[_Candidate],
    Dict[str, int],
    Dict[str, set[str]],
    Dict[str, set[str]],
]:
    candidates = []
    reference_counts: Counter[str] = Counter()
    reference_recipes: Dict[str, set[str]] = defaultdict(set)
    reference_datasets: Dict[str, set[str]] = defaultdict(set)
    for path in paths:
        for record in iter_csv_records(path):
            reference = normalize_reference(record.reference)
            reaction_text_id = _reaction_text_id(record.reaction_smiles)
            recipe_id = _raw_recipe_id(
                record.catalyst_cas,
                record.reagent_cas,
                record.solvent_cas,
            )
            candidate_id = (
                f"{record.source_path}\0{record.source_row_number}"
            )
            candidates.append(
                _Candidate(
                    candidate_id=candidate_id,
                    source_dataset=record.source_dataset,
                    source_path=record.source_path,
                    source_row_number=record.source_row_number,
                    reaction_id=record.reaction_id,
                    reference_id=reference.reference_id,
                    reaction_text_id=reaction_text_id,
                    raw_recipe_id=recipe_id,
                    stage_class=_stage_class(record.stages),
                    condition_class=_condition_class(
                        record.catalyst_cas,
                        record.reagent_cas,
                        record.solvent_cas,
                    ),
                    yield_class=_yield_class(record.yield_pct),
                    procedure_class=(
                        "procedure:present"
                        if record.experimental_procedure
                        else "procedure:missing"
                    ),
                    priority=_priority(
                        seed,
                        record.source_dataset,
                        record.reaction_id,
                        record.reaction_smiles,
                        record.source_row_number,
                    ),
                )
            )
            if reference.reference_id:
                reference_counts[reference.reference_id] += 1
                reference_recipes[reference.reference_id].add(recipe_id)
                reference_datasets[reference.reference_id].add(
                    record.source_dataset
                )
    return (
        candidates,
        dict(reference_counts),
        reference_recipes,
        reference_datasets,
    )


def _assign_connected_groups(
    candidates: list[_Candidate],
    *,
    reference_counts: Mapping[str, int],
    reference_recipes: Mapping[str, set[str]],
    reference_datasets: Mapping[str, set[str]],
) -> None:
    parents = list(range(len(candidates)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    owners: Dict[Tuple[str, str], int] = {}
    tokens_by_candidate: list[Tuple[str, ...]] = []
    for index, candidate in enumerate(candidates):
        tokens = []
        if candidate.reference_id:
            tokens.append(f"reference:{candidate.reference_id}")
        if candidate.reaction_text_id:
            tokens.append(f"reaction_text:{candidate.reaction_text_id}")
        if not tokens:
            tokens.append(f"observation:{candidate.candidate_id}")
        tokens_by_candidate.append(tuple(tokens))
        for token in tokens:
            key = tuple(token.split(":", maxsplit=1))
            previous = owners.setdefault(key, index)
            union(index, previous)

    tokens_by_root: Dict[int, set[str]] = defaultdict(set)
    for index, tokens in enumerate(tokens_by_candidate):
        tokens_by_root[find(index)].update(tokens)
    group_ids = {
        root: _sha256_text("SG1:", "\0".join(sorted(tokens)))
        for root, tokens in tokens_by_root.items()
    }
    for index, candidate in enumerate(candidates):
        candidate.group_id = group_ids[find(index)]
        if not candidate.reference_id:
            candidate.reference_class = "reference:missing"
            continue
        count = int(reference_counts.get(candidate.reference_id, 0))
        recipes = len(reference_recipes.get(candidate.reference_id, set()))
        datasets = len(reference_datasets.get(candidate.reference_id, set()))
        if datasets > 1:
            candidate.reference_class = "reference:cross_file"
        elif count <= 1:
            candidate.reference_class = "reference:singleton"
        elif recipes <= 1:
            candidate.reference_class = "reference:repeated_single_recipe"
        else:
            candidate.reference_class = "reference:repeated_multi_recipe"


def _assign_partitions(
    candidates: list[_Candidate],
    *,
    seed: int,
    development_size: int,
    validation_size: int,
    test_size: int,
) -> None:
    total = development_size + validation_size + test_size
    if total < 1:
        raise ValueError("At least one primary sample target must be positive")
    development_threshold = development_size / total
    validation_threshold = (development_size + validation_size) / total
    group_partitions = {}
    for group_id in sorted({candidate.group_id for candidate in candidates}):
        value = int(_priority(seed, "partition", group_id), 16) / (2**256)
        partition = (
            "development"
            if value < development_threshold
            else "validation"
            if value < validation_threshold
            else "untouched_test"
        )
        group_partitions[group_id] = partition
    for candidate in candidates:
        candidate.partition = group_partitions[candidate.group_id]


def _allocate_source_quotas(
    available: Mapping[str, int],
    target: int,
) -> Dict[str, int]:
    sources = tuple(sorted(source for source, count in available.items() if count))
    if target <= 0 or not sources:
        return {source: 0 for source in sources}
    target = min(target, sum(available.values()))
    minimum = 1 if target >= len(sources) else 0
    quotas = {
        source: min(minimum, int(available[source])) for source in sources
    }
    remaining = target - sum(quotas.values())
    capacities = {
        source: int(available[source]) - quotas[source] for source in sources
    }
    while remaining:
        active = tuple(source for source in sources if capacities[source] > 0)
        if not active:
            break
        weights = {source: math.sqrt(capacities[source]) for source in active}
        total_weight = sum(weights.values())
        shares = {
            source: remaining * weights[source] / total_weight for source in active
        }
        additions = {
            source: min(capacities[source], int(shares[source]))
            for source in active
        }
        added = sum(additions.values())
        if added == 0:
            source = sorted(
                active,
                key=lambda item: (
                    -(shares[item] - int(shares[item])),
                    item,
                ),
            )[0]
            additions[source] = 1
            added = 1
        for source, count in additions.items():
            quotas[source] += count
            capacities[source] -= count
        remaining -= added
    return quotas


def _round_robin_source(
    candidates: Iterable[_Candidate],
    *,
    limit: int,
    group_counts: Counter[str],
    max_rows_per_group: int,
) -> list[_Candidate]:
    buckets: Dict[Tuple[str, ...], deque[_Candidate]] = {}
    grouped: Dict[Tuple[str, ...], list[_Candidate]] = defaultdict(list)
    for candidate in candidates:
        grouped[candidate.strata].append(candidate)
    for stratum, values in grouped.items():
        buckets[stratum] = deque(sorted(values, key=lambda item: item.priority))
    selected = []
    while len(selected) < limit and buckets:
        progressed = False
        for stratum in tuple(sorted(buckets)):
            values = buckets[stratum]
            while values:
                candidate = values.popleft()
                if group_counts[candidate.group_id] >= max_rows_per_group:
                    continue
                selected.append(candidate)
                group_counts[candidate.group_id] += 1
                progressed = True
                break
            if not values:
                del buckets[stratum]
            if len(selected) >= limit:
                break
        if not progressed:
            break
    return selected


def _balanced_select(
    candidates: Iterable[_Candidate],
    *,
    target: int,
    max_rows_per_group: int,
) -> Tuple[_Candidate, ...]:
    values = tuple(candidates)
    if target <= 0 or not values:
        return ()
    by_source: Dict[str, list[_Candidate]] = defaultdict(list)
    for candidate in values:
        by_source[candidate.source_dataset].append(candidate)
    quotas = _allocate_source_quotas(
        {source: len(rows) for source, rows in by_source.items()},
        target,
    )
    selected = []
    selected_ids = set()
    group_counts: Counter[str] = Counter()
    for source in sorted(by_source):
        rows = _round_robin_source(
            by_source[source],
            limit=quotas.get(source, 0),
            group_counts=group_counts,
            max_rows_per_group=max_rows_per_group,
        )
        selected.extend(rows)
        selected_ids.update(row.candidate_id for row in rows)
    if len(selected) < min(target, len(values)):
        remaining = sorted(
            (
                candidate
                for candidate in values
                if candidate.candidate_id not in selected_ids
            ),
            key=lambda item: item.priority,
        )
        for enforce_cap in (True, False):
            for candidate in remaining:
                if candidate.candidate_id in selected_ids:
                    continue
                if (
                    enforce_cap
                    and group_counts[candidate.group_id] >= max_rows_per_group
                ):
                    continue
                selected.append(candidate)
                selected_ids.add(candidate.candidate_id)
                group_counts[candidate.group_id] += 1
                if len(selected) >= min(target, len(values)):
                    break
            if len(selected) >= min(target, len(values)):
                break
    return tuple(sorted(selected, key=lambda item: item.candidate_id))


def _source_headers(paths: Tuple[Path, ...]) -> Tuple[str, ...]:
    headers = []
    seen = set(_METADATA_FIELDS)
    for path in paths:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            for field in next(reader, ()):
                name = str(field)
                if name not in seen:
                    seen.add(name)
                    headers.append(name)
    return _METADATA_FIELDS + tuple(headers)


def _write_sample_csvs(
    paths: Tuple[Path, ...],
    destination: Path,
    *,
    selected: Mapping[str, Tuple[_Candidate, ...]],
) -> None:
    fieldnames = _source_headers(paths)
    selected_lookup: Dict[Tuple[str, int], _Candidate] = {}
    smoke_ids = {candidate.candidate_id for candidate in selected["smoke"]}
    for partition in _PRIMARY_PARTITIONS:
        for candidate in selected[partition]:
            selected_lookup[(candidate.source_path, candidate.source_row_number)] = (
                candidate
            )
    handles = {
        name: (destination / f"{name}.csv").open(
            "w", encoding="utf-8-sig", newline=""
        )
        for name in ("smoke", *_PRIMARY_PARTITIONS)
    }
    try:
        writers = {
            name: csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
            for name, handle in handles.items()
        }
        for writer in writers.values():
            writer.writeheader()
        for path in paths:
            with path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                for row_number, row in enumerate(reader, start=2):
                    candidate = selected_lookup.get((str(path), row_number))
                    if candidate is None:
                        continue
                    output = dict(row)
                    output.update(
                        {
                            "source_dataset": candidate.source_dataset,
                            "source_path": candidate.source_path,
                            "source_row_number": candidate.source_row_number,
                            "sample_partition": candidate.partition,
                            "sample_is_smoke": str(
                                candidate.candidate_id in smoke_ids
                            ).lower(),
                            "sampling_group_id": candidate.group_id,
                            "reference_id": candidate.reference_id,
                            "reaction_text_id": candidate.reaction_text_id,
                            "sampling_strata": json.dumps(candidate.strata),
                            "sampling_priority": candidate.priority,
                        }
                    )
                    writers[candidate.partition].writerow(output)
                    if candidate.candidate_id in smoke_ids:
                        writers["smoke"].writerow(output)
    finally:
        for handle in handles.values():
            handle.close()


def _partition_leakage(
    selected: Mapping[str, Tuple[_Candidate, ...]],
    field: str,
) -> int:
    owners: Dict[str, set[str]] = defaultdict(set)
    for partition in _PRIMARY_PARTITIONS:
        for candidate in selected[partition]:
            value = str(getattr(candidate, field))
            if value:
                owners[value].add(partition)
    return sum(len(partitions) > 1 for partitions in owners.values())


def build_reference_safe_samples(
    dataset_path: str | Path,
    output_dir: str | Path,
    *,
    smoke_size: int = 500,
    development_size: int = 5000,
    validation_size: int = 2000,
    test_size: int = 2000,
    seed: int = 17,
    max_rows_per_group: int = 20,
) -> Dict[str, Any]:
    """Build deterministic source-balanced samples without reference leakage."""
    sizes = {
        "smoke": smoke_size,
        "development": development_size,
        "validation": validation_size,
        "untouched_test": test_size,
    }
    if any(value < 0 for value in sizes.values()):
        raise ValueError("Sample sizes cannot be negative")
    if smoke_size > development_size:
        raise ValueError("Smoke sample must be a subset of development")
    if max_rows_per_group < 1:
        raise ValueError("max_rows_per_group must be positive")
    paths = discover_csv_datasets(dataset_path)
    if not paths:
        raise ValueError(f"No CSV datasets found at {dataset_path}")

    (
        candidates,
        reference_counts,
        reference_recipes,
        reference_datasets,
    ) = _scan_candidates(paths, seed=seed)
    _assign_connected_groups(
        candidates,
        reference_counts=reference_counts,
        reference_recipes=reference_recipes,
        reference_datasets=reference_datasets,
    )
    _assign_partitions(
        candidates,
        seed=seed,
        development_size=development_size,
        validation_size=validation_size,
        test_size=test_size,
    )
    selected = {
        partition: _balanced_select(
            (
                candidate
                for candidate in candidates
                if candidate.partition == partition
            ),
            target=sizes[partition],
            max_rows_per_group=max_rows_per_group,
        )
        for partition in _PRIMARY_PARTITIONS
    }
    selected["smoke"] = _balanced_select(
        selected["development"],
        target=smoke_size,
        max_rows_per_group=max_rows_per_group,
    )

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    _write_sample_csvs(paths, destination, selected=selected)

    smoke_ids = {candidate.candidate_id for candidate in selected["smoke"]}
    entries = []
    for partition in _PRIMARY_PARTITIONS:
        for candidate in selected[partition]:
            entries.append(
                {
                    "source_dataset": candidate.source_dataset,
                    "source_path": candidate.source_path,
                    "source_row_number": candidate.source_row_number,
                    "reaction_id": candidate.reaction_id,
                    "reference_id": candidate.reference_id,
                    "reaction_text_id": candidate.reaction_text_id,
                    "sampling_group_id": candidate.group_id,
                    "partition": partition,
                    "smoke": candidate.candidate_id in smoke_ids,
                    "strata": candidate.strata,
                    "priority": candidate.priority,
                }
            )
    entries.sort(
        key=lambda item: (
            item["partition"],
            item["source_dataset"],
            item["source_row_number"],
        )
    )
    actual_sizes = {
        name: len(selected[name]) for name in ("smoke", *_PRIMARY_PARTITIONS)
    }
    source_counts = {
        partition: dict(
            sorted(
                Counter(
                    candidate.source_dataset
                    for candidate in selected[partition]
                ).items()
            )
        )
        for partition in ("smoke", *_PRIMARY_PARTITIONS)
    }
    manifest = {
        "schema_version": SAMPLER_SCHEMA_VERSION,
        "sampler_version": SAMPLER_VERSION,
        "dataset_path": str(Path(dataset_path)),
        "seed": seed,
        "targets": sizes,
        "actual_sizes": actual_sizes,
        "max_rows_per_group": max_rows_per_group,
        "input_files": [
            {
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": _sha256_file(path),
            }
            for path in paths
        ],
        "source_counts": source_counts,
        "reference_leakage_count": _partition_leakage(
            selected, "reference_id"
        ),
        "reaction_text_leakage_count": _partition_leakage(
            selected, "reaction_text_id"
        ),
        "entries": entries,
    }
    report = {
        key: value for key, value in manifest.items() if key != "entries"
    }
    (destination / "sample_manifest.v1.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (destination / "sampling_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


__all__ = [
    "SAMPLER_SCHEMA_VERSION",
    "SAMPLER_VERSION",
    "build_reference_safe_samples",
]
