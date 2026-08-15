"""Deterministic corpus conversion for observed route-action replay."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import random
import tempfile
import time
from typing import Any, Iterable, Iterator, Mapping, Optional, Sequence

from .generic_library import load_generic_library
from .generic_models import GenericTemplateLibrary
from .hierarchical_ranking import build_completion_prior_index
from .route_action_evaluation import (
    ROUTE_ACTION_EVALUATION_ALGORITHM_VERSION,
    ROUTE_ACTION_EVALUATION_SCHEMA_VERSION,
    RouteActionEvaluation,
    RouteActionEvaluationConfig,
    evaluate_route_actions,
)
from .route_conversion import iter_route_trees
from .observed_route_action import load_observed_route_action_label_policy


ROUTE_ACTION_CONVERTER_VERSION = "3.1"
DEFAULT_ROUTE_ACTION_SAMPLE_SEED = 20_260_814
_WORKER_LIBRARY: Any = None
_WORKER_PRIOR_INDEX: Any = None
_WORKER_CONFIG: Optional[RouteActionEvaluationConfig] = None


def _sha256_file(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def stable_route_shard(route_id: str, shard_count: int) -> int:
    """Return the stable content-derived shard for a route identity."""

    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    if not route_id:
        raise ValueError("route_id cannot be empty")
    token = hashlib.sha256(route_id.encode("utf-8")).digest()
    return int.from_bytes(token[:8], byteorder="big") % shard_count


def _atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        os.replace(temporary, path)
    except BaseException:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise


def _initialize_worker(library_path: str, config_values: dict[str, Any]) -> None:
    global _WORKER_LIBRARY, _WORKER_PRIOR_INDEX, _WORKER_CONFIG
    _WORKER_CONFIG = RouteActionEvaluationConfig(**config_values)
    _WORKER_LIBRARY = (
        load_generic_library(library_path)
        if _WORKER_CONFIG.run_search
        else GenericTemplateLibrary(
            templates=(),
            source_row_count=0,
            accepted_observation_count=0,
            rejection_counts={},
            definition={},
        )
    )
    _WORKER_PRIOR_INDEX = (
        build_completion_prior_index(_WORKER_LIBRARY)
        if _WORKER_CONFIG.run_search and _WORKER_CONFIG.use_hierarchical_ranking
        else None
    )


def _evaluation_task(
    tree: Any,
) -> tuple[str, Optional[RouteActionEvaluation], str, str]:
    route_id = tree.source_route_id or tree.tree_id
    try:
        if _WORKER_LIBRARY is None or _WORKER_CONFIG is None:
            raise RuntimeError("route-action worker is not initialized")
        evaluation = evaluate_route_actions(
            tree,
            _WORKER_LIBRARY,
            config=_WORKER_CONFIG,
            completion_prior_index=_WORKER_PRIOR_INDEX,
        )
        return route_id, evaluation, "", ""
    except (RuntimeError, TypeError, ValueError) as exc:
        return route_id, None, type(exc).__name__, str(exc)


def _selected_trees(
    source: Path,
    *,
    sample_size: Optional[int],
    seed: int,
    shard_count: int,
    shard_index: int,
) -> Iterator[Any]:
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    if shard_index < 0 or shard_index >= shard_count:
        raise ValueError("shard_index must be in [0, shard_count)")
    if sample_size is None:
        trees: Iterable[Any] = iter_route_trees(source)
    else:
        ordered = tuple(
            sorted(
                iter_route_trees(source),
                key=lambda tree: tree.source_route_id or tree.tree_id,
            )
        )
        if sample_size < 1 or sample_size > len(ordered):
            raise ValueError(f"sample_size must be between 1 and {len(ordered)}")
        trees = random.Random(seed).sample(ordered, sample_size)
    for tree in trees:
        route_id = tree.source_route_id or tree.tree_id
        if stable_route_shard(route_id, shard_count) == shard_index:
            yield tree


def _evaluation_results(
    source: Path,
    library_path: Path,
    *,
    config: RouteActionEvaluationConfig,
    sample_size: Optional[int],
    seed: int,
    workers: int,
    shard_count: int,
    shard_index: int,
) -> Iterator[tuple[str, Optional[RouteActionEvaluation], str, str]]:
    trees = _selected_trees(
        source,
        sample_size=sample_size,
        seed=seed,
        shard_count=shard_count,
        shard_index=shard_index,
    )
    config_values = asdict(config)
    if workers == 1:
        _initialize_worker(str(library_path), config_values)
        for tree in trees:
            yield _evaluation_task(tree)
        return
    materialized = tuple(trees)
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_initialize_worker,
        initargs=(str(library_path), config_values),
    ) as executor:
        yield from executor.map(_evaluation_task, materialized, chunksize=2)


def _new_metrics() -> dict[str, Any]:
    return {
        "route_count": 0,
        "step_count": 0,
        "search_eligible_step_count": 0,
        "searched_step_count": 0,
        "steps_with_candidates": 0,
        "candidate_count": 0,
        "source_patent_overlap_count": 0,
        "eligibility_facets": Counter(),
        "core_quality": Counter(),
        "limitations": Counter(),
        "operator_admission": Counter(),
        "departing_edit_descriptors": Counter(),
        "named_annotations": Counter(),
        "promoted": Counter(),
        "outcomes": Counter(),
        "candidate_supervision": Counter(),
        "diagnostics": Counter(),
        "rank_denominators": Counter(),
        "ranks": {
            name: Counter()
            for name in ("exact_precursor", "site", "operator", "synthon", "strategy")
        },
    }


def _accumulate(metrics: dict[str, Any], evaluation: RouteActionEvaluation) -> None:
    metrics["route_count"] += 1
    metrics["step_count"] += len(evaluation.steps)
    for step in evaluation.steps:
        observed = step.observed_action
        facets = {
            "product_site": observed.product_site_verified,
            "retained_edits": observed.retained_edits_verified,
            "synthon_partition": observed.synthon_partition_verified,
            "exact_precursors": observed.exact_precursors_verified,
            "strategy": observed.strategy_verified,
            "realization": observed.realization_verified,
            "operator_roundtrip": observed.operator_roundtrip_verified,
            "search": observed.search_eligible,
        }
        for name, available in facets.items():
            metrics["eligibility_facets"][name] += int(available)
        metrics["core_quality"][observed.core_quality_status] += 1
        metrics["limitations"].update(observed.limitations)
        metrics["operator_admission"][
            observed.operator_admission_reason or "accepted"
        ] += 1
        metrics["departing_edit_descriptors"].update(
            observed.departing_edit_descriptors
        )
        metrics["named_annotations"][
            observed.named_annotation or "unannotated"
        ] += 1
        if observed.retained_edits_verified and not observed.operator_roundtrip_verified:
            metrics["promoted"]["all"] += 1
            metrics["promoted"][
                "strategy" if observed.strategy_verified else "retained_only"
            ] += 1
        if not observed.search_eligible:
            continue
        metrics["search_eligible_step_count"] += 1
        if step.search_status != "searched":
            continue
        metrics["searched_step_count"] += 1
        metrics["steps_with_candidates"] += int(bool(step.candidates))
        metrics["candidate_count"] += len(step.candidates)
        metrics["candidate_supervision"].update(
            candidate.supervision_label for candidate in step.candidates
        )
        metrics["source_patent_overlap_count"] += int(
            step.source_patent_precedent_overlap
        )
        metrics["outcomes"][step.outcome] += 1
        rank_values = {
            "exact_precursor": (
                step.exact_precursor_rank,
                observed.exact_precursors_verified,
            ),
            "site": (step.site_rank, observed.product_site_verified),
            "operator": (step.operator_rank, observed.retained_edits_verified),
            "synthon": (step.synthon_rank, observed.synthon_partition_verified),
            "strategy": (step.strategy_rank, observed.strategy_verified),
        }
        for name, (rank, available) in rank_values.items():
            if not available:
                continue
            metrics["rank_denominators"][name] += 1
            for limit in (1, 5, 10, 25):
                metrics["ranks"][name][limit] += int(
                    rank is not None and rank <= limit
                )
        diagnostics = step.search_diagnostics
        for name in (
            "proposed_action_count",
            "validation_attempt_count",
            "valid_action_count",
        ):
            metrics["diagnostics"][name] += int(diagnostics.get(name) or 0)


def _finalize_metrics(
    metrics: dict[str, Any],
    *,
    top_k: int,
) -> dict[str, Any]:
    search_eligible = int(metrics["search_eligible_step_count"])
    searched = int(metrics["searched_step_count"])
    search_denominator = max(1, searched)
    step_count = int(metrics["step_count"])
    facet_counts = {
        name: int(value)
        for name, value in sorted(metrics["eligibility_facets"].items())
    }
    recall_limits = tuple(limit for limit in (1, 5, 10, 25) if limit <= top_k)
    return {
        "route_count": int(metrics["route_count"]),
        "step_count": step_count,
        "eligibility_facet_counts": facet_counts,
        "eligibility_facet_fractions": {
            name: count / max(1, step_count)
            for name, count in facet_counts.items()
        },
        "search_eligible_step_count": search_eligible,
        "searched_step_count": searched,
        "steps_with_candidates": int(metrics["steps_with_candidates"]),
        "candidate_coverage": (
            int(metrics["steps_with_candidates"]) / search_denominator
            if searched
            else None
        ),
        "candidate_count": int(metrics["candidate_count"]),
        "mean_candidates_per_search_eligible_step": (
            int(metrics["candidate_count"]) / search_denominator
            if searched
            else None
        ),
        "source_patent_precedent_overlap_count": int(
            metrics["source_patent_overlap_count"]
        ),
        "core_quality_status_counts": dict(sorted(metrics["core_quality"].items())),
        "limitation_counts": dict(sorted(metrics["limitations"].items())),
        "operator_admission_counts": dict(
            sorted(metrics["operator_admission"].items())
        ),
        "departing_edit_descriptor_counts": dict(
            sorted(metrics["departing_edit_descriptors"].items())
        ),
        "named_annotation_counts": dict(
            sorted(metrics["named_annotations"].items())
        ),
        "promoted_step_count": int(metrics["promoted"]["all"]),
        "promoted_strategy_step_count": int(metrics["promoted"]["strategy"]),
        "promoted_retained_only_step_count": int(
            metrics["promoted"]["retained_only"]
        ),
        "outcome_counts": dict(sorted(metrics["outcomes"].items())),
        "candidate_supervision_label_counts": dict(
            sorted(metrics["candidate_supervision"].items())
        ),
        "search_diagnostics": dict(sorted(metrics["diagnostics"].items())),
        "recall": {
            name: {
                "denominator": int(metrics["rank_denominators"][name]),
                **{
                    f"top{limit}": counts[limit]
                    / max(1, int(metrics["rank_denominators"][name]))
                    if int(metrics["rank_denominators"][name])
                    else None
                    for limit in recall_limits
                },
            }
            for name, counts in metrics["ranks"].items()
        },
    }


def convert_route_action_corpus(
    source_trees: str | Path,
    operator_library: str | Path,
    output_jsonl: str | Path,
    *,
    manifest_path: Optional[str | Path] = None,
    config: RouteActionEvaluationConfig = RouteActionEvaluationConfig(),
    sample_size: Optional[int] = None,
    seed: int = DEFAULT_ROUTE_ACTION_SAMPLE_SEED,
    overwrite: bool = False,
    strict: bool = True,
    workers: int = 1,
    shard_count: int = 1,
    shard_index: int = 0,
    resume: bool = False,
) -> dict[str, Any]:
    """Replay a route-tree corpus and write deterministic route-level JSONL."""

    source = Path(source_trees)
    library_path = Path(operator_library)
    output = Path(output_jsonl)
    manifest = (
        Path(manifest_path)
        if manifest_path is not None
        else Path(f"{output}.manifest.json")
    )
    for required in (source, library_path):
        if not required.is_file():
            raise FileNotFoundError(required)
    if workers < 1:
        raise ValueError("workers must be positive")
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    if shard_index < 0 or shard_index >= shard_count:
        raise ValueError("shard_index must be in [0, shard_count)")
    if overwrite and resume:
        raise ValueError("overwrite and resume are mutually exclusive")
    if resume and (output.exists() or manifest.exists()):
        if not output.is_file() or not manifest.is_file():
            raise ValueError("resume requires both output and manifest")
        report = json.loads(manifest.read_text(encoding="utf-8"))
        expected_selection = {
            "sample_size": sample_size,
            "seed": seed if sample_size is not None else None,
            "shard_count": shard_count,
            "shard_index": shard_index,
        }
        actual_selection = dict(report.get("selection") or {})
        comparable_selection = {
            key: actual_selection.get(key) for key in expected_selection
        }
        checks = {
            "schema": report.get("route_action_schema_version")
            == ROUTE_ACTION_EVALUATION_SCHEMA_VERSION,
            "algorithm": report.get("route_action_algorithm_version")
            == ROUTE_ACTION_EVALUATION_ALGORITHM_VERSION,
            "converter": report.get("converter_version")
            == ROUTE_ACTION_CONVERTER_VERSION,
            "label_policy": report.get("label_policy_definition_id")
            == load_observed_route_action_label_policy().definition_id,
            "source": (report.get("source") or {}).get("sha256")
            == _sha256_file(source),
            "library": (report.get("operator_library") or {}).get("sha256")
            == _sha256_file(library_path),
            "config": report.get("search_config") == config.to_dict(),
            "selection": comparable_selection == expected_selection,
            "output": (report.get("output") or {}).get("sha256")
            == _sha256_file(output),
        }
        failed = sorted(name for name, passed in checks.items() if not passed)
        if failed:
            raise ValueError(
                "existing shard is not resumable; mismatched " + ", ".join(failed)
            )
        resumed_count = 0
        for evaluation in iter_route_action_evaluations(output):
            route_id = evaluation.source_route_id or evaluation.tree_id
            if stable_route_shard(route_id, shard_count) != shard_index:
                raise ValueError("existing shard contains a misassigned route")
            if evaluation.search_config_id != config.config_id:
                raise ValueError("existing shard contains a mismatched search config")
            resumed_count += 1
        if resumed_count != int((report.get("output") or {}).get("record_count") or 0):
            raise ValueError("existing shard record count contradicts its manifest")
        return report
    if not overwrite:
        for candidate in (output, manifest):
            if candidate.exists():
                raise FileExistsError(candidate)

    metrics = _new_metrics()
    split_metrics: dict[str, dict[str, Any]] = {}
    rejections: Counter[str] = Counter()
    rejection_examples: dict[str, list[str]] = {}
    started = time.perf_counter()
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    )
    os.close(descriptor)
    temporary_output = Path(temporary_name)
    try:
        with temporary_output.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_handle, mtime=0
            ) as gzip_handle:
                with io.TextIOWrapper(gzip_handle, encoding="utf-8") as handle:
                    for route_id, evaluation, error_type, detail in _evaluation_results(
                        source,
                        library_path,
                        config=config,
                        sample_size=sample_size,
                        seed=seed,
                        workers=workers,
                        shard_count=shard_count,
                        shard_index=shard_index,
                    ):
                        if evaluation is None:
                            rejections[error_type] += 1
                            examples = rejection_examples.setdefault(error_type, [])
                            if len(examples) < 5:
                                examples.append(f"{route_id}: {detail}")
                            if strict:
                                raise RuntimeError(
                                    f"Route-action replay failed for {route_id}: {detail}"
                                )
                            continue
                        handle.write(
                            json.dumps(
                                evaluation.to_dict(),
                                sort_keys=True,
                                separators=(",", ":"),
                            )
                        )
                        handle.write("\n")
                        _accumulate(metrics, evaluation)
                        split_key = evaluation.split or "unknown"
                        scoped = split_metrics.setdefault(split_key, _new_metrics())
                        _accumulate(scoped, evaluation)
        os.replace(temporary_output, output)
    except BaseException:
        try:
            temporary_output.unlink()
        except FileNotFoundError:
            pass
        raise

    elapsed = time.perf_counter() - started
    report = {
        "route_action_schema_version": ROUTE_ACTION_EVALUATION_SCHEMA_VERSION,
        "route_action_algorithm_version": (
            ROUTE_ACTION_EVALUATION_ALGORITHM_VERSION
        ),
        "converter_version": ROUTE_ACTION_CONVERTER_VERSION,
        "label_policy_definition_id": (
            load_observed_route_action_label_policy().definition_id
        ),
        "source": {
            "path": str(source.resolve()),
            "sha256": _sha256_file(source),
        },
        "operator_library": {
            "path": str(library_path.resolve()),
            "sha256": _sha256_file(library_path),
        },
        "selection": {
            "sample_size": sample_size,
            "seed": seed if sample_size is not None else None,
            "workers": workers,
            "shard_count": shard_count,
            "shard_index": shard_index,
        },
        "search_config": config.to_dict(),
        "metrics": _finalize_metrics(metrics, top_k=config.top_k),
        "metrics_by_split": {
            key: _finalize_metrics(value, top_k=config.top_k)
            for key, value in sorted(split_metrics.items())
        },
        "conversion": {
            "rejected_route_count": sum(rejections.values()),
            "rejection_counts": dict(sorted(rejections.items())),
            "rejection_examples": rejection_examples,
            "elapsed_seconds": round(elapsed, 3),
            "processed_steps_per_second": round(
                int(metrics["step_count"]) / max(elapsed, 1e-9),
                3,
            ),
        },
        "output": {
            "path": str(output.resolve()),
            "sha256": _sha256_file(output),
            "record_count": int(metrics["route_count"]),
        },
        "warnings": [
            "Eligibility is reported independently for each supervision task.",
            "Candidate generation and hard graph validation remain deterministic.",
            "Review-grade route labels do not relax operator-library admission.",
            "Unchosen alternatives are not asserted to be chemically negative.",
            "Patent overlap flags cover visible retained precedents only.",
            "A leakage-controlled library is required for release-quality metrics.",
            "Elapsed time is manifest metadata and is not part of record identity.",
        ],
    }
    _atomic_write_json(manifest, report)
    return report


def merge_route_action_shards(
    shard_outputs: Sequence[str | Path],
    output_jsonl: str | Path,
    *,
    manifest_path: Optional[str | Path] = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Validate and deterministically merge a complete route-action shard set."""

    if not shard_outputs:
        raise ValueError("at least one shard output is required")
    output = Path(output_jsonl)
    manifest = (
        Path(manifest_path)
        if manifest_path is not None
        else Path(f"{output}.manifest.json")
    )
    if not overwrite:
        for candidate in (output, manifest):
            if candidate.exists():
                raise FileExistsError(candidate)

    entries = []
    common: Optional[dict[str, Any]] = None
    expected_count: Optional[int] = None
    for raw_path in shard_outputs:
        shard_path = Path(raw_path)
        shard_manifest_path = Path(f"{shard_path}.manifest.json")
        if not shard_path.is_file() or not shard_manifest_path.is_file():
            raise FileNotFoundError(
                shard_path if not shard_path.is_file() else shard_manifest_path
            )
        shard_manifest = json.loads(
            shard_manifest_path.read_text(encoding="utf-8")
        )
        if shard_manifest.get("converter_version") != ROUTE_ACTION_CONVERTER_VERSION:
            raise ValueError(
                f"shard converter version mismatch: {shard_manifest_path}"
            )
        if (
            shard_manifest.get("route_action_schema_version")
            != ROUTE_ACTION_EVALUATION_SCHEMA_VERSION
            or shard_manifest.get("route_action_algorithm_version")
            != ROUTE_ACTION_EVALUATION_ALGORITHM_VERSION
            or shard_manifest.get("label_policy_definition_id")
            != load_observed_route_action_label_policy().definition_id
        ):
            raise ValueError(f"shard contract mismatch: {shard_manifest_path}")
        selection = dict(shard_manifest.get("selection") or {})
        shard_count = int(selection.get("shard_count") or 0)
        raw_shard_index = selection.get("shard_index")
        shard_index = int(raw_shard_index) if raw_shard_index is not None else -1
        if shard_count < 1 or shard_index < 0 or shard_index >= shard_count:
            raise ValueError(f"invalid shard metadata: {shard_manifest_path}")
        if expected_count is None:
            expected_count = shard_count
        elif shard_count != expected_count:
            raise ValueError("shards disagree on shard_count")
        comparable = {
            "route_action_schema_version": shard_manifest.get(
                "route_action_schema_version"
            ),
            "route_action_algorithm_version": shard_manifest.get(
                "route_action_algorithm_version"
            ),
            "converter_version": shard_manifest.get("converter_version"),
            "label_policy_definition_id": shard_manifest.get(
                "label_policy_definition_id"
            ),
            "source_sha256": (shard_manifest.get("source") or {}).get("sha256"),
            "library_sha256": (shard_manifest.get("operator_library") or {}).get(
                "sha256"
            ),
            "search_config": shard_manifest.get("search_config"),
            "sample_size": selection.get("sample_size"),
            "seed": selection.get("seed"),
            "shard_count": shard_count,
        }
        if common is None:
            common = comparable
        elif comparable != common:
            raise ValueError("shard manifests describe incompatible conversions")
        if (shard_manifest.get("output") or {}).get("sha256") != _sha256_file(
            shard_path
        ):
            raise ValueError(f"shard output checksum mismatch: {shard_path}")
        entries.append((shard_index, shard_path, shard_manifest))

    assert expected_count is not None and common is not None
    indices = [index for index, _, _ in entries]
    if len(set(indices)) != len(indices):
        raise ValueError("duplicate shard index")
    if set(indices) != set(range(expected_count)):
        missing = sorted(set(range(expected_count)) - set(indices))
        raise ValueError(f"incomplete shard set; missing indices {missing}")

    records: dict[str, RouteActionEvaluation] = {}
    for shard_index, shard_path, shard_manifest in sorted(entries):
        shard_record_count = 0
        for evaluation in iter_route_action_evaluations(shard_path):
            route_id = evaluation.source_route_id or evaluation.tree_id
            if stable_route_shard(route_id, expected_count) != shard_index:
                raise ValueError(
                    f"route {route_id} is stored in the wrong shard"
                )
            if route_id in records:
                raise ValueError(f"duplicate route identity across shards: {route_id}")
            if evaluation.search_config_id != dict(common["search_config"]).get(
                "config_id"
            ):
                raise ValueError(f"route {route_id} has a mismatched search config")
            records[route_id] = evaluation
            shard_record_count += 1
        if shard_record_count != int(
            (shard_manifest.get("output") or {}).get("record_count") or 0
        ):
            raise ValueError(
                f"shard record count contradicts manifest: {shard_path}"
            )

    metrics = _new_metrics()
    split_metrics: dict[str, dict[str, Any]] = {}
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    )
    os.close(descriptor)
    temporary_output = Path(temporary_name)
    try:
        with temporary_output.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_handle, mtime=0
            ) as gzip_handle:
                with io.TextIOWrapper(gzip_handle, encoding="utf-8") as handle:
                    for route_id in sorted(records):
                        evaluation = records[route_id]
                        handle.write(
                            json.dumps(
                                evaluation.to_dict(),
                                sort_keys=True,
                                separators=(",", ":"),
                            )
                        )
                        handle.write("\n")
                        _accumulate(metrics, evaluation)
                        split_key = evaluation.split or "unknown"
                        _accumulate(
                            split_metrics.setdefault(split_key, _new_metrics()),
                            evaluation,
                        )
        os.replace(temporary_output, output)
    except BaseException:
        try:
            temporary_output.unlink()
        except FileNotFoundError:
            pass
        raise

    top_k = int(dict(common["search_config"]).get("top_k") or 25)
    elapsed = sum(
        float((item[2].get("conversion") or {}).get("elapsed_seconds") or 0)
        for item in entries
    )
    rejection_counts: Counter[str] = Counter()
    rejection_examples: dict[str, list[str]] = {}
    for _, _, item in entries:
        conversion = dict(item.get("conversion") or {})
        rejection_counts.update(conversion.get("rejection_counts") or {})
        for reason, examples in (conversion.get("rejection_examples") or {}).items():
            retained = rejection_examples.setdefault(str(reason), [])
            retained.extend(str(example) for example in examples[: 5 - len(retained)])
    first_manifest = sorted(entries)[0][2]
    report = {
        "route_action_schema_version": common["route_action_schema_version"],
        "route_action_algorithm_version": common[
            "route_action_algorithm_version"
        ],
        "converter_version": ROUTE_ACTION_CONVERTER_VERSION,
        "label_policy_definition_id": common["label_policy_definition_id"],
        "source": first_manifest["source"],
        "operator_library": first_manifest["operator_library"],
        "selection": {
            "sample_size": common["sample_size"],
            "seed": common["seed"],
            "shard_count": expected_count,
            "merged_shard_indices": list(range(expected_count)),
        },
        "search_config": common["search_config"],
        "metrics": _finalize_metrics(metrics, top_k=top_k),
        "metrics_by_split": {
            key: _finalize_metrics(value, top_k=top_k)
            for key, value in sorted(split_metrics.items())
        },
        "conversion": {
            "rejected_route_count": sum(rejection_counts.values()),
            "rejection_counts": dict(sorted(rejection_counts.items())),
            "rejection_examples": rejection_examples,
            "summed_shard_elapsed_seconds": round(elapsed, 3),
        },
        "shards": [
            {
                "shard_index": index,
                "path": str(path.resolve()),
                "sha256": (item.get("output") or {}).get("sha256"),
                "record_count": (item.get("output") or {}).get("record_count"),
            }
            for index, path, item in sorted(entries)
        ],
        "output": {
            "path": str(output.resolve()),
            "sha256": _sha256_file(output),
            "record_count": len(records),
        },
        "warnings": list(first_manifest.get("warnings") or ()),
    }
    _atomic_write_json(manifest, report)
    return report


def iter_route_action_evaluations(
    source: str | Path,
) -> Iterator[RouteActionEvaluation]:
    """Yield validated route-action evaluations from JSONL or gzip JSONL."""

    path = Path(source)
    if path.suffix.lower() == ".gz":
        handle = gzip.open(path, "rt", encoding="utf-8")
    else:
        handle = path.open("r", encoding="utf-8")
    with handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(
                    f"Route-action record {line_number} must be an object"
                )
            yield RouteActionEvaluation.from_dict(value)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_trees")
    parser.add_argument("operator_library")
    parser.add_argument("output_jsonl")
    parser.add_argument("--manifest")
    parser.add_argument("--sample-size", type=int)
    parser.add_argument("--seed", type=int, default=DEFAULT_ROUTE_ACTION_SAMPLE_SEED)
    parser.add_argument("--top-k", type=int, default=25)
    parser.add_argument("--max-templates-to-apply", type=int, default=500)
    parser.add_argument("--max-candidates-to-validate", type=int, default=100)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--labels-only", action="store_true")
    parser.add_argument("--allow-rejections", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run route-action replay conversion as a standalone command."""

    arguments = _parser().parse_args(argv)
    report = convert_route_action_corpus(
        arguments.source_trees,
        arguments.operator_library,
        arguments.output_jsonl,
        manifest_path=arguments.manifest,
        config=RouteActionEvaluationConfig(
            top_k=arguments.top_k,
            max_templates_to_apply=arguments.max_templates_to_apply,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
            run_search=not arguments.labels_only,
        ),
        sample_size=arguments.sample_size,
        seed=arguments.seed,
        overwrite=arguments.overwrite,
        strict=not arguments.allow_rejections,
        workers=arguments.workers,
        shard_count=arguments.shard_count,
        shard_index=arguments.shard_index,
        resume=arguments.resume,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "DEFAULT_ROUTE_ACTION_SAMPLE_SEED",
    "ROUTE_ACTION_CONVERTER_VERSION",
    "convert_route_action_corpus",
    "iter_route_action_evaluations",
    "merge_route_action_shards",
    "stable_route_shard",
]
