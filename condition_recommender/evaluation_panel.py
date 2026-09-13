"""Freeze reproducible evaluation inputs before recommendation outcomes exist."""

from __future__ import annotations

import hashlib
import random
from pathlib import Path
from typing import Any, Iterable

from .conversion.atomic import atomic_json
from .evaluation import grouped_holdout_split, select_evaluation_queries, _leakage_summary
from .generic_indexing import GenericIndexedReaction, build_generic_index_from_rows, load_generic_index
from .sqlite_indexing import save_sqlite_generic_index
from .support import mapping_equivalence_key


def _tokens(row: GenericIndexedReaction) -> set[str]:
    return {f"{key}:{value}" for key, value in (
        ("observation", row.observation_id), ("reference", row.reference_id),
        ("reaction", row.canonical_reaction_id), ("mapping", mapping_equivalence_key(row)),
    ) if value}


def freeze_evaluation_panel(
    source_path: str | Path,
    output_dir: str | Path,
    *,
    size: int,
    seed: int,
    exclusion_paths: Iterable[str | Path] = (),
    test_fraction: float = 0.2,
    max_queries: int = 500,
    split_modes: tuple[str, ...] = ("grouped_random", "scaffold_disjoint"),
) -> dict[str, Any]:
    """Sample source observations and freeze all split/query IDs and code hashes.

    Prior observations, publications, canonical reactions and map-independent
    exact core identities are excluded by direct token overlap. This does not
    assert transitive exclusion through source rows outside the sampled panel.
    Neither retrieval results nor condition outcomes influence selection.
    """
    if size < 2 or max_queries < 1:
        raise ValueError("Panel size must be at least two and max_queries positive")
    destination = Path(output_dir)
    manifest_path = destination / "frozen_panel_manifest.json"
    panel_path = destination / "frozen_panel.sqlite"
    if manifest_path.exists() or panel_path.exists():
        raise FileExistsError("Refusing to overwrite frozen evaluation inputs")
    source = load_generic_index(source_path)
    if size > len(source.rows):
        raise ValueError("Requested panel exceeds source size")
    excluded = set()
    exclusions = []
    for path in exclusion_paths:
        prior = load_generic_index(path)
        for row in prior.rows:
            excluded.update(_tokens(row))
        exclusions.append({"path": str(path), "row_count": len(prior.rows)})
    order = random.Random(seed).sample(range(len(source.rows)), len(source.rows))
    selected = []
    positions = []
    skipped = 0
    for start in range(0, len(order), 500):
        chunk = order[start:start + 500]
        for position, row in zip(chunk, source.select(chunk)):
            if _tokens(row) & excluded:
                skipped += 1
                continue
            selected.append(row)
            positions.append(position)
            if len(selected) == size:
                break
        if len(selected) == size:
            break
    if len(selected) < size:
        raise ValueError("Insufficient observations after prior-evidence exclusions")
    index = build_generic_index_from_rows(selected)
    splits = {}
    for mode in split_modes:
        split = grouped_holdout_split(index.rows, seed=seed, test_fraction=test_fraction, split_mode=mode)
        queries = select_evaluation_queries(split.test_rows, seed=seed, max_queries=max_queries)
        splits[mode] = {
            "train_ids": [row.observation_id for row in split.train_rows],
            "test_ids": [row.observation_id for row in split.test_rows],
            "query_ids": [row.observation_id for row in queries],
            "train_groups": split.train_group_ids, "test_groups": split.test_group_ids,
            **_leakage_summary(split),
        }
    destination.mkdir(parents=True, exist_ok=True)
    save_sqlite_generic_index(index, panel_path)
    root = Path(__file__).resolve().parents[1]
    hashes = {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for package in ("reactive_taxonomy", "condition_registry", "condition_recommender")
        for path in sorted((root / package).rglob("*"))
        if path.suffix in {".py", ".json"} and "__pycache__" not in path.parts
    }
    manifest = {
        "schema_version": "1.0", "status": "frozen_before_retrieval_outcomes",
        "source_path": str(source_path), "source_count": len(source.rows),
        "source_identity": getattr(source.rows, "artifact_identity", None),
        "row_count": size, "positions_in_draw_order": positions,
        "prior_exclusions": exclusions, "skipped_prior_overlap_count": skipped,
        "prior_overlap_policy": "direct observation/reference/reaction/mapping token exclusion",
        "seed": seed, "test_fraction": test_fraction, "max_queries": max_queries,
        "protocol": {"top_k": 5, "minimum_pool_size": 2,
                     "engines": ["baseline", "shared_core_v2"],
                     "no_tuning_after_outcomes": True,
                     "maximum_coverage_loss": 0.05,
                     "maximum_top5_recipe_recovery_loss": 0.02,
                     "independent_precision_review": "pending"},
        "splits": splits, "code_definition_sha256": hashes,
        "panel_sha256": hashlib.sha256(panel_path.read_bytes()).hexdigest(),
    }
    atomic_json(manifest_path, manifest)
    return manifest
