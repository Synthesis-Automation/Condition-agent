"""Deterministic corpus conversion for route-core projections."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import random
import tempfile
from typing import Any, Iterator, Optional, Sequence

from .route_conversion import iter_route_trees
from .route_core import (
    ROUTE_CORE_ALGORITHM_VERSION,
    ROUTE_CORE_SCHEMA_VERSION,
    RouteCoreProjection,
    build_route_core_projection,
)


DEFAULT_ROUTE_CORE_SAMPLE_SEED = 20_260_814
ROUTE_CORE_CONVERTER_VERSION = "1.0"


def _projection_task(tree: Any) -> tuple[str, RouteCoreProjection | None, str, str]:
    """Build one projection in a process-safe, context-preserving envelope."""

    route_id = tree.source_route_id or tree.tree_id
    try:
        return route_id, build_route_core_projection(tree), "", ""
    except (RuntimeError, TypeError, ValueError) as exc:
        return route_id, None, type(exc).__name__, str(exc)


def _sha256_file(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def _selected_trees(
    source: Path,
    *,
    sample_size: Optional[int],
    seed: int,
) -> Iterator[Any]:
    if sample_size is None:
        yield from iter_route_trees(source)
        return
    trees = tuple(
        sorted(
            iter_route_trees(source),
            key=lambda tree: tree.source_route_id or tree.tree_id,
        )
    )
    if sample_size < 1 or sample_size > len(trees):
        raise ValueError(f"sample_size must be between 1 and {len(trees)}")
    yield from random.Random(seed).sample(trees, sample_size)


def _projection_results(
    source: Path,
    *,
    sample_size: Optional[int],
    seed: int,
    workers: int,
) -> Iterator[tuple[str, RouteCoreProjection | None, str, str]]:
    trees = _selected_trees(source, sample_size=sample_size, seed=seed)
    if workers == 1:
        for tree in trees:
            yield _projection_task(tree)
        return
    materialized = tuple(trees)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        yield from executor.map(_projection_task, materialized, chunksize=8)


def convert_route_core_corpus(
    source_trees: str | Path,
    output_jsonl: str | Path,
    *,
    manifest_path: Optional[str | Path] = None,
    sample_size: Optional[int] = None,
    seed: int = DEFAULT_ROUTE_CORE_SAMPLE_SEED,
    overwrite: bool = False,
    strict: bool = True,
    workers: int = 1,
) -> dict[str, Any]:
    """Convert canonical route trees to deterministic route-core JSONL."""

    source = Path(source_trees)
    output = Path(output_jsonl)
    manifest = (
        Path(manifest_path)
        if manifest_path is not None
        else Path(f"{output}.manifest.json")
    )
    if not source.is_file():
        raise FileNotFoundError(source)
    if workers < 1:
        raise ValueError("workers must be positive")
    if not overwrite:
        for candidate in (output, manifest):
            if candidate.exists():
                raise FileExistsError(candidate)

    route_count = 0
    reaction_count = 0
    signature_count = 0
    core_count = 0
    display_count = 0
    fully_chemistry_resolved_count = 0
    fully_lineage_connected_count = 0
    quality_counts: Counter[str] = Counter()
    lineage_counts: Counter[str] = Counter()
    split_counts: Counter[str] = Counter()
    depth_counts: Counter[int] = Counter()
    typed_motif_counts: Counter[str] = Counter()
    shape_motif_counts: Counter[str] = Counter()
    motif_definitions: dict[str, dict[str, Any]] = {}
    rejection_counts: Counter[str] = Counter()
    rejection_examples: dict[str, list[str]] = {}

    output.parent.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
    )
    os.close(descriptor)
    temporary_output = Path(temporary_name)
    try:
        with temporary_output.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_handle,
                mtime=0,
            ) as gzip_handle:
                with io.TextIOWrapper(gzip_handle, encoding="utf-8") as handle:
                    for route_id, projection, error_type, error_detail in _projection_results(
                        source,
                        sample_size=sample_size,
                        seed=seed,
                        workers=workers,
                    ):
                        if projection is None:
                            rejection_counts[error_type] += 1
                            examples = rejection_examples.setdefault(error_type, [])
                            if len(examples) < 5:
                                examples.append(route_id)
                            if strict:
                                raise RuntimeError(
                                    "Route-core conversion failed for "
                                    f"{route_id}: {error_detail}"
                                )
                            continue
                        handle.write(
                            json.dumps(
                                projection.to_dict(),
                                sort_keys=True,
                                separators=(",", ":"),
                            )
                        )
                        handle.write("\n")
                        route_count += 1
                        reaction_count += projection.reaction_count
                        split_counts[projection.split or "unknown"] += 1
                        depth_counts[projection.maximum_depth] += 1
                        fully_chemistry_resolved_count += int(
                            projection.fully_chemistry_resolved
                        )
                        fully_lineage_connected_count += int(
                            projection.fully_lineage_connected
                        )
                        typed_motif_counts.update(projection.typed_motif_keys)
                        shape_motif_counts.update(projection.shape_motif_keys)
                        for motif in projection.motifs:
                            motif_definitions.setdefault(
                                motif.motif_id,
                                motif.to_dict(),
                            )
                        for step in projection.steps:
                            quality_counts[step.quality_status] += 1
                            signature_count += int(
                                step.reaction_signature_id is not None
                            )
                            core_count += int(step.reaction_core_id is not None)
                            display_count += int(
                                step.render_reaction_smiles is not None
                            )
                        lineage_counts.update(
                            link.status for link in projection.lineage_links
                        )
        os.replace(temporary_output, output)
    except BaseException:
        try:
            temporary_output.unlink()
        except FileNotFoundError:
            pass
        raise

    report: dict[str, Any] = {
        "route_core_schema_version": ROUTE_CORE_SCHEMA_VERSION,
        "route_core_algorithm_version": ROUTE_CORE_ALGORITHM_VERSION,
        "converter_version": ROUTE_CORE_CONVERTER_VERSION,
        "source": {
            "path": str(source.resolve()),
            "sha256": _sha256_file(source),
        },
        "selection": {
            "sample_size": sample_size,
            "seed": seed if sample_size is not None else None,
            "workers": workers,
        },
        "conversion": {
            "route_count": route_count,
            "reaction_count": reaction_count,
            "reaction_signature_count": signature_count,
            "reaction_core_count": core_count,
            "display_projection_count": display_count,
            "fully_chemistry_resolved_route_count": (
                fully_chemistry_resolved_count
            ),
            "fully_lineage_connected_route_count": (
                fully_lineage_connected_count
            ),
            "quality_status_counts": dict(sorted(quality_counts.items())),
            "lineage_status_counts": dict(sorted(lineage_counts.items())),
            "split_counts": dict(sorted(split_counts.items())),
            "maximum_depth_counts": {
                str(key): depth_counts[key] for key in sorted(depth_counts)
            },
            "rejected_count": sum(rejection_counts.values()),
            "rejection_counts": dict(sorted(rejection_counts.items())),
            "rejection_examples": rejection_examples,
        },
        "motifs": {
            "unique_typed_motif_count": len(typed_motif_counts),
            "unique_shape_motif_count": len(shape_motif_counts),
            "top_typed_motifs": [
                {
                    **motif_definitions[key],
                    "route_count": count,
                }
                for key, count in typed_motif_counts.most_common(25)
            ],
            "top_shape_motifs": [
                {
                    **motif_definitions[key],
                    "route_count": count,
                }
                for key, count in shape_motif_counts.most_common(25)
            ],
        },
        "output": {
            "path": str(output.resolve()),
            "sha256": _sha256_file(output),
            "record_count": route_count,
        },
        "warnings": [
            "Route-core projections are not executable multistep templates.",
            "Cross-step symmetry ambiguity is retained as candidate lineage.",
            "Algorithm-generated higher-level reactions are weak annotations only.",
        ],
    }
    manifest.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return report


def iter_route_core_projections(
    source: str | Path,
) -> Iterator[RouteCoreProjection]:
    """Yield validated typed route-core projections from JSONL or gzip JSONL."""

    path = Path(source)
    opener = gzip.open if path.suffix.lower() == ".gz" else Path.open
    if path.suffix.lower() == ".gz":
        handle = opener(path, "rt", encoding="utf-8")
    else:
        handle = opener(path, "r", encoding="utf-8")
    with handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(
                    f"Route-core record {line_number} is not an object"
                )
            yield RouteCoreProjection.from_dict(value)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build chemistry-rich route-core projections"
    )
    parser.add_argument("source_trees")
    parser.add_argument("output_jsonl")
    parser.add_argument("--manifest")
    parser.add_argument("--sample-size", type=int)
    parser.add_argument("--seed", type=int, default=DEFAULT_ROUTE_CORE_SAMPLE_SEED)
    parser.add_argument("--allow-rejections", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run route-core corpus conversion."""

    arguments = _parser().parse_args(argv)
    report = convert_route_core_corpus(
        arguments.source_trees,
        arguments.output_jsonl,
        manifest_path=arguments.manifest,
        sample_size=arguments.sample_size,
        seed=arguments.seed,
        overwrite=arguments.overwrite,
        strict=not arguments.allow_rejections,
        workers=arguments.workers,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "DEFAULT_ROUTE_CORE_SAMPLE_SEED",
    "ROUTE_CORE_CONVERTER_VERSION",
    "convert_route_core_corpus",
    "iter_route_core_projections",
]
