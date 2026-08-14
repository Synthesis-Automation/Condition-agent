"""Conversion of curated observed routes into the shared route-tree contract."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import random
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Optional, Sequence

from .chemistry import digest
from .route_contract import (
    ROUTE_TREE_SCHEMA_VERSION,
    MoleculeOccurrenceNode,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    assert_valid_route_tree,
    iter_molecule_occurrences,
)
from .row_io import iter_rows


DEFAULT_OBSERVED_ROUTE_DATASET_ID = "higher_level_retrosynthesis_figshare_v2"
DEFAULT_OBSERVED_ROUTE_SAMPLE_SEED = 20_260_814
ROUTE_CONVERTER_VERSION = "1.0"


class ObservedRouteConversionError(ValueError):
    """A stable conversion rejection with route-level context."""

    def __init__(self, reason: str, detail: str = "") -> None:
        self.reason = reason
        self.detail = detail
        message = reason if not detail else f"{reason}: {detail}"
        super().__init__(message)


def _required_string(
    source: Mapping[str, Any],
    key: str,
    *,
    allow_empty: bool = False,
) -> str:
    value = source.get(key)
    if not isinstance(value, str) or (not value and not allow_empty):
        raise ObservedRouteConversionError("invalid_string_field", key)
    return value


def _required_string_list(source: Mapping[str, Any], key: str) -> tuple[str, ...]:
    value = source.get(key)
    if not isinstance(value, list) or not all(
        isinstance(item, str) and item for item in value
    ):
        raise ObservedRouteConversionError("invalid_string_list_field", key)
    return tuple(sorted(value))


def build_observed_route_tree(
    record: Mapping[str, Any],
    *,
    source_dataset_id: str = DEFAULT_OBSERVED_ROUTE_DATASET_ID,
) -> ReactionRouteTree:
    """Build one observed route tree without trusting source step order."""

    route_id = _required_string(record, "route_id")
    patent_id = _required_string(record, "patent_id")
    target_smiles = _required_string(record, "target_smiles")
    split = _required_string(record, "split")
    steps_value = record.get("steps")
    if not isinstance(steps_value, list) or not steps_value:
        raise ObservedRouteConversionError("missing_route_steps", route_id)

    step_by_product: dict[str, Mapping[str, Any]] = {}
    source_reaction_ids: set[str] = set()
    for step in steps_value:
        if not isinstance(step, dict):
            raise ObservedRouteConversionError("invalid_step_record", route_id)
        product = _required_string(step, "product_smiles")
        source_reaction_id = _required_string(step, "source_reaction_id")
        if product in step_by_product:
            raise ObservedRouteConversionError(
                "duplicate_step_product",
                product,
            )
        if source_reaction_id in source_reaction_ids:
            raise ObservedRouteConversionError(
                "duplicate_source_reaction_id",
                source_reaction_id,
            )
        step_by_product[product] = step
        source_reaction_ids.add(source_reaction_id)
    if target_smiles not in step_by_product:
        raise ObservedRouteConversionError(
            "target_has_no_recorded_reaction",
            target_smiles,
        )

    visited_products: set[str] = set()
    active_products: set[str] = set()
    fingerprint_tokens: list[str] = []

    def molecule_node(
        smiles: str,
        path: str,
        depth: int,
    ) -> MoleculeOccurrenceNode:
        occurrence_id = digest("ROUTEMOL2", route_id, path, smiles)
        step = step_by_product.get(smiles)
        if step is None:
            return MoleculeOccurrenceNode(
                occurrence_id=occurrence_id,
                smiles=smiles,
                depth=depth,
                terminal=True,
                terminal_evidence="observed_route_leaf",
                unresolved_reason=None,
            )
        if smiles in active_products:
            raise ObservedRouteConversionError("cyclic_route", smiles)
        if smiles in visited_products:
            raise ObservedRouteConversionError(
                "step_product_reused",
                smiles,
            )
        active_products.add(smiles)
        visited_products.add(smiles)
        source_reaction_id = _required_string(step, "source_reaction_id")
        reaction_smiles = _required_string(step, "reaction_smiles")
        if reaction_smiles.count(">") != 2:
            raise ObservedRouteConversionError(
                "invalid_reaction_syntax",
                source_reaction_id,
            )
        precursors = _required_string_list(step, "precursor_smiles")
        children = tuple(
            molecule_node(value, f"{path}.{index}", depth + 1)
            for index, value in enumerate(precursors)
        )
        active_products.remove(smiles)
        abstracted = _required_string(
            step,
            "abstracted_reaction_smiles",
            allow_empty=True,
        )
        reaction_depth = depth + 1
        fingerprint_tokens.append(
            f"reaction:{reaction_depth}:{source_reaction_id}:"
            f"{smiles}:{'.'.join(precursors)}"
        )
        evidence = RouteStepEvidence(
            evidence_kind="observed",
            source_dataset_id=source_dataset_id,
            source_route_id=route_id,
            source_reaction_id=source_reaction_id,
            patent_id=patent_id,
            connectivity_method=(
                "reconstructed_from_canonical_molecule_identity"
            ),
            confidence=None,
            reactants_smiles=_required_string(step, "reactants_smiles"),
            reagents_smiles=_required_string(
                step,
                "reagents_smiles",
                allow_empty=True,
            ),
            product_smiles_mapped=_required_string(
                step,
                "product_smiles_mapped",
            ),
            abstracted_reaction_smiles=abstracted,
            abstraction_status=(
                "algorithm_generated"
                if abstracted
                else "absorbed_by_algorithmic_abstraction"
            ),
            warnings=(
                "Route connectivity is inferred from patent-local molecular identity.",
                "Higher-level abstraction is an algorithm-generated weak label.",
            ),
        )
        reaction = RouteReactionNode(
            reaction_node_id=digest(
                "ROUTERXN2",
                route_id,
                path,
                source_reaction_id,
            ),
            step_id=f"{route_id}:{source_reaction_id}",
            depth=reaction_depth,
            reaction_smiles=reaction_smiles,
            evidence=evidence,
            children=children,
        )
        return MoleculeOccurrenceNode(
            occurrence_id=occurrence_id,
            smiles=smiles,
            depth=depth,
            terminal=False,
            terminal_evidence="expanded_observed_intermediate",
            unresolved_reason=None,
            reaction=reaction,
        )

    root = molecule_node(target_smiles, "0", 0)
    if visited_products != set(step_by_product):
        unreachable = sorted(set(step_by_product) - visited_products)
        raise ObservedRouteConversionError(
            "disconnected_route_steps",
            ",".join(unreachable[:5]),
        )
    expected_count = record.get("original_reaction_count")
    if expected_count != len(step_by_product):
        raise ObservedRouteConversionError(
            "source_reaction_count_mismatch",
            f"{expected_count!r}!={len(step_by_product)}",
        )
    maximum_depth = max(
        node.depth for node in iter_molecule_occurrences(
            ReactionRouteTree(
                tree_id="pending",
                route_kind="observed",
                target_smiles=target_smiles,
                root=root,
                reaction_count=len(step_by_product),
                maximum_depth=0,
                fingerprint_tokens=tuple(sorted(fingerprint_tokens)),
            )
        )
    )
    tokens = tuple(sorted(fingerprint_tokens))
    tree = ReactionRouteTree(
        tree_id=digest(
            "ROUTETREE2",
            source_dataset_id,
            route_id,
            target_smiles,
            *tokens,
        ),
        route_kind="observed",
        target_smiles=target_smiles,
        root=root,
        reaction_count=len(step_by_product),
        maximum_depth=maximum_depth,
        fingerprint_tokens=tokens,
        source_dataset_id=source_dataset_id,
        source_route_id=route_id,
        patent_id=patent_id,
        split=split,
        source_record_schema_version=str(record.get("schema_version") or ""),
        connectivity_method="reconstructed_from_canonical_molecule_identity",
        confidence=None,
        higher_level_reaction_count=record.get("higher_level_reaction_count"),
        higher_level_depth=record.get("higher_level_depth"),
        abstraction_reduction=record.get("abstraction_reduction"),
        warnings=(
            "Source reactions are observed, but route connectivity is inferred.",
            "Higher-level abstractions are algorithm-generated weak labels.",
            "One-pot and tandem source records remain atomic reaction steps.",
        ),
    )
    assert_valid_route_tree(tree)
    return tree


def _sha256_file(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def _write_deterministic_gzip(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            mtime=0,
        ) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8") as handle:
                for row in rows:
                    handle.write(
                        json.dumps(
                            row,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                    )
                    handle.write("\n")


def _selected_records(
    records: Sequence[Mapping[str, Any]],
    *,
    sample_size: Optional[int],
    seed: int,
) -> tuple[Mapping[str, Any], ...]:
    ordered = tuple(
        sorted(records, key=lambda row: str(row.get("route_id") or ""))
    )
    if sample_size is None:
        return ordered
    if sample_size < 1 or sample_size > len(ordered):
        raise ValueError(
            f"sample_size must be between 1 and {len(ordered)}"
        )
    return tuple(random.Random(seed).sample(ordered, sample_size))


def convert_observed_route_corpus(
    source_routes: str | Path,
    output_jsonl: str | Path,
    *,
    manifest_path: Optional[str | Path] = None,
    source_dataset_id: str = DEFAULT_OBSERVED_ROUTE_DATASET_ID,
    sample_size: Optional[int] = None,
    seed: int = DEFAULT_OBSERVED_ROUTE_SAMPLE_SEED,
    overwrite: bool = False,
    strict: bool = True,
) -> dict[str, Any]:
    """Convert curated route rows into deterministic nested route-tree JSONL."""

    source = Path(source_routes)
    output = Path(output_jsonl)
    manifest = (
        Path(manifest_path)
        if manifest_path is not None
        else Path(f"{output}.manifest.json")
    )
    if not source.is_file():
        raise FileNotFoundError(source)
    if not overwrite:
        for candidate in (output, manifest):
            if candidate.exists():
                raise FileExistsError(candidate)
    records = tuple(iter_rows(source))
    selected = _selected_records(records, sample_size=sample_size, seed=seed)
    converted: list[ReactionRouteTree] = []
    rejection_counts: Counter[str] = Counter()
    rejection_examples: dict[str, list[str]] = {}
    for record in selected:
        try:
            converted.append(
                build_observed_route_tree(
                    record,
                    source_dataset_id=source_dataset_id,
                )
            )
        except ObservedRouteConversionError as exc:
            rejection_counts[exc.reason] += 1
            examples = rejection_examples.setdefault(exc.reason, [])
            if len(examples) < 5:
                examples.append(str(record.get("route_id") or "unknown"))
    if strict and rejection_counts:
        summary = ", ".join(
            f"{key}={value}" for key, value in sorted(rejection_counts.items())
        )
        raise RuntimeError(f"Observed route conversion rejected records: {summary}")

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
        _write_deterministic_gzip(
            temporary_output,
            (tree.to_dict() for tree in converted),
        )
        os.replace(temporary_output, output)
    except BaseException:
        try:
            temporary_output.unlink()
        except FileNotFoundError:
            pass
        raise

    route_depths = Counter(tree.maximum_depth for tree in converted)
    split_counts = Counter(tree.split or "unknown" for tree in converted)
    report: dict[str, Any] = {
        "route_tree_schema_version": ROUTE_TREE_SCHEMA_VERSION,
        "converter_version": ROUTE_CONVERTER_VERSION,
        "source_dataset_id": source_dataset_id,
        "source": {
            "path": str(source.resolve()),
            "sha256": _sha256_file(source),
            "record_count": len(records),
        },
        "selection": {
            "sample_size": sample_size,
            "seed": seed if sample_size is not None else None,
            "selected_count": len(selected),
        },
        "conversion": {
            "converted_count": len(converted),
            "rejected_count": sum(rejection_counts.values()),
            "rejection_counts": dict(sorted(rejection_counts.items())),
            "rejection_examples": rejection_examples,
            "split_counts": dict(sorted(split_counts.items())),
            "maximum_depth_counts": {
                str(key): route_depths[key] for key in sorted(route_depths)
            },
        },
        "output": {
            "path": str(output.resolve()),
            "sha256": _sha256_file(output),
            "record_count": len(converted),
        },
        "warnings": [
            "Route connectivity is inferred from patent-local molecular identity.",
            "Higher-level abstractions are algorithm-generated weak labels.",
        ],
    }
    manifest.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return report


def iter_route_trees(source: str | Path) -> Iterator[ReactionRouteTree]:
    """Yield validated typed trees from canonical route-tree JSONL."""

    for row in iter_rows(source):
        yield ReactionRouteTree.from_dict(row)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert curated observations into canonical route trees"
    )
    parser.add_argument("source_routes")
    parser.add_argument("output_jsonl")
    parser.add_argument("--manifest")
    parser.add_argument("--source-dataset-id", default=DEFAULT_OBSERVED_ROUTE_DATASET_ID)
    parser.add_argument("--sample-size", type=int)
    parser.add_argument("--seed", type=int, default=DEFAULT_OBSERVED_ROUTE_SAMPLE_SEED)
    parser.add_argument("--allow-rejections", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run observed route-tree conversion."""

    arguments = _parser().parse_args(argv)
    report = convert_observed_route_corpus(
        arguments.source_routes,
        arguments.output_jsonl,
        manifest_path=arguments.manifest,
        source_dataset_id=arguments.source_dataset_id,
        sample_size=arguments.sample_size,
        seed=arguments.seed,
        overwrite=arguments.overwrite,
        strict=not arguments.allow_rejections,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "DEFAULT_OBSERVED_ROUTE_DATASET_ID",
    "DEFAULT_OBSERVED_ROUTE_SAMPLE_SEED",
    "ROUTE_CONVERTER_VERSION",
    "ObservedRouteConversionError",
    "build_observed_route_tree",
    "convert_observed_route_corpus",
    "iter_route_trees",
]
