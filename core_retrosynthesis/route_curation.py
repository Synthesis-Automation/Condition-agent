"""Deterministic curation of external multistep retrosynthesis routes.

The raw higher-level retrosynthesis release stores reactions in an unordered
list.  This module validates molecular evidence, reconstructs a linear route
from molecule identity, and writes a small patent-disjoint corpus without
modifying the source snapshot.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import heapq
import io
import json
import os
import tempfile
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Optional, Sequence, TextIO

from rdkit import Chem, rdBase


ROUTE_CORPUS_SCHEMA_VERSION = "1.0"
ROUTE_CURATION_POLICY_VERSION = "1.0"
DEFAULT_SELECTION_SEED = "higher-level-route-corpus-v1"
DEFAULT_SPLIT_SEED = "higher-level-route-splits-v1"


@dataclass(frozen=True)
class RouteSubsetPolicy:
    """Versioned, chemistry-agnostic policy for a compact route corpus."""

    minimum_steps: int = 3
    maximum_steps: int = 6
    maximum_routes: int = 5_000
    require_linear: bool = True
    require_single_subtree: bool = True
    require_abstraction_reduction: bool = True
    train_percent: int = 80
    validation_percent: int = 10
    selection_seed: str = DEFAULT_SELECTION_SEED
    split_seed: str = DEFAULT_SPLIT_SEED
    policy_version: str = ROUTE_CURATION_POLICY_VERSION

    def __post_init__(self) -> None:
        if self.minimum_steps < 2:
            raise ValueError("minimum_steps must be at least two")
        if self.maximum_steps < self.minimum_steps:
            raise ValueError("maximum_steps must not be below minimum_steps")
        if self.maximum_routes < 1:
            raise ValueError("maximum_routes must be positive")
        if not 0 < self.train_percent < 100:
            raise ValueError("train_percent must be between zero and 100")
        if not 0 <= self.validation_percent < 100:
            raise ValueError(
                "validation_percent must be between zero and 99"
            )
        if self.train_percent + self.validation_percent >= 100:
            raise ValueError("train and validation percentages must total below 100")


class RouteQualityError(ValueError):
    """A typed route rejection with a stable machine-readable reason."""

    def __init__(self, reason: str, detail: str = "") -> None:
        self.reason = reason
        self.detail = detail
        message = reason if not detail else f"{reason}: {detail}"
        super().__init__(message)


@dataclass(frozen=True)
class _ParsedReaction:
    source_reaction_id: str
    reaction_smiles: str
    abstracted_reaction_smiles: str
    reactants_smiles: str
    reagents_smiles: str
    product_smiles_mapped: str
    product_smiles: str
    precursor_smiles: tuple[str, ...]
    reaction_identity: str


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _iter_jsonl(path: Path) -> Iterator[tuple[int, dict[str, Any]]]:
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise RouteQualityError(
                    "invalid_json",
                    f"line {line_number}: {exc}",
                ) from exc
            if not isinstance(value, dict):
                raise RouteQualityError(
                    "record_not_object",
                    f"line {line_number}",
                )
            yield line_number, value


def _stable_integer(seed: str, namespace: str, value: str) -> int:
    payload = f"{seed}\0{namespace}\0{value}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest(), "big")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _patent_id(route_id: object) -> str:
    if not isinstance(route_id, str) or "_" not in route_id:
        raise RouteQualityError("invalid_route_id", repr(route_id))
    patent_id, suffix = route_id.rsplit("_", 1)
    if not patent_id or not suffix:
        raise RouteQualityError("invalid_route_id", route_id)
    return patent_id


def _source_reaction_id(full_id: object, subtree_id: object) -> str:
    prefix = f"{subtree_id}_"
    if not isinstance(full_id, str) or not full_id.startswith(prefix):
        raise RouteQualityError(
            "reaction_id_prefix_mismatch",
            repr(full_id),
        )
    value = full_id[len(prefix) :]
    if not value:
        raise RouteQualityError("empty_source_reaction_id")
    return value


def _canonical_unmapped_molecule(molecule: Chem.Mol) -> str:
    value = Chem.Mol(molecule)
    for atom in value.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(value, canonical=True, isomericSmiles=True)


def _canonical_unmapped_smiles(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise RouteQualityError("molecule_parse_failed", smiles[:160])
    return _canonical_unmapped_molecule(molecule)


def _validate_complete_mapping(
    reactants: Chem.Mol,
    product: Chem.Mol,
) -> None:
    reactant_maps = [atom.GetAtomMapNum() for atom in reactants.GetAtoms()]
    product_maps = [atom.GetAtomMapNum() for atom in product.GetAtoms()]
    if (
        not reactant_maps
        or not product_maps
        or any(value <= 0 for value in reactant_maps + product_maps)
    ):
        raise RouteQualityError("incomplete_atom_mapping")
    if len(reactant_maps) != len(set(reactant_maps)):
        raise RouteQualityError("duplicate_reactant_atom_map")
    if len(product_maps) != len(set(product_maps)):
        raise RouteQualityError("duplicate_product_atom_map")
    if not set(product_maps).issubset(reactant_maps):
        raise RouteQualityError("product_atom_without_reactant_source")


def _validate_abstracted_reaction(reaction_smiles: str) -> None:
    if reaction_smiles.count(">") != 2:
        raise RouteQualityError("invalid_abstracted_reaction_syntax")
    reactants, reagents, products = reaction_smiles.split(">")
    if (
        not reactants
        or not products
        or Chem.MolFromSmiles(reactants) is None
        or Chem.MolFromSmiles(products) is None
        or (reagents and Chem.MolFromSmiles(reagents) is None)
    ):
        raise RouteQualityError("abstracted_reaction_parse_failed")


def _parse_reaction(
    reaction: Mapping[str, Any],
    subtree_id: str,
) -> _ParsedReaction:
    reaction_smiles = reaction.get("reaction_smiles")
    if not isinstance(reaction_smiles, str) or reaction_smiles.count(">") != 2:
        raise RouteQualityError("invalid_reaction_syntax")
    reactants, reagents, product = reaction_smiles.split(">")
    if not reactants or not product:
        raise RouteQualityError("empty_reactant_or_product")
    reactant_molecule = Chem.MolFromSmiles(reactants)
    product_molecule = Chem.MolFromSmiles(product)
    if reactant_molecule is None or product_molecule is None:
        raise RouteQualityError("reaction_parse_failed")
    if len(Chem.GetMolFrags(product_molecule)) != 1:
        raise RouteQualityError("multiple_product_components")
    if reagents and Chem.MolFromSmiles(reagents) is None:
        raise RouteQualityError("reagent_parse_failed")
    _validate_complete_mapping(reactant_molecule, product_molecule)

    precursor_values = tuple(
        sorted(
            _canonical_unmapped_smiles(component)
            for component in reactants.split(".")
        )
    )
    product_value = _canonical_unmapped_molecule(product_molecule)
    if product_value in precursor_values:
        raise RouteQualityError("unchanged_product_as_reactant")

    abstracted = reaction.get("abstracted_reaction_smiles", "")
    if not isinstance(abstracted, str):
        raise RouteQualityError("invalid_abstracted_reaction_type")
    if abstracted:
        _validate_abstracted_reaction(abstracted)

    source_id = _source_reaction_id(reaction.get("_id"), subtree_id)
    reaction_identity = f"{'.'.join(precursor_values)}>>{product_value}"
    return _ParsedReaction(
        source_reaction_id=source_id,
        reaction_smiles=reaction_smiles,
        abstracted_reaction_smiles=abstracted,
        reactants_smiles=reactants,
        reagents_smiles=reagents,
        product_smiles_mapped=product,
        product_smiles=product_value,
        precursor_smiles=precursor_values,
        reaction_identity=reaction_identity,
    )


def _cheap_rejection_reason(
    route: Mapping[str, Any],
    policy: RouteSubsetPolicy,
) -> Optional[str]:
    try:
        _patent_id(route.get("route_id"))
        original_tree = route["original_tree"]
        reaction_count = original_tree["num_reactions"]
        depth = original_tree["depth"]
        reaction_ids = original_tree["reaction_ids"]
        subtrees = route["subtrees"]
    except (KeyError, TypeError):
        return "invalid_route_schema"
    if not isinstance(reaction_count, int) or not (
        policy.minimum_steps <= reaction_count <= policy.maximum_steps
    ):
        return "outside_step_range"
    if policy.require_linear and depth != reaction_count:
        return "nonlinear_route"
    if not isinstance(reaction_ids, list) or len(reaction_ids) != reaction_count:
        return "reaction_id_count_mismatch"
    if len(reaction_ids) != len(set(reaction_ids)):
        return "duplicate_source_reaction_id"
    if not isinstance(subtrees, list) or not subtrees:
        return "missing_subtree"
    if policy.require_single_subtree and len(subtrees) != 1:
        return "multiple_subtrees"
    subtree = subtrees[0]
    try:
        subtree_counts = subtree["num_reactions"]
        reactions = subtree["reactions"]
    except (KeyError, TypeError):
        return "invalid_subtree_schema"
    if (
        not isinstance(subtree_counts, list)
        or len(subtree_counts) != 2
        or subtree_counts[0] != reaction_count
        or not isinstance(reactions, list)
        or len(reactions) != reaction_count
    ):
        return "incomplete_subtree"
    if (
        policy.require_abstraction_reduction
        and subtree_counts[1] >= subtree_counts[0]
    ):
        return "no_abstraction_reduction"
    for reaction in reactions:
        if not isinstance(reaction, dict):
            return "invalid_reaction_record"
        reaction_smiles = reaction.get("reaction_smiles")
        if not isinstance(reaction_smiles, str) or reaction_smiles.count(">") != 2:
            return "invalid_reaction_syntax"
    return None


def normalize_linear_route_record(
    route: Mapping[str, Any],
    *,
    split: str,
    excluded_reaction_identities: frozenset[str] = frozenset(),
    excluded_target_smiles: frozenset[str] = frozenset(),
) -> dict[str, Any]:
    """Validate and normalize one source route into retrosynthetic order.

    Array order from the source is ignored.  Step order is reconstructed from
    canonical, atom-map-free molecule identity while the mapped reaction and
    reagent fields are retained as source evidence.
    """

    route_id = route.get("route_id")
    patent_id = _patent_id(route_id)
    original_tree = route.get("original_tree")
    subtrees = route.get("subtrees")
    if not isinstance(original_tree, dict) or not isinstance(subtrees, list):
        raise RouteQualityError("invalid_route_schema")
    if len(subtrees) != 1 or not isinstance(subtrees[0], dict):
        raise RouteQualityError("route_not_single_subtree")
    subtree = subtrees[0]
    subtree_id = subtree.get("subtree_id")
    if not isinstance(subtree_id, str):
        raise RouteQualityError("invalid_subtree_id")
    raw_reactions = subtree.get("reactions")
    if not isinstance(raw_reactions, list) or not raw_reactions:
        raise RouteQualityError("missing_reactions")

    parsed = tuple(
        _parse_reaction(reaction, subtree_id) for reaction in raw_reactions
    )
    source_ids = {reaction.source_reaction_id for reaction in parsed}
    expected_source_ids = original_tree.get("reaction_ids")
    if not isinstance(expected_source_ids, list) or source_ids != set(
        expected_source_ids
    ):
        raise RouteQualityError("source_reaction_ids_do_not_match")

    product_to_reaction: dict[str, _ParsedReaction] = {}
    for reaction in parsed:
        if reaction.product_smiles in product_to_reaction:
            raise RouteQualityError("duplicate_step_product")
        product_to_reaction[reaction.product_smiles] = reaction
    product_set = set(product_to_reaction)
    precursor_occurrences = Counter(
        precursor
        for reaction in parsed
        for precursor in reaction.precursor_smiles
    )
    root_products = tuple(
        product for product in product_set if precursor_occurrences[product] == 0
    )
    if len(root_products) != 1:
        raise RouteQualityError(
            "ambiguous_route_target",
            f"found {len(root_products)} roots",
        )
    root_product = root_products[0]
    if root_product in excluded_target_smiles:
        raise RouteQualityError("excluded_test_target")

    ordered: list[tuple[_ParsedReaction, tuple[str, ...], tuple[str, ...]]] = []
    visited: set[str] = set()
    current_product = root_product
    while current_product in product_to_reaction:
        if current_product in visited:
            raise RouteQualityError("cyclic_route")
        visited.add(current_product)
        reaction = product_to_reaction[current_product]
        if reaction.reaction_identity in excluded_reaction_identities:
            raise RouteQualityError("excluded_test_reaction")
        internal = tuple(
            value for value in reaction.precursor_smiles if value in product_set
        )
        terminal = tuple(
            value for value in reaction.precursor_smiles if value not in product_set
        )
        ordered.append((reaction, internal, terminal))
        if not internal:
            break
        if len(internal) != 1:
            raise RouteQualityError("route_not_linear_after_reconstruction")
        if precursor_occurrences[internal[0]] != 1:
            raise RouteQualityError("intermediate_reused")
        current_product = internal[0]

    if len(ordered) != len(parsed) or len(visited) != len(product_set):
        raise RouteQualityError("disconnected_route")
    for position, (_, internal, _) in enumerate(ordered):
        if position < len(ordered) - 1 and len(internal) != 1:
            raise RouteQualityError("broken_linear_chain")
        if position == len(ordered) - 1 and internal:
            raise RouteQualityError("unterminated_linear_chain")

    subtree_counts = subtree.get("num_reactions")
    subtree_depths = subtree.get("depth")
    if (
        not isinstance(subtree_counts, list)
        or len(subtree_counts) != 2
        or not isinstance(subtree_depths, list)
        or len(subtree_depths) != 2
    ):
        raise RouteQualityError("invalid_subtree_metrics")

    steps = []
    for position, (reaction, internal, terminal) in enumerate(ordered):
        steps.append(
            {
                "retrosynthetic_position": position,
                "source_reaction_id": reaction.source_reaction_id,
                "reaction_smiles": reaction.reaction_smiles,
                "abstracted_reaction_smiles": (
                    reaction.abstracted_reaction_smiles
                ),
                "reactants_smiles": reaction.reactants_smiles,
                "reagents_smiles": reaction.reagents_smiles,
                "product_smiles_mapped": reaction.product_smiles_mapped,
                "product_smiles": reaction.product_smiles,
                "precursor_smiles": list(reaction.precursor_smiles),
                "internal_precursor_smiles": list(internal),
                "terminal_precursor_smiles": list(terminal),
            }
        )

    return {
        "schema_version": ROUTE_CORPUS_SCHEMA_VERSION,
        "route_id": route_id,
        "patent_id": patent_id,
        "split": split,
        "target_smiles": root_product,
        "source_subtree_id": subtree_id,
        "original_reaction_count": len(parsed),
        "original_depth": original_tree.get("depth"),
        "higher_level_reaction_count": subtree_counts[1],
        "higher_level_depth": subtree_depths[1],
        "abstraction_reduction": len(parsed) - subtree_counts[1],
        "quality": {
            "atom_mapping": "complete_unique_product_subset",
            "topology": "connected_linear_unique_target",
            "source_reaction_ids_verified": True,
            "source_array_order_used": False,
        },
        "steps": steps,
    }


def _canonical_reaction_identity(reaction_smiles: str) -> tuple[str, str]:
    if reaction_smiles.count(">") != 2:
        raise RouteQualityError("invalid_test_reaction_syntax")
    reactants, _, product = reaction_smiles.split(">")
    precursor_values = sorted(
        _canonical_unmapped_smiles(value) for value in reactants.split(".")
    )
    product_value = _canonical_unmapped_smiles(product)
    return f"{'.'.join(precursor_values)}>>{product_value}", product_value


def load_reaction_exclusions(
    path: Optional[str | Path],
) -> tuple[frozenset[str], frozenset[str]]:
    """Load canonical reaction and product exclusions from a SMILES file."""

    if path is None:
        return frozenset(), frozenset()
    reaction_ids: set[str] = set()
    products: set[str] = set()
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            value = line.strip()
            if not value:
                continue
            try:
                reaction_id, product = _canonical_reaction_identity(value)
            except RouteQualityError as exc:
                raise RouteQualityError(
                    exc.reason,
                    f"test-set line {line_number}: {exc.detail}",
                ) from exc
            reaction_ids.add(reaction_id)
            products.add(product)
    return frozenset(reaction_ids), frozenset(products)


def _split_for_patent(patent_id: str, policy: RouteSubsetPolicy) -> str:
    bucket = _stable_integer(policy.split_seed, "patent", patent_id) % 100
    if bucket < policy.train_percent:
        return "train"
    if bucket < policy.train_percent + policy.validation_percent:
        return "validation"
    return "test"


def _write_deterministic_jsonl_gzip(
    path: Path,
    records: Iterable[Mapping[str, Any]],
) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            mtime=0,
        ) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8") as text_handle:
                for record in records:
                    text_handle.write(
                        json.dumps(
                            record,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                    )
                    text_handle.write("\n")


def _atomic_json_write(path: Path, value: Mapping[str, Any]) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def curate_route_subset(
    source_routes: str | Path,
    output_jsonl: str | Path,
    *,
    manifest_path: Optional[str | Path] = None,
    testset_path: Optional[str | Path] = None,
    policy: RouteSubsetPolicy = RouteSubsetPolicy(),
    overwrite: bool = False,
) -> dict[str, Any]:
    """Create a deterministic, patent-disjoint, normalized route subset."""

    source = Path(source_routes)
    output = Path(output_jsonl)
    manifest = (
        Path(manifest_path)
        if manifest_path is not None
        else Path(f"{output}.manifest.json")
    )
    if not source.is_file():
        raise FileNotFoundError(source)
    if testset_path is not None and not Path(testset_path).is_file():
        raise FileNotFoundError(testset_path)
    if not overwrite:
        existing = [path for path in (output, manifest) if path.exists()]
        if existing:
            raise FileExistsError(existing[0])
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    scan_counts: Counter[str] = Counter()
    best_route_by_patent: dict[str, tuple[int, str]] = {}
    for _, route in _iter_jsonl(source):
        scan_counts["source_routes"] += 1
        reason = _cheap_rejection_reason(route, policy)
        if reason is not None:
            scan_counts[f"rejected_{reason}"] += 1
            continue
        scan_counts["cheap_eligible_routes"] += 1
        route_id = str(route["route_id"])
        patent_id = _patent_id(route_id)
        rank = _stable_integer(policy.selection_seed, "route", route_id)
        previous = best_route_by_patent.get(patent_id)
        if previous is None or (rank, route_id) < previous:
            best_route_by_patent[patent_id] = (rank, route_id)

    selected_route_ids = {
        route_id: patent_id
        for patent_id, (_, route_id) in best_route_by_patent.items()
    }
    excluded_reactions, excluded_targets = load_reaction_exclusions(testset_path)
    validation_counts: Counter[str] = Counter()
    retained: list[tuple[int, str, dict[str, Any]]] = []
    for _, route in _iter_jsonl(source):
        route_id = route.get("route_id")
        patent_id = selected_route_ids.get(route_id)
        if patent_id is None:
            continue
        validation_counts["candidate_patent_routes"] += 1
        split = _split_for_patent(patent_id, policy)
        try:
            normalized = normalize_linear_route_record(
                route,
                split=split,
                excluded_reaction_identities=excluded_reactions,
                excluded_target_smiles=excluded_targets,
            )
        except RouteQualityError as exc:
            validation_counts[f"rejected_{exc.reason}"] += 1
            continue
        validation_counts["fully_valid_routes"] += 1
        patent_rank = _stable_integer(
            policy.selection_seed,
            "patent",
            patent_id,
        )
        entry = (-patent_rank, str(route_id), normalized)
        if len(retained) < policy.maximum_routes:
            heapq.heappush(retained, entry)
        elif patent_rank < -retained[0][0]:
            heapq.heapreplace(retained, entry)

    if len(retained) < policy.maximum_routes:
        raise RuntimeError(
            "Only "
            f"{len(retained)} routes passed validation; "
            f"{policy.maximum_routes} requested"
        )
    ordered = sorted(
        ((-rank, route_id, record) for rank, route_id, record in retained),
        key=lambda value: (value[0], value[1]),
    )
    records = [record for _, _, record in ordered]

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
    )
    os.close(descriptor)
    temporary_output = Path(temporary_name)
    try:
        _write_deterministic_jsonl_gzip(temporary_output, records)
        os.replace(temporary_output, output)
    except BaseException:
        try:
            temporary_output.unlink()
        except FileNotFoundError:
            pass
        raise

    split_counts = Counter(record["split"] for record in records)
    step_counts = Counter(record["original_reaction_count"] for record in records)
    reduction_counts = Counter(record["abstraction_reduction"] for record in records)
    report: dict[str, Any] = {
        "schema_version": ROUTE_CORPUS_SCHEMA_VERSION,
        "policy": asdict(policy),
        "source": {
            "route_file": str(source.resolve()),
            "route_file_sha256": _sha256_file(source),
            "testset_file": (
                str(Path(testset_path).resolve())
                if testset_path is not None
                else None
            ),
            "testset_file_sha256": (
                _sha256_file(Path(testset_path))
                if testset_path is not None
                else None
            ),
            "excluded_reaction_count": len(excluded_reactions),
            "excluded_target_count": len(excluded_targets),
        },
        "toolkit": {"rdkit_version": rdBase.rdkitVersion},
        "counts": {
            "scan": dict(sorted(scan_counts.items())),
            "cheap_eligible_patents": len(best_route_by_patent),
            "validation": dict(sorted(validation_counts.items())),
            "selected_routes": len(records),
            "selected_patents": len({record["patent_id"] for record in records}),
            "splits": dict(sorted(split_counts.items())),
            "steps": {str(key): step_counts[key] for key in sorted(step_counts)},
            "abstraction_reduction": {
                str(key): reduction_counts[key] for key in sorted(reduction_counts)
            },
        },
        "output": {
            "route_file": str(output.resolve()),
            "route_file_sha256": _sha256_file(output),
            "record_count": len(records),
        },
        "interpretation_warning": (
            "Higher-level abstractions are algorithm-generated labels from the "
            "source publication, not direct observations of chemist intent."
        ),
    }
    _atomic_json_write(manifest, report)
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Curate a deterministic higher-level route POC corpus"
    )
    parser.add_argument("source_routes")
    parser.add_argument("output_jsonl")
    parser.add_argument("--manifest")
    parser.add_argument("--testset")
    parser.add_argument("--minimum-steps", type=int, default=3)
    parser.add_argument("--maximum-steps", type=int, default=6)
    parser.add_argument("--maximum-routes", type=int, default=5_000)
    parser.add_argument("--allow-no-abstraction-reduction", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run route curation as a standalone module."""

    arguments = _parser().parse_args(argv)
    policy = RouteSubsetPolicy(
        minimum_steps=arguments.minimum_steps,
        maximum_steps=arguments.maximum_steps,
        maximum_routes=arguments.maximum_routes,
        require_abstraction_reduction=(
            not arguments.allow_no_abstraction_reduction
        ),
    )
    report = curate_route_subset(
        arguments.source_routes,
        arguments.output_jsonl,
        manifest_path=arguments.manifest,
        testset_path=arguments.testset,
        policy=policy,
        overwrite=arguments.overwrite,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "ROUTE_CORPUS_SCHEMA_VERSION",
    "ROUTE_CURATION_POLICY_VERSION",
    "RouteQualityError",
    "RouteSubsetPolicy",
    "curate_route_subset",
    "load_reaction_exclusions",
    "normalize_linear_route_record",
]
