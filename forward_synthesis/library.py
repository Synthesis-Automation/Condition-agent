"""Build, validate, persist, and retrieve forward reaction operators."""

from __future__ import annotations

import gzip
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from reactive_taxonomy import (
    BidirectionalReactionOperator,
    ReactionOperatorPrecedent,
    apply_forward_operator,
    canonical_molecule_collection,
    reverse_recovers_precursors,
)
from reactive_taxonomy.reaction_operators import clear_operator_application_cache

from .models import (
    FORWARD_LIBRARY_SCHEMA_VERSION,
    ForwardOperatorLibrary,
    ForwardPrecursorIndex,
)


def _value(source: Any, name: str, default: Any = None) -> Any:
    if isinstance(source, Mapping):
        return source.get(name, default)
    return getattr(source, name, default)


def _without_stereo(smiles: str) -> Optional[str]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    Chem.RemoveStereochemistry(molecule)
    return canonical_molecule_collection(
        Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    )


def _product_matches(expected: str, observed: str, stereo_policy: str) -> bool:
    expected_canonical = canonical_molecule_collection(expected)
    if expected_canonical is None:
        return False
    if expected_canonical == observed:
        return True
    return bool(
        stereo_policy == "relaxed"
        and _without_stereo(expected_canonical) == _without_stereo(observed)
    )


def _as_precedent(value: Any) -> ReactionOperatorPrecedent:
    return ReactionOperatorPrecedent(
        reaction_id=str(_value(value, "reaction_id", "") or ""),
        reference_id=str(_value(value, "reference_id", "") or ""),
        precursor_smiles=str(_value(value, "precursor_smiles", "") or ""),
        product_smiles=str(_value(value, "product_smiles", "") or ""),
    )


def _as_operator(template: Any) -> BidirectionalReactionOperator:
    precursor_smarts = str(_value(template, "precursor_smarts", "") or "")
    product_smarts = str(_value(template, "product_smarts", "") or "")
    return BidirectionalReactionOperator(
        operator_id=str(_value(template, "operator_id", "") or ""),
        realization_id=str(_value(template, "realization_id", "") or ""),
        template_id=str(_value(template, "template_id", "") or ""),
        abstraction_level=str(_value(template, "abstraction_level", "") or ""),
        forward_smarts=f"{precursor_smarts}>>{product_smarts}",
        reverse_smarts=f"{product_smarts}>>{precursor_smarts}",
        precursor_smarts=precursor_smarts,
        product_smarts=product_smarts,
        edit_tokens=tuple(_value(template, "edit_tokens", ()) or ()),
        operator_signature=str(_value(template, "operator_signature", "") or ""),
        stereo_policy=str(_value(template, "stereo_policy", "exact") or "exact"),
        observation_support=int(_value(template, "observation_support", 0) or 0),
        independent_reference_support=int(
            _value(template, "independent_reference_support", 0) or 0
        ),
        precedents=tuple(
            _as_precedent(value) for value in (_value(template, "precedents", ()) or ())
        ),
        named_annotations=tuple(_value(template, "named_annotations", ()) or ()),
    )


def _source_round_trip(operator: BidirectionalReactionOperator) -> bool:
    if not operator.precedents:
        return False
    for precedent in operator.precedents:
        outcomes = apply_forward_operator(
            operator,
            precedent.precursor_smiles,
            max_assignments=64,
            max_outcomes=256,
        )
        matching = tuple(
            outcome
            for outcome in outcomes
            if _product_matches(
                precedent.product_smiles,
                outcome.product_smiles,
                operator.stereo_policy,
            )
        )
        if not matching or not any(
            reverse_recovers_precursors(operator, outcome) for outcome in matching
        ):
            return False
    return True


def _operator_source_contract(operator: BidirectionalReactionOperator) -> tuple[Any, ...]:
    """Fields retained from the first admitted source, before support merging."""

    return (
        operator.operator_id, operator.realization_id, operator.template_id,
        operator.abstraction_level, operator.forward_smarts, operator.reverse_smarts,
        operator.precursor_smarts, operator.product_smarts, operator.edit_tokens,
        operator.operator_signature, operator.stereo_policy, operator.schema_version,
    )


def validate_forward_library_source(
    forward_library: ForwardOperatorLibrary,
    generic_library: Any,
) -> None:
    """Reject a prepared library whose admitted content contradicts its source.

    This inexpensive metadata check accepts the same structural source contract
    as :func:`build_forward_library` and performs no reaction execution, SMARTS
    compilation, or source re-admission. Source chemistry was checked offline.
    The existing artifact records a source definition and template count, not a
    digest of rejected templates; this check therefore establishes compatibility
    of admitted operators rather than byte identity of the whole source artifact.

    Duplicate directional operators may merge support, annotations and source
    precedents. Some duplicates can independently fail source admission, so the
    saved projection must match a nonempty subset of their source projections.
    """

    if (
        forward_library.schema_version != FORWARD_LIBRARY_SCHEMA_VERSION
        or forward_library.definition_id != "forward_operator_library.v1"
    ):
        raise ValueError("unsupported prepared forward-library schema or definition")
    source_definition = _value(generic_library, "definition", {}) or {}
    definition_id = str(_value(source_definition, "definition_id", "") or "")
    if (
        not definition_id
        or forward_library.source_library_definition_id != definition_id
    ):
        raise ValueError("prepared forward-library source definition does not match")
    templates = tuple(_value(generic_library, "templates", ()) or ())
    if forward_library.source_template_count != len(templates):
        raise ValueError("prepared forward-library source template count does not match")

    sources: dict[str, list[BidirectionalReactionOperator]] = {}
    for template in templates:
        try:
            source = _as_operator(template)
        except (TypeError, ValueError):
            # Invalid source contracts can legitimately have been rejected at build.
            continue
        sources.setdefault(source.forward_operator_id, []).append(source)

    for operator in forward_library.operators:
        precedents = frozenset(operator.precedents)
        annotations = frozenset(operator.named_annotations)
        compatible = tuple(
            source for source in sources.get(operator.forward_operator_id, ())
            if frozenset(source.precedents) <= precedents
            and frozenset(source.named_annotations) <= annotations
            and source.observation_support <= operator.observation_support
            and source.independent_reference_support <= operator.independent_reference_support
        )
        if not compatible or not any(
            _operator_source_contract(source) == _operator_source_contract(operator)
            for source in compatible
        ):
            raise ValueError(
                "prepared forward operator does not match its source contract: "
                + operator.forward_operator_id
            )
        if (
            frozenset(item for source in compatible for item in source.precedents)
            != precedents
            or frozenset(item for source in compatible for item in source.named_annotations)
            != annotations
            or max(source.observation_support for source in compatible)
            != operator.observation_support
            or max(source.independent_reference_support for source in compatible)
            != operator.independent_reference_support
        ):
            raise ValueError(
                "prepared forward operator provenance or support does not match its source: "
                + operator.forward_operator_id
            )


def _source_round_trip_passes(operator: BidirectionalReactionOperator) -> bool:
    """Apply the same source-admission gate in serial or worker processes."""

    try:
        return _source_round_trip(operator)
    except Exception:
        return False


def _worker_source_round_trip(operator: BidirectionalReactionOperator) -> bool:
    """Bound compiled-query memory in an independent source-check worker."""

    try:
        return _source_round_trip_passes(operator)
    finally:
        clear_operator_application_cache()


def _required_atomic_numbers(
    operator: BidirectionalReactionOperator,
) -> tuple[int, ...]:
    reaction = rdChemReactions.ReactionFromSmarts(operator.forward_smarts)
    if reaction is None:
        return ()
    values = []
    for index in range(reaction.GetNumReactantTemplates()):
        for atom in reaction.GetReactantTemplate(index).GetAtoms():
            atomic_number = int(atom.GetAtomicNum())
            if atomic_number > 0:
                values.append(atomic_number)
    return tuple(sorted(values))


def build_forward_precursor_index(
    operators: Iterable[BidirectionalReactionOperator],
) -> ForwardPrecursorIndex:
    """Build a conservative element and component-count retrieval index."""

    by_count: dict[int, list[str]] = {}
    requirements = {}
    for operator in operators:
        reaction = rdChemReactions.ReactionFromSmarts(operator.forward_smarts)
        if reaction is None:
            continue
        operator_id = operator.forward_operator_id
        by_count.setdefault(int(reaction.GetNumReactantTemplates()), []).append(
            operator_id
        )
        requirements[operator_id] = _required_atomic_numbers(operator)
    return ForwardPrecursorIndex(
        component_count_to_operator_ids={
            count: tuple(sorted(set(operator_ids)))
            for count, operator_ids in sorted(by_count.items())
        },
        operator_required_atomic_numbers=dict(sorted(requirements.items())),
    )


def build_forward_library(
    generic_library: Any,
    *,
    require_source_round_trip: bool = True,
    workers: int = 1,
) -> ForwardOperatorLibrary:
    """Project a generic library into independently forward-admitted operators.

    The input is deliberately structural: either an object exposing ``templates``
    or its serialized mapping.  This package therefore does not import or own a
    retrosynthesis model. Offline builds can distribute source checks across
    workers; result ordering and admission rules are unchanged.
    """

    if workers < 1:
        raise ValueError("workers must be positive")
    templates = tuple(_value(generic_library, "templates", ()) or ())
    rejection_counts: Counter[str] = Counter()
    admitted: dict[str, BidirectionalReactionOperator] = {}
    candidates = []
    for template in templates:
        try:
            operator = _as_operator(template)
            rdChemReactions.ReactionFromSmarts(operator.forward_smarts)
        except Exception:
            rejection_counts["invalid_operator_contract"] += 1
            continue
        candidates.append(operator)

    def admit(operator: BidirectionalReactionOperator) -> None:
        current = admitted.get(operator.forward_operator_id)
        if current is None:
            admitted[operator.forward_operator_id] = operator
            return
        precedents = {
            (
                item.reaction_id,
                item.reference_id,
                item.precursor_smiles,
                item.product_smiles,
            ): item
            for item in (*current.precedents, *operator.precedents)
        }
        admitted[operator.forward_operator_id] = replace(
            current,
            observation_support=max(
                current.observation_support,
                operator.observation_support,
            ),
            independent_reference_support=max(
                current.independent_reference_support,
                operator.independent_reference_support,
            ),
            precedents=tuple(precedents[key] for key in sorted(precedents)),
            named_annotations=tuple(
                sorted(set(current.named_annotations + operator.named_annotations))
            ),
        )
    parallel = workers > 1 and require_source_round_trip
    executor = ProcessPoolExecutor(max_workers=workers) if parallel else nullcontext()
    with executor as pool:
        if not require_source_round_trip:
            checks = (True for _ in candidates)
        elif pool is not None:
            checks = pool.map(_worker_source_round_trip, candidates, chunksize=16)
        else:
            checks = map(_source_round_trip_passes, candidates)
        for operator, accepted in zip(candidates, checks):
            if accepted:
                admit(operator)
            else:
                rejection_counts["source_forward_round_trip_failed"] += 1
    operators = tuple(
        sorted(
            admitted.values(),
            key=lambda item: (
                item.operator_id,
                item.abstraction_level,
                item.forward_operator_id,
            ),
        )
    )
    source_definition = _value(generic_library, "definition", {}) or {}
    return ForwardOperatorLibrary(
        operators=operators,
        source_template_count=len(templates),
        admitted_operator_count=len(operators),
        rejection_counts=dict(sorted(rejection_counts.items())),
        precursor_index=build_forward_precursor_index(operators),
        source_library_definition_id=str(
            _value(source_definition, "definition_id", "")
            if not isinstance(source_definition, Mapping)
            else source_definition.get("definition_id") or ""
        ),
    )


def _library_from_dict(value: Mapping[str, Any]) -> ForwardOperatorLibrary:
    raw_index = value.get("precursor_index") or {}
    index = ForwardPrecursorIndex(
        component_count_to_operator_ids={
            int(count): tuple(str(item) for item in operator_ids or ())
            for count, operator_ids in (
                raw_index.get("component_count_to_operator_ids") or {}
            ).items()
        },
        operator_required_atomic_numbers={
            str(operator_id): tuple(int(item) for item in atomic_numbers or ())
            for operator_id, atomic_numbers in (
                raw_index.get("operator_required_atomic_numbers") or {}
            ).items()
        },
        definition_id=str(
            raw_index.get("definition_id") or "forward_precursor_index.v1"
        ),
    )
    operators = tuple(
        BidirectionalReactionOperator.from_dict(item)
        for item in value.get("operators") or ()
    )
    current_ids = {operator.forward_operator_id for operator in operators}
    indexed_ids = set(index.operator_required_atomic_numbers)
    if indexed_ids != current_ids:
        # Directional IDs include the application-engine version. Rebuild this
        # derived index when loading an artifact created by an older engine;
        # the admitted chemistry and its source-round-trip evidence are intact.
        index = build_forward_precursor_index(operators)
    return ForwardOperatorLibrary(
        operators=operators,
        source_template_count=int(value.get("source_template_count") or 0),
        admitted_operator_count=int(value.get("admitted_operator_count") or 0),
        rejection_counts={
            str(key): int(count)
            for key, count in (value.get("rejection_counts") or {}).items()
        },
        precursor_index=index,
        definition_id=str(value.get("definition_id") or "forward_operator_library.v1"),
        source_library_definition_id=str(
            value.get("source_library_definition_id") or ""
        ),
        schema_version=str(value.get("schema_version") or "1.0"),
    )


def save_forward_library(
    library: ForwardOperatorLibrary,
    path: str | Path,
) -> None:
    """Write a deterministic JSON or gzip-compressed JSON library."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        library.to_dict(),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    if destination.suffix.casefold() == ".gz":
        with gzip.open(destination, "wt", encoding="utf-8", newline="\n") as handle:
            handle.write(payload)
            handle.write("\n")
    else:
        destination.write_text(payload + "\n", encoding="utf-8")


def load_forward_library(path: str | Path) -> ForwardOperatorLibrary:
    """Load and validate a JSON or gzip-compressed JSON forward library."""

    source = Path(path)
    if source.suffix.casefold() == ".gz":
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            value = json.load(handle)
    else:
        value = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError("forward library must contain one JSON object")
    return _library_from_dict(value)


def indexed_forward_operators(
    starting_materials: str,
    library: ForwardOperatorLibrary,
    *,
    allow_self_reaction: bool = False,
) -> tuple[BidirectionalReactionOperator, ...]:
    """Retrieve operators using only conservative precursor-observable facts."""

    canonical = canonical_molecule_collection(starting_materials)
    if canonical is None:
        return ()
    molecule = Chem.MolFromSmiles(canonical)
    if molecule is None:
        return ()
    input_components = len(canonical.split("."))
    atomic_counts = Counter(
        int(atom.GetAtomicNum())
        for atom in molecule.GetAtoms()
        if int(atom.GetAtomicNum()) > 0
    )
    eligible_ids = set()
    for (
        required_count,
        operator_ids,
    ) in library.precursor_index.component_count_to_operator_ids.items():
        if required_count <= input_components or allow_self_reaction:
            eligible_ids.update(operator_ids)
    selected = []
    for operator in library.operators:
        operator_id = operator.forward_operator_id
        if operator_id not in eligible_ids:
            continue
        required_counts = Counter(
            library.precursor_index.operator_required_atomic_numbers.get(
                operator_id,
                (),
            )
        )
        sufficient_inventory = not any(
            atomic_counts[number] < count
            for number, count in required_counts.items()
        )
        if not sufficient_inventory:
            if not allow_self_reaction or any(
                atomic_counts[number] < 1 for number in required_counts
            ):
                continue
        selected.append(operator)
    return tuple(selected)


__all__ = [
    "build_forward_library",
    "build_forward_precursor_index",
    "indexed_forward_operators",
    "load_forward_library",
    "save_forward_library",
    "validate_forward_library_source",
]
