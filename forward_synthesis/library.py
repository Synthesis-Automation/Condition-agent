"""Build, validate, persist, and retrieve forward reaction operators."""

from __future__ import annotations

import gzip
import json
from collections import Counter
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

from .models import ForwardOperatorLibrary, ForwardPrecursorIndex


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
) -> ForwardOperatorLibrary:
    """Project a generic library into independently forward-admitted operators.

    The input is deliberately structural: either an object exposing ``templates``
    or its serialized mapping.  This package therefore does not import or own a
    retrosynthesis model.
    """

    templates = tuple(_value(generic_library, "templates", ()) or ())
    rejection_counts: Counter[str] = Counter()
    admitted: dict[str, BidirectionalReactionOperator] = {}
    for template in templates:
        try:
            operator = _as_operator(template)
            rdChemReactions.ReactionFromSmarts(operator.forward_smarts)
        except Exception:
            rejection_counts["invalid_operator_contract"] += 1
            continue
        if require_source_round_trip:
            try:
                accepted = _source_round_trip(operator)
            except Exception:
                accepted = False
            if not accepted:
                rejection_counts["source_forward_round_trip_failed"] += 1
                continue
        current = admitted.get(operator.forward_operator_id)
        if current is None:
            admitted[operator.forward_operator_id] = operator
            continue
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
    return ForwardOperatorLibrary(
        operators=tuple(
            BidirectionalReactionOperator.from_dict(item)
            for item in value.get("operators") or ()
        ),
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
        if required_count <= input_components:
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
        if any(
            atomic_counts[number] < count for number, count in required_counts.items()
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
]
