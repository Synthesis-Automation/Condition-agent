"""Necessary product-feature index for data-derived graph operators."""

from __future__ import annotations

import json
from collections import defaultdict
from rdkit import Chem

from retrosynthesis_poc.chemistry import canonical_smiles

from .generic_models import (
    GenericRetrievalIndex,
    GenericTemplateLibrary,
)


def _element_pair(first: str, second: str) -> tuple[str, str]:
    return tuple(sorted((first, second)))  # type: ignore[return-value]


def operator_required_product_tokens(operator_signature: str) -> tuple[str, ...]:
    """Extract necessary, product-observable tokens from an edit signature."""

    try:
        edits = json.loads(operator_signature)
    except (TypeError, ValueError):
        return ()
    bond_tokens = set()
    atom_tokens = set()
    for edit in edits if isinstance(edits, list) else ():
        if not isinstance(edit, list) or len(edit) != 4:
            continue
        edit_type, raw_atoms, _, new_order = edit
        atoms = raw_atoms if isinstance(raw_atoms, list) else ()
        elements = [
            str(atom[0])
            for atom in atoms
            if isinstance(atom, list) and atom and str(atom[0]) != "H"
        ]
        atom_tokens.update(f"A:{element}" for element in elements)
        if (
            edit_type in {"formed", "order_changed"}
            and str(new_order) != "NONE"
            and len(elements) == 2
        ):
            first, second = _element_pair(elements[0], elements[1])
            bond_tokens.add(f"B:{first}:{second}:{new_order}")
    return tuple(sorted(bond_tokens or atom_tokens))


def target_product_tokens(target_smiles: str) -> frozenset[str]:
    """Return coarse atom and bond tokens observable in one target product."""

    canonical = canonical_smiles(target_smiles)
    if canonical is None or "." in canonical:
        raise ValueError("target must be one valid molecule")
    molecule = Chem.MolFromSmiles(canonical)
    if molecule is None:
        raise ValueError("target could not be parsed")
    values = {f"A:{atom.GetSymbol()}" for atom in molecule.GetAtoms()}
    for bond in molecule.GetBonds():
        first, second = _element_pair(
            bond.GetBeginAtom().GetSymbol(),
            bond.GetEndAtom().GetSymbol(),
        )
        values.add(f"B:{first}:{second}:{bond.GetBondType().name}")
    return frozenset(values)


def build_generic_retrieval_index(
    library: GenericTemplateLibrary,
) -> GenericRetrievalIndex:
    """Build a deterministic inverted index without excluding valid matches."""

    inverted: dict[str, set[str]] = defaultdict(set)
    requirements = {}
    fallback = []
    for template in library.templates:
        required = operator_required_product_tokens(template.operator_signature)
        requirements[template.template_id] = required
        if not required:
            fallback.append(template.template_id)
            continue
        for token in required:
            inverted[token].add(template.template_id)
    return GenericRetrievalIndex(
        token_to_template_ids={
            token: tuple(sorted(template_ids))
            for token, template_ids in sorted(inverted.items())
        },
        template_required_tokens={
            template_id: tuple(tokens)
            for template_id, tokens in sorted(requirements.items())
        },
        fallback_template_ids=tuple(sorted(fallback)),
    )


def indexed_template_ids(
    target_smiles: str,
    index: GenericRetrievalIndex,
) -> frozenset[str]:
    """Return templates whose necessary product tokens occur in the target."""

    target_tokens = target_product_tokens(target_smiles)
    candidates = set(index.fallback_template_ids)
    for token in target_tokens:
        candidates.update(index.token_to_template_ids.get(token, ()))
    return frozenset(
        template_id
        for template_id in candidates
        if set(index.template_required_tokens.get(template_id, ())).issubset(
            target_tokens
        )
    )


__all__ = [
    "build_generic_retrieval_index",
    "indexed_template_ids",
    "operator_required_product_tokens",
    "target_product_tokens",
]
