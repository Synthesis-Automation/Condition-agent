"""Context-aware local before/after labels for normalized reaction edits."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .context import classify_context
from .labels import load_rendering_taxonomy, render_context
from .reaction_models import ReactionComponent, ReactionEdit


@dataclass(frozen=True)
class ContextualTransformationLabel:
    """Local reactant and product motifs derived from one validated edit."""

    before: str
    after: str
    concise: str
    evidence: str
    confidence: float


def _atom_key(edit: ReactionEdit, first: bool) -> Tuple[int, int]:
    atom = edit.atom_1 if first else edit.atom_2
    if atom is None:
        raise ValueError("Contextual bond labels require two endpoints")
    return int(atom.component_index), int(atom.atom_index)


def _hydrogen_delta(
    edits: Sequence[ReactionEdit], component_index: int, atom_index: int
) -> int:
    delta = 0
    for edit in edits:
        if edit.edit_type != "hydrogen_change":
            continue
        atom = edit.atom_1
        if (
            atom.component_index != component_index
            or atom.atom_index != atom_index
        ):
            continue
        delta += 1 if edit.new_order and not edit.old_order else -1
    return delta


def _external_contexts(
    mol: object, atom_index: int, opposite_index: int, *, style: str
) -> Tuple[str, ...]:
    atom = mol.GetAtomWithIdx(atom_index)
    labels = []
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() == opposite_index or neighbor.GetAtomicNum() <= 1:
            continue
        context = classify_context(
            mol,
            neighbor.GetIdx(),
            excluded=(atom_index, opposite_index),
        )
        labels.append(render_context(context.token, style=style))
    return tuple(sorted(labels))


def _hydrogen_text(count: int) -> str:
    if count <= 0:
        return ""
    return "H" if count == 1 else f"H{count}"


def _endpoint_label(
    *, element: str, hydrogens: int, contexts: Tuple[str, ...], side: str, bond: str
) -> str:
    hydrogen = _hydrogen_text(hydrogens)
    if element == "C":
        center = f"C{hydrogen}" if side == "right" or contexts else f"{hydrogen}C"
    elif side == "right":
        center = f"{element}{hydrogen}"
    else:
        center = f"{hydrogen}{element}"
    if not contexts:
        return center
    primary, *additional = contexts
    suffix = "".join(f"({context})" for context in additional)
    if side == "left":
        if element == "C":
            center = f"C{hydrogen}"
        return f"{primary}{bond}{center}{suffix}"
    return f"{center}{suffix}{bond}{primary}"


def _bond_symbol(order: Optional[str], *, style: str) -> str:
    styling = load_rendering_taxonomy()["styles"][style]
    return {
        "SINGLE": styling["bond"],
        "DOUBLE": styling["double"],
        "TRIPLE": styling["triple"],
        "AROMATIC": ":",
    }.get(str(order or "SINGLE").upper(), styling["bond"])


def build_contextual_transformation_label(
    reactants: Tuple[ReactionComponent, ...],
    edits: Sequence[ReactionEdit],
    *,
    style: str = "unicode",
) -> Optional[ContextualTransformationLabel]:
    """Render one bond-order transformation with retained-neighbor context."""
    order_edits = [edit for edit in edits if edit.edit_type == "order_changed"]
    if len(order_edits) != 1 or any(
        edit.edit_type not in {"order_changed", "hydrogen_change"}
        for edit in edits
    ):
        return None
    changed = order_edits[0]
    if changed.atom_2 is None or changed.atom_1.component_index != changed.atom_2.component_index:
        return None
    component = next(
        (
            item
            for item in reactants
            if item.component_index == changed.atom_1.component_index
        ),
        None,
    )
    if component is None:
        return None
    mol = parse_smiles(component.input_smiles)
    if mol is None:
        return None
    left_key = _atom_key(changed, True)
    right_key = _atom_key(changed, False)
    endpoint_data = []
    for key, atom, opposite_index in (
        (left_key, changed.atom_1, right_key[1]),
        (right_key, changed.atom_2, left_key[1]),
    ):
        source_atom = mol.GetAtomWithIdx(key[1])
        before_h = int(source_atom.GetTotalNumHs(includeNeighbors=True))
        contexts = _external_contexts(
            mol, key[1], opposite_index, style=style
        )
        endpoint_data.append(
            {
                "key": key,
                "element": atom.element,
                "before_h": before_h,
                "after_h": before_h + _hydrogen_delta(edits, *key),
                "contexts": contexts,
            }
        )
    if any(item["after_h"] < 0 for item in endpoint_data):
        return None
    endpoint_data.sort(
        key=lambda item: (
            -len(item["contexts"]),
            item["contexts"],
            item["element"] != "C",
            item["key"],
        )
    )
    left, right = endpoint_data
    styling = load_rendering_taxonomy()["styles"][style]
    single = styling["bond"]
    before = (
        _endpoint_label(
            element=left["element"],
            hydrogens=left["before_h"],
            contexts=left["contexts"],
            side="left",
            bond=single,
        )
        + _bond_symbol(changed.old_order, style=style)
        + _endpoint_label(
            element=right["element"],
            hydrogens=right["before_h"],
            contexts=right["contexts"],
            side="right",
            bond=single,
        )
    )
    after = (
        _endpoint_label(
            element=left["element"],
            hydrogens=left["after_h"],
            contexts=left["contexts"],
            side="left",
            bond=single,
        )
        + _bond_symbol(changed.new_order, style=style)
        + _endpoint_label(
            element=right["element"],
            hydrogens=right["after_h"],
            contexts=right["contexts"],
            side="right",
            bond=single,
        )
    )
    arrow = "→" if style == "unicode" else "->"
    return ContextualTransformationLabel(
        before=before,
        after=after,
        concise=f"{before} {arrow} {after}",
        evidence=changed.evidence,
        confidence=min(edit.confidence for edit in edits),
    )


__all__ = [
    "ContextualTransformationLabel",
    "build_contextual_transformation_label",
]
