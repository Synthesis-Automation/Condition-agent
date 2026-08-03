"""Single terminal renderer for chemist-facing reaction labels.

All upstream contracts remain structured chemistry.  This module is the only
place where an entire reaction is collapsed into display text.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence

from .chemistry.rdkit_utils import parse_smiles
from .context import classify_context
from .labels import render_context
from .notation import format_chemist_text, notation_style
from .reaction_models import (
    PartialProductTransformation,
    ReactionComponent,
    ReactionEdit,
    ReactionObservation,
    RenderedReactionLabel,
)
from .reaction_render_context import ReactionRenderContext

_REACTION_LABEL_DEFINITION_FILES = (
    "chemist_notation.v1.json",
    "context_facets.v2.json",
    "rendering.v1.json",
)

_SITE_PRIORITY = {
    "leaving_group": 90,
    "transfer_group": 85,
    "pronucleophile_XH": 80,
    "nucleophile_anion": 75,
    "electrophilic_center": 70,
    "aromatic_CH": 65,
    "unsaturated_bond": 60,
    "heteroatom_bond": 55,
    "dipolar_group": 50,
    "addition_donor": 45,
    "eliminable_pair": 40,
}
_SITE_DISPLAY_ORDER = {
    "leaving_group": 0,
    "electrophilic_center": 1,
    "unsaturated_bond": 2,
    "dipolar_group": 3,
    "transfer_group": 4,
    "pronucleophile_XH": 5,
    "nucleophile_anion": 6,
    "aromatic_CH": 7,
    "heteroatom_bond": 8,
    "addition_donor": 9,
    "eliminable_pair": 10,
}


def reaction_label_definition_versions() -> dict[str, str]:
    """Return content-addressed versions of the complete label vocabulary."""
    root = Path(__file__).with_name("definitions")
    versions = {}
    for filename in _REACTION_LABEL_DEFINITION_FILES:
        raw = (root / filename).read_bytes()
        payload = json.loads(raw)
        version = str(payload.get("schema_version") or "unknown")
        versions[filename] = (
            f"{version}@sha256:{hashlib.sha256(raw).hexdigest()[:16]}"
        )
    return dict(sorted(versions.items()))


def _bond(order: str | None, style: str) -> str:
    glyphs = notation_style(style)
    return {
        "DOUBLE": glyphs["double"],
        "TRIPLE": glyphs["triple"],
        "AROMATIC": glyphs["aromatic"],
    }.get(str(order or "SINGLE").upper(), glyphs["single"])


def _participating_atoms(edits: Sequence[ReactionEdit]) -> dict[int, set[int]]:
    atoms: dict[int, set[int]] = {}
    for edit in edits:
        for endpoint in (edit.atom_1, edit.atom_2):
            if endpoint is None or endpoint.side != "reactant":
                continue
            atoms.setdefault(endpoint.component_index, set()).add(endpoint.atom_index)
    return atoms


def _site_terms(
    components: Sequence[ReactionComponent], edits: Sequence[ReactionEdit]
) -> tuple[str, ...]:
    """Select the most specific site label for each participating component."""
    participating = _participating_atoms(edits)
    terms: list[tuple[int, str, int]] = []
    for component in components:
        active = participating.get(component.component_index)
        if not active:
            continue
        candidates = []
        for site in component.molecule_analysis.reactive_site_hypotheses:
            overlap = len(active.intersection(site.atom_indices))
            if not overlap:
                continue
            candidates.append(
                (
                    -overlap,
                    -_SITE_PRIORITY.get(site.site_type, 0),
                    len(site.atom_indices),
                    site.hypothesis_id,
                    site.chemist_label,
                )
            )
        if candidates:
            selected = min(candidates)
            label = selected[-1]
            site_type = next(
                site.site_type
                for site in component.molecule_analysis.reactive_site_hypotheses
                if site.hypothesis_id == selected[-2]
            )
            order = _SITE_DISPLAY_ORDER.get(site_type, 99)
        else:
            label = _component_atom_term(component, active)
            order = 99
        terms.append((order, label, component.component_index))
    return tuple(label for _, label, _ in sorted(terms))


def _component_atom_term(component: ReactionComponent, atom_indices: Iterable[int]) -> str:
    molecule = parse_smiles(component.input_smiles)
    if molecule is None:
        return "?"
    symbols = [molecule.GetAtomWithIdx(index).GetSymbol() for index in sorted(atom_indices)]
    return "".join(symbols) or "?"


def _retained_heavy_degree(
    component: ReactionComponent, atom_index: int, opposite_index: int
) -> int:
    molecule = parse_smiles(component.input_smiles)
    if molecule is None:
        return 0
    return sum(
        neighbor.GetAtomicNum() > 1 and int(neighbor.GetIdx()) != opposite_index
        for neighbor in molecule.GetAtomWithIdx(atom_index).GetNeighbors()
    )


def _hydrogen_delta(edits: Sequence[ReactionEdit], component: int, atom: int) -> int:
    delta = 0
    for edit in edits:
        if edit.edit_type != "hydrogen_change":
            continue
        if (edit.atom_1.component_index, edit.atom_1.atom_index) != (component, atom):
            continue
        delta += 1 if edit.new_order and not edit.old_order else -1
    return delta


def _external_contexts(
    component: ReactionComponent,
    atom_index: int,
    excluded: set[int],
    *,
    style: str,
) -> tuple[str, ...]:
    molecule = parse_smiles(component.input_smiles)
    if molecule is None:
        return ()
    atom = molecule.GetAtomWithIdx(atom_index)
    labels = []
    for neighbor in atom.GetNeighbors():
        neighbor_index = int(neighbor.GetIdx())
        if neighbor_index in excluded or neighbor.GetAtomicNum() <= 1:
            continue
        context = classify_context(
            molecule,
            neighbor_index,
            excluded=tuple(sorted(excluded | {atom_index})),
        )
        labels.append(render_context(context.token, style=style))
    return tuple(sorted(labels))


def _local_endpoint(
    component: ReactionComponent,
    atom_index: int,
    excluded: set[int],
    edits: Sequence[ReactionEdit],
    *,
    product_side: bool,
    style: str,
) -> str:
    molecule = parse_smiles(component.input_smiles)
    if molecule is None:
        return "?"
    atom = molecule.GetAtomWithIdx(atom_index)
    element = atom.GetSymbol()
    if element == "C" and atom.GetIsAromatic():
        aromatic_context = classify_context(molecule, atom_index)
        return render_context(aromatic_context.token, style=style)
    carbonyl_oxygen_indices = {
        int(neighbor.GetIdx())
        for neighbor in atom.GetNeighbors()
        if element == "C"
        and neighbor.GetSymbol() == "O"
        and str(molecule.GetBondBetweenAtoms(atom_index, int(neighbor.GetIdx())).GetBondType()).upper()
        == "DOUBLE"
    }
    carbonyl_oxygen_indices -= excluded
    if product_side:
        changed_from_carbonyl = {
            opposite.atom_index
            for edit in edits
            if edit.edit_type == "order_changed"
            and str(edit.old_order or "").upper() == "DOUBLE"
            and str(edit.new_order or "").upper() != "DOUBLE"
            for endpoint, opposite in (
                (edit.atom_1, edit.atom_2),
                (edit.atom_2, edit.atom_1),
            )
            if endpoint is not None
            and opposite is not None
            and endpoint.component_index == component.component_index
            and endpoint.atom_index == atom_index
            and opposite.element == "O"
        }
        carbonyl_oxygen_indices -= changed_from_carbonyl
    contexts = _external_contexts(
        component,
        atom_index,
        excluded | carbonyl_oxygen_indices,
        style=style,
    )
    hydrogen_count = int(atom.GetTotalNumHs(includeNeighbors=True))
    if product_side:
        hydrogen_count += _hydrogen_delta(
            edits, component.component_index, atom_index
        )
    hydrogen = "" if hydrogen_count <= 0 else "H" if hydrogen_count == 1 else f"H{hydrogen_count}"
    center = f"{element}{hydrogen}"
    single = notation_style(style)["single"]
    if carbonyl_oxygen_indices:
        center = "C(=O)"
    if not contexts:
        return center
    if element in {"N", "O", "S", "P"}:
        return single.join((center, *contexts))
    return single.join((*contexts, center))


def _context_equation(
    components: Sequence[ReactionComponent],
    edits: Sequence[ReactionEdit],
    *,
    style: str,
) -> str | None:
    """Render a local before/after equation directly from edited atoms."""
    component_by_index = {item.component_index: item for item in components}
    bond_edits = [
        edit
        for edit in edits
        if edit.atom_2 is not None
        and edit.atom_1.side == "reactant"
        and edit.atom_2.side == "reactant"
    ]
    departed: dict[tuple[int, int], set[int]] = {}
    for edit in bond_edits:
        if edit.edit_type != "broken" or edit.atom_2 is None:
            continue
        if edit.atom_1.component_index == edit.atom_2.component_index:
            departed.setdefault(
                (edit.atom_1.component_index, edit.atom_1.atom_index), set()
            ).add(edit.atom_2.atom_index)
            departed.setdefault(
                (edit.atom_2.component_index, edit.atom_2.atom_index), set()
            ).add(edit.atom_1.atom_index)
    if not bond_edits:
        return None
    before_terms: list[str] = []
    after_terms: list[str] = []
    for edit in bond_edits:
        assert edit.atom_2 is not None
        left_atom, right_atom = edit.atom_1, edit.atom_2
        if (
            edit.edit_type == "formed"
            and left_atom.element != "C"
            and right_atom.element == "C"
        ):
            left_atom, right_atom = right_atom, left_atom
        elif (
            edit.edit_type == "order_changed"
            and left_atom.component_index == right_atom.component_index
        ):
            component = component_by_index.get(left_atom.component_index)
            if component is not None and _retained_heavy_degree(
                component, right_atom.atom_index, left_atom.atom_index
            ) > _retained_heavy_degree(
                component, left_atom.atom_index, right_atom.atom_index
            ):
                left_atom, right_atom = right_atom, left_atom
        first = component_by_index.get(left_atom.component_index)
        second = component_by_index.get(right_atom.component_index)
        if first is None or second is None:
            return None
        same_component_excluded_left = (
            {right_atom.atom_index} if first is second else set()
        )
        same_component_excluded_right = (
            {left_atom.atom_index} if first is second else set()
        )
        left_before = _local_endpoint(
            first,
            left_atom.atom_index,
            same_component_excluded_left,
            edits,
            product_side=False,
            style=style,
        )
        right_before = _local_endpoint(
            second,
            right_atom.atom_index,
            same_component_excluded_right,
            edits,
            product_side=False,
            style=style,
        )
        left_after = _local_endpoint(
            first,
            left_atom.atom_index,
            same_component_excluded_left
            | departed.get((left_atom.component_index, left_atom.atom_index), set()),
            edits,
            product_side=True,
            style=style,
        )
        right_after = _local_endpoint(
            second,
            right_atom.atom_index,
            same_component_excluded_right
            | departed.get((right_atom.component_index, right_atom.atom_index), set()),
            edits,
            product_side=True,
            style=style,
        )
        if edit.old_order:
            before_terms.append(left_before + _bond(edit.old_order, style) + right_before)
        if edit.new_order:
            after_terms.append(left_after + _bond(edit.new_order, style) + right_after)
    if not before_terms and not after_terms:
        return None
    plus = " + "
    arrow = notation_style(style)["arrow"]
    return f"{plus.join(before_terms) or '∅'} {arrow} {plus.join(after_terms) or '∅'}"


def _site_equation(
    components: Sequence[ReactionComponent],
    edits: Sequence[ReactionEdit],
    *,
    style: str,
) -> str | None:
    terms = _site_terms(components, edits)
    formed = [edit for edit in edits if edit.edit_type == "formed" and edit.atom_2]
    if not terms or not formed:
        return None
    local = _context_equation(components, edits, style=style)
    if local is None or notation_style(style)["arrow"] not in local:
        return None
    product = local.split(notation_style(style)["arrow"], 1)[1].strip()
    return f"{' + '.join(terms)} {notation_style(style)['arrow']} {product}"


def _partial_equation(
    transformation: PartialProductTransformation,
    *,
    style: str,
) -> str:
    single = notation_style(style)["single"]
    arrow = notation_style(style)["arrow"]
    old = transformation.removed_attachment.element
    installed = transformation.installed_fragment.canonical_fragment_smiles or transformation.added_attachment.element
    installed = installed.replace("#", notation_style(style)["triple"])
    if transformation.transformation_class == "acyl_heteroatom_substitution":
        before = f"R{single}C(=O){single}{old}"
        after = f"R{single}C(=O){single}{installed}"
    else:
        center = transformation.reactant_center.element
        before = f"{center}{single}{old}"
        after = f"{center}{single}{installed}"
    return f"{before} {arrow} {after}"


def _ring_equation(
    observation: ReactionObservation,
    reactants: Sequence[ReactionComponent],
    *,
    style: str,
) -> str | None:
    topology = observation.topology
    if topology is None or not topology.ring_changes:
        return None
    products = []
    for change in topology.ring_changes:
        counts = Counter(change.element_sequence)
        formula = "".join(
            element + (str(counts[element]) if counts[element] > 1 else "")
            for element in sorted(counts, key=lambda value: (value != "C", value))
        )
        aromaticity = "aromatic" if change.aromatic_after else "nonaromatic"
        products.append(f"{aromaticity} {change.ring_size}-membered {formula} ring")
    terms = _site_terms(reactants, observation.edits)
    left = " + ".join(terms) if terms else "reactive fragments"
    return f"{left} {notation_style(style)['arrow']} {' + '.join(products)}"


def render_reaction(
    context: ReactionRenderContext,
) -> RenderedReactionLabel:
    """Render the sole reaction label from completed upstream contracts."""
    observation = context.observation
    reactants = context.reactants
    signature = context.signature
    partial_transformation = context.partial_transformation
    style = context.style
    notation_style(style)
    warnings = tuple(sorted(set(observation.warnings)))
    versions = reaction_label_definition_versions()
    definition_version = ";".join(f"{key}={value}" for key, value in versions.items())
    conflict = observation.evidence_quality in {
        "conflicting_edit_evidence",
        "conflicting_stereochemical_evidence",
    }
    if conflict:
        return RenderedReactionLabel(
            text="Unavailable",
            status="conflicting_evidence",
            basis="unavailable",
            evidence=observation.evidence_quality,
            confidence=observation.evidence_confidence,
            warnings=warnings,
            style=style,
            definition_version=definition_version,
        )
    if partial_transformation is not None and not observation.edits:
        text = _partial_equation(partial_transformation, style=style)
        return RenderedReactionLabel(
            text=format_chemist_text(text, style=style),
            status="partial_product_correspondence",
            basis="partial_product_correspondence",
            evidence=partial_transformation.evidence,
            confidence=partial_transformation.confidence,
            warnings=tuple(sorted(set(warnings + partial_transformation.warnings))),
            style=style,
            definition_version=definition_version,
        )
    if not observation.edits:
        return RenderedReactionLabel(
            text="Unavailable",
            status="unavailable",
            basis="unavailable",
            evidence=observation.evidence_quality,
            confidence=observation.evidence_confidence,
            warnings=warnings,
            style=style,
            definition_version=definition_version,
        )
    event_count = len(signature.events) if signature is not None else (
        observation.core.event_count if observation.core is not None else 1
    )
    equation = _ring_equation(observation, reactants, style=style)
    basis = "ring_topology" if equation is not None else "reaction_sites"
    if equation is None:
        equation = _site_equation(reactants, observation.edits, style=style)
    if equation is None:
        equation = _context_equation(reactants, observation.edits, style=style)
        basis = "local_context"
    if equation is None:
        fragments = []
        for edit in observation.edits:
            endpoint = edit.atom_1.element
            if edit.atom_2 is not None:
                endpoint += _bond(edit.new_order or edit.old_order, style) + edit.atom_2.element
            fragments.append(f"{edit.edit_type}:{endpoint}")
        equation = "; ".join(fragments)
        basis = "literal_edits"
    status = "multi_event" if event_count > 1 else "structural_equation"
    if observation.topology and observation.topology.ring_changes:
        status = "ring_formation"
        basis = "ring_topology"
    return RenderedReactionLabel(
        text=format_chemist_text(re.sub(r"\s+", " ", equation).strip(), style=style),
        status=status,
        basis=basis,  # type: ignore[arg-type]
        evidence=observation.evidence_quality,
        confidence=observation.evidence_confidence,
        warnings=warnings,
        style=style,
        definition_version=definition_version,
        event_count=event_count,
    )


__all__ = ["reaction_label_definition_versions", "render_reaction"]
