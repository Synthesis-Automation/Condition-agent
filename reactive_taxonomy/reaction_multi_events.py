"""Exact reconstruction of balanced, unmapped multi-event reactions."""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, Iterable, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .connectivity_rewrite import connectivity_rewrite_for_grammar
from .reaction_models import ReactionComponent, ReactionSiteReference
from .reaction_site_interfaces import normalize_reaction_assignment


RawCandidate = Tuple[Dict[str, Any], Dict[str, ReactionSiteReference]]


def _event_sites(
    candidate: RawCandidate,
    reactants: Tuple[ReactionComponent, ...],
) -> tuple[
    ReactionSiteReference,
    ReactionSiteReference,
    int,
    int,
    int,
] | None:
    grammar, assignment = candidate
    rewrite = connectivity_rewrite_for_grammar(str(grammar.get("id") or ""))
    if rewrite is None:
        return None
    bindings = rewrite.grammar_role_bindings.get(str(grammar["id"]))
    if not bindings or {
        "leaving_source",
        "joining_partner",
    } - set(bindings):
        return None
    source = assignment.get(bindings["leaving_source"])
    partner = assignment.get(bindings["joining_partner"])
    if source is None or partner is None:
        return None
    normalized = normalize_reaction_assignment(
        {"source": source, "partner": partner},
        reactants,
    )
    if (
        len(normalized["source"].reactive_links) != 1
        or len(normalized["partner"].connection_endpoints) != 1
    ):
        return None
    link = normalized["source"].reactive_links[0]
    endpoint = normalized["partner"].connection_endpoints[0].endpoint
    if (
        link.source_kind != "explicit_bond"
        or link.endpoint_a.atom_index is None
        or link.endpoint_b.atom_index is None
        or endpoint.atom_index is None
    ):
        return None
    return (
        source,
        partner,
        link.endpoint_a.atom_index,
        link.endpoint_b.atom_index,
        endpoint.atom_index,
    )


def _operation_key(candidate: RawCandidate) -> Tuple[Any, ...]:
    grammar, assignment = candidate
    rewrite = connectivity_rewrite_for_grammar(str(grammar.get("id") or ""))
    bindings = (
        rewrite.grammar_role_bindings.get(str(grammar["id"]))
        if rewrite is not None
        else None
    )
    if not bindings:
        return (str(grammar.get("id") or ""),)
    electrophile = assignment[bindings["leaving_source"]]
    partner = assignment[bindings["joining_partner"]]
    anchor_role = (
        "anchor" if "anchor" in electrophile.atom_roles else "center"
    )
    return (
        grammar["id"],
        electrophile.component_index,
        int(electrophile.atom_roles[anchor_role][0]),
        partner.component_index,
        int(partner.atom_roles["center"][0]),
    )


def _component_by_index(
    components: Tuple[ReactionComponent, ...],
    component_index: int,
) -> ReactionComponent:
    return next(
        component
        for component in components
        if component.component_index == component_index
    )


def _fragment_to_remove(
    components: Tuple[ReactionComponent, ...],
    component_index: int,
    retained_atom: int,
    discarded_atom: int,
) -> set[int]:
    from rdkit import Chem

    molecule = parse_smiles(
        _component_by_index(components, component_index).input_smiles
    )
    if molecule is None:
        return set()
    rw = Chem.RWMol(molecule)
    rw.RemoveBond(retained_atom, discarded_atom)
    return next(
        (
            set(int(index) for index in fragment)
            for fragment in Chem.GetMolFrags(
                rw.GetMol(),
                asMols=False,
                sanitizeFrags=False,
            )
            if discarded_atom in fragment
        ),
        set(),
    )


def apply_rewrite_sequence(
    operations: Sequence[RawCandidate],
    components: Tuple[ReactionComponent, ...],
) -> str | None:
    """Apply compatible release-and-connect programs as one edit graph."""
    from rdkit import Chem

    if len(operations) < 2:
        return None
    resolved = [
        _event_sites(candidate, components) for candidate in operations
    ]
    if any(event is None for event in resolved):
        return None
    events = [event for event in resolved if event is not None]
    participants: dict[int, ReactionSiteReference] = {}
    removals: Dict[int, set[int]] = {}
    joins: list[Tuple[Tuple[int, int], Tuple[int, int]]] = []
    used_sources: set[Tuple[int, int]] = set()
    used_partners: set[Tuple[int, int]] = set()
    for source, partner, retained, discarded, joining in events:
        source_key = (source.component_index, retained)
        partner_key = (partner.component_index, joining)
        if source_key in used_sources or partner_key in used_partners:
            return None
        used_sources.add(source_key)
        used_partners.add(partner_key)
        fragment = _fragment_to_remove(
            components,
            source.component_index,
            retained,
            discarded,
        )
        if not fragment:
            return None
        removals.setdefault(source.component_index, set()).update(fragment)
        joins.append((source_key, partner_key))
        participants[source.component_index] = source
        participants[partner.component_index] = partner

    used_indices = sorted(participants)
    molecules = []
    offsets: Dict[int, int] = {}
    total = 0
    for component_index in used_indices:
        molecule = parse_smiles(
            _component_by_index(
                components, component_index
            ).input_smiles
        )
        if molecule is None:
            return None
        offsets[component_index] = total
        total += molecule.GetNumAtoms()
        molecules.append(molecule)
    combined = molecules[0]
    for molecule in molecules[1:]:
        combined = Chem.CombineMols(combined, molecule)
    removed_global = {
        offsets[component_index] + atom_index
        for component_index, atom_indices in removals.items()
        for atom_index in atom_indices
    }
    global_joins = [
        (
            offsets[left_component] + left_atom,
            offsets[right_component] + right_atom,
        )
        for (left_component, left_atom), (
            right_component,
            right_atom,
        ) in joins
    ]
    if any(
        left in removed_global or right in removed_global
        for left, right in global_joins
    ):
        return None
    rw = Chem.RWMol(combined)
    for atom_index in sorted(removed_global, reverse=True):
        rw.RemoveAtom(atom_index)

    removed_ascending = tuple(sorted(removed_global))

    def shifted(index: int) -> int:
        return index - sum(
            removed < index for removed in removed_ascending
        )

    for left_global, right_global in global_joins:
        left, right = shifted(left_global), shifted(right_global)
        if (
            left == right
            or rw.GetBondBetweenAtoms(left, right) is not None
        ):
            return None
        rw.AddBond(left, right, Chem.BondType.SINGLE)
    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        return Chem.MolToSmiles(
            product, canonical=True, isomericSmiles=True
        )
    except Exception:
        return None


def _interpretation_key(candidates: Sequence[RawCandidate]) -> Tuple[Any, ...]:
    return tuple(
        sorted(
            (
                grammar["id"],
                tuple(
                    sorted(
                        site.canonical_signature for site in assignment.values()
                    )
                ),
            )
            for grammar, assignment in candidates
        )
    )


def exact_multi_event_reconstructions(
    raw_candidates: Iterable[RawCandidate],
    reactants: Tuple[ReactionComponent, ...],
    observed_products: set[str],
    *,
    max_events: int = 4,
    max_combinations: int = 5000,
) -> Tuple[Tuple[RawCandidate, ...], ...]:
    """Return composite rewrite sets that exactly reconstruct one product."""
    eligible = tuple(
        sorted(
            (
                candidate
                for candidate in raw_candidates
                if _event_sites(candidate, reactants) is not None
            ),
            key=_operation_key,
        )
    )
    if len(eligible) < 2 or not observed_products:
        return ()
    exact = []
    attempted = 0
    upper = min(max_events, len(eligible))
    for event_count in range(2, upper + 1):
        for selected in combinations(eligible, event_count):
            attempted += 1
            if attempted > max_combinations:
                return tuple(exact)
            predicted = apply_rewrite_sequence(selected, reactants)
            if predicted in observed_products:
                exact.append(tuple(selected))
        if exact:
            break
    if not exact:
        return ()
    return tuple(
        sorted(
            exact,
            key=lambda selected: (
                _interpretation_key(selected),
                tuple(_operation_key(candidate) for candidate in selected),
            ),
        )
    )


def equivalent_multi_event_interpretations(
    reconstructions: Sequence[Sequence[RawCandidate]],
) -> bool:
    """Return whether exact alternatives differ only by equivalent sites."""
    return len({_interpretation_key(selected) for selected in reconstructions}) <= 1


__all__ = [
    "RawCandidate",
    "apply_rewrite_sequence",
    "equivalent_multi_event_interpretations",
    "exact_multi_event_reconstructions",
]
