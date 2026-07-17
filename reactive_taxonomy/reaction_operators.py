"""Generic graph-edit operators for v1 reaction grammars."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Tuple

from .chemistry.rdkit_utils import parse_smiles

from .reaction_models import BondChange, ReactionComponent, ReactionSiteReference


def _component_by_index(
    components: Tuple[ReactionComponent, ...], index: int
) -> ReactionComponent:
    return next(
        component for component in components if component.component_index == index
    )


def _fragment_to_remove(
    components: Tuple[ReactionComponent, ...],
    component_index: int,
    anchor_index: int,
    handle_index: int,
) -> set[int]:
    """Return the complete handle-side fragment after cutting its anchor bond."""
    from rdkit import Chem

    mol = parse_smiles(_component_by_index(components, component_index).input_smiles)
    if mol is None or mol.GetBondBetweenAtoms(anchor_index, handle_index) is None:
        return {handle_index}
    rw = Chem.RWMol(mol)
    rw.RemoveBond(anchor_index, handle_index)
    for fragment in Chem.GetMolFrags(rw.GetMol(), asMols=False, sanitizeFrags=False):
        atoms = set(int(index) for index in fragment)
        if handle_index in atoms and anchor_index not in atoms:
            return atoms
    return {handle_index}


def _bonded_role_atom(
    mol: Any,
    site: ReactionSiteReference,
    anchor_role: str,
    candidate_roles: Iterable[str],
) -> int:
    anchor = site.atom_roles[anchor_role][0]
    for role in candidate_roles:
        for atom_index in site.atom_roles.get(role, ()):
            if mol.GetBondBetweenAtoms(anchor, atom_index) is not None:
                return atom_index
    raise ValueError("No handle atom bonded to reactive anchor")


def _build_product(
    components: Tuple[ReactionComponent, ...],
    participants: List[ReactionSiteReference],
    removals: Dict[int, set[int]],
    join: Tuple[Tuple[int, int], Tuple[int, int]],
    *,
    add_hydrogen_on_left: bool = False,
    consume_hydrogen_on_right: bool = False,
) -> str | None:
    from rdkit import Chem

    used_indices = sorted({site.component_index for site in participants})
    mols = []
    offsets: Dict[int, int] = {}
    total = 0
    for component_index in used_indices:
        mol = parse_smiles(
            _component_by_index(components, component_index).input_smiles
        )
        if mol is None:
            return None
        offsets[component_index] = total
        total += mol.GetNumAtoms()
        mols.append(mol)
    combined = mols[0]
    for mol in mols[1:]:
        combined = Chem.CombineMols(combined, mol)
    remove_global = sorted(
        {offsets[ci] + ai for ci, atoms in removals.items() for ai in atoms},
        reverse=True,
    )
    join_global = [offsets[ci] + ai for ci, ai in join]
    left_hydrogen_count = int(combined.GetAtomWithIdx(join_global[0]).GetTotalNumHs())
    rw = Chem.RWMol(combined)
    for atom_index in remove_global:
        rw.RemoveAtom(atom_index)

    def shifted(index: int) -> int:
        return index - sum(removed < index for removed in remove_global)

    left, right = shifted(join_global[0]), shifted(join_global[1])
    if left == right or rw.GetBondBetweenAtoms(left, right) is not None:
        return None
    if add_hydrogen_on_left:
        product_atom = rw.GetAtomWithIdx(left)
        product_atom.SetNumExplicitHs(left_hydrogen_count + 1)
        product_atom.SetNoImplicit(True)
    if consume_hydrogen_on_right:
        partner_atom = rw.GetAtomWithIdx(right)
        explicit_h = int(partner_atom.GetNumExplicitHs())
        if explicit_h > 0:
            partner_atom.SetNumExplicitHs(explicit_h - 1)
    rw.AddBond(left, right, Chem.BondType.SINGLE)
    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        return Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def apply_operator(
    grammar: Dict[str, Any],
    assignment: Dict[str, ReactionSiteReference],
    components: Tuple[ReactionComponent, ...],
) -> Tuple[str | None, Tuple[BondChange, ...]]:
    operator = grammar["operator"]
    if operator["id"] == "join_two_anchors":
        left, right = assignment["electrophile"], assignment["transfer_partner"]
        left_mol = parse_smiles(
            _component_by_index(components, left.component_index).input_smiles
        )
        right_mol = parse_smiles(
            _component_by_index(components, right.component_index).input_smiles
        )
        if left_mol is None or right_mol is None:
            return None, ()
        left_handle = _bonded_role_atom(
            left_mol, left, "anchor", ("connector", "handle", "center")
        )
        right_handle = _bonded_role_atom(
            right_mol, right, "anchor", ("center", "handle")
        )
        removals: Dict[int, set[int]] = {}
        removals.setdefault(left.component_index, set()).update(
            _fragment_to_remove(
                components,
                left.component_index,
                left.atom_roles["anchor"][0],
                left_handle,
            )
        )
        removals.setdefault(right.component_index, set()).update(
            _fragment_to_remove(
                components,
                right.component_index,
                right.atom_roles["anchor"][0],
                right_handle,
            )
        )
        predicted = _build_product(
            components,
            [left, right],
            removals,
            (
                (left.component_index, left.atom_roles["anchor"][0]),
                (right.component_index, right.atom_roles["anchor"][0]),
            ),
        )
        changes = (
            BondChange(
                "broken",
                "electrophile.anchor",
                "electrophile.handle",
                "SINGLE",
                None,
                "grammar_operator",
            ),
            BondChange(
                "broken",
                "transfer_partner.anchor",
                "transfer_partner.center",
                "SINGLE",
                None,
                "grammar_operator",
            ),
            BondChange(
                "formed",
                "electrophile.anchor",
                "transfer_partner.anchor",
                None,
                "SINGLE",
                "grammar_operator",
            ),
        )
        return predicted, changes
    if operator["id"] == "replace_handle_with_alkene_endpoint":
        electrophile_role = str(operator["electrophile_role"])
        alkene_role = str(operator["alkene_role"])
        electrophile = assignment[electrophile_role]
        alkene = assignment[alkene_role]
        electrophile_mol = parse_smiles(
            _component_by_index(
                components, electrophile.component_index
            ).input_smiles
        )
        if electrophile_mol is None:
            return None, ()
        leaving_atom = _bonded_role_atom(
            electrophile_mol,
            electrophile,
            "anchor",
            ("connector", "handle", "center"),
        )
        endpoint_roles = ("endpoint_a", "endpoint_b")
        endpoint_h_counts = tuple(
            int(value) for value in alkene.details.get("endpoint_h_counts") or ()
        )
        if len(endpoint_h_counts) != 2 or max(endpoint_h_counts) < 1:
            return None, ()
        endpoint_role = min(
            endpoint_roles,
            key=lambda role: (
                -endpoint_h_counts[endpoint_roles.index(role)],
                int(alkene.atom_roles[role][0]),
            ),
        )
        predicted = _build_product(
            components,
            [electrophile, alkene],
            {
                electrophile.component_index: _fragment_to_remove(
                    components,
                    electrophile.component_index,
                    int(electrophile.atom_roles["anchor"][0]),
                    leaving_atom,
                )
            },
            (
                (
                    electrophile.component_index,
                    int(electrophile.atom_roles["anchor"][0]),
                ),
                (alkene.component_index, int(alkene.atom_roles[endpoint_role][0])),
            ),
            consume_hydrogen_on_right=True,
        )
        changes = (
            BondChange(
                "broken",
                f"{electrophile_role}.anchor",
                f"{electrophile_role}.handle",
                "SINGLE",
                None,
                "grammar_operator",
            ),
            BondChange(
                "formed",
                f"{electrophile_role}.anchor",
                f"{alkene_role}.{endpoint_role}",
                None,
                "SINGLE",
                "grammar_operator",
            ),
            BondChange(
                "hydrogen_change",
                f"{alkene_role}.{endpoint_role}",
                None,
                "SINGLE",
                None,
                "grammar_operator",
            ),
        )
        return predicted, changes
    if operator["id"] == "replace_handle_with_center":
        e_role, p_role = operator["electrophile_role"], operator["partner_role"]
        electrophile, partner = assignment[e_role], assignment[p_role]
        e_anchor_role = "anchor" if "anchor" in electrophile.atom_roles else "center"
        leaving_role = (
            "handle"
            if "handle" in electrophile.atom_roles
            else "leaving_or_activatable"
        )
        electrophile_mol = parse_smiles(
            _component_by_index(components, electrophile.component_index).input_smiles
        )
        if electrophile_mol is None:
            return None, ()
        leaving_atom = _bonded_role_atom(
            electrophile_mol,
            electrophile,
            e_anchor_role,
            ("connector", leaving_role, "center"),
        )
        removals = {
            electrophile.component_index: _fragment_to_remove(
                components,
                electrophile.component_index,
                electrophile.atom_roles[e_anchor_role][0],
                leaving_atom,
            )
        }
        predicted = _build_product(
            components,
            [electrophile, partner],
            removals,
            (
                (
                    electrophile.component_index,
                    electrophile.atom_roles[e_anchor_role][0],
                ),
                (partner.component_index, partner.atom_roles["center"][0]),
            ),
            consume_hydrogen_on_right=True,
        )
        changes = (
            BondChange(
                "broken",
                f"{e_role}.{e_anchor_role}",
                f"{e_role}.{leaving_role}",
                "SINGLE",
                None,
                "grammar_operator",
            ),
            BondChange(
                "formed",
                f"{e_role}.{e_anchor_role}",
                f"{p_role}.center",
                None,
                "SINGLE",
                "grammar_operator",
            ),
            BondChange(
                "hydrogen_change",
                f"{p_role}.center",
                None,
                "SINGLE",
                None,
                "grammar_operator",
            ),
        )
        return predicted, changes
    if operator["id"] == "replace_carbonyl_oxygen_with_amine":
        carbonyl_role = str(operator["carbonyl_role"])
        amine_role = str(operator["amine_role"])
        carbonyl = assignment[carbonyl_role]
        amine = assignment[amine_role]
        carbonyl_center = int(carbonyl.atom_roles["center"][0])
        oxygen = int(carbonyl.atom_roles["heteroatom"][0])
        nitrogen = int(amine.atom_roles["center"][0])
        predicted = _build_product(
            components,
            [carbonyl, amine],
            {carbonyl.component_index: {oxygen}},
            (
                (carbonyl.component_index, carbonyl_center),
                (amine.component_index, nitrogen),
            ),
            add_hydrogen_on_left=True,
            consume_hydrogen_on_right=True,
        )
        changes = (
            BondChange(
                "broken",
                f"{carbonyl_role}.center",
                f"{carbonyl_role}.heteroatom",
                "DOUBLE",
                None,
                "grammar_operator",
            ),
            BondChange(
                "formed",
                f"{carbonyl_role}.center",
                f"{amine_role}.center",
                None,
                "SINGLE",
                "grammar_operator",
            ),
            BondChange(
                "hydrogen_change",
                f"{carbonyl_role}.center",
                None,
                None,
                "SINGLE",
                "grammar_operator",
            ),
            BondChange(
                "hydrogen_change",
                f"{amine_role}.center",
                None,
                "SINGLE",
                None,
                "grammar_operator",
            ),
        )
        return predicted, changes
    return None, ()


__all__ = ["apply_operator"]
