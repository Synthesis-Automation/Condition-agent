"""Generic graph-edit operators for v1 reaction grammars."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Tuple

from chemtools.core.rdkit import parse_smiles

from .reaction_models import BondChange, ReactionComponent, ReactionSiteReference


def _component_by_index(components: Tuple[ReactionComponent, ...], index: int) -> ReactionComponent:
    return next(component for component in components if component.component_index == index)


def _build_product(
    components: Tuple[ReactionComponent, ...],
    participants: List[ReactionSiteReference],
    removals: Dict[int, set[int]],
    join: Tuple[Tuple[int, int], Tuple[int, int]],
    *,
    consume_hydrogen_on_right: bool = False,
) -> str | None:
    from rdkit import Chem
    used_indices = sorted({site.component_index for site in participants})
    mols = []
    offsets: Dict[int, int] = {}
    total = 0
    for component_index in used_indices:
        mol = parse_smiles(_component_by_index(components, component_index).input_smiles)
        if mol is None: return None
        offsets[component_index] = total; total += mol.GetNumAtoms(); mols.append(mol)
    combined = mols[0]
    for mol in mols[1:]: combined = Chem.CombineMols(combined, mol)
    remove_global = sorted({offsets[ci] + ai for ci, atoms in removals.items() for ai in atoms}, reverse=True)
    join_global = [offsets[ci] + ai for ci, ai in join]
    rw = Chem.RWMol(combined)
    for atom_index in remove_global: rw.RemoveAtom(atom_index)
    def shifted(index: int) -> int:
        return index - sum(removed < index for removed in remove_global)
    left, right = shifted(join_global[0]), shifted(join_global[1])
    if left == right or rw.GetBondBetweenAtoms(left, right) is not None: return None
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
        removals = {
            left.component_index: set(left.atom_roles["handle"]),
            right.component_index: set(right.atom_roles["handle"]),
        }
        predicted = _build_product(components, [left, right], removals, ((left.component_index, left.atom_roles["anchor"][0]), (right.component_index, right.atom_roles["anchor"][0])))
        changes = (
            BondChange("broken", "electrophile.anchor", "electrophile.handle", "SINGLE", None, "grammar_operator"),
            BondChange("broken", "transfer_partner.anchor", "transfer_partner.center", "SINGLE", None, "grammar_operator"),
            BondChange("formed", "electrophile.anchor", "transfer_partner.anchor", None, "SINGLE", "grammar_operator"),
        )
        return predicted, changes
    if operator["id"] == "replace_handle_with_center":
        e_role, p_role = operator["electrophile_role"], operator["partner_role"]
        electrophile, partner = assignment[e_role], assignment[p_role]
        e_anchor_role = "anchor" if "anchor" in electrophile.atom_roles else "center"
        leaving_role = "handle" if "handle" in electrophile.atom_roles else "leaving_or_activatable"
        removals = {electrophile.component_index: set(electrophile.atom_roles[leaving_role])}
        predicted = _build_product(
            components, [electrophile, partner], removals,
            ((electrophile.component_index, electrophile.atom_roles[e_anchor_role][0]), (partner.component_index, partner.atom_roles["center"][0])),
            consume_hydrogen_on_right=True,
        )
        changes = (
            BondChange("broken", f"{e_role}.{e_anchor_role}", f"{e_role}.{leaving_role}", "SINGLE", None, "grammar_operator"),
            BondChange("formed", f"{e_role}.{e_anchor_role}", f"{p_role}.center", None, "SINGLE", "grammar_operator"),
        )
        return predicted, changes
    return None, ()


__all__ = ["apply_operator"]
