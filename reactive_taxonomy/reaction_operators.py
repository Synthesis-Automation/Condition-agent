"""Generic graph-edit operators for v1 reaction grammars."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Tuple

from .chemistry.rdkit_utils import parse_smiles

from .reaction_models import BondChange, ReactionComponent, ReactionSiteReference


def _component_by_index(components: Tuple[ReactionComponent, ...], index: int) -> ReactionComponent:
    return next(component for component in components if component.component_index == index)


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


def _bonded_role_atom(mol: Any, site: ReactionSiteReference, anchor_role: str, candidate_roles: Iterable[str]) -> int:
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
        left_mol = parse_smiles(_component_by_index(components, left.component_index).input_smiles)
        right_mol = parse_smiles(_component_by_index(components, right.component_index).input_smiles)
        if left_mol is None or right_mol is None: return None, ()
        left_handle = _bonded_role_atom(left_mol, left, "anchor", ("connector", "handle", "center"))
        right_handle = _bonded_role_atom(right_mol, right, "anchor", ("center", "handle"))
        removals = {
            left.component_index: _fragment_to_remove(components, left.component_index, left.atom_roles["anchor"][0], left_handle),
            right.component_index: _fragment_to_remove(components, right.component_index, right.atom_roles["anchor"][0], right_handle),
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
        electrophile_mol = parse_smiles(_component_by_index(components, electrophile.component_index).input_smiles)
        if electrophile_mol is None: return None, ()
        leaving_atom = _bonded_role_atom(electrophile_mol, electrophile, e_anchor_role, ("connector", leaving_role, "center"))
        removals = {electrophile.component_index: _fragment_to_remove(components, electrophile.component_index, electrophile.atom_roles[e_anchor_role][0], leaving_atom)}
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
