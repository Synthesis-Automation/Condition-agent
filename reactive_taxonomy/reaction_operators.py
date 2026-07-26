"""Generic graph-edit operators for v1 reaction grammars."""

from __future__ import annotations

from typing import Any, Callable, Dict, Iterable, List, Sequence, Tuple

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


def _captured_join_stereochemistry(
    molecule: Any,
    *,
    removed_indices: set[int],
    join_indices: Sequence[int],
) -> tuple[tuple[int, int, int, int, Any], ...]:
    """Capture defined alkene stereo whose leaving substituent is replaced."""
    from rdkit import Chem

    if len(join_indices) != 2:
        return ()
    replacement_pairs = (
        (int(join_indices[0]), int(join_indices[1])),
        (int(join_indices[1]), int(join_indices[0])),
    )
    captured = []
    for bond in molecule.GetBonds():
        stereo = bond.GetStereo()
        stereo_atoms = tuple(int(index) for index in bond.GetStereoAtoms())
        if (
            bond.GetBondType() != Chem.BondType.DOUBLE
            or stereo == Chem.BondStereo.STEREONONE
            or len(stereo_atoms) != 2
            or not any(index in removed_indices for index in stereo_atoms)
        ):
            continue
        endpoints = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        remapped_stereo_atoms = []
        valid = True
        for stereo_atom in stereo_atoms:
            if stereo_atom not in removed_indices:
                remapped_stereo_atoms.append(stereo_atom)
                continue
            replacements = [
                replacement
                for endpoint, replacement in replacement_pairs
                if endpoint in endpoints
                and molecule.GetBondBetweenAtoms(endpoint, stereo_atom) is not None
                and replacement not in removed_indices
            ]
            if len(replacements) != 1:
                valid = False
                break
            remapped_stereo_atoms.append(replacements[0])
        if valid and not any(index in removed_indices for index in endpoints):
            captured.append(
                (
                    int(endpoints[0]),
                    int(endpoints[1]),
                    int(remapped_stereo_atoms[0]),
                    int(remapped_stereo_atoms[1]),
                    stereo,
                )
            )
    return tuple(captured)


def _restore_join_stereochemistry(
    molecule: Any,
    *,
    captured: Sequence[tuple[int, int, int, int, Any]],
    shifted: Callable[[int], int],
) -> None:
    """Restore captured alkene stereo after atom deletion and bond formation."""
    from rdkit import Chem

    restored = False
    for endpoint_1, endpoint_2, stereo_atom_1, stereo_atom_2, stereo in captured:
        bond = molecule.GetBondBetweenAtoms(
            shifted(endpoint_1), shifted(endpoint_2)
        )
        if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        bond.SetStereoAtoms(
            shifted(stereo_atom_1),
            shifted(stereo_atom_2),
        )
        bond.SetStereo(stereo)
        restored = True
    if restored:
        Chem.SetDoubleBondNeighborDirections(molecule)
        Chem.AssignStereochemistry(molecule, cleanIt=False, force=True)


def _build_product(
    components: Tuple[ReactionComponent, ...],
    participants: List[ReactionSiteReference],
    removals: Dict[int, set[int]],
    join: Tuple[Tuple[int, int], Tuple[int, int]],
    *,
    add_hydrogen_on_left: bool = False,
    consume_hydrogen_on_right: bool = False,
    neutralize_negative_charge_on_right: bool = False,
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
    # A broadly matched candidate can place one join endpoint inside the
    # leaving fragment selected by the other site. Removing that fragment
    # makes the endpoint impossible to remap and can otherwise send an
    # out-of-range atom index into RDKit.
    if any(endpoint in remove_global for endpoint in join_global):
        return None
    captured_stereochemistry = _captured_join_stereochemistry(
        combined,
        removed_indices=set(remove_global),
        join_indices=join_global,
    )
    left_hydrogen_count = int(combined.GetAtomWithIdx(join_global[0]).GetTotalNumHs())
    rw = Chem.RWMol(combined)
    for atom_index in remove_global:
        rw.RemoveAtom(atom_index)

    def shifted(index: int) -> int:
        return index - sum(removed < index for removed in remove_global)

    left, right = shifted(join_global[0]), shifted(join_global[1])
    atom_count = rw.GetNumAtoms()
    if (
        not 0 <= left < atom_count
        or not 0 <= right < atom_count
        or left == right
        or rw.GetBondBetweenAtoms(left, right) is not None
    ):
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
    if neutralize_negative_charge_on_right:
        partner_atom = rw.GetAtomWithIdx(right)
        if partner_atom.GetFormalCharge() < 0:
            partner_atom.SetFormalCharge(partner_atom.GetFormalCharge() + 1)
    rw.AddBond(left, right, Chem.BondType.SINGLE)
    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        _restore_join_stereochemistry(
            product,
            captured=captured_stereochemistry,
            shifted=shifted,
        )
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
        partner_h_count = int(partner.details.get("h_count", 0))
        partner_charge = int(partner.details.get("formal_charge", 0))
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
            consume_hydrogen_on_right=partner_h_count > 0,
            neutralize_negative_charge_on_right=partner_charge < 0,
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
        )
        if partner_h_count > 0:
            changes += (
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


def apply_operator_sequence(
    operations: Sequence[
        Tuple[Dict[str, Any], Dict[str, ReactionSiteReference]]
    ],
    components: Tuple[ReactionComponent, ...],
) -> str | None:
    """Apply distinct handle-replacement events as one composite graph edit."""
    from rdkit import Chem

    if len(operations) < 2 or any(
        grammar.get("operator", {}).get("id") != "replace_handle_with_center"
        for grammar, _ in operations
    ):
        return None
    participants: dict[int, ReactionSiteReference] = {}
    removals: Dict[int, set[int]] = {}
    joins: list[Tuple[Tuple[int, int], Tuple[int, int]]] = []
    used_electrophiles: set[Tuple[int, int]] = set()
    used_partners: set[Tuple[int, int]] = set()
    for grammar, assignment in operations:
        operator = grammar["operator"]
        e_role = str(operator["electrophile_role"])
        p_role = str(operator["partner_role"])
        electrophile = assignment[e_role]
        partner = assignment[p_role]
        e_anchor_role = (
            "anchor" if "anchor" in electrophile.atom_roles else "center"
        )
        leaving_role = (
            "handle"
            if "handle" in electrophile.atom_roles
            else "leaving_or_activatable"
        )
        anchor_index = int(electrophile.atom_roles[e_anchor_role][0])
        partner_index = int(partner.atom_roles["center"][0])
        electrophile_key = (electrophile.component_index, anchor_index)
        partner_key = (partner.component_index, partner_index)
        if electrophile_key in used_electrophiles or partner_key in used_partners:
            return None
        used_electrophiles.add(electrophile_key)
        used_partners.add(partner_key)
        electrophile_mol = parse_smiles(
            _component_by_index(
                components, electrophile.component_index
            ).input_smiles
        )
        if electrophile_mol is None:
            return None
        leaving_atom = _bonded_role_atom(
            electrophile_mol,
            electrophile,
            e_anchor_role,
            ("connector", leaving_role, "center"),
        )
        removals.setdefault(electrophile.component_index, set()).update(
            _fragment_to_remove(
                components,
                electrophile.component_index,
                anchor_index,
                leaving_atom,
            )
        )
        joins.append((electrophile_key, partner_key))
        participants[electrophile.component_index] = electrophile
        participants[partner.component_index] = partner

    used_indices = sorted(participants)
    molecules = []
    offsets: Dict[int, int] = {}
    total = 0
    for component_index in used_indices:
        molecule = parse_smiles(
            _component_by_index(components, component_index).input_smiles
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
    remove_global = sorted(removed_global, reverse=True)
    global_joins = [
        (
            offsets[left_component] + left_atom,
            offsets[right_component] + right_atom,
        )
        for (left_component, left_atom), (right_component, right_atom) in joins
    ]
    # Candidate events can overlap even when their anchors are individually
    # distinct.  In particular, one event's leaving fragment may contain the
    # partner center (or anchor) used by another event.  Such a sequence is not
    # a coherent simultaneous graph edit and its atom indices cannot be
    # remapped after removal.
    if any(
        left in removed_global or right in removed_global
        for left, right in global_joins
    ):
        return None
    rw = Chem.RWMol(combined)
    for atom_index in remove_global:
        rw.RemoveAtom(atom_index)

    def shifted(index: int) -> int:
        return index - sum(removed < index for removed in remove_global)

    for left_global, right_global in global_joins:
        left = shifted(left_global)
        right = shifted(right_global)
        atom_count = rw.GetNumAtoms()
        if (
            not 0 <= left < atom_count
            or not 0 <= right < atom_count
            or left == right
            or rw.GetBondBetweenAtoms(left, right) is not None
        ):
            return None
        rw.AddBond(left, right, Chem.BondType.SINGLE)
    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        return Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
    except Exception:
        return None


__all__ = ["apply_operator", "apply_operator_sequence"]
