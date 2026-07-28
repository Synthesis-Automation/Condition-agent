"""Generic graph-edit operators for v1 reaction grammars."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_graph_editing import (
    bond_type as _bond_type,
    capture_join_stereochemistry as _captured_join_stereochemistry,
    restore_join_stereochemistry as _restore_join_stereochemistry,
    set_total_hydrogens as _set_total_hydrogens,
)

from .reaction_models import (
    BondChange,
    OperatorOutcome,
    ReactionComponent,
    ReactionSiteReference,
)


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


def _change_site_bond_order(
    components: Tuple[ReactionComponent, ...],
    site: ReactionSiteReference,
    *,
    atom_role_1: str,
    atom_role_2: str,
    old_order: str,
    new_order: str,
    hydrogen_changes: Sequence[Dict[str, Any]] = (),
) -> str | None:
    """Apply one atom-provenanced bond-order edit to a reactant component."""
    from rdkit import Chem

    component = _component_by_index(components, site.component_index)
    molecule = parse_smiles(component.input_smiles)
    if molecule is None:
        return None
    indices_1 = site.atom_roles.get(atom_role_1) or ()
    indices_2 = site.atom_roles.get(atom_role_2) or ()
    if len(indices_1) != 1 or len(indices_2) != 1:
        return None
    atom_1, atom_2 = int(indices_1[0]), int(indices_2[0])
    bond = molecule.GetBondBetweenAtoms(atom_1, atom_2)
    expected = _bond_type(old_order)
    replacement = _bond_type(new_order)
    if bond is None or bond.GetBondType() != expected:
        return None
    rw = Chem.RWMol(molecule)
    edited_bond = rw.GetBondBetweenAtoms(atom_1, atom_2)
    if edited_bond is None:
        return None
    edited_bond.SetBondType(replacement)
    edited_bond.SetStereo(Chem.BondStereo.STEREONONE)
    edited_bond.SetBondDir(Chem.BondDir.NONE)
    if replacement in {Chem.BondType.DOUBLE, Chem.BondType.TRIPLE}:
        rw.GetAtomWithIdx(atom_1).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        rw.GetAtomWithIdx(atom_2).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    hydrogen_deltas: Dict[int, int] = {}
    for change in hydrogen_changes:
        indices = site.atom_roles.get(str(change["atom_role"])) or ()
        if len(indices) != 1:
            return None
        atom_index = int(indices[0])
        direction = str(change["direction"])
        hydrogen_deltas[atom_index] = hydrogen_deltas.get(atom_index, 0) + (
            1 if direction == "added" else -1
        )
    for atom_index, delta in hydrogen_deltas.items():
        source_atom = molecule.GetAtomWithIdx(atom_index)
        target_count = int(
            source_atom.GetTotalNumHs(includeNeighbors=True)
        ) + delta
        if target_count < 0:
            return None
        product_atom = rw.GetAtomWithIdx(atom_index)
        product_atom.SetNumExplicitHs(target_count)
        product_atom.SetNoImplicit(True)
    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        Chem.AssignStereochemistry(product, cleanIt=True, force=True)
        return Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def _combined_participants(
    components: Tuple[ReactionComponent, ...],
    participants: Sequence[ReactionSiteReference],
) -> tuple[Any, Dict[int, int]] | None:
    """Return one editable molecule and component offsets for retained participants."""
    from rdkit import Chem

    component_indices = sorted({site.component_index for site in participants})
    molecules = []
    offsets: Dict[int, int] = {}
    total = 0
    for component_index in component_indices:
        molecule = parse_smiles(
            _component_by_index(components, component_index).input_smiles
        )
        if molecule is None:
            return None
        offsets[component_index] = total
        total += molecule.GetNumAtoms()
        molecules.append(molecule)
    if not molecules:
        return None
    combined = molecules[0]
    for molecule in molecules[1:]:
        combined = Chem.CombineMols(combined, molecule)
    return combined, offsets


def _pair_addition_outcomes(
    grammar: Dict[str, Any],
    assignment: Dict[str, ReactionSiteReference],
    components: Tuple[ReactionComponent, ...],
) -> Tuple[OperatorOutcome, ...]:
    """Enumerate constitutional A-B orientations across one unsaturated bond."""
    from rdkit import Chem

    operator = grammar["operator"]
    acceptor_role = str(operator.get("acceptor_role") or "acceptor")
    donor_role = str(operator.get("donor_role") or "donor")
    acceptor = assignment[acceptor_role]
    donor = assignment[donor_role]
    endpoint_roles = tuple(
        str(value)
        for value in operator.get(
            "acceptor_endpoint_roles", ("endpoint_a", "endpoint_b")
        )
    )
    if len(endpoint_roles) != 2:
        return ()
    endpoint_a_values = acceptor.atom_roles.get(endpoint_roles[0]) or ()
    endpoint_b_values = acceptor.atom_roles.get(endpoint_roles[1]) or ()
    if len(endpoint_a_values) != 1 or len(endpoint_b_values) != 1:
        return ()
    source_kind = (
        "implicit_hydrogen"
        if donor.site_type == "pronucleophile_XH"
        else str(donor.details.get("source_kind") or "")
    )
    addend_a_role = (
        "center"
        if donor.site_type == "pronucleophile_XH"
        else str(operator.get("donor_addend_a_role") or "addend_a")
    )
    hydrogen_carrier_role = (
        "center"
        if donor.site_type == "pronucleophile_XH"
        else str(operator.get("hydrogen_carrier_role") or "hydrogen_carrier")
    )
    addend_a_values = donor.atom_roles.get(addend_a_role) or ()
    if len(addend_a_values) != 1:
        return ()
    addend_b_role = str(operator.get("donor_addend_b_role") or "addend_b")
    addend_b_values = donor.atom_roles.get(addend_b_role) or ()
    if source_kind == "explicit_bond" and len(addend_b_values) != 1:
        return ()
    if source_kind not in {"explicit_bond", "implicit_hydrogen"}:
        return ()

    endpoint_assignments = (
        (
            endpoint_roles[0],
            endpoint_roles[1],
            int(endpoint_a_values[0]),
            int(endpoint_b_values[0]),
        ),
        (
            endpoint_roles[1],
            endpoint_roles[0],
            int(endpoint_b_values[0]),
            int(endpoint_a_values[0]),
        ),
    )
    outcomes: list[OperatorOutcome] = []
    seen_products: set[str] = set()
    old_order = str(operator["old_order"]).upper()
    new_order = str(operator["new_order"]).upper()
    for addend_endpoint_role, other_endpoint_role, addend_endpoint, other_endpoint in (
        endpoint_assignments
    ):
        combined_result = _combined_participants(
            components, (acceptor, donor)
        )
        if combined_result is None:
            continue
        combined, offsets = combined_result
        rw = Chem.RWMol(combined)
        acceptor_offset = offsets[acceptor.component_index]
        donor_offset = offsets[donor.component_index]
        endpoint_global = acceptor_offset + addend_endpoint
        other_endpoint_global = acceptor_offset + other_endpoint
        addend_a_global = donor_offset + int(addend_a_values[0])
        acceptor_bond = rw.GetBondBetweenAtoms(
            acceptor_offset + int(endpoint_a_values[0]),
            acceptor_offset + int(endpoint_b_values[0]),
        )
        if acceptor_bond is None or acceptor_bond.GetBondType() != _bond_type(
            old_order
        ):
            continue
        acceptor_bond.SetBondType(_bond_type(new_order))
        acceptor_bond.SetStereo(Chem.BondStereo.STEREONONE)
        acceptor_bond.SetBondDir(Chem.BondDir.NONE)

        changes: Tuple[BondChange, ...] = (
            BondChange(
                "order_changed",
                f"{acceptor_role}.{endpoint_roles[0]}",
                f"{acceptor_role}.{endpoint_roles[1]}",
                old_order,
                new_order,
                "grammar_operator",
            ),
        )
        if source_kind == "explicit_bond":
            addend_b_global = donor_offset + int(addend_b_values[0])
            source_bond = rw.GetBondBetweenAtoms(addend_a_global, addend_b_global)
            source_order = str(
                donor.details.get("source_bond_order") or "SINGLE"
            ).upper()
            if source_bond is None or source_bond.GetBondType() != _bond_type(
                source_order
            ):
                continue
            rw.RemoveBond(addend_a_global, addend_b_global)
            if (
                rw.GetBondBetweenAtoms(endpoint_global, addend_a_global) is not None
                or rw.GetBondBetweenAtoms(other_endpoint_global, addend_b_global)
                is not None
            ):
                continue
            rw.AddBond(endpoint_global, addend_a_global, Chem.BondType.SINGLE)
            rw.AddBond(other_endpoint_global, addend_b_global, Chem.BondType.SINGLE)
            changes += (
                BondChange(
                    "broken",
                    f"{donor_role}.{addend_a_role}",
                    f"{donor_role}.{addend_b_role}",
                    source_order,
                    None,
                    "grammar_operator",
                ),
                BondChange(
                    "formed",
                    f"{acceptor_role}.{addend_endpoint_role}",
                    f"{donor_role}.{addend_a_role}",
                    None,
                    "SINGLE",
                    "grammar_operator",
                ),
                BondChange(
                    "formed",
                    f"{acceptor_role}.{other_endpoint_role}",
                    f"{donor_role}.{addend_b_role}",
                    None,
                    "SINGLE",
                    "grammar_operator",
                ),
            )
        else:
            carrier_values = donor.atom_roles.get(hydrogen_carrier_role) or ()
            if len(carrier_values) != 1:
                continue
            if rw.GetBondBetweenAtoms(endpoint_global, addend_a_global) is not None:
                continue
            rw.AddBond(endpoint_global, addend_a_global, Chem.BondType.SINGLE)
            if not _set_total_hydrogens(
                combined, rw, donor_offset + int(carrier_values[0]), -1
            ):
                continue
            if not _set_total_hydrogens(combined, rw, other_endpoint_global, 1):
                continue
            changes += (
                BondChange(
                    "formed",
                    f"{acceptor_role}.{addend_endpoint_role}",
                    f"{donor_role}.{addend_a_role}",
                    None,
                    "SINGLE",
                    "grammar_operator",
                ),
                BondChange(
                    "hydrogen_change",
                    f"{donor_role}.{hydrogen_carrier_role}",
                    None,
                    "SINGLE",
                    None,
                    "grammar_operator",
                ),
                BondChange(
                    "hydrogen_change",
                    f"{acceptor_role}.{other_endpoint_role}",
                    None,
                    None,
                    "SINGLE",
                    "grammar_operator",
                ),
            )
        product = rw.GetMol()
        try:
            product.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(product)
            Chem.AssignStereochemistry(product, cleanIt=True, force=True)
            smiles = Chem.MolToSmiles(
                product, canonical=True, isomericSmiles=True
            )
        except Exception:
            continue
        if smiles in seen_products:
            continue
        seen_products.add(smiles)
        outcomes.append(
            OperatorOutcome(
                outcome_id=(
                    f"{addend_endpoint_role}_addend_a__"
                    f"{other_endpoint_role}_addend_b"
                ),
                predicted_product_smiles=smiles,
                predicted_bond_changes=changes,
            )
        )
    return tuple(sorted(outcomes, key=lambda outcome: outcome.outcome_id))


def _pair_elimination_outcomes(
    grammar: Dict[str, Any],
    assignment: Dict[str, ReactionSiteReference],
    components: Tuple[ReactionComponent, ...],
) -> Tuple[OperatorOutcome, ...]:
    """Apply one conservative vicinal heavy-group/H elimination."""
    from rdkit import Chem

    operator = grammar["operator"]
    substrate_role = str(operator.get("substrate_role") or "substrate")
    substrate = assignment[substrate_role]
    endpoint_a_role = str(operator.get("endpoint_a_role") or "endpoint_a")
    endpoint_b_role = str(operator.get("endpoint_b_role") or "endpoint_b")
    departing_role = str(operator.get("departing_a_role") or "departing_a")
    hydrogen_role = str(
        operator.get("hydrogen_carrier_b_role") or "hydrogen_carrier_b"
    )
    role_values = {
        role: substrate.atom_roles.get(role) or ()
        for role in (
            endpoint_a_role,
            endpoint_b_role,
            departing_role,
            hydrogen_role,
        )
    }
    if any(len(values) != 1 for values in role_values.values()):
        return ()
    endpoint_a = int(role_values[endpoint_a_role][0])
    endpoint_b = int(role_values[endpoint_b_role][0])
    departing = int(role_values[departing_role][0])
    component = _component_by_index(components, substrate.component_index)
    molecule = parse_smiles(component.input_smiles)
    if molecule is None:
        return ()
    backbone = molecule.GetBondBetweenAtoms(endpoint_a, endpoint_b)
    leaving_bond = molecule.GetBondBetweenAtoms(endpoint_a, departing)
    old_order = str(operator["old_order"]).upper()
    new_order = str(operator["new_order"]).upper()
    if (
        backbone is None
        or backbone.GetBondType() != _bond_type(old_order)
        or leaving_bond is None
    ):
        return ()
    removals = sorted(
        _fragment_to_remove(
            components,
            substrate.component_index,
            endpoint_a,
            departing,
        ),
        reverse=True,
    )
    if endpoint_a in removals or endpoint_b in removals:
        return ()
    rw = Chem.RWMol(molecule)
    for atom_index in removals:
        rw.RemoveAtom(atom_index)

    def shifted(index: int) -> int:
        return index - sum(removed < index for removed in removals)

    shifted_a = shifted(endpoint_a)
    shifted_b = shifted(endpoint_b)
    edited_backbone = rw.GetBondBetweenAtoms(shifted_a, shifted_b)
    if edited_backbone is None:
        return ()
    edited_backbone.SetBondType(_bond_type(new_order))
    edited_backbone.SetStereo(Chem.BondStereo.STEREONONE)
    edited_backbone.SetBondDir(Chem.BondDir.NONE)
    source_hydrogen_count = int(
        molecule.GetAtomWithIdx(endpoint_b).GetTotalNumHs(includeNeighbors=True)
    )
    if source_hydrogen_count < 1:
        return ()
    product_hydrogen_atom = rw.GetAtomWithIdx(shifted_b)
    product_hydrogen_atom.SetNumExplicitHs(source_hydrogen_count - 1)
    product_hydrogen_atom.SetNoImplicit(True)
    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        Chem.AssignStereochemistry(product, cleanIt=True, force=True)
        smiles = Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
    except Exception:
        return ()
    leaving_order = {
        1: "SINGLE",
        2: "DOUBLE",
        3: "TRIPLE",
    }.get(int(round(float(leaving_bond.GetBondTypeAsDouble()))), "SINGLE")
    changes = (
        BondChange(
            "broken",
            f"{substrate_role}.{endpoint_a_role}",
            f"{substrate_role}.{departing_role}",
            leaving_order,
            None,
            "grammar_operator",
        ),
        BondChange(
            "hydrogen_change",
            f"{substrate_role}.{hydrogen_role}",
            None,
            "SINGLE",
            None,
            "grammar_operator",
        ),
        BondChange(
            "order_changed",
            f"{substrate_role}.{endpoint_a_role}",
            f"{substrate_role}.{endpoint_b_role}",
            old_order,
            new_order,
            "grammar_operator",
        ),
    )
    return (
        OperatorOutcome(
            outcome_id="vicinal_pair",
            predicted_product_smiles=smiles,
            predicted_bond_changes=changes,
        ),
    )


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
    if operator["id"] == "center_replacement":
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
    if operator["id"] == "change_bond_order":
        site_role = str(operator["site_role"])
        atom_role_1 = str(operator["atom_role_1"])
        atom_role_2 = str(operator["atom_role_2"])
        old_order = str(operator["old_order"]).upper()
        new_order = str(operator["new_order"]).upper()
        hydrogen_changes = tuple(operator.get("hydrogen_changes") or ())
        site = assignment[site_role]
        predicted = _change_site_bond_order(
            components,
            site,
            atom_role_1=atom_role_1,
            atom_role_2=atom_role_2,
            old_order=old_order,
            new_order=new_order,
            hydrogen_changes=hydrogen_changes,
        )
        changes = (
            BondChange(
                "order_changed",
                f"{site_role}.{atom_role_1}",
                f"{site_role}.{atom_role_2}",
                old_order,
                new_order,
                "grammar_operator",
            ),
        )
        for hydrogen_change in hydrogen_changes:
            direction = str(hydrogen_change["direction"])
            changes += (
                BondChange(
                    "hydrogen_change",
                    f"{site_role}.{hydrogen_change['atom_role']}",
                    None,
                    "SINGLE" if direction == "removed" else None,
                    "SINGLE" if direction == "added" else None,
                    "grammar_operator",
                ),
            )
        return predicted, changes
    return None, ()


def enumerate_operator_outcomes(
    grammar: Dict[str, Any],
    assignment: Dict[str, ReactionSiteReference],
    components: Tuple[ReactionComponent, ...],
) -> Tuple[OperatorOutcome, ...]:
    """Return every distinct constitutional outcome for one grammar assignment."""
    operator_id = str(grammar.get("operator", {}).get("id") or "")
    if operator_id == "pair_addition":
        return _pair_addition_outcomes(grammar, assignment, components)
    if operator_id == "pair_elimination":
        return _pair_elimination_outcomes(grammar, assignment, components)
    predicted, changes = apply_operator(grammar, assignment, components)
    return (
        OperatorOutcome(
            outcome_id="default",
            predicted_product_smiles=predicted,
            predicted_bond_changes=changes,
        ),
    )


def apply_operator_sequence(
    operations: Sequence[
        Tuple[Dict[str, Any], Dict[str, ReactionSiteReference]]
    ],
    components: Tuple[ReactionComponent, ...],
) -> str | None:
    """Apply distinct handle-replacement events as one composite graph edit."""
    from rdkit import Chem

    if len(operations) < 2 or any(
        grammar.get("operator", {}).get("id") != "center_replacement"
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


__all__ = [
    "apply_operator",
    "apply_operator_sequence",
    "enumerate_operator_outcomes",
]
