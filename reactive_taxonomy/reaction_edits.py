"""Normalize mapped and inferred graph edits into typed contracts."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, replace
from typing import Any, Dict, Optional, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import (
    AtomStateTransition,
    ConnectivityEditGraph,
    HydrogenDelta,
    ReactionAtomReference,
    ReactionComponent,
    ReactionEdit,
    ReactionEditHypothesis,
    ReactionEvidenceCandidate,
    ReactionStereoChange,
)
from .reaction_connectivity import (
    atom_provenance_key,
    bond_state,
    build_connectivity_edit_graph,
    connectivity_graph_from_reaction_edits,
    endpoint_absent_state,
    make_bond_transition,
    unknown_bond_state,
)
from .reaction_correspondence import (
    infer_fragmented_scaffold_correspondence_candidates,
    infer_global_correspondence_candidates,
    infer_scaffold_correspondence_candidates,
)


@dataclass(frozen=True)
class EditNormalizationResult:
    """Normalized edits plus validation and reconciliation evidence."""

    edits: Tuple[ReactionEdit, ...]
    evidence: str
    confidence: float
    warnings: Tuple[str, ...] = ()
    valid: bool = True
    stereo_changes: Tuple[ReactionStereoChange, ...] = ()
    connectivity_edit_graph: Optional[ConnectivityEditGraph] = None
    evidence_candidates: Tuple[ReactionEvidenceCandidate, ...] = ()
    edit_hypotheses: Tuple[ReactionEditHypothesis, ...] = ()
    atom_correspondence: Tuple[Tuple[int, int, int, int], ...] = ()


@dataclass(frozen=True)
class _MappedSide:
    atoms: Dict[int, ReactionAtomReference]
    bonds: Dict[Tuple[int, int], str]
    hydrogen_counts: Dict[int, int]
    warnings: Tuple[str, ...]
    mapped_atom_count: int
    heavy_atom_count: int
    mapped_heavy_atom_count: int


def _environment_id(mol: Any, atom_index: int) -> str:
    atom = mol.GetAtomWithIdx(atom_index)
    payload = {
        "aromatic": bool(atom.GetIsAromatic()),
        "charge": int(atom.GetFormalCharge()),
        "element": atom.GetSymbol(),
        "hybridization": str(atom.GetHybridization()),
        "neighbors": sorted(
            (
                neighbor.GetSymbol(),
                str(mol.GetBondBetweenAtoms(atom_index, neighbor.GetIdx()).GetBondType()),
                bool(neighbor.GetIsAromatic()),
                int(neighbor.GetFormalCharge()),
            )
            for neighbor in atom.GetNeighbors()
        ),
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return "AE1:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:20]


def _atom_reference(
    component: ReactionComponent,
    atom_index: int,
    *,
    side: Optional[str] = None,
) -> ReactionAtomReference:
    from rdkit import Chem

    mol = parse_smiles(component.input_smiles)
    if mol is None:
        raise ValueError(f"Cannot parse component {component.component_index}")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    atom = mol.GetAtomWithIdx(int(atom_index))
    map_number = int(atom.GetAtomMapNum()) or None
    chiral_tag = (
        str(atom.GetChiralTag())
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        else None
    )
    return ReactionAtomReference(
        side=side or component.side,
        component_index=component.component_index,
        atom_index=int(atom_index),
        atom_map_number=map_number,
        element=atom.GetSymbol(),
        formal_charge=int(atom.GetFormalCharge()),
        aromatic=bool(atom.GetIsAromatic()),
        hybridization=str(atom.GetHybridization()),
        local_environment_id=_environment_id(mol, int(atom_index)),
        chiral_tag=chiral_tag,
        cip_code=(
            str(atom.GetProp("_CIPCode"))
            if atom.HasProp("_CIPCode")
            else None
        ),
    )


def reaction_atom_reference(
    component: ReactionComponent,
    atom_index: int,
    *,
    side: Optional[str] = None,
) -> ReactionAtomReference:
    """Build a normalized atom reference for an explicitly identified atom."""
    return _atom_reference(component, atom_index, side=side)


def _mapped_side(components: Tuple[ReactionComponent, ...]) -> _MappedSide:
    atoms: Dict[int, ReactionAtomReference] = {}
    bonds: Dict[Tuple[int, int], str] = {}
    hydrogen_counts: Dict[int, int] = {}
    warnings = []
    mapped_atom_count = 0
    heavy_atom_count = 0
    mapped_heavy_atom_count = 0
    for component in components:
        mol = parse_smiles(component.input_smiles)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() > 1:
                heavy_atom_count += 1
            map_number = int(atom.GetAtomMapNum())
            if not map_number:
                continue
            mapped_atom_count += 1
            if atom.GetAtomicNum() > 1:
                mapped_heavy_atom_count += 1
            if map_number in atoms:
                warnings.append(
                    f"DUPLICATE_ATOM_MAP:{component.side}:{map_number}"
                )
                continue
            atoms[map_number] = _atom_reference(component, atom.GetIdx())
            if atom.GetSymbol() != "H":
                hydrogen_counts[map_number] = int(
                    atom.GetTotalNumHs(includeNeighbors=True)
                )
        for bond in mol.GetBonds():
            left = int(bond.GetBeginAtom().GetAtomMapNum())
            right = int(bond.GetEndAtom().GetAtomMapNum())
            if not left or not right:
                continue
            pair = tuple(sorted((left, right)))
            order = str(bond.GetBondType()).upper()
            if pair in bonds and bonds[pair] != order:
                warnings.append(
                    f"CONTRADICTORY_MAPPED_BOND:{component.side}:{left}:{right}"
                )
            bonds[pair] = order
    return _MappedSide(
        atoms,
        bonds,
        hydrogen_counts,
        tuple(sorted(set(warnings))),
        mapped_atom_count,
        heavy_atom_count,
        mapped_heavy_atom_count,
    )


def _mapped_connectivity_graph(
    left: _MappedSide,
    right: _MappedSide,
) -> ConnectivityEditGraph:
    """Preserve exact versus projected states before compatibility collapse."""
    transitions = []
    hydrogens = []
    atom_states = []
    warnings = []
    right_heavy_mapping_complete = (
        right.mapped_heavy_atom_count == right.heavy_atom_count
    )
    for pair in sorted(set(left.bonds).union(right.bonds)):
        old_order = left.bonds.get(pair)
        new_order = right.bonds.get(pair)
        if old_order == new_order:
            continue
        atom_1 = left.atoms.get(pair[0]) or right.atoms.get(pair[0])
        atom_2 = left.atoms.get(pair[1]) or right.atoms.get(pair[1])
        if atom_1 is None or atom_2 is None:
            warnings.append(
                f"CONNECTIVITY_TRANSITION_ENDPOINT_MISSING:{pair[0]}:{pair[1]}"
            )
            continue
        left_has_both = pair[0] in left.atoms and pair[1] in left.atoms
        right_has_both = pair[0] in right.atoms and pair[1] in right.atoms
        before_state = (
            bond_state(old_order)
            if old_order is not None
            else bond_state(None)
            if left_has_both
            else unknown_bond_state()
        )
        if new_order is not None:
            after_state = bond_state(new_order)
        elif right_has_both:
            after_state = bond_state(None)
        elif old_order is not None and right_heavy_mapping_complete:
            after_state = endpoint_absent_state()
        else:
            after_state = unknown_bond_state()
        if (
            before_state.state_kind in {"bond", "no_bond"}
            and after_state.state_kind in {"bond", "no_bond"}
            and left_has_both
            and right_has_both
        ):
            scope = "observed_product"
            confidence = 1.0
        elif after_state.state_kind == "endpoint_absent":
            scope = "main_product_projection"
            confidence = 1.0
            warnings.append(
                f"PROJECTED_ATTACHMENT_NOT_FULLY_OBSERVED:{pair[0]}:{pair[1]}"
            )
        else:
            scope = "unresolved"
            confidence = 0.0
            warnings.append(
                f"PRODUCT_ENDPOINT_WITHOUT_REACTANT_PROVENANCE:{pair[0]}:{pair[1]}"
            )
        transitions.append(
            make_bond_transition(
                atom_1=atom_1,
                atom_2=atom_2,
                before_state=before_state,
                after_state=after_state,
                observation_scope=scope,
                evidence="supplied_atom_mapping",
                confidence=confidence,
            )
        )
    for map_number in sorted(
        set(left.hydrogen_counts).intersection(right.hydrogen_counts)
    ):
        before_count = left.hydrogen_counts[map_number]
        after_count = right.hydrogen_counts[map_number]
        if before_count == after_count:
            continue
        hydrogens.append(
            HydrogenDelta(
                atom=left.atoms[map_number],
                before_count=before_count,
                after_count=after_count,
                delta_count=after_count - before_count,
                observation_scope="observed_product",
                evidence="supplied_atom_mapping",
                confidence=1.0,
            )
        )
    for map_number in sorted(set(left.atoms).intersection(right.atoms)):
        before = left.atoms[map_number]
        after = right.atoms[map_number]
        if before.formal_charge == after.formal_charge:
            continue
        atom_states.append(
            AtomStateTransition(
                reactant_atom=before,
                product_atom=after,
                before_formal_charge=before.formal_charge,
                after_formal_charge=after.formal_charge,
                before_radical_electrons=None,
                after_radical_electrons=None,
                before_isotope=None,
                after_isotope=None,
                observation_scope="observed_product",
                evidence="supplied_atom_mapping",
                confidence=1.0,
            )
        )
    return build_connectivity_edit_graph(
        bond_transitions=transitions,
        hydrogen_deltas=hydrogens,
        atom_state_transitions=atom_states,
        evidence="validated_atom_mapping",
        confidence=1.0 if transitions or hydrogens or atom_states else 0.0,
        warnings=warnings,
    )


def _reactant_hydrogen_counts(
    components: Tuple[ReactionComponent, ...],
) -> Dict[Tuple[object, ...], int]:
    counts: Dict[Tuple[object, ...], int] = {}
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            if atom.GetSymbol() == "H":
                continue
            reference = _atom_reference(component, int(atom.GetIdx()))
            counts[atom_provenance_key(reference)] = int(
                atom.GetTotalNumHs(includeNeighbors=True)
            )
    return counts


def normalize_mapped_edits(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
) -> EditNormalizationResult:
    """Validate supplied atom maps and extract typed observed bond edits."""
    left = _mapped_side(reactants)
    right = _mapped_side(products)
    warnings = list(left.warnings + right.warnings)
    if not left.mapped_atom_count and not right.mapped_atom_count:
        return EditNormalizationResult((), "no_atom_mapping", 0.0, valid=False)
    if not left.mapped_atom_count or not right.mapped_atom_count:
        warnings.append("INCOMPLETE_REACTION_ATOM_MAPPING")
    if (
        left.mapped_heavy_atom_count < left.heavy_atom_count
        or right.mapped_heavy_atom_count < right.heavy_atom_count
    ):
        warnings.append("PARTIAL_ATOM_MAPPING")
    if set(right.atoms) - set(left.atoms):
        warnings.append("PRODUCT_MAPS_MISSING_FROM_REACTANTS")
    if set(left.atoms) - set(right.atoms):
        warnings.append("REACTANT_MAPS_MISSING_FROM_PRODUCTS")
    for map_number in sorted(set(left.atoms).intersection(right.atoms)):
        before = left.atoms[map_number]
        after = right.atoms[map_number]
        if before.element != after.element:
            warnings.append(
                f"ATOM_MAP_ELEMENT_MISMATCH:{map_number}:{before.element}:{after.element}"
            )
    fatal = any(
        warning.startswith(
            ("DUPLICATE_ATOM_MAP", "CONTRADICTORY_MAPPED_BOND", "ATOM_MAP_ELEMENT_MISMATCH")
        )
        for warning in warnings
    )
    if fatal:
        return EditNormalizationResult(
            (), "invalid_atom_mapping", 0.0, tuple(sorted(set(warnings))), False
        )
    left_only_maps = set(left.atoms) - set(right.atoms)
    ignored_reactant_internal_bonds = 0
    edits = []
    for pair in sorted(set(left.bonds).union(right.bonds)):
        old_order = left.bonds.get(pair)
        new_order = right.bonds.get(pair)
        if old_order == new_order:
            continue
        # A main-product reaction record commonly omits a complete departing
        # fragment. Bonds wholly inside that one-sided fragment are inventory
        # differences, not reaction-center edits. Retain the boundary bond
        # whenever one endpoint survives, so C-X and C-O leaving-group
        # departures remain explicit structural evidence.
        if (
            old_order is not None
            and new_order is None
            and set(pair).issubset(left_only_maps)
        ):
            ignored_reactant_internal_bonds += 1
            continue
        if old_order is None:
            edit_type = "formed"
        elif new_order is None:
            edit_type = "broken"
        else:
            edit_type = "order_changed"
        atom_1 = left.atoms.get(pair[0]) or right.atoms.get(pair[0])
        atom_2 = left.atoms.get(pair[1]) or right.atoms.get(pair[1])
        if atom_1 is None or atom_2 is None:
            warnings.append(f"MAPPED_EDIT_ENDPOINT_MISSING:{pair[0]}:{pair[1]}")
            continue
        edits.append(
            ReactionEdit(
                edit_type=edit_type,
                atom_1=atom_1,
                atom_2=atom_2,
                old_order=old_order,
                new_order=new_order,
                evidence="supplied_atom_mapping",
                confidence=1.0,
            )
        )
    if ignored_reactant_internal_bonds:
        warnings.append(
            "REACTANT_ONLY_FRAGMENT_INTERNAL_BONDS_IGNORED:"
            f"{ignored_reactant_internal_bonds}"
        )
    for map_number in sorted(
        set(left.hydrogen_counts).intersection(right.hydrogen_counts)
    ):
        old_count = left.hydrogen_counts[map_number]
        new_count = right.hydrogen_counts[map_number]
        delta = new_count - old_count
        if not delta:
            continue
        center = left.atoms[map_number]
        for _ in range(abs(delta)):
            edits.append(
                ReactionEdit(
                    edit_type="hydrogen_change",
                    atom_1=center,
                    atom_2=None,
                    old_order="SINGLE" if delta < 0 else None,
                    new_order="SINGLE" if delta > 0 else None,
                    evidence="supplied_atom_mapping",
                    confidence=1.0,
                )
            )
    if not edits:
        warnings.append("NO_MAPPED_BOND_EDITS")
    atom_correspondence = tuple(
        (
            left.atoms[map_number].component_index,
            left.atoms[map_number].atom_index,
            right.atoms[map_number].component_index,
            right.atoms[map_number].atom_index,
        )
        for map_number in sorted(set(left.atoms).intersection(right.atoms))
    )
    stereo_changes = _correspondence_stereo_changes(
        atom_correspondence,
        reactants,
        products,
        evidence="supplied_atom_mapping",
        confidence=1.0,
    )
    connectivity_graph = _mapped_connectivity_graph(left, right)
    return EditNormalizationResult(
        tuple(edits),
        "validated_atom_mapping" if edits else "atom_mapping_without_edits",
        1.0 if edits else 0.0,
        tuple(sorted(set(warnings))),
        bool(edits),
        stereo_changes,
        connectivity_graph,
    )


def _correspondence_edits(
    mapping: Tuple[Tuple[int, int, int, int], ...],
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    evidence: str = "unique_scaffold_correspondence",
    confidence: float = 0.85,
) -> Tuple[ReactionEdit, ...]:
    reactant_components = {
        component.component_index: component for component in reactants
    }
    product_components = {
        component.component_index: component for component in products
    }
    forward = {
        (reactant_component, reactant_atom): (product_component, product_atom)
        for reactant_component, reactant_atom, product_component, product_atom
        in mapping
    }
    reverse = {product: reactant for reactant, product in forward.items()}
    edits = []
    for component in reactants:
        mol = parse_smiles(component.input_smiles)
        if mol is None:
            continue
        for bond in mol.GetBonds():
            left = (component.component_index, bond.GetBeginAtomIdx())
            right = (component.component_index, bond.GetEndAtomIdx())
            mapped_left = forward.get(left)
            mapped_right = forward.get(right)
            old_order = str(bond.GetBondType()).upper()
            if mapped_left is not None and mapped_right is not None:
                if mapped_left[0] != mapped_right[0]:
                    continue
                product_component = product_components[mapped_left[0]]
                product_mol = parse_smiles(product_component.input_smiles)
                if product_mol is None:
                    continue
                product_bond = product_mol.GetBondBetweenAtoms(
                    mapped_left[1], mapped_right[1]
                )
                new_order = (
                    str(product_bond.GetBondType()).upper()
                    if product_bond is not None
                    else None
                )
                if new_order == old_order:
                    continue
                edit_type = "broken" if new_order is None else "order_changed"
            elif (mapped_left is None) == (mapped_right is None):
                continue
            else:
                unmapped_index = right[1] if mapped_left is not None else left[1]
                if mol.GetAtomWithIdx(unmapped_index).GetAtomicNum() <= 1:
                    continue
                new_order = None
                edit_type = "broken"
            edits.append(
                ReactionEdit(
                    edit_type=edit_type,
                    atom_1=_atom_reference(component, left[1]),
                    atom_2=_atom_reference(component, right[1]),
                    old_order=old_order,
                    new_order=new_order,
                    evidence=evidence,
                    confidence=confidence,
                )
            )
    for component in products:
        mol = parse_smiles(component.input_smiles)
        if mol is None:
            continue
        for bond in mol.GetBonds():
            left = (component.component_index, bond.GetBeginAtomIdx())
            right = (component.component_index, bond.GetEndAtomIdx())
            reactant_left = reverse.get(left)
            reactant_right = reverse.get(right)
            if reactant_left is None or reactant_right is None:
                continue
            if reactant_left[0] != reactant_right[0]:
                old_bond = None
            else:
                reactant_component = reactant_components[reactant_left[0]]
                reactant_mol = parse_smiles(reactant_component.input_smiles)
                old_bond = (
                    reactant_mol.GetBondBetweenAtoms(
                        reactant_left[1], reactant_right[1]
                    )
                    if reactant_mol is not None
                    else None
                )
            if old_bond is not None:
                continue
            reactant_component = reactant_components[reactant_left[0]]
            edits.append(
                ReactionEdit(
                    edit_type="formed",
                    atom_1=_atom_reference(reactant_component, reactant_left[1]),
                    atom_2=_atom_reference(
                        reactant_components[reactant_right[0]], reactant_right[1]
                    ),
                    old_order=None,
                    new_order=str(bond.GetBondType()).upper(),
                    evidence=evidence,
                    confidence=confidence,
                )
            )
    for reactant_key, product_key in sorted(forward.items()):
        reactant_component = reactant_components[reactant_key[0]]
        product_component = product_components[product_key[0]]
        reactant_mol = parse_smiles(reactant_component.input_smiles)
        product_mol = parse_smiles(product_component.input_smiles)
        if reactant_mol is None or product_mol is None:
            continue
        before = reactant_mol.GetAtomWithIdx(reactant_key[1])
        after = product_mol.GetAtomWithIdx(product_key[1])
        old_count = int(before.GetTotalNumHs(includeNeighbors=True))
        new_count = int(after.GetTotalNumHs(includeNeighbors=True))
        delta = new_count - old_count
        for _ in range(abs(delta)):
            edits.append(
                ReactionEdit(
                    edit_type="hydrogen_change",
                    atom_1=_atom_reference(reactant_component, reactant_key[1]),
                    atom_2=None,
                    old_order="SINGLE" if delta < 0 else None,
                    new_order="SINGLE" if delta > 0 else None,
                    evidence=evidence,
                    confidence=confidence,
                )
            )
    return tuple(edits)


def _correspondence_connectivity_graph(
    mapping: Tuple[Tuple[int, int, int, int], ...],
    edits: Tuple[ReactionEdit, ...],
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    evidence: str,
    confidence: float,
) -> ConnectivityEditGraph:
    """Build inferred transitions while retaining projected attachment loss."""
    mapped_reactant_atoms = {
        (reactant_component, reactant_atom)
        for reactant_component, reactant_atom, _, _ in mapping
    }
    transitions = []
    warnings = []
    for edit in edits:
        if edit.edit_type == "hydrogen_change":
            continue
        if edit.atom_2 is None:
            warnings.append("CONNECTIVITY_EDIT_ENDPOINT_MISSING")
            continue
        left_position = (edit.atom_1.component_index, edit.atom_1.atom_index)
        right_position = (edit.atom_2.component_index, edit.atom_2.atom_index)
        projected = (
            edit.edit_type == "broken"
            and (
                left_position not in mapped_reactant_atoms
                or right_position not in mapped_reactant_atoms
            )
        )
        transitions.append(
            make_bond_transition(
                atom_1=edit.atom_1,
                atom_2=edit.atom_2,
                before_state=bond_state(edit.old_order),
                after_state=(
                    endpoint_absent_state()
                    if projected
                    else bond_state(edit.new_order)
                ),
                observation_scope=(
                    "main_product_projection"
                    if projected
                    else "correspondence_inference"
                ),
                evidence=edit.evidence,
                confidence=edit.confidence,
            )
        )
        if projected:
            warnings.append("PROJECTED_ATTACHMENT_NOT_FULLY_OBSERVED")
    hydrogen_graph = connectivity_graph_from_reaction_edits(
        tuple(edit for edit in edits if edit.edit_type == "hydrogen_change"),
        observation_scope="correspondence_inference",
        evidence=evidence,
        confidence=confidence,
        hydrogen_before_counts=_reactant_hydrogen_counts(reactants),
    )
    reactant_components = {
        component.component_index: component for component in reactants
    }
    product_components = {
        component.component_index: component for component in products
    }
    atom_states = []
    for (
        reactant_component_index,
        reactant_atom_index,
        product_component_index,
        product_atom_index,
    ) in mapping:
        reactant_component = reactant_components[reactant_component_index]
        product_component = product_components[product_component_index]
        before = _atom_reference(reactant_component, reactant_atom_index)
        after = _atom_reference(product_component, product_atom_index)
        if before.formal_charge == after.formal_charge:
            continue
        atom_states.append(
            AtomStateTransition(
                reactant_atom=before,
                product_atom=after,
                before_formal_charge=before.formal_charge,
                after_formal_charge=after.formal_charge,
                before_radical_electrons=None,
                after_radical_electrons=None,
                before_isotope=None,
                after_isotope=None,
                observation_scope="correspondence_inference",
                evidence=evidence,
                confidence=confidence,
            )
        )
    return build_connectivity_edit_graph(
        bond_transitions=transitions,
        hydrogen_deltas=hydrogen_graph.hydrogen_deltas,
        atom_state_transitions=atom_states,
        evidence=evidence,
        confidence=confidence,
        warnings=tuple(warnings) + hydrogen_graph.warnings,
    )


def _atom_stereo_descriptor(molecule: Any, atom_index: int) -> Optional[str]:
    from rdkit import Chem

    atom = molecule.GetAtomWithIdx(int(atom_index))
    if atom.HasProp("_CIPCode"):
        return str(atom.GetProp("_CIPCode"))
    tag = atom.GetChiralTag()
    if tag == Chem.ChiralType.CHI_UNSPECIFIED:
        return None
    return str(tag).replace("CHI_", "")


def _bond_stereo_descriptor(bond: Any) -> Optional[str]:
    from rdkit import Chem

    if bond is None:
        return None
    stereo = bond.GetStereo()
    if stereo in {Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY}:
        return None
    value = str(stereo).replace("STEREO", "")
    return {"CIS": "Z", "TRANS": "E"}.get(value, value)


def _stereo_change_type(
    old_descriptor: Optional[str],
    new_descriptor: Optional[str],
) -> str:
    if old_descriptor == new_descriptor:
        return "retained"
    if old_descriptor is None:
        return "created"
    if new_descriptor is None:
        return "destroyed"
    return "descriptor_changed"


def _correspondence_stereo_changes(
    mapping: Tuple[Tuple[int, int, int, int], ...],
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    evidence: str,
    confidence: float,
) -> Tuple[ReactionStereoChange, ...]:
    """Extract explicit atom and E/Z descriptors across one correspondence."""
    reactant_components = {
        component.component_index: component for component in reactants
    }
    product_components = {
        component.component_index: component for component in products
    }
    reactant_molecules = {
        index: parse_smiles(component.input_smiles)
        for index, component in reactant_components.items()
    }
    product_molecules = {
        index: parse_smiles(component.input_smiles)
        for index, component in product_components.items()
    }
    molecules = tuple(
        molecule
        for molecule in (
            tuple(reactant_molecules.values())
            + tuple(product_molecules.values())
        )
        if molecule is not None
    )
    if not any(
        _has_explicit_stereochemistry(molecule) for molecule in molecules
    ):
        return ()
    from rdkit import Chem

    for molecule in molecules:
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    forward = {
        (reactant_component, reactant_atom): (product_component, product_atom)
        for reactant_component, reactant_atom, product_component, product_atom
        in mapping
    }
    reverse = {product: reactant for reactant, product in forward.items()}
    changes = []
    for reactant_key, product_key in sorted(forward.items()):
        reactant_molecule = reactant_molecules.get(reactant_key[0])
        product_molecule = product_molecules.get(product_key[0])
        if reactant_molecule is None or product_molecule is None:
            continue
        old_descriptor = _atom_stereo_descriptor(
            reactant_molecule, reactant_key[1]
        )
        new_descriptor = _atom_stereo_descriptor(
            product_molecule, product_key[1]
        )
        if old_descriptor is None and new_descriptor is None:
            continue
        changes.append(
            ReactionStereoChange(
                stereo_type="atom",
                atom_1=_atom_reference(
                    reactant_components[reactant_key[0]], reactant_key[1]
                ),
                atom_2=None,
                old_descriptor=old_descriptor,
                new_descriptor=new_descriptor,
                change_type=_stereo_change_type(
                    old_descriptor, new_descriptor
                ),
                evidence=evidence,
                confidence=confidence,
            )
        )
    for product_component_index, product_molecule in product_molecules.items():
        if product_molecule is None:
            continue
        for product_bond in product_molecule.GetBonds():
            product_left = (
                product_component_index,
                int(product_bond.GetBeginAtomIdx()),
            )
            product_right = (
                product_component_index,
                int(product_bond.GetEndAtomIdx()),
            )
            reactant_left = reverse.get(product_left)
            reactant_right = reverse.get(product_right)
            if reactant_left is None or reactant_right is None:
                continue
            old_bond = None
            if reactant_left[0] == reactant_right[0]:
                reactant_molecule = reactant_molecules.get(reactant_left[0])
                if reactant_molecule is not None:
                    old_bond = reactant_molecule.GetBondBetweenAtoms(
                        reactant_left[1], reactant_right[1]
                    )
            old_descriptor = _bond_stereo_descriptor(old_bond)
            new_descriptor = _bond_stereo_descriptor(product_bond)
            if old_descriptor is None and new_descriptor is None:
                continue
            left_component = reactant_components[reactant_left[0]]
            right_component = reactant_components[reactant_right[0]]
            changes.append(
                ReactionStereoChange(
                    stereo_type="bond",
                    atom_1=_atom_reference(
                        left_component, reactant_left[1]
                    ),
                    atom_2=_atom_reference(
                        right_component, reactant_right[1]
                    ),
                    old_descriptor=old_descriptor,
                    new_descriptor=new_descriptor,
                    change_type=_stereo_change_type(
                        old_descriptor, new_descriptor
                    ),
                    evidence=evidence,
                    confidence=confidence,
                )
            )
    return tuple(sorted(changes, key=_stereo_comparison_key))


def _stereo_comparison_key(change: ReactionStereoChange) -> Tuple[Any, ...]:
    endpoints = tuple(
        sorted(
            (
                atom.component_index,
                atom.atom_index,
                atom.element,
                atom.local_environment_id,
            )
            for atom in (change.atom_1, change.atom_2)
            if atom is not None
        )
    )
    return (
        change.stereo_type,
        endpoints,
        change.old_descriptor or "NONE",
        change.new_descriptor or "NONE",
        change.change_type,
    )


def _chemistry_stereo_key(change: ReactionStereoChange) -> Tuple[Any, ...]:
    endpoints = tuple(
        sorted(
            (
                atom.element,
                atom.formal_charge,
                atom.aromatic,
                atom.hybridization,
                atom.local_environment_id,
            )
            for atom in (change.atom_1, change.atom_2)
            if atom is not None
        )
    )
    return (
        change.stereo_type,
        endpoints,
        change.old_descriptor or "NONE",
        change.new_descriptor or "NONE",
        change.change_type,
    )


def _chemistry_edit_key(edit: ReactionEdit) -> Tuple[Any, ...]:
    endpoints = tuple(
        sorted(
            (
                atom.element,
                atom.formal_charge,
                atom.aromatic,
                atom.hybridization,
                atom.local_environment_id,
            )
            for atom in (edit.atom_1, edit.atom_2)
            if atom is not None
        )
    )
    return (
        edit.edit_type,
        endpoints,
        edit.old_order or "NONE",
        edit.new_order or "NONE",
    )


def _correspondence_edit_cost(
    edits: Tuple[ReactionEdit, ...],
) -> Tuple[int, int, int]:
    """Rank inferred mappings by the smallest chemically explicit edit set."""
    heavy_edits = tuple(
        edit for edit in edits if edit.edit_type != "hydrogen_change"
    )
    weighted = sum(
        2 if edit.edit_type in {"formed", "broken"} else 1
        for edit in heavy_edits
    )
    hydrogen_edits = len(edits) - len(heavy_edits)
    return weighted, len(heavy_edits), hydrogen_edits


def _correspondence_provider(evidence: str) -> str:
    return {
        "unique_scaffold_correspondence": "scaffold_correspondence",
        "global_atom_correspondence": "global_correspondence",
        "fragmented_scaffold_correspondence": (
            "fragmented_scaffold_correspondence"
        ),
    }.get(evidence, "correspondence")


def _edit_hypotheses(
    candidate_results: Tuple[
        Tuple[
            Tuple[Tuple[int, int, int, int], ...],
            Tuple[ReactionEdit, ...],
            Tuple[ReactionStereoChange, ...],
        ],
        ...,
    ],
    *,
    evidence: str,
    confidence: float,
) -> Tuple[ReactionEditHypothesis, ...]:
    """Collapse mapping multiplicity while retaining distinct edit alternatives."""
    grouped: Dict[
        Tuple[Tuple[Any, ...], Tuple[Any, ...]],
        list[
            Tuple[
                Tuple[Tuple[int, int, int, int], ...],
                Tuple[ReactionEdit, ...],
                Tuple[ReactionStereoChange, ...],
            ]
        ],
    ] = {}
    for result in candidate_results:
        _, edits, stereo = result
        if not edits:
            continue
        key = (
            tuple(sorted(_chemistry_edit_key(edit) for edit in edits)),
            tuple(
                sorted(_chemistry_stereo_key(change) for change in stereo)
            ),
        )
        grouped.setdefault(key, []).append(result)

    hypotheses = []
    provider = _correspondence_provider(evidence)
    for key, values in sorted(
        grouped.items(),
        key=lambda item: json.dumps(
            item[0],
            sort_keys=True,
            separators=(",", ":"),
        ),
    ):
        selected = min(
            values,
            key=lambda result: (
                tuple(sorted(_comparison_key(edit) for edit in result[1])),
                tuple(
                    sorted(
                        _stereo_comparison_key(change)
                        for change in result[2]
                    )
                ),
            ),
        )
        canonical = json.dumps(
            {
                "provider": provider,
                "evidence": evidence,
                "edits": key[0],
                "stereo": key[1],
                "schema_version": "1.0",
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        hypothesis_id = (
            "REH1:"
            + hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:24]
        )
        hypotheses.append(
            ReactionEditHypothesis(
                hypothesis_id=hypothesis_id,
                provider=provider,
                evidence=evidence,
                confidence=confidence,
                edits=selected[1],
                stereo_changes=selected[2],
                correspondence_count=len(values),
                edit_cost=_correspondence_edit_cost(selected[1]),
                atom_correspondence=selected[0],
                warnings=("UNVERIFIED_EDIT_HYPOTHESIS",),
            )
        )
    return tuple(hypotheses)


def normalize_inferred_scaffold_edits(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
) -> EditNormalizationResult:
    """Infer edits only when all best scaffold mappings imply one chemistry."""
    correspondence = infer_scaffold_correspondence_candidates(reactants, products)
    correspondence_warnings = set(correspondence.warnings)
    evidence = "unique_scaffold_correspondence"
    confidence = 0.85
    inferred_warning = "INFERRED_ATOM_CORRESPONDENCE"
    if not correspondence.valid:
        correspondence = infer_global_correspondence_candidates(
            reactants, products
        )
        correspondence_warnings.update(correspondence.warnings)
        evidence = "global_atom_correspondence"
        confidence = 0.8
        inferred_warning = "INFERRED_GLOBAL_ATOM_CORRESPONDENCE"
    if not correspondence.valid:
        correspondence = infer_fragmented_scaffold_correspondence_candidates(
            reactants, products
        )
        correspondence_warnings.update(correspondence.warnings)
        evidence = "fragmented_scaffold_correspondence"
        confidence = 0.75
        inferred_warning = "INFERRED_FRAGMENTED_SCAFFOLD_CORRESPONDENCE"
    if not correspondence.valid:
        unresolved_candidate = ReactionEvidenceCandidate(
            provider="correspondence",
            status="unresolved",
            evidence="unresolved",
            confidence=0.0,
            warnings=tuple(sorted(correspondence_warnings)),
        )
        return EditNormalizationResult(
            (),
            "unresolved",
            0.0,
            tuple(sorted(correspondence_warnings)),
            False,
            evidence_candidates=(unresolved_candidate,),
        )
    all_candidate_results = tuple(
        (
            mapping,
            _correspondence_edits(
                mapping,
                reactants,
                products,
                evidence=evidence,
                confidence=confidence,
            ),
            _correspondence_stereo_changes(
                mapping,
                reactants,
                products,
                evidence=evidence,
                confidence=confidence,
            ),
        )
        for mapping in correspondence.candidates
    )
    if evidence in {
        "global_atom_correspondence",
        "fragmented_scaffold_correspondence",
    }:
        nonempty_costs = tuple(
            _correspondence_edit_cost(edits)
            for _, edits, _ in all_candidate_results
            if edits
            and (
                evidence != "fragmented_scaffold_correspondence"
                or sum(
                    edit.edit_type != "hydrogen_change"
                    for edit in edits
                )
                <= 12
            )
        )
        best_cost = min(nonempty_costs) if nonempty_costs else None
        candidate_results = tuple(
            result
            for result in all_candidate_results
            for _, edits, _ in (result,)
            if edits
            and (
                evidence != "fragmented_scaffold_correspondence"
                or sum(
                    edit.edit_type != "hydrogen_change"
                    for edit in edits
                )
                <= 12
            )
            and (
                _correspondence_edit_cost(edits) == best_cost
                if evidence != "fragmented_scaffold_correspondence"
                else best_cost is not None
                and _correspondence_edit_cost(edits)[0] <= best_cost[0] + 3
                and _correspondence_edit_cost(edits)[1] <= best_cost[1] + 2
            )
        )
    else:
        candidate_results = all_candidate_results
    nonempty = tuple(result for result in candidate_results if result[1])
    if not nonempty:
        unresolved_candidate = ReactionEvidenceCandidate(
            provider=_correspondence_provider(evidence),
            status="unresolved",
            evidence="unresolved",
            confidence=0.0,
            warnings=("SCAFFOLD_CORRESPONDENCE_WITHOUT_EDITS",),
        )
        return EditNormalizationResult(
            (),
            "unresolved",
            0.0,
            ("SCAFFOLD_CORRESPONDENCE_WITHOUT_EDITS",),
            False,
            evidence_candidates=(unresolved_candidate,),
        )
    hypotheses = _edit_hypotheses(
        nonempty,
        evidence=evidence,
        confidence=confidence,
    )
    if len(hypotheses) != 1 or len(nonempty) != len(candidate_results):
        ambiguous_candidate = ReactionEvidenceCandidate(
            provider=_correspondence_provider(evidence),
            status="ambiguous",
            evidence="ambiguous_atom_correspondence",
            confidence=0.0,
            edit_hypotheses=hypotheses,
            warnings=(
                f"AMBIGUOUS_SCAFFOLD_CORRESPONDENCE:{len(hypotheses)}",
            ),
        )
        return EditNormalizationResult(
            (),
            "ambiguous_atom_correspondence",
            0.0,
            (f"AMBIGUOUS_SCAFFOLD_CORRESPONDENCE:{len(hypotheses)}",),
            False,
            evidence_candidates=(ambiguous_candidate,),
            edit_hypotheses=hypotheses,
        )
    selected = min(
        nonempty,
        key=lambda result: (
            tuple(sorted(_comparison_key(edit) for edit in result[1])),
            tuple(
                sorted(
                    _stereo_comparison_key(change) for change in result[2]
                )
            ),
        ),
    )
    connectivity_graph = _correspondence_connectivity_graph(
        selected[0],
        selected[1],
        reactants,
        products,
        evidence=evidence,
        confidence=confidence,
    )
    success_warnings = tuple(
        sorted(
            {inferred_warning}.union(
                warning
                for warning in correspondence.warnings
                if warning.startswith("INFERRED_")
            )
        )
    )
    verified_candidate = ReactionEvidenceCandidate(
        provider=_correspondence_provider(evidence),
        status="verified",
        evidence=evidence,
        confidence=confidence,
        edits=selected[1],
        stereo_changes=selected[2],
        warnings=success_warnings,
    )
    return EditNormalizationResult(
        selected[1],
        evidence,
        confidence,
        success_warnings,
        True,
        selected[2],
        connectivity_graph,
        evidence_candidates=(verified_candidate,),
        atom_correspondence=selected[0],
    )


def _comparison_key(edit: ReactionEdit) -> Tuple[Any, ...]:
    endpoints = []
    for atom in (edit.atom_1, edit.atom_2):
        if atom is None:
            endpoints.append(("H",))
        elif atom.atom_map_number is not None:
            endpoints.append(("map", atom.atom_map_number))
        else:
            endpoints.append(
                ("atom", atom.component_index, atom.atom_index, atom.element)
            )
    return (
        edit.edit_type,
        tuple(sorted(endpoints)),
        edit.old_order or "NONE",
        edit.new_order or "NONE",
    )


def _has_explicit_stereochemistry(molecule: Any) -> bool:
    from rdkit import Chem

    return any(
        atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        for atom in molecule.GetAtoms()
    ) or any(
        bond.GetStereo() not in {
            Chem.BondStereo.STEREONONE,
            Chem.BondStereo.STEREOANY,
        }
        for bond in molecule.GetBonds()
    )


def _provider_evidence_candidate(
    provider: str,
    result: EditNormalizationResult,
) -> Optional[ReactionEvidenceCandidate]:
    """Expose one attempted provider without changing evidence authority."""
    if result.evidence in {
        "no_atom_mapping",
        "no_predicted_edits",
        "no_predicted_multi_event_edits",
    }:
        return None
    status = (
        "verified"
        if result.valid and bool(result.edits)
        else "invalid"
        if result.evidence == "invalid_atom_mapping"
        else "ambiguous"
        if result.edit_hypotheses
        else "unresolved"
    )
    return ReactionEvidenceCandidate(
        provider=provider,
        status=status,
        evidence=result.evidence,
        confidence=result.confidence,
        edits=result.edits if status == "verified" else (),
        stereo_changes=result.stereo_changes if status == "verified" else (),
        edit_hypotheses=result.edit_hypotheses,
        warnings=result.warnings,
    )


def _provider_candidates(
    values: Tuple[Tuple[str, EditNormalizationResult], ...],
) -> Tuple[ReactionEvidenceCandidate, ...]:
    return tuple(
        candidate
        for provider, result in values
        for candidate in (_provider_evidence_candidate(provider, result),)
        if candidate is not None
    )


def resolve_structural_evidence(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    mapped_override: Optional[EditNormalizationResult] = None,
    mapped_provider: str = "supplied_atom_mapping",
    additional_warnings: Tuple[str, ...] = (),
) -> EditNormalizationResult:
    """Resolve only observed mapping and graph-correspondence evidence.

    This is the canonical low-level path.  It deliberately does not evaluate
    reaction operators, reconstruction candidates, or named chemistry patterns.
    """
    mapped = (
        mapped_override
        if mapped_override is not None
        else normalize_mapped_edits(reactants, products)
    )
    mapped_candidates = _provider_candidates(((mapped_provider, mapped),))
    if mapped.evidence == "invalid_atom_mapping":
        return EditNormalizationResult(
            (),
            "unresolved",
            0.0,
            tuple(sorted(set(mapped.warnings + additional_warnings))),
            False,
            evidence_candidates=mapped_candidates,
        )
    if mapped.valid:
        return replace(
            mapped,
            warnings=tuple(sorted(set(mapped.warnings + additional_warnings))),
            evidence_candidates=mapped_candidates,
        )

    inferred = normalize_inferred_scaffold_edits(reactants, products)
    combined_candidates = mapped_candidates + inferred.evidence_candidates
    if inferred.valid:
        return replace(
            inferred,
            warnings=tuple(sorted(set(inferred.warnings + additional_warnings))),
            evidence_candidates=combined_candidates,
        )
    return EditNormalizationResult(
        (),
        inferred.evidence,
        0.0,
        tuple(
            sorted(
                set(mapped.warnings + inferred.warnings + additional_warnings)
            )
        ),
        False,
        connectivity_edit_graph=(
            mapped.connectivity_edit_graph
            if mapped.evidence == "atom_mapping_without_edits"
            else inferred.connectivity_edit_graph
        ),
        evidence_candidates=combined_candidates,
        edit_hypotheses=inferred.edit_hypotheses,
    )


__all__ = [
    "EditNormalizationResult",
    "normalize_mapped_edits",
    "normalize_inferred_scaffold_edits",
    "reaction_atom_reference",
    "resolve_structural_evidence",
]
