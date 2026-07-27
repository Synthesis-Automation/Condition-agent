"""Normalize mapped and operator-predicted graph edits into typed contracts."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import (
    ReactionAtomReference,
    ReactionCandidate,
    ReactionComponent,
    ReactionEdit,
    ReactionSiteReference,
    ReactionStereoChange,
)
from .reaction_correspondence import (
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
    edits = []
    for pair in sorted(set(left.bonds).union(right.bonds)):
        old_order = left.bonds.get(pair)
        new_order = right.bonds.get(pair)
        if old_order == new_order:
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
    stereo_changes = _correspondence_stereo_changes(
        tuple(
            (
                left.atoms[map_number].component_index,
                left.atoms[map_number].atom_index,
                right.atoms[map_number].component_index,
                right.atoms[map_number].atom_index,
            )
            for map_number in sorted(set(left.atoms).intersection(right.atoms))
        ),
        reactants,
        products,
        evidence="supplied_atom_mapping",
        confidence=1.0,
    )
    return EditNormalizationResult(
        tuple(edits),
        "validated_atom_mapping" if edits else "atom_mapping_without_edits",
        1.0 if edits else 0.0,
        tuple(sorted(set(warnings))),
        bool(edits),
        stereo_changes,
    )


def _component(
    components: Tuple[ReactionComponent, ...], component_index: int
) -> ReactionComponent:
    return next(
        component
        for component in components
        if component.component_index == component_index
    )


def _role_atom(
    role_path: str,
    assignments: Dict[str, ReactionSiteReference],
    components: Tuple[ReactionComponent, ...],
) -> ReactionAtomReference:
    partner_role, atom_role = role_path.split(".", 1)
    site = assignments[partner_role]
    indices = site.atom_roles.get(atom_role)
    if not indices and atom_role == "handle":
        indices = site.atom_roles.get("leaving_or_activatable")
    if not indices:
        raise ValueError(f"Missing atom role {role_path}")
    return _atom_reference(
        _component(components, site.component_index), int(indices[0])
    )


def normalize_predicted_edits(
    selected: Optional[ReactionCandidate],
    reactants: Tuple[ReactionComponent, ...],
) -> EditNormalizationResult:
    """Convert operator changes for an exact selected candidate to typed edits."""
    if selected is None or selected.verification not in {
        "exact_product_reconstruction",
        "exact_multi_event_reconstruction",
    }:
        return EditNormalizationResult((), "no_exact_reconstruction", 0.0, valid=False)
    edits = []
    warnings = []
    for change in selected.predicted_bond_changes:
        try:
            atom_1 = _role_atom(
                change.atom_1_role, selected.role_assignments, reactants
            )
            atom_2 = (
                _role_atom(change.atom_2_role, selected.role_assignments, reactants)
                if change.atom_2_role is not None
                else None
            )
        except (KeyError, StopIteration, ValueError) as exc:
            warnings.append(f"PREDICTED_EDIT_PROVENANCE_ERROR:{exc}")
            continue
        edits.append(
            ReactionEdit(
                edit_type=change.change_type,
                atom_1=atom_1,
                atom_2=atom_2,
                old_order=change.old_order.upper() if change.old_order else None,
                new_order=change.new_order.upper() if change.new_order else None,
                evidence=change.evidence,
                confidence=1.0,
            )
        )
    return EditNormalizationResult(
        tuple(edits),
        "exact_product_reconstruction" if edits else "edit_normalization_failed",
        1.0 if edits else 0.0,
        tuple(sorted(set(warnings))),
        bool(edits),
    )


def normalize_predicted_multi_event_edits(
    selected_events: Tuple[ReactionCandidate, ...],
    reactants: Tuple[ReactionComponent, ...],
) -> EditNormalizationResult:
    """Normalize an exactly reconstructed collection of reaction events."""
    if len(selected_events) < 2:
        return EditNormalizationResult((), "no_exact_multi_event_reconstruction", 0.0, valid=False)
    normalized = tuple(
        normalize_predicted_edits(candidate, reactants)
        for candidate in selected_events
    )
    warnings = tuple(
        sorted({warning for result in normalized for warning in result.warnings})
    )
    if not all(result.valid for result in normalized):
        return EditNormalizationResult(
            (), "multi_event_edit_normalization_failed", 0.0, warnings, False
        )
    edits = tuple(edit for result in normalized for edit in result.edits)
    return EditNormalizationResult(
        edits,
        "exact_multi_event_reconstruction",
        1.0,
        warnings,
        True,
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


def normalize_inferred_scaffold_edits(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
) -> EditNormalizationResult:
    """Infer edits only when all best scaffold mappings imply one chemistry."""
    correspondence = infer_scaffold_correspondence_candidates(reactants, products)
    evidence = "unique_scaffold_correspondence"
    confidence = 0.85
    inferred_warning = "INFERRED_ATOM_CORRESPONDENCE"
    if not correspondence.valid:
        correspondence = infer_global_correspondence_candidates(
            reactants, products
        )
        evidence = "global_atom_correspondence"
        confidence = 0.8
        inferred_warning = "INFERRED_GLOBAL_ATOM_CORRESPONDENCE"
    if not correspondence.valid:
        return EditNormalizationResult(
            (), "unresolved", 0.0, correspondence.warnings, False
        )
    all_candidate_results = tuple(
        (
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
    if evidence == "global_atom_correspondence":
        nonempty_costs = tuple(
            _correspondence_edit_cost(edits)
            for edits, _ in all_candidate_results
            if edits
        )
        best_cost = min(nonempty_costs) if nonempty_costs else None
        candidate_results = tuple(
            result
            for result in all_candidate_results
            for edits, _ in (result,)
            if edits and _correspondence_edit_cost(edits) == best_cost
        )
    else:
        candidate_results = all_candidate_results
    nonempty = tuple(result for result in candidate_results if result[0])
    if not nonempty:
        return EditNormalizationResult(
            (),
            "unresolved",
            0.0,
            ("SCAFFOLD_CORRESPONDENCE_WITHOUT_EDITS",),
            False,
        )
    edit_sets = {
        (
            tuple(sorted(_chemistry_edit_key(edit) for edit in edits)),
            tuple(
                sorted(_chemistry_stereo_key(change) for change in stereo)
            ),
        )
        for edits, stereo in nonempty
    }
    if len(edit_sets) != 1 or len(nonempty) != len(candidate_results):
        return EditNormalizationResult(
            (),
            "ambiguous_atom_correspondence",
            0.0,
            (f"AMBIGUOUS_SCAFFOLD_CORRESPONDENCE:{len(edit_sets)}",),
            False,
        )
    selected = min(
        nonempty,
        key=lambda result: (
            tuple(sorted(_comparison_key(edit) for edit in result[0])),
            tuple(
                sorted(
                    _stereo_comparison_key(change) for change in result[1]
                )
            ),
        ),
    )
    return EditNormalizationResult(
        selected[0],
        evidence,
        confidence,
        (inferred_warning,),
        True,
        selected[1],
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


def _canonical_stereo_pair(molecule: Any) -> Tuple[str, str]:
    from rdkit import Chem

    copy = Chem.Mol(molecule)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return (
        Chem.MolToSmiles(copy, canonical=True, isomericSmiles=False),
        Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True),
    )


def _stereochemical_reconstruction_conflict(
    candidates: Tuple[ReactionCandidate, ...],
    products: Tuple[ReactionComponent, ...],
) -> bool:
    """Detect a structurally matching but explicitly opposite prediction."""
    product_molecules = tuple(
        molecule
        for component in products
        for molecule in (parse_smiles(component.input_smiles),)
        if molecule is not None and molecule.GetNumHeavyAtoms() > 0
    )
    if len(product_molecules) != 1:
        return False
    observed = product_molecules[0]
    if not _has_explicit_stereochemistry(observed):
        return False
    observed_nonisomeric, observed_isomeric = _canonical_stereo_pair(observed)
    for candidate in candidates:
        predicted = parse_smiles(candidate.predicted_product_smiles or "")
        if predicted is None or not _has_explicit_stereochemistry(predicted):
            continue
        predicted_nonisomeric, predicted_isomeric = _canonical_stereo_pair(
            predicted
        )
        if (
            predicted_nonisomeric == observed_nonisomeric
            and predicted_isomeric != observed_isomeric
        ):
            return True
    return False


def normalize_reaction_edits(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    selected: Optional[ReactionCandidate],
    selected_events: Tuple[ReactionCandidate, ...] = (),
    candidates: Tuple[ReactionCandidate, ...] = (),
) -> EditNormalizationResult:
    """Choose observed edits when valid and reconcile exact operator evidence."""
    mapped = normalize_mapped_edits(reactants, products)
    predicted = normalize_predicted_edits(selected, reactants)
    predicted_multi = normalize_predicted_multi_event_edits(
        selected_events, reactants
    )
    warnings = tuple(
        sorted(
            set(mapped.warnings + predicted.warnings + predicted_multi.warnings)
        )
    )
    if mapped.evidence == "invalid_atom_mapping":
        return EditNormalizationResult(
            (),
            "unresolved",
            0.0,
            warnings,
            False,
        )
    if mapped.valid and predicted.valid:
        mapped_keys = {_comparison_key(edit) for edit in mapped.edits}
        predicted_keys = {_comparison_key(edit) for edit in predicted.edits}
        if mapped_keys == predicted_keys:
            return EditNormalizationResult(
                mapped.edits,
                "validated_mapping_and_exact_reconstruction",
                1.0,
                warnings,
                True,
                mapped.stereo_changes,
            )
        return EditNormalizationResult(
            mapped.edits,
            "conflicting_edit_evidence",
            0.5,
            tuple(sorted(set(warnings + ("MAPPING_RECONSTRUCTION_CONFLICT",)))),
            True,
            mapped.stereo_changes,
        )
    if mapped.valid and predicted_multi.valid:
        mapped_keys = {_comparison_key(edit) for edit in mapped.edits}
        predicted_keys = {_comparison_key(edit) for edit in predicted_multi.edits}
        if mapped_keys == predicted_keys:
            return EditNormalizationResult(
                mapped.edits,
                "validated_mapping_and_exact_multi_event_reconstruction",
                1.0,
                warnings,
                True,
                mapped.stereo_changes,
            )
        return EditNormalizationResult(
            mapped.edits,
            "conflicting_edit_evidence",
            0.5,
            tuple(sorted(set(warnings + ("MAPPING_RECONSTRUCTION_CONFLICT",)))),
            True,
            mapped.stereo_changes,
        )
    if mapped.valid:
        if _stereochemical_reconstruction_conflict(candidates, products):
            return EditNormalizationResult(
                mapped.edits,
                "conflicting_stereochemical_evidence",
                0.5,
                tuple(
                    sorted(
                        set(
                            warnings
                            + ("STEREOCHEMICAL_RECONSTRUCTION_CONFLICT",)
                        )
                    )
                ),
                True,
                mapped.stereo_changes,
            )
        return EditNormalizationResult(
            mapped.edits,
            mapped.evidence,
            mapped.confidence,
            warnings,
            True,
            mapped.stereo_changes,
        )
    if predicted.valid:
        return EditNormalizationResult(
            predicted.edits, predicted.evidence, predicted.confidence, warnings
        )
    if predicted_multi.valid:
        return EditNormalizationResult(
            predicted_multi.edits,
            predicted_multi.evidence,
            predicted_multi.confidence,
            warnings,
        )
    inferred = normalize_inferred_scaffold_edits(reactants, products)
    if inferred.valid:
        if _stereochemical_reconstruction_conflict(candidates, products):
            return EditNormalizationResult(
                inferred.edits,
                "conflicting_stereochemical_evidence",
                0.5,
                tuple(
                    sorted(
                        set(
                            inferred.warnings
                            + ("STEREOCHEMICAL_RECONSTRUCTION_CONFLICT",)
                        )
                    )
                ),
                True,
                inferred.stereo_changes,
            )
        return inferred
    return EditNormalizationResult(
        (),
        inferred.evidence,
        0.0,
        tuple(sorted(set(warnings + inferred.warnings))),
        False,
    )


__all__ = [
    "EditNormalizationResult",
    "normalize_mapped_edits",
    "normalize_inferred_scaffold_edits",
    "normalize_predicted_edits",
    "normalize_predicted_multi_event_edits",
    "normalize_reaction_edits",
    "reaction_atom_reference",
]
