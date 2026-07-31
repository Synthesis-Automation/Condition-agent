"""Optional structural annotations layered over the generic core.

Annotations may improve chemist-facing explanations, but the generic edit,
atom-state, and remote-subgraph observations remain usable without them.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

from ..reaction_models import ReactionEdit, ReactionTopology
from .common import Location as _Location, atom_identity as _atom_identity
from .keys import digest as _digest
from .models import ReactionCoreAbstraction


def _aromatic_center_class(molecule: Any, atom_index: int) -> str:
    """Classify the connected aromatic system around one active carbon."""
    remaining = {int(atom_index)}
    visited = set()
    while remaining:
        current = remaining.pop()
        if current in visited:
            continue
        atom = molecule.GetAtomWithIdx(current)
        if not atom.GetIsAromatic():
            continue
        visited.add(current)
        remaining.update(
            int(neighbor.GetIdx())
            for neighbor in atom.GetNeighbors()
            if neighbor.GetIsAromatic()
        )
    return (
        "heteroaryl"
        if any(
            molecule.GetAtomWithIdx(index).GetAtomicNum() != 6
            for index in visited
        )
        else "aryl"
    )


def _transfer_center_limiter(
    molecule: Any,
    atom_index: int,
    carboxyl_atom_index: int,
) -> tuple[str, str]:
    """Return a typed token and readable limiter for the transferred group."""
    atom = molecule.GetAtomWithIdx(atom_index)
    if atom.GetIsAromatic():
        center_class = _aromatic_center_class(molecule, atom_index)
        label = (
            "HetAr– (heteroaryl)"
            if center_class == "heteroaryl"
            else "Ar– (aryl)"
        )
        return center_class, label
    if any(
        neighbor.GetAtomicNum() == 8
        and str(
            molecule.GetBondBetweenAtoms(
                atom_index,
                int(neighbor.GetIdx()),
            ).GetBondType()
        ).upper()
        == "DOUBLE"
        for neighbor in atom.GetNeighbors()
    ):
        return "acyl", "R′–C(=O)– (acyl)"
    hybridization = str(atom.GetHybridization()).upper()
    if hybridization == "SP":
        return "alkynyl", "R′–C≡C– (alkynyl)"
    if hybridization == "SP2":
        return "alkenyl", "R′–CH=CH– (alkenyl)"
    hydrogens = int(atom.GetTotalNumHs(includeNeighbors=True))
    carbon_neighbors = sum(
        neighbor.GetAtomicNum() == 6
        and int(neighbor.GetIdx()) != int(carboxyl_atom_index)
        for neighbor in atom.GetNeighbors()
    )
    if hydrogens >= 3 and carbon_neighbors == 0:
        return "methyl", "CH₃– (methyl)"
    if hydrogens == 2 and carbon_neighbors == 1:
        return "primary_alkyl", "R′–CH₂– (primary alkyl)"
    if hydrogens == 1 and carbon_neighbors == 2:
        return "secondary_alkyl", "R′R″CH– (secondary alkyl)"
    if hydrogens == 0 and carbon_neighbors >= 3:
        return "tertiary_alkyl", "R′R″R‴C– (tertiary alkyl)"
    return "alkyl_other", "alkyl group"


def decarboxylative_abstraction(
    *,
    edits: Sequence[ReactionEdit],
    reactant_by_map: Mapping[int, _Location],
    product_by_map: Mapping[int, _Location],
    topology: Optional[ReactionTopology],
) -> Optional[ReactionCoreAbstraction]:
    """Recognize C–C formation coupled to loss of a carboxylic-acid carbon."""
    formed = tuple(
        edit
        for edit in edits
        if edit.edit_type == "formed" and edit.atom_2 is not None
    )
    for broken in edits:
        if broken.edit_type != "broken" or broken.atom_2 is None:
            continue
        endpoints = (broken.atom_1, broken.atom_2)
        for carboxyl, transfer in (endpoints, endpoints[::-1]):
            if carboxyl.element != "C" or transfer.element != "C":
                continue
            carboxyl_identity = _atom_identity(carboxyl)
            transfer_identity = _atom_identity(transfer)
            if carboxyl_identity[0] != "map" or transfer_identity[0] != "map":
                continue
            carboxyl_map = int(carboxyl_identity[1])
            transfer_map = int(transfer_identity[1])
            carboxyl_location = reactant_by_map.get(carboxyl_map)
            transfer_location = reactant_by_map.get(transfer_map)
            if (
                carboxyl_location is None
                or transfer_location is None
                or carboxyl_map in product_by_map
                or transfer_map not in product_by_map
            ):
                continue
            component, molecule, carboxyl_index = carboxyl_location
            if not any(
                str(group.group_id) == "carboxylic_acid"
                and int(carboxyl_index) in set(group.atom_indices)
                for group in component.compound_analysis.functional_groups
            ):
                continue
            partner_reference = None
            for edit in formed:
                formed_endpoints = (edit.atom_1, edit.atom_2)
                formed_identities = tuple(
                    _atom_identity(endpoint) for endpoint in formed_endpoints
                )
                if transfer_identity not in formed_identities:
                    continue
                partner_reference = formed_endpoints[
                    1 if formed_identities[0] == transfer_identity else 0
                ]
                if partner_reference.element == "C":
                    break
                partner_reference = None
            if partner_reference is None:
                continue
            transfer_component, transfer_molecule, transfer_index = (
                transfer_location
            )
            if transfer_component.component_index != component.component_index:
                continue
            transfer_token, transfer_label = _transfer_center_limiter(
                transfer_molecule,
                transfer_index,
                carboxyl_index,
            )
            partner_identity = _atom_identity(partner_reference)
            partner_location = (
                reactant_by_map.get(int(partner_identity[1]))
                if partner_identity[0] == "map"
                else None
            )
            partner_token = "carbon"
            partner_label = "C"
            partner_component_index = None
            if partner_location is not None:
                partner_component, partner_molecule, partner_index = (
                    partner_location
                )
                partner_component_index = partner_component.component_index
                partner_atom = partner_molecule.GetAtomWithIdx(partner_index)
                if partner_atom.GetIsAromatic():
                    partner_token = _aromatic_center_class(
                        partner_molecule,
                        partner_index,
                    )
                    partner_label = (
                        "HetAr" if partner_token == "heteroaryl" else "Ar"
                    )
            same_component = (
                partner_component_index is not None
                and partner_component_index
                == transfer_component.component_index
            )
            ring_sizes = (
                tuple(topology.formed_ring_sizes)
                if topology is not None
                else ()
            )
            is_cyclization = bool(
                same_component
                and topology is not None
                and topology.reaction_scope in {"intramolecular", "mixed"}
                and (
                    ring_sizes
                    or (
                        topology.ring_count_delta is not None
                        and topology.ring_count_delta > 0
                    )
                )
            )
            if is_cyclization and len(ring_sizes) == 1:
                general_label = (
                    "R–C(=O)OH + Ar–H → R–Ar; intramolecular, "
                    f"{ring_sizes[0]}-membered ring"
                )
            elif is_cyclization:
                general_label = (
                    "R–C(=O)OH + Ar–H → R–Ar; "
                    "intramolecular cyclization"
                )
            elif same_component:
                general_label = (
                    "R–C(=O)OH + Ar–H → R–Ar; intramolecular"
                )
            else:
                general_label = "R–C(=O)OH + Ar–H → R–Ar"
            motif_tokens = (
                "bond_formed:C-C",
                "departing_handle:carboxylic_acid",
                "motif:decarboxylative_coupling",
            )
            limiter_tokens = tuple(
                sorted(
                    (
                        f"partner_center:{partner_token}",
                        f"transfer_center:{transfer_token}",
                        (
                            "topology:intramolecular_cyclization"
                            if is_cyclization
                            else "topology:intramolecular"
                            if same_component
                            else "topology:intermolecular"
                        ),
                        *(
                            (f"ring_size:{ring_sizes[0]}",)
                            if is_cyclization and len(ring_sizes) == 1
                            else ()
                        ),
                    )
                )
            )
            return ReactionCoreAbstraction(
                motif_id="decarboxylative_c_c_coupling",
                motif_key=_digest("RCM1", motif_tokens),
                general_label=general_label,
                limiter_label=(
                    (
                        f"acyl center = {transfer_label}; "
                        f"partner center = {partner_label}"
                    )
                    if same_component
                    else f"R = {transfer_label}; Ar = {partner_label}"
                ),
                motif_tokens=motif_tokens,
                limiter_tokens=limiter_tokens,
            )
    return None


__all__ = ["decarboxylative_abstraction"]

