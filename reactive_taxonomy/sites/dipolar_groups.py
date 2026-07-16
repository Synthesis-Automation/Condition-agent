"""Multi-atom dipolar reactive-group detection."""

from __future__ import annotations

from typing import Any, List

from ..models import SiteCandidate
from ..patterns import MatchIndex


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    """Return carbon- or silicon-attached organic azides with typed N roles."""
    candidate_proximal = match_index.role_atoms(
        "dipolar_group", "proximal_nitrogen"
    )
    candidate_central = match_index.role_atoms("dipolar_group", "central_nitrogen")
    candidate_terminal = match_index.role_atoms(
        "dipolar_group", "terminal_nitrogen"
    )
    candidate_attachments = match_index.role_atoms("dipolar_group", "attachment")
    sites: List[SiteCandidate] = []
    for central in mol.GetAtoms():
        if central.GetIdx() not in candidate_central:
            continue
        nitrogen_neighbors = [
            neighbor for neighbor in central.GetNeighbors() if neighbor.GetAtomicNum() == 7
        ]
        if len(nitrogen_neighbors) != 2:
            continue
        for proximal in nitrogen_neighbors:
            if proximal.GetIdx() not in candidate_proximal:
                continue
            attachments = [
                neighbor
                for neighbor in proximal.GetNeighbors()
                if neighbor.GetIdx() != central.GetIdx()
                and neighbor.GetIdx() in candidate_attachments
                and neighbor.GetAtomicNum() in {6, 14}
            ]
            terminal = next(
                (
                    neighbor
                    for neighbor in nitrogen_neighbors
                    if neighbor.GetIdx() != proximal.GetIdx()
                    and neighbor.GetIdx() in candidate_terminal
                ),
                None,
            )
            if len(attachments) != 1 or terminal is None:
                continue
            attachment = attachments[0]
            atoms = (
                attachment.GetIdx(),
                proximal.GetIdx(),
                central.GetIdx(),
                terminal.GetIdx(),
            )
            bonds = tuple(
                int(mol.GetBondBetweenAtoms(left, right).GetIdx())
                for left, right in zip(atoms, atoms[1:])
            )
            patterns = tuple(
                definition["id"]
                for definition in match_index.patterns_for_atom(
                    "dipolar_group", "central_nitrogen", central.GetIdx()
                )
            )
            sites.append(
                SiteCandidate(
                    site_type="dipolar_group",
                    topology="edge",
                    atom_roles={
                        "attachment": (attachment.GetIdx(),),
                        "proximal_nitrogen": (proximal.GetIdx(),),
                        "central_nitrogen": (central.GetIdx(),),
                        "terminal_nitrogen": (terminal.GetIdx(),),
                        "group": (
                            proximal.GetIdx(),
                            central.GetIdx(),
                            terminal.GetIdx(),
                        ),
                    },
                    atom_indices=atoms,
                    bond_indices=bonds,
                    canonical_signature="DG|Azide|Organic",
                    render_kind="named_handle",
                    render_data={"template_id": "organic_azide"},
                    matched_patterns=patterns,
                    details={
                        "handle_token": "Azide",
                        "attachment_element": attachment.GetSymbol(),
                        "net_group_charge": sum(
                            atom.GetFormalCharge()
                            for atom in (proximal, central, terminal)
                        ),
                        "reaction_modes": ["cycloaddition", "reduction"],
                        "resonance_normalized": True,
                    },
                    availability="available",
                )
            )
            break
    return sites


__all__ = ["detect"]
