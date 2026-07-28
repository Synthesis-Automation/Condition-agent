"""Curated two-addend donors for addition across unsaturated bonds."""

from __future__ import annotations

from typing import Any, List

from ..models import SiteCandidate
from ..patterns import MatchIndex


def _bond_order_name(bond: Any) -> str:
    value = int(round(float(bond.GetBondTypeAsDouble())))
    return {1: "SINGLE", 2: "DOUBLE", 3: "TRIPLE"}.get(value, "UNKNOWN")


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    """Return curated explicit A-B and implicit A-H addition donors."""
    sites: List[SiteCandidate] = []
    seen: set[tuple[object, ...]] = set()
    candidate_atoms = match_index.role_atoms("addition_donor", "addend_a")
    for atom_index in sorted(candidate_atoms):
        for definition, roles in match_index.role_matches_for_atom(
            "addition_donor", "addend_a", atom_index
        ):
            source_kind = str(definition.get("source_kind") or "")
            donor_class = str(definition.get("donor_class") or "addition_donor")
            pattern_id = str(definition["id"])
            addend_a_indices = roles.get("addend_a") or ()
            if len(addend_a_indices) != 1:
                continue
            addend_a = int(addend_a_indices[0])
            atom_a = mol.GetAtomWithIdx(addend_a)
            if source_kind == "implicit_hydrogen":
                hydrogen_count = int(atom_a.GetTotalNumHs(includeNeighbors=True))
                if hydrogen_count < 1:
                    continue
                key = (pattern_id, addend_a, "H")
                if key in seen:
                    continue
                seen.add(key)
                token = str((definition.get("tokens") or ["XH"])[0])
                sites.append(
                    SiteCandidate(
                        site_type="addition_donor",
                        topology="atom",
                        atom_roles={
                            "addend_a": (addend_a,),
                            "hydrogen_carrier": (addend_a,),
                        },
                        atom_indices=(addend_a,),
                        bond_indices=(),
                        canonical_signature=(
                            f"AD|ATOM_HYDROGEN|{atom_a.GetSymbol()}-H|H{hydrogen_count}"
                        ),
                        render_kind="named_handle",
                        render_data={"template_id": "addition_donor"},
                        matched_patterns=(pattern_id,),
                        details={
                            "handle_token": token,
                            "source_kind": source_kind,
                            "donor_class": donor_class,
                            "addend_elements": [atom_a.GetSymbol(), "H"],
                            "source_bond_order": "SINGLE",
                            "hydrogen_count": hydrogen_count,
                        },
                        availability="transferable",
                    )
                )
                continue

            addend_b_indices = roles.get("addend_b") or ()
            if source_kind != "explicit_bond" or len(addend_b_indices) != 1:
                continue
            addend_b = int(addend_b_indices[0])
            bond = mol.GetBondBetweenAtoms(addend_a, addend_b)
            if bond is None:
                continue
            elements = (atom_a.GetSymbol(), mol.GetAtomWithIdx(addend_b).GetSymbol())
            key = (pattern_id, min(addend_a, addend_b), max(addend_a, addend_b))
            if key in seen:
                continue
            seen.add(key)
            token = str((definition.get("tokens") or ["AddendPair"])[0])
            sites.append(
                SiteCandidate(
                    site_type="addition_donor",
                    topology="bond",
                    atom_roles={
                        "addend_a": (addend_a,),
                        "addend_b": (addend_b,),
                    },
                    atom_indices=(addend_a, addend_b),
                    bond_indices=(int(bond.GetIdx()),),
                    canonical_signature=(
                        f"AD|ATOM_ATOM|{elements[0]}-{elements[1]}|"
                        f"{_bond_order_name(bond)}"
                    ),
                    render_kind="named_handle",
                    render_data={"template_id": "addition_donor"},
                    matched_patterns=(pattern_id,),
                    details={
                        "handle_token": token,
                        "source_kind": source_kind,
                        "donor_class": donor_class,
                        "addend_elements": list(elements),
                        "source_bond_order": _bond_order_name(bond),
                    },
                    availability="transferable",
                )
            )
    return sites


__all__ = ["detect"]
