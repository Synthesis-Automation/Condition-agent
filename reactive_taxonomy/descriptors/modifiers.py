"""Orthogonal handle and liability modifiers for reactivity profiles."""

from __future__ import annotations

from typing import Any, Iterable, Tuple

from .models import DescriptorEvidence, ReactivityModifier


def _role_indices(site: Any, role: str) -> Tuple[int, ...]:
    roles = site.details.get("atom_roles") or {}
    value = roles.get(role) or ()
    if isinstance(value, int):
        return (int(value),)
    return tuple(int(index) for index in value)


def build_reactivity_modifiers(
    mol: Any,
    site: Any,
    center: int,
    *,
    beta_hydrogen_count: int | None = None,
    ring_sizes: Iterable[int] = (),
    nearby_groups: Iterable[dict[str, Any]] = (),
) -> Tuple[ReactivityModifier, ...]:
    """Build deterministic site-handle and broadly valid liability modifiers."""
    values = []
    details = site.details or {}
    handle_token = str(
        details.get("handle_token")
        or details.get("leaving_group")
        or details.get("carrier")
        or ""
    )
    if site.site_type == "leaving_group":
        handle_atoms = _role_indices(site, "handle")
        handle_elements = tuple(
            sorted(mol.GetAtomWithIdx(index).GetSymbol() for index in handle_atoms)
        )
        modifier_id = handle_token or "-".join(handle_elements) or "leaving_group"
        values.append(
            ReactivityModifier(
                modifier_type="leaving_group",
                modifier_id=modifier_id,
                class_name=modifier_id,
                attributes=tuple(
                    sorted(
                        (
                            ("element_set", ",".join(handle_elements)),
                            ("site_type", str(site.site_type)),
                        )
                    )
                ),
                evidence=DescriptorEvidence(
                    source="reactive_site",
                    method="site_handle_modifier_v1",
                    confidence=float(site.confidence),
                    contributing_atom_indices=handle_atoms,
                ),
            )
        )
    if site.site_type in {"transfer_group", "addition_donor"}:
        handle_atoms = _role_indices(site, "handle") or _role_indices(
            site, "center"
        )
        elements = tuple(
            sorted(mol.GetAtomWithIdx(index).GetSymbol() for index in handle_atoms)
        )
        modifier_id = handle_token or "-".join(elements) or str(site.site_type)
        values.append(
            ReactivityModifier(
                modifier_type="transfer_carrier",
                modifier_id=modifier_id,
                class_name=str(site.site_type),
                attributes=tuple(
                    sorted((("element_set", ",".join(elements)),))
                ),
                evidence=DescriptorEvidence(
                    source="reactive_site",
                    method="site_transfer_modifier_v1",
                    confidence=float(site.confidence),
                    contributing_atom_indices=handle_atoms,
                ),
            )
        )
    if site.site_type in {"pronucleophile_XH", "aromatic_CH"}:
        values.append(
            ReactivityModifier(
                modifier_type="removable_hydrogen",
                modifier_id=f"{mol.GetAtomWithIdx(center).GetSymbol()}-H",
                class_name="explicit_or_schema_level_hydrogen",
                attributes=(("hydrogen_count", str(
                    mol.GetAtomWithIdx(center).GetTotalNumHs()
                )),),
                evidence=DescriptorEvidence(
                    source="reactive_site",
                    method="site_hydrogen_modifier_v1",
                    confidence=float(site.confidence),
                    contributing_atom_indices=(center,),
                ),
            )
        )
    if beta_hydrogen_count is not None:
        values.append(
            ReactivityModifier(
                modifier_type="elimination",
                modifier_id="beta_hydrogen",
                class_name=(
                    "beta_h_present" if beta_hydrogen_count > 0 else "no_beta_h"
                ),
                attributes=(("count", str(beta_hydrogen_count)),),
                evidence=DescriptorEvidence(
                    source="molecular_graph",
                    method="beta_hydrogen_count_v1",
                    confidence=1.0,
                    contributing_atom_indices=(center,),
                ),
            )
        )
    small_rings = tuple(sorted(size for size in ring_sizes if int(size) <= 4))
    if small_rings:
        values.append(
            ReactivityModifier(
                modifier_type="strain",
                modifier_id="small_ring",
                class_name="ring_strain",
                attributes=(("ring_sizes", ",".join(map(str, small_rings))),),
                evidence=DescriptorEvidence(
                    source="molecular_graph",
                    method="ring_size_liability_v1",
                    confidence=1.0,
                    contributing_atom_indices=(center,),
                ),
            )
        )
    coordinating = tuple(
        sorted(
            {
                str(group.get("group_id") or "")
                for group in nearby_groups
                if set(group.get("tags") or ())
                & {"strong_metal_binding", "palladium_poisoning"}
            }
            - {""}
        )
    )
    if coordinating:
        values.append(
            ReactivityModifier(
                modifier_type="coordination",
                modifier_id="nearby_coordination",
                class_name="potential_catalyst_coordination",
                attributes=(("group_ids", ",".join(coordinating)),),
                evidence=DescriptorEvidence(
                    source="functional_group_detection",
                    method="nearby_group_modifier_v1",
                    confidence=1.0,
                    contributing_atom_indices=(center,),
                ),
            )
        )
    return tuple(
        sorted(
            values,
            key=lambda item: (
                item.modifier_type,
                item.modifier_id,
                item.class_name,
            ),
        )
    )


__all__ = ["build_reactivity_modifiers"]
