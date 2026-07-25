"""Verified product-connection construction from selected reactant roles."""

from __future__ import annotations

from .reaction_labels import load_reaction_rendering, render_product_label
from .reaction_models import (
    ProductConnection,
    ProductConnectionEndpoint,
    ReactionCandidate,
)


def build_product_connection(
    selected: ReactionCandidate | None,
    evidence_quality: str,
    *,
    style: str = "unicode",
) -> ProductConnection | None:
    """Build a role-preserving connection for verified two-anchor C–C joins."""
    if selected is None or evidence_quality != "exact_product_reconstruction":
        return None
    rule = load_reaction_rendering().get(selected.grammar_id) or {}
    kind = rule.get("product_kind")
    if kind == "join_contexts":
        left_role, right_role = str(rule["left_role"]), str(rule["right_role"])
        right_atom_role = "anchor"
        right_context_key = "anchor_context"
        connection_type = "C_C"
    elif kind == "nitrogen_substitution":
        left_role, right_role = str(rule["anchor_role"]), str(rule["nitrogen_role"])
        right_atom_role = "center"
        right_context_key = "center_token"
        connection_type = "C_N"
    elif kind == "heteroatom_substitution":
        left_role, right_role = str(rule["anchor_role"]), str(rule["partner_role"])
        right_atom_role = "center"
        right_context_key = "center_token"
        connection_type = f"C_{rule['element']}"
    elif kind == "activated_carbon_substitution":
        left_role, right_role = str(rule["anchor_role"]), str(rule["carbon_role"])
        right_atom_role = "center"
        right_context_key = "center_token"
        connection_type = "C_C"
    else:
        return None
    left, right = (
        selected.role_assignments[left_role],
        selected.role_assignments[right_role],
    )
    left_context = str(left.details.get("anchor_context") or "Other")
    right_context = str(right.details.get(right_context_key) or "N")
    concise_label = render_product_label(
        {"id": selected.grammar_id},
        selected.role_assignments,
        style=style,
    )
    return ProductConnection(
        endpoint_1=ProductConnectionEndpoint(
            role=left_role,
            component_index=left.component_index,
            atom_index=int(left.atom_roles["anchor"][0]),
            context=left_context,
            source_site_id=left.site_id,
        ),
        endpoint_2=ProductConnectionEndpoint(
            role=right_role,
            component_index=right.component_index,
            atom_index=int(right.atom_roles[right_atom_role][0]),
            context=right_context,
            source_site_id=right.site_id,
        ),
        bond_order="single",
        connection_type=connection_type,
        concise_label=concise_label,
        evidence=evidence_quality,
    )


__all__ = ["build_product_connection"]
