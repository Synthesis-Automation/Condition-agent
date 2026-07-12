"""Verified product-connection construction from selected reactant roles."""

from __future__ import annotations

from .reaction_labels import load_reaction_rendering, render_context_connection
from .reaction_models import ProductConnection, ProductConnectionEndpoint, ReactionCandidate


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
    if rule.get("product_kind") != "join_contexts":
        return None
    left_role, right_role = str(rule["left_role"]), str(rule["right_role"])
    left, right = selected.role_assignments[left_role], selected.role_assignments[right_role]
    left_context = str(left.details.get("anchor_context") or "Other")
    right_context = str(right.details.get("anchor_context") or "Other")
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
            atom_index=int(right.atom_roles["anchor"][0]),
            context=right_context,
            source_site_id=right.site_id,
        ),
        bond_order="single",
        concise_label=render_context_connection(left_context, right_context, style=style),
        evidence=evidence_quality,
    )


__all__ = ["build_product_connection"]
