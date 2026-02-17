from __future__ import annotations

from chemtools.taxonomy.motif_scope_index import build_motif_scope_index


def test_build_motif_scope_index_generates_expected_parent_child_edges() -> None:
    payload = build_motif_scope_index()
    scope_map = payload.get("scope_map", {}) or {}

    alkyl_children = set(scope_map.get("Alkyl-B(OH)2") or [])
    assert {"RCH2-B(OH)2", "R2CH-B(OH)2", "R3C-B(OH)2"} <= alkyl_children

    ar_children = set(scope_map.get("Ar-Ar") or [])
    assert "HeteroAr-Ar" in ar_children
