"""Regression tests for the compact deterministic target audit."""

from reactive_taxonomy import audit_target


def test_target_audit_preserves_assigned_and_unassigned_stereochemistry() -> None:
    assigned = audit_target("C[C@H](O)Cl")
    unassigned = audit_target("CC(O)Cl")

    assert assigned.valid is True
    assert assigned.canonical_smiles
    assert len(assigned.stereocenters) == 1
    assert assigned.stereocenters[0].assigned is True
    assert assigned.stereocenters[0].assignment in {"R", "S"}
    assert len(unassigned.stereocenters) == 1
    assert unassigned.stereocenters[0].assigned is False


def test_target_audit_exposes_graph_derived_motifs_and_sites() -> None:
    audit = audit_target("O=C(O)CCBr")

    assert audit.valid is True
    assert audit.component_count == 1
    assert audit.heavy_atom_count == 6
    assert any(item.site_type == "leaving_group" for item in audit.reactive_sites)
    assert all(item.hypothesis_id for item in audit.reactive_sites)


def test_target_audit_rejects_invalid_input_without_guessing() -> None:
    audit = audit_target("not-a-smiles")

    assert audit.valid is False
    assert audit.canonical_smiles is None
    assert audit.error == "INVALID_SMILES"
    assert audit.stereocenters == ()
