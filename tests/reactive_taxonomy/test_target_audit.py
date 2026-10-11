"""Regression tests for the compact deterministic target audit."""

from reactive_taxonomy import audit_target
import pytest
from rdkit import Chem

from reactive_taxonomy.molecule_inspection import inspect_reactive_sites
from reactive_taxonomy.validation import validate_taxonomy

TARGET = "Cn(cc(c1ncc(CC)cc1)c2c3CCC2)c3=O"


@pytest.mark.parametrize("smiles", [TARGET, "O=c1cccc[nH]1", "Cn1ccccc1=O", "CN1C=CC=CC1=O"])
def test_pyridone_recognition_preserves_graph_evidence_without_reactivity_tags(smiles: str) -> None:
    audit = audit_target(smiles)
    motifs = {item.motif_id: item for item in audit.motifs}
    assert {"pyridin_2_one", "aromatic_lactam"} <= motifs.keys()
    assert motifs["pyridin_2_one"].tags == motifs["aromatic_lactam"].tags == ()
    assert "amide" not in motifs
    assert audit.carbonyl_nitrogen_observations[0].cn_bond_in_ring


@pytest.mark.parametrize("smiles", ["n1ccccc1", "Oc1ccccn1", "O=c1cc[nH]cc1", "O=c1cccco1", "CC(=O)NC", "O=C1CCCN1"])
def test_pyridone_subtype_does_not_absorb_other_heterocycles_or_amides(smiles: str) -> None:
    audit = audit_target(smiles)
    assert "pyridin_2_one" not in {m.motif_id for m in audit.motifs}
    if smiles in {"CC(=O)NC", "O=C1CCCN1"}:
        assert "amide" in {m.motif_id for m in audit.motifs}
        assert audit.carbonyl_nitrogen_observations


def test_tautomers_are_not_merged_and_equivalent_serializations_are_consistent() -> None:
    aromatic, kekule, enol = (audit_target(s) for s in ("O=c1cccc[nH]1", "O=C1C=CC=CN1", "Oc1ccccn1"))
    assert aromatic.canonical_smiles == kekule.canonical_smiles
    assert aromatic.motifs == kekule.motifs
    assert aromatic.canonical_smiles != enol.canonical_smiles
    assert not enol.carbonyl_nitrogen_observations


@pytest.mark.parametrize("smiles", [TARGET, "[Na+].O=C([O-])CCBr", "[CH3:7][C@H:2](O)Cl", "BrCC(=O)O.[Na+]"])
def test_audit_atom_projection_matches_graph_components_and_focused_inspection(smiles: str) -> None:
    audit = audit_target(smiles)
    original, canonical = Chem.MolFromSmiles(smiles), Chem.MolFromSmiles(audit.canonical_smiles)
    assert len(audit.canonical_atom_to_input_atom) == canonical.GetNumAtoms()
    for index, source in enumerate(audit.canonical_atom_to_input_atom):
        a, b = canonical.GetAtomWithIdx(index), original.GetAtomWithIdx(source)
        assert (a.GetSymbol(), a.GetFormalCharge(), a.GetIsotope()) == (b.GetSymbol(), b.GetFormalCharge(), b.GetIsotope())
        assert audit.input_atom_map_numbers[source] == b.GetAtomMapNum()
    for bond in canonical.GetBonds():
        a, b = (audit.canonical_atom_to_input_atom[i] for i in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        assert original.GetBondBetweenAtoms(a, b).GetBondType() == bond.GetBondType()
    for motif in audit.motifs:
        component = audit.components[motif.component_index]
        ids = [component.canonical_atom_indices[i] for i in motif.atom_indices]
        assert all(0 <= i < canonical.GetNumAtoms() for i in ids)
    for site in audit.reactive_sites:
        component = audit.components[site.component_index]
        for index in site.bond_indices:
            bond = canonical.GetBondWithIdx(component.canonical_bond_indices[index])
            assert bond.GetBeginAtomIdx() in component.canonical_atom_indices
    if audit.component_count == 1:
        inspection = inspect_reactive_sites(smiles)
        assert inspection.molecule.canonical_smiles == audit.canonical_smiles
        assert {(m.motif_id, m.atom_indices) for m in inspection.motifs} == {(m.motif_id, m.atom_indices) for m in audit.motifs}


def test_planning_audit_reports_coverage_versions_and_raw_observations_even_without_motifs(monkeypatch) -> None:
    from reactive_taxonomy import api
    monkeypatch.setattr(api, "detect_molecular_motifs", lambda *a, **kw: [])
    audit = audit_target(TARGET)
    assert audit.molecular_formula == "C16H18N2O"
    assert audit.carbonyl_nitrogen_observations
    assert audit.rings and audit.definition_versions and audit.rdkit_version
    assert audit.coverage == "registered_patterns_only_not_exhaustive"
    assert not validate_taxonomy()


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
