"""Deterministic structural-core observation regressions."""

from __future__ import annotations

from copy import deepcopy
import json

from reactive_taxonomy import (
    load_structural_core_matching_definition,
    load_structural_core_observation_definition,
    observe_structural_cores,
    render_structural_core_review_html,
    validate_structural_core_matching_definition,
    validate_structural_core_observation_definition,
    validate_taxonomy,
    write_structural_core_review_html,
)
from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles


def test_structural_core_definitions_are_versioned_and_valid() -> None:
    observation = load_structural_core_observation_definition()
    matching = load_structural_core_matching_definition()

    assert observation["definition_id"] == "structural_core_observations.v1"
    assert observation["definition_version"] == "1.1.0"
    assert matching["definition_id"] == "structural_core_matching.v1"
    assert observation["maximum_observations"] == 5
    assert not validate_structural_core_observation_definition(observation)
    assert not validate_structural_core_matching_definition(matching)
    assert not validate_taxonomy()


def test_structural_core_definition_validation_rejects_policy_drift() -> None:
    observation = deepcopy(load_structural_core_observation_definition())
    matching = deepcopy(load_structural_core_matching_definition())
    observation["maximum_observations"] = 0
    matching["shape_atom_fields"] = ["element"]

    assert (
        "structural_core_observations:invalid_maximum_observations"
        in validate_structural_core_observation_definition(observation)
    )
    assert (
        "structural_core_matching:invalid_shape_atom_fields"
        in validate_structural_core_matching_definition(matching)
    )


def test_simple_ring_alone_abstains_instead_of_claiming_a_useful_core() -> None:
    analysis = observe_structural_cores("c1ccccc1")

    assert analysis.valid is True
    assert analysis.observations == ()
    assert "NO_STRUCTURAL_CORE_OBSERVATION" in analysis.warnings


def test_biphenyl_produces_focused_bridge_and_linker_observations() -> None:
    analysis = observe_structural_cores("c1ccc(-c2ccccc2)cc1")

    assert analysis.valid is True
    assert [item.kind for item in analysis.observations] == [
        "bridge_interface",
        "linker_region",
    ]
    bridge, linker = analysis.observations
    assert bridge.focus_atom_indices == linker.focus_atom_indices
    assert bridge.focus_bond_keys == linker.focus_bond_keys
    assert len(linker.atom_indices) == 2
    assert "DIRECT_RING_TO_RING_BOND" in linker.evidence_codes
    assert dict(bridge.descriptor_values)["primary_component_fraction"] == 0.5


def test_murcko_scaffold_linker_and_carbon_framework_are_diverse() -> None:
    analysis = observe_structural_cores("Cc1cc(Oc2nccc(CCC)c2)ccc1")

    kinds = [item.kind for item in analysis.observations]

    assert kinds == [
        "scaffold_backbone",
        "bridge_interface",
        "linker_region",
        "carbon_framework",
        "ring_system",
    ]
    assert (
        "BEMIS_MURCKO_RING_LINKER_SCAFFOLD"
        in analysis.observations[0].evidence_codes
    )
    linker = next(item for item in analysis.observations if item.kind == "linker_region")
    assert "HETEROATOM_LINKER" in linker.evidence_codes
    assert len(linker.focus_bond_keys) == 2


def test_continuous_carbon_framework_is_observed_separately() -> None:
    analysis = observe_structural_cores("OCC(C)(C)CCO")

    framework = next(
        item for item in analysis.observations if item.kind == "carbon_framework"
    )

    assert "CONTINUOUS_CARBON_FRAMEWORK" in framework.evidence_codes
    assert "BRANCHED_CARBON_FRAMEWORK" in framework.evidence_codes
    assert dict(framework.descriptor_values)["carbon_atom_count"] == 6.0
    assert dict(framework.descriptor_values)["heteroatom_count"] == 0.0


def test_stereochemical_region_is_local_and_focuses_the_stereocenter() -> None:
    analysis = observe_structural_cores("CCCC[C@H](O)CCCCC")

    assert analysis.valid is True
    stereo = next(
        item
        for item in analysis.observations
        if item.kind == "stereo_backbone_region"
    )
    assert len(stereo.focus_atom_indices) == 1
    assert stereo.focus_atom_indices[0] in stereo.atom_indices
    assert len(stereo.atom_indices) == 4
    assert "LOCAL_STEREOCHEMICAL_BACKBONE" in stereo.evidence_codes
    assert dict(stereo.descriptor_values)["stereocenter_count"] == 1.0


def test_shape_key_relaxes_atom_type_but_typed_key_does_not() -> None:
    alcohol = observe_structural_cores("CC(C)(C)O").observations[0]
    amine = observe_structural_cores("CC(C)(C)N").observations[0]

    assert alcohol.structural_shape_key == amine.structural_shape_key
    assert alcohol.structural_typed_key != amine.structural_typed_key
    assert alcohol.structural_exact_key != amine.structural_exact_key


def test_serialization_and_component_order_do_not_change_observations() -> None:
    left = observe_structural_cores("c1ccc(-c2ccccc2)cc1.[Na+]")
    right = observe_structural_cores("[Na+].c1ccc(cc1)c2ccccc2")

    assert left.canonical_smiles == right.canonical_smiles
    assert left.molecule_id == right.molecule_id
    assert left.observations == right.observations
    assert "MULTICOMPONENT_TARGET_PRIMARY_COMPONENT_SELECTED" in left.warnings


def test_atom_references_resolve_against_the_normalized_target() -> None:
    analysis = observe_structural_cores("c1ccc(-c2ccccc2)cc1.[Na+]")
    molecule = parse_smiles(analysis.canonical_smiles or "")

    assert molecule is not None
    for observation in analysis.observations:
        for reference in observation.atom_references:
            atom = molecule.GetAtomWithIdx(reference.atom_index)
            assert atom.GetSymbol() == reference.element
            assert atom.GetFormalCharge() == reference.formal_charge
            assert atom.GetIsAromatic() is reference.aromatic
        assert set(observation.focus_atom_indices) <= set(observation.atom_indices)
        for key in observation.focus_bond_keys:
            _, left, right, *_ = key.split(":")
            assert molecule.GetBondBetweenAtoms(int(left), int(right)) is not None


def test_coverage_and_kind_quotas_bound_every_observation() -> None:
    definition = load_structural_core_observation_definition()
    analysis = observe_structural_cores("Cc1cc(Oc2nccc(CCC)c2)ccc1")

    kind_counts: dict[str, int] = {}
    for observation in analysis.observations:
        kind_counts[observation.kind] = kind_counts.get(observation.kind, 0) + 1
        coverage = dict(observation.descriptor_values)["primary_component_fraction"]
        assert coverage <= definition["maximum_coverage_by_kind"][observation.kind]
    assert all(
        count <= definition["maximum_per_kind"][kind]
        for kind, count in kind_counts.items()
    )


def test_unstructured_and_invalid_targets_abstain_cleanly() -> None:
    linear = observe_structural_cores("CCCC")
    invalid = observe_structural_cores("not smiles")

    assert linear.valid is True
    assert linear.observations == ()
    assert "NO_STRUCTURAL_CORE_OBSERVATION" in linear.warnings
    assert invalid.valid is False
    assert invalid.observations == ()
    assert invalid.error == "INVALID_TARGET_SMILES"


def test_analysis_serialization_is_json_compatible() -> None:
    analysis = observe_structural_cores("c1ccc(-c2ccccc2)cc1")

    encoded = json.dumps(analysis.to_dict(), sort_keys=True)

    assert analysis.molecule_id in encoded
    assert analysis.observations[0].core_observation_id in encoded


def test_review_renderer_highlights_observations_and_writes_packet(tmp_path) -> None:
    analyses = (
        observe_structural_cores("c1ccc(-c2ccccc2)cc1"),
        observe_structural_cores("not smiles"),
    )

    rendered = render_structural_core_review_html(
        analyses,
        title="Blind core review",
    )
    destination = write_structural_core_review_html(
        analyses,
        tmp_path / "structural-core-review.html",
        title="Blind core review",
    )

    assert "Blind core review" in rendered
    assert "<svg" in rendered
    assert analyses[0].observations[0].core_observation_id in rendered
    assert "not reaction atom maps" in rendered
    assert "Focus bonds:" in rendered
    assert rendered.count("<svg") == 1
    assert destination.read_text(encoding="utf-8") == rendered
