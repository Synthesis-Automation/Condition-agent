"""Role-neutral synthetic partition and landscape regressions."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import core_retrosynthesis.cli as cli_module
from core_retrosynthesis import (
    GenericDisconnectionCandidate,
    ModuleAnnotation,
    MoleculeOccurrenceNode,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    SYNTHETIC_PARTITION_POLICY_PATH,
    SyntheticPartition,
    SyntheticPartitionLandscape,
    analyze_partition_target,
    build_module_id,
    build_operator_partition_landscape,
    create_synthetic_partition,
    load_synthetic_partition_policy,
    project_reaction_to_target,
    project_route_partitions,
    render_partition_landscape_html,
    validate_synthetic_partition_policy,
    write_partition_landscape_review,
)


TARGET = "CNC(C)=O"
PROJECT_ROOT = Path(__file__).resolve().parents[2]
AMIDE_FORMATION = (
    "[CH3:1][C:2](=[O:3])[OH:6].[NH2:4][CH3:5]"
    ">>[CH3:1][C:2](=[O:3])[NH:4][CH3:5]"
)
N_METHYLATION = (
    "[CH3:1][C:2](=[O:3])[NH2:4].[CH3:5][Br:6]"
    ">>[CH3:1][C:2](=[O:3])[NH:4][CH3:5]"
)


def _candidate(
    name: str,
    reaction: str,
    *,
    score: float = 0.9,
    support: int = 2,
    context: float = 0.7,
    status: str = "verified_signature",
) -> GenericDisconnectionCandidate:
    precursors, _ = reaction.split(">>")
    return GenericDisconnectionCandidate(
        target_smiles=TARGET,
        precursor_smiles=precursors,
        proposed_reaction_smiles=reaction,
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="test",
        template_id=f"template:{name}",
        score=score,
        context_similarity=context,
        product_similarity=score,
        precursor_similarity=score,
        template_specificity=1.0,
        independent_reference_support=support,
        forward_validation_status=status,
        center_transition_key=f"center:{name}",
        disconnection_site_key=f"site:{name}",
        precedent_reaction_ids=(f"precedent:{name}",),
        operator_id=f"operator:{name}",
        realization_id=f"realization:{name}",
        operator_signature=f"signature:{name}",
        synthon_signature=f"synthon:{name}",
        condition_query_reaction_smiles=reaction,
    )


def _observed_two_step_tree() -> ReactionRouteTree:
    ammonia = MoleculeOccurrenceNode(
        occurrence_id="ammonia",
        smiles="N",
        depth=2,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    methyl_bromide = MoleculeOccurrenceNode(
        occurrence_id="methyl_bromide",
        smiles="CBr",
        depth=2,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    methylamine_reaction = RouteReactionNode(
        reaction_node_id="reaction:methylamine",
        step_id="step:methylamine",
        depth=2,
        reaction_smiles="[NH3:11].[CH3:12][Br:13]>>[NH2:11][CH3:12]",
        evidence=RouteStepEvidence(evidence_kind="observed"),
        children=(ammonia, methyl_bromide),
    )
    methylamine = MoleculeOccurrenceNode(
        occurrence_id="methylamine",
        smiles="CN",
        depth=1,
        terminal=False,
        terminal_evidence="expanded",
        unresolved_reason=None,
        reaction=methylamine_reaction,
    )
    acetic_acid = MoleculeOccurrenceNode(
        occurrence_id="acetic_acid",
        smiles="CC(=O)O",
        depth=1,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    root_reaction = RouteReactionNode(
        reaction_node_id="reaction:amide",
        step_id="step:amide",
        depth=1,
        reaction_smiles=AMIDE_FORMATION,
        evidence=RouteStepEvidence(evidence_kind="observed"),
        children=(acetic_acid, methylamine),
    )
    root = MoleculeOccurrenceNode(
        occurrence_id="target",
        smiles=TARGET,
        depth=0,
        terminal=False,
        terminal_evidence="expanded",
        unresolved_reason=None,
        reaction=root_reaction,
    )
    return ReactionRouteTree(
        tree_id="tree:partition-test",
        route_kind="observed",
        target_smiles=TARGET,
        root=root,
        reaction_count=2,
        maximum_depth=2,
        fingerprint_tokens=("amide", "methylamine"),
    )


def test_partition_identity_is_symmetric_and_role_neutral() -> None:
    canonical, target_id, atoms, _ = analyze_partition_target(TARGET)
    target_maps = tuple(atom.atom_map for atom in atoms)
    left = target_maps[:2]
    right = target_maps[2:]
    plain = create_synthetic_partition(
        canonical,
        (left, right),
        source_kind="structural_baseline",
        evidence_level="E0",
    )
    annotated = create_synthetic_partition(
        canonical,
        (right, left),
        source_kind="structural_baseline",
        evidence_level="E0",
        annotations=(
            ModuleAnnotation(
                module_id=build_module_id(target_id, left),
                label="optional label",
                proposed_role="complex_fragment",
                confidence=0.6,
            ),
        ),
    )

    assert plain.target_atom_maps == target_maps
    assert plain.partition_id == annotated.partition_id
    assert tuple(module.target_atom_maps for module in plain.modules) == (
        left,
        right,
    )
    assert SyntheticPartition.from_dict(annotated.to_dict()) == annotated


def test_partition_rejects_duplicate_or_incomplete_atom_ownership() -> None:
    with pytest.raises(ValueError, match="duplicates target atom ownership"):
        create_synthetic_partition(
            TARGET,
            ((1, 2, 3), (3, 4, 5)),
            source_kind="structural_baseline",
            evidence_level="E0",
        )

    with pytest.raises(ValueError, match="cover every target heavy atom"):
        create_synthetic_partition(
            TARGET,
            ((1, 2), (3, 4)),
            source_kind="structural_baseline",
            evidence_level="E0",
        )


def test_mapped_reaction_projection_separates_target_and_tactical_atoms() -> None:
    projection = project_reaction_to_target(TARGET, AMIDE_FORMATION)

    assert projection.k == 2
    assert sorted(map(len, projection.module_atom_sets)) == [2, 3]
    assert len(projection.target_boundary_bonds) == 1
    assert 6 not in {
        atom_map
        for module in projection.module_atom_sets
        for atom_map in module
    }
    assert projection.mapping_evidence == "supplied_atom_mapping"


def test_reaction_projection_reports_inferred_mapping_and_rejects_missing_source() -> None:
    inferred = project_reaction_to_target(
        TARGET,
        "CC(=O)O.CN>>CNC(C)=O",
    )

    assert inferred.k == 2
    assert inferred.mapping_evidence == "global_atom_correspondence"
    assert inferred.mapping_confidence < 1.0
    with pytest.raises(ValueError, match="do not cover current target ownership"):
        project_reaction_to_target(
            "CO",
            "[CH3:1]>>[CH3:1][OH:2]",
        )


def test_observed_route_projection_returns_k1_k2_and_k3_frontiers() -> None:
    projection = project_route_partitions(_observed_two_step_tree())

    assert projection.unresolved_occurrence_ids == ()
    assert [frontier.depth for frontier in projection.frontiers] == [0, 1, 2]
    assert [frontier.partition.k for frontier in projection.frontiers] == [1, 2, 3]
    assert all(
        frontier.partition.evidence_level == "E4"
        for frontier in projection.frontiers
    )
    assert len(projection.frontiers[-1].latent_states) == 3
    assert projection.to_landscape().searched_k_values == (1, 2, 3)


def test_operator_landscape_combines_independent_views_without_realizing_them() -> None:
    amide = _candidate("amide", AMIDE_FORMATION)
    methylation = _candidate("methylation", N_METHYLATION, score=0.85)
    tactical_variant = _candidate("amide-variant", AMIDE_FORMATION, score=0.8)

    landscape = build_operator_partition_landscape(
        TARGET,
        (amide, methylation, tactical_variant),
    )

    assert landscape.abstained is False
    assert {partition.k for partition in landscape.partitions} == {1, 2, 3}
    combined = next(partition for partition in landscape.partitions if partition.k == 3)
    assert combined.source_kind == "operator_combination_unrealized"
    assert combined.evidence_level == "E1"
    assert "INTERFACE_COMBINATION_NOT_JOINTLY_REALIZED" in combined.warnings
    amide_partition = next(
        partition
        for partition in landscape.partitions
        if any(
            "operator:amide" in interface.candidate_operator_ids
            for interface in partition.interfaces
        )
        and partition.k == 2
    )
    assert any(
        set(interface.candidate_operator_ids)
        == {"operator:amide", "operator:amide-variant"}
        for interface in amide_partition.interfaces
    )
    assert landscape.diagnostics.accepted_seed_count == 2


def test_landscape_abstains_without_a_partition_changing_interface() -> None:
    invalid = _candidate(
        "unverified",
        AMIDE_FORMATION,
        status="unresolved",
    )

    landscape = build_operator_partition_landscape(TARGET, (invalid,))

    assert landscape.abstained is True
    assert landscape.abstention_reasons == (
        "NO_VALIDATED_PARTITION_CHANGING_INTERFACE",
    )
    assert [partition.k for partition in landscape.partitions] == [1]
    assert landscape.diagnostics.rejection_counts == (
        ("candidate_not_forward_verified", 1),
    )


def test_policy_and_static_review_outputs_are_valid(tmp_path) -> None:
    policy = load_synthetic_partition_policy()
    assert policy.minimum_k == 1
    raw_policy = json.loads(
        SYNTHETIC_PARTITION_POLICY_PATH.read_text(encoding="utf-8")
    )
    assert validate_synthetic_partition_policy(raw_policy) == []
    raw_policy["default_k_range"]["minimum"] = 0
    assert "invalid_k_range" in validate_synthetic_partition_policy(raw_policy)

    landscape = build_operator_partition_landscape(
        TARGET,
        (_candidate("amide", AMIDE_FORMATION),),
    )
    assert SyntheticPartitionLandscape.from_dict(landscape.to_dict()) == landscape
    html = render_partition_landscape_html(landscape)
    assert "Synthetic partition landscape" in html
    assert "Role-neutral modules" in html
    json_output = tmp_path / "landscape.json"
    html_output = tmp_path / "landscape.html"
    write_partition_landscape_review(landscape, json_output, html_output)
    assert json.loads(json_output.read_text(encoding="utf-8"))["target_id"]
    assert "<svg" in html_output.read_text(encoding="utf-8")


def test_phase0_baseline_keeps_panels_hash_bound() -> None:
    baseline = json.loads(
        (
            PROJECT_ROOT
            / "docs/new/synthetic_partition_phase0_baseline.v1.json"
        ).read_text(encoding="utf-8")
    )

    assert baseline["test_baseline"]["passed"] == 1193
    assert baseline["test_baseline"]["failed"] == 0
    for panel in baseline["panels"]:
        payload = (PROJECT_ROOT / panel["path"]).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == panel["sha256"]


def test_partition_landscape_cli_writes_review_artifacts(
    monkeypatch,
    tmp_path,
    capsys,
) -> None:
    monkeypatch.setattr(cli_module, "load_generic_library", lambda _: object())
    monkeypatch.setattr(
        cli_module,
        "disconnect_operator_ladder",
        lambda *args, **kwargs: (_candidate("amide", AMIDE_FORMATION),),
    )
    json_output = tmp_path / "cli-landscape.json"
    html_output = tmp_path / "cli-landscape.html"

    status = cli_module.main(
        (
            "partition-landscape",
            "library.json.gz",
            TARGET,
            str(json_output),
            str(html_output),
            "--maximum-k",
            "3",
        )
    )

    assert status == 0
    assert json_output.is_file()
    assert html_output.is_file()
    summary = json.loads(capsys.readouterr().out)
    assert summary["partition_count"] == 2
    assert summary["abstained"] is False
