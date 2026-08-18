"""Executable core-admission and source round-trip regressions."""

from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path

import pytest
from reactive_taxonomy import featurize_reaction

from core_retrosynthesis.cli import _parser
from core_retrosynthesis.core_admission import (
    assess_generic_core_admission,
    load_generic_core_admission_policy,
)
from core_retrosynthesis.full_scale import FullScaleBuildConfig
from core_retrosynthesis.generic_compiler import compile_generic_templates
from core_retrosynthesis.generic_library import build_generic_library
from core_retrosynthesis.generic_search import disconnect_generic_target_detailed
from core_retrosynthesis.mapping import materialize_atom_mapping
from core_retrosynthesis.round_trip_audit import audit_generic_round_trips


NITRATION = (
    "[CH3:1][O:2][C:3](=[O:4])[c:5]1[cH:6][cH:10][c:11]2"
    "[c:12]([cH:13]1)[O:14][CH2:15][CH2:16][O:17]2."
    "[N+:7]([O-:8])(=[O:9])[OH:18]>>"
    "[CH3:1][O:2][C:3](=[O:4])[c:5]1[c:6]([N+:7]([O-:8])=[O:9])"
    "[cH:10][c:11]2[c:12]([cH:13]1)[O:14][CH2:15][CH2:16][O:17]2"
)
NITRO_REDUCTION = (
    "[CH3:1][O:2][C:3](=[O:4])[c:5]1[c:6]([N+:7]([O-:16])=[O:17])"
    "[cH:8][c:9]2[c:10]([cH:11]1)[O:12][CH2:13][CH2:14][O:15]2>>"
    "[CH3:1][O:2][C:3](=[O:4])[c:5]1[c:6]([NH2:7])[cH:8][c:9]2"
    "[c:10]([cH:11]1)[O:12][CH2:13][CH2:14][O:15]2"
)
ROUND_TRIP_FAILURE = (
    "[C:1]1(=[O:2])[CH:3]=[CH:4][C:5](=[O:7])[O:6]1."
    "[cH:8]1[cH:9][c:10]([Cl:11])[c:12]([O:13][CH3:14])"
    "[cH:15][cH:16]1>>"
    "[C:1](=[O:2])(/[CH:3]=[CH:4]/[C:5]([OH:6])=[O:7])"
    "[c:8]1[cH:9][c:10]([Cl:11])[c:12]([O:13][CH3:14])"
    "[cH:15][cH:16]1"
)


def _analysis_dict(reaction_smiles: str) -> dict:
    materialized = materialize_atom_mapping(reaction_smiles)
    assert materialized is not None
    analysis = featurize_reaction(materialized.reaction_smiles)
    assert analysis.valid and analysis.reaction_core is not None
    return analysis.to_dict()


def test_validated_departure_policy_is_explicit_and_versioned() -> None:
    strict = load_generic_core_admission_policy("pass_only")
    departure = load_generic_core_admission_policy("validated_departures")

    assert strict.allowed_core_statuses == ("pass",)
    assert departure.definition_id.endswith(":validated_departures")
    assert departure.allowed_review_reasons == (
        "not_all_edits_graph_checked",
    )
    assert departure.require_validated_departures_for_review


@pytest.mark.parametrize(
    ("reaction_smiles", "expected_levels"),
    (
        (NITRATION, {"L0", "L1", "L2"}),
        (NITRO_REDUCTION, {"L0"}),
    ),
)
def test_opt_in_departure_review_round_trips_real_aryl_sequence(
    reaction_smiles: str,
    expected_levels: set[str],
) -> None:
    row = {"reaction_id": "real-route-step", "reaction_smiles": reaction_smiles}

    strict = compile_generic_templates(
        row,
        levels=("L0", "L1", "L2"),
        admission_mode="data_driven",
    )
    experimental = compile_generic_templates(
        row,
        levels=("L0", "L1", "L2"),
        admission_mode="data_driven",
        core_admission_policy="validated_departures",
    )

    assert strict.rejection_reason == "materialized_core_not_verified"
    assert experimental.rejection_reason is None
    assert {item.abstraction_level for item in experimental.templates} == (
        expected_levels
    )
    assert all(item.stereo_policy == "exact" for item in experimental.templates)
    assert experimental.diagnostics is not None
    assert experimental.diagnostics["core_admission_reason"] == (
        "validated_departure_review"
    )
    level_status = {
        item["level"]: item["status"]
        for item in experimental.diagnostics["round_trip_levels"]
    }
    assert {level for level, status in level_status.items() if status == "exact"} == (
        expected_levels
    )


@pytest.mark.parametrize(
    ("mutation", "expected_issue"),
    (
        ("partial_mapping", "active_atom_mapping_below_policy"),
        ("blocking_reason", "core_has_blocking_reasons"),
        ("missing_completeness", "product_completeness_not_verified"),
        ("unsafe_unchecked_edit", "unchecked_edit_is_not_validated_departure"),
    ),
)
def test_departure_review_rejects_unsafe_evidence(
    mutation: str,
    expected_issue: str,
) -> None:
    value = _analysis_dict(NITRATION)
    core = copy.deepcopy(value["reaction_core"])
    observation = copy.deepcopy(value["observation"])
    if mutation == "partial_mapping":
        core["quality"]["active_atom_mapping_coverage"] = 0.5
    elif mutation == "blocking_reason":
        core["quality"]["blocking_reasons"] = ["formed_bond_state_inconsistent"]
    elif mutation == "missing_completeness":
        observation["completeness"]["status"] = "review"
    elif mutation == "unsafe_unchecked_edit":
        departure = next(
            item
            for item in observation["edits"]
            if item["edit_type"] == "broken"
        )
        departure["edit_type"] = "formed"
        departure["new_order"] = "SINGLE"
        departure["old_order"] = None

    decision = assess_generic_core_admission(
        core,
        observation,
        policy_name="validated_departures",
    )

    assert not decision.admitted
    assert expected_issue in decision.issues


def test_round_trip_failure_remains_a_hard_rejection() -> None:
    result = compile_generic_templates(
        {"reaction_id": "negative", "reaction_smiles": ROUND_TRIP_FAILURE},
        levels=("L0", "L1", "L2"),
        admission_mode="data_driven",
        core_admission_policy="validated_departures",
    )

    assert result.rejection_reason == "source_round_trip_failed"
    assert result.rejection_stage == "round_trip"
    assert result.diagnostics is not None
    assert all(
        item["status"] == "failed"
        for item in result.diagnostics["round_trip_levels"]
    )


def test_round_tripped_aryl_sequence_produces_forward_verified_candidates() -> None:
    rows = (
        {
            "reaction_id": "nitration",
            "reference_id": "US04011323",
            "reaction_smiles": NITRATION,
        },
        {
            "reaction_id": "reduction",
            "reference_id": "US04011323",
            "reaction_smiles": NITRO_REDUCTION,
        },
    )
    library = build_generic_library(
        rows,
        levels=("L0", "L1", "L2"),
        admission_mode="data_driven",
        core_admission_policy="validated_departures",
    )

    assert library.accepted_observation_count == 2
    assert len(library.templates) == 4
    assert len(library.operators) == 2
    for row in rows:
        target = row["reaction_smiles"].split(">>", 1)[1]
        candidates, diagnostics = disconnect_generic_target_detailed(
            target,
            library,
            top_k=10,
            max_templates_to_apply=20,
            max_candidates_to_validate=20,
            use_context=False,
        )
        assert diagnostics.invalid_forward_count == 0
        assert len(candidates) == 1
        assert candidates[0].forward_validation_status == "verified_signature"


def test_round_trip_audit_writes_level_evidence(tmp_path: Path) -> None:
    source = tmp_path / "records.jsonl.gz"
    with gzip.open(source, "wt", encoding="utf-8") as stream:
        stream.write(
            json.dumps({"reaction_id": "nitration", "reaction_smiles": NITRATION})
            + "\n"
        )
        stream.write(
            json.dumps(
                {"reaction_id": "reduction", "reaction_smiles": NITRO_REDUCTION}
            )
            + "\n"
        )
    output = tmp_path / "audit.json"

    report = audit_generic_round_trips(
        source,
        output,
        reaction_ids=("nitration", "reduction"),
        core_admission_policy="validated_departures",
    )

    assert report["accepted_reaction_count"] == 2
    assert output.is_file()
    stored = json.loads(output.read_text(encoding="utf-8"))
    assert stored["results"][0]["status"] == "accepted"
    assert stored["results"][0]["departure_edit_descriptors"] == [
        "broken:N-O:SINGLE:7:18"
    ]


def test_full_build_and_cli_expose_opt_in_policy() -> None:
    strict = FullScaleBuildConfig()
    experimental = FullScaleBuildConfig(
        core_admission_policy="validated_departures"
    )
    arguments = _parser().parse_args(
        [
            "audit-operator-round-trips",
            "records.jsonl.gz",
            "audit.json",
            "--reaction-id",
            "nitration",
            "--core-admission-policy",
            "validated_departures",
        ]
    )

    assert strict.config_id != experimental.config_id
    assert arguments.core_admission_policy == "validated_departures"
