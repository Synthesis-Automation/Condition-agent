import csv
import json
from pathlib import Path

import pytest
from reactive_taxonomy import (
    REACTION_SIGNATURE_SCHEMA_VERSION,
    reaction_signature_definition_versions,
)

from condition_recommender.evaluation import (
    compare_generic_baselines,
    evaluate_generic_index,
    grouped_holdout_split,
)
from condition_recommender.evaluation_features import reaction_scaffold_tokens
from condition_recommender.calibration import (
    calibrate_generic_ranking,
    load_generic_calibration_rules,
    validate_generic_calibration_rules,
)
from condition_recommender.chemist_review import generate_chemist_review_packet
from condition_recommender.chemist_review_adjudication import (
    adjudicate_chemist_review,
    valid_completed_chemist_review,
)
from condition_recommender.generic_indexing import (
    build_generic_index,
    load_generic_index,
    load_persisted_generic_index,
    save_generic_index,
    validate_generic_index_artifact,
)
from condition_recommender.release_validation import (
    build_release_readiness_report,
)


def _signature(index: int) -> dict:
    return {
        "schema_version": REACTION_SIGNATURE_SCHEMA_VERSION,
        "definition_versions": reaction_signature_definition_versions(),
        "signature_id": f"RS3:{index}",
        "exact_signature_key": f"L0:{index}",
        "handle_signature_key": "L1:shared",
        "transformation_signature_key": "L2:shared",
        "bond_edit_signature_key": "L3:shared",
        "environment_signature_key": f"L4:{index % 2}",
        "formed_bond_types": ["C-N:SINGLE"],
        "broken_bond_types": ["B-C:SINGLE"],
        "order_changes": [],
        "partners": [
            {
                "component_index": 0,
                "handle_tokens": ["B(OH)2", "N-H"],
                "anchor_contexts": ["Ar"],
                "role": "transfer_partner",
                "reactivity_profile": {
                    "schema_version": "1.0",
                    "context_kind": "aromatic",
                    "context": {
                        "context_kind": "aromatic",
                        "ring_family": "benzene",
                        "ring_sizes": [6],
                        "ortho_occupancy_count": 0,
                        "ortho_capacity": 2,
                        "ortho_burden_class": "none",
                        "heteroatoms": [],
                    },
                    "steric": {
                        "accessibility_class": "open",
                        "accessibility_score": 0.0,
                    },
                    "electronic": {
                        "activation_axis": "electronic_demand",
                        "activation_class": "balanced",
                        "activation_score": 0.0,
                    },
                    "reactive_center": {},
                    "modifiers": [],
                },
                "nearby_groups": [],
            }
        ],
        "spectator_groups": [],
        "transformation_class": "generic_c_n_coupling",
        "named_family": None,
        "family_confidence": 0.0,
        "topology": {
            "reaction_scope": "intermolecular",
            "formed_bond_scopes": ["intermolecular"],
            "formed_ring_sizes": [],
            "ring_count_delta": 0,
        },
    }


def _record(index: int, *, canonical_group: str | None = None) -> dict:
    recipe_id = f"RCR1:{index % 2}"
    recipe_core_id = f"RCORE1:{index % 2}"
    return {
        "schema_version": "3.1",
        "converter_definition_version": "generic_conversion.v2.1",
        "admission_tier": "verified",
        "index_eligibility": "eligible",
        "chemistry_status": "verified",
        "condition_status": "resolved_complete",
        "condition_stage_status": "single_stage",
        "outcome_status": "usable",
        "reaction_id": f"reaction-{index}",
        "observation_id": f"observation-{index}",
        "canonical_reaction_id": canonical_group or f"CRX1:{index}",
        "reaction_smiles": "C.N>>CN",
        "yield_pct": 50.0 + index,
        "source_dataset": f"dataset-{index % 2}",
        "reference_id": f"REF1:{index}",
        "reference_condition_series_id": f"RCS1:{index}",
        "reaction_signature": _signature(index),
        "resolved_recipe_id": recipe_id,
        "resolved_recipe_core_id": recipe_core_id,
        "resolved_recipe": {
            "recipe_id": recipe_id,
            "recipe_core_id": recipe_core_id,
        },
        "condition_resolution": {"has_uncertainty": False},
    }


def test_persisted_index_round_trip_is_deterministic(tmp_path: Path) -> None:
    index = build_generic_index([_record(1), _record(2)])
    first_path = tmp_path / "first.json"
    second_path = tmp_path / "second.json"
    first = save_generic_index(index, first_path)
    second = save_generic_index(index, second_path)
    loaded = load_persisted_generic_index(first_path)
    assert first["index_id"] == second["index_id"]
    assert first_path.read_bytes() == second_path.read_bytes()
    assert loaded == index
    assert load_generic_index(first_path) == index
    payload = json.loads(first_path.read_text(encoding="utf-8"))
    assert payload["schema_version"] == "2.0"
    assert payload["reaction_signature_schema_version"] == "3.0"
    assert payload["record_schema_versions"] == ["3.1"]
    assert payload["maps"]["environment_features"]
    integrity = validate_generic_index_artifact(first_path)
    assert integrity["valid"]
    assert integrity["row_count"] == 2


def test_compressed_persisted_index_round_trip(tmp_path: Path) -> None:
    index = build_generic_index([_record(1), _record(2)])
    path = tmp_path / "generic_index.json.gz"

    report = save_generic_index(index, path)

    assert path.read_bytes().startswith(b"\x1f\x8b")
    assert report["row_count"] == 2
    assert load_persisted_generic_index(path) == index
    assert load_generic_index(path) == index
    assert validate_generic_index_artifact(path)["valid"]


def test_index_rejects_stale_converted_records() -> None:
    stale = _record(1)
    stale["reaction_signature"] = {
        **stale["reaction_signature"],
        "schema_version": "1.0",
    }
    with pytest.raises(ValueError, match="regenerate converted records"):
        build_generic_index([stale])


def test_loader_rejects_stale_index_schema(tmp_path: Path) -> None:
    path = tmp_path / "index.json"
    save_generic_index(build_generic_index([_record(1)]), path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["schema_version"] = "1.0"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="rebuild the index"):
        load_persisted_generic_index(path)


def test_persisted_index_detects_tampering(tmp_path: Path) -> None:
    path = tmp_path / "index.json"
    save_generic_index(build_generic_index([_record(1)]), path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["rows"][0]["yield_pct"] = 99.0
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="integrity"):
        load_persisted_generic_index(path)


def test_grouped_split_keeps_duplicate_reactions_together() -> None:
    index = build_generic_index(
        [
            _record(1, canonical_group="CRX1:duplicate"),
            _record(2, canonical_group="CRX1:duplicate"),
            _record(3),
            _record(4),
        ]
    )
    first = grouped_holdout_split(index.rows, test_fraction=0.5, seed=7)
    second = grouped_holdout_split(index.rows, test_fraction=0.5, seed=7)
    assert first == second
    assert not set(first.train_group_ids).intersection(first.test_group_ids)
    duplicate_locations = {
        "train" if row in first.train_rows else "test"
        for row in index.rows
        if row.canonical_reaction_id == "CRX1:duplicate"
    }
    assert len(duplicate_locations) == 1


def test_grouped_split_keeps_publication_references_together() -> None:
    records = [_record(index) for index in range(1, 5)]
    records[0]["reference_id"] = "REF1:shared"
    records[1]["reference_id"] = "REF1:shared"
    index = build_generic_index(records)

    split = grouped_holdout_split(index.rows, test_fraction=0.5, seed=7)
    locations = {
        "train" if row in split.train_rows else "test"
        for row in index.rows
        if row.reference_id == "REF1:shared"
    }

    assert len(locations) == 1


def test_grouped_evaluation_writes_leakage_safe_metrics(tmp_path: Path) -> None:
    index_path = tmp_path / "index.json"
    output = tmp_path / "evaluation"
    save_generic_index(
        build_generic_index([_record(index) for index in range(10)]),
        index_path,
    )
    report = evaluate_generic_index(
        index_path,
        output,
        test_fraction=0.3,
        seed=11,
        top_k=2,
        minimum_pool_size=1,
    )
    assert report["split"]["leakage_group_count"] == 0
    assert report["split"]["reference_overlap_count"] == 0
    assert report["split"]["canonical_reaction_overlap_count"] == 0
    assert report["split"]["test_group_count"] == 3
    assert report["metrics"]["query_count"] == 3
    assert report["metrics"]["coverage_rate"] == 1.0
    assert report["metrics"]["hard_incompatible_recommendation_count"] == 0
    assert report["schema_version"] == "1.4"
    assert report["definition_versions"] == {
        "compatibility": "1.1",
        "retrieval": "1.5",
        "similarity": "1.0",
        "ranking": "1.0",
    }
    assert (output / "evaluation_report.json").is_file()
    cases = [
        json.loads(line)
        for line in (output / "evaluation_cases.jsonl").read_text().splitlines()
    ]
    assert len(cases) == 3
    assert cases[0]["top_recommendation_score_trace"] is not None
    assert report["metrics"]["explanation_complete_rate"] == 1.0
    assert "generic_c_n_coupling" in report["by_transformation_class"]


def test_scaffold_disjoint_split_has_no_reactive_scaffold_overlap() -> None:
    reactions = (
        "Brc1ccccc1.N>>Nc1ccccc1",
        "BrC1CCCCC1.N>>NC1CCCCC1",
        "Brc1cccc2ccccc12.N>>Nc1cccc2ccccc12",
        "BrC1CCCC1.N>>NC1CCCC1",
    )
    records = []
    for index, reaction_smiles in enumerate(reactions, start=1):
        record = _record(index)
        record["reaction_smiles"] = reaction_smiles
        records.append(record)
    split = grouped_holdout_split(
        build_generic_index(records).rows,
        test_fraction=0.5,
        seed=7,
        split_mode="scaffold_disjoint",
    )

    train_tokens = {
        token for row in split.train_rows for token in row.scaffold_tokens
    }
    test_tokens = {
        token for row in split.test_rows for token in row.scaffold_tokens
    }
    assert split.split_mode == "scaffold_disjoint"
    assert not train_tokens.intersection(test_tokens)


def test_murcko_scaffold_clears_orphaned_double_bond_stereo() -> None:
    reaction = (
        "Ic1ccccc1."
        "CCCC/C=C/C(=C(/C=C\\c1ccccc1)[Si](C)(C)O)c1ccccc1"
        ">>CCCC/C=C/C(=C(/C=C\\c1ccccc1)c1ccccc1)c1ccccc1"
    )

    tokens = reaction_scaffold_tokens(
        reaction,
        {"partners": [{"component_index": 1}]},
    )

    assert tokens
    assert all(token.startswith("BM:") for token in tokens)


def test_forward_time_split_uses_latest_publication_groups() -> None:
    records = []
    for index, year in enumerate((2001, 2005, 2010, 2020), start=1):
        record = _record(index)
        record["reference_identity"] = {"publication_year": year}
        records.append(record)
    split = grouped_holdout_split(
        build_generic_index(records).rows,
        test_fraction=0.5,
        split_mode="forward_time",
    )

    train_years = {
        row.publication_year for row in split.train_rows
    }
    test_years = {
        row.publication_year for row in split.test_rows
    }
    assert split.split_mode == "forward_time"
    assert max(train_years) <= min(test_years)
    assert split.cutoff_year == min(test_years)


def test_source_disjoint_split_has_no_dataset_overlap() -> None:
    index = build_generic_index([_record(index) for index in range(8)])
    split = grouped_holdout_split(
        index.rows,
        test_fraction=0.5,
        seed=7,
        split_mode="source_disjoint",
    )

    assert not {
        row.source_dataset for row in split.train_rows
    }.intersection(row.source_dataset for row in split.test_rows)


def test_grouped_split_validates_inputs() -> None:
    one_group = build_generic_index([_record(1)]).rows
    with pytest.raises(ValueError, match="two canonical"):
        grouped_holdout_split(one_group)
    with pytest.raises(ValueError, match="between zero and one"):
        grouped_holdout_split(one_group, test_fraction=1.0)


def test_baseline_comparison_uses_one_split_and_writes_report(
    tmp_path: Path,
) -> None:
    index_path = tmp_path / "index.json"
    output = tmp_path / "baselines"
    save_generic_index(
        build_generic_index([_record(index) for index in range(10)]),
        index_path,
    )

    comparison = compare_generic_baselines(
        index_path,
        output,
        test_fraction=0.3,
        seed=11,
        top_k=2,
        minimum_pool_size=1,
    )

    assert set(comparison["strategies"]) == {
        "family_only",
        "generic_only",
        "hybrid",
        "transformation_prior",
        "legacy_pilot",
    }
    assert comparison["strategies"]["family_only"]["coverage_rate"] == 0.0
    assert comparison["strategies"]["generic_only"]["coverage_rate"] == 1.0
    assert (
        comparison["split"]
        == json.loads(
            (output / "hybrid" / "evaluation_report.json").read_text()
        )["split"]
    )
    assert (output / "baseline_comparison.json").is_file()
    assert (output / "baseline_comparison.md").is_file()


def test_calibration_selects_on_development_and_gates_on_validation(
    tmp_path: Path,
) -> None:
    development_path = tmp_path / "development.json"
    validation_path = tmp_path / "validation.json"
    output = tmp_path / "calibration"
    records = [_record(index) for index in range(30)]
    save_generic_index(build_generic_index(records), development_path)
    save_generic_index(build_generic_index(records), validation_path)

    report = calibrate_generic_ranking(
        development_path,
        validation_path,
        output,
    )

    assert report["selected_candidate_id"] is not None
    assert report["promoted"]
    assert not report["promotion_reasons"]
    assert (output / "calibration_report.json").is_file()
    assert (output / "recommended_ranking_definition.json").is_file()
    assert (output / "recommended_retrieval_definition.json").is_file()


def test_calibration_definition_rejects_stale_schema() -> None:
    rules = {**load_generic_calibration_rules(), "schema_version": "stale"}

    with pytest.raises(ValueError, match="calibration schema"):
        validate_generic_calibration_rules(rules)


def test_chemist_review_packet_randomizes_candidates_and_separates_key(
    tmp_path: Path,
) -> None:
    index_path = tmp_path / "index.json"
    records = [_record(index) for index in range(12)]
    for record in records[8:]:
        record["reaction_signature"] = {
            **record["reaction_signature"],
            "bond_edit_signature_key": "L3:alternate",
        }
    save_generic_index(build_generic_index(records), index_path)
    output = tmp_path / "review"

    report = generate_chemist_review_packet(
        index_path,
        output,
        max_cases=4,
        seed=17,
        top_k=1,
        minimum_pool_size=1,
    )

    assert report["case_count"] == 2
    assert report["negative_control_count"] >= 1
    assert (output / "review_packet.jsonl").is_file()
    assert (output / "review_packet.html").is_file()
    assert (output / "answer_key.jsonl").is_file()
    assert (output / "review_form.csv").is_file()
    assert len(tuple((output / "structures").glob("*.svg"))) == 2
    packet_text = (output / "review_packet.jsonl").read_text()
    assert "negative_control" not in packet_text
    assert "Control candidate" not in packet_text
    assert "negative_control" in (output / "answer_key.jsonl").read_text()


def test_completed_chemist_review_is_unblinded_with_bound_provenance(
    tmp_path: Path,
) -> None:
    index_path = tmp_path / "index.json"
    records = [_record(index) for index in range(12)]
    for record in records[8:]:
        record["reaction_signature"] = {
            **record["reaction_signature"],
            "bond_edit_signature_key": "L3:alternate",
        }
    save_generic_index(build_generic_index(records), index_path)
    review_dir = tmp_path / "review"
    generate_chemist_review_packet(
        index_path,
        review_dir,
        max_cases=4,
        top_k=1,
        minimum_pool_size=1,
    )
    form_path = review_dir / "review_form.csv"
    with form_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
        fieldnames = tuple(rows[0])
    for row in rows:
        row["chemist_decision"] = "incompatible"
        row["chemist_confidence"] = "high"
    with form_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary = adjudicate_chemist_review(
        review_dir,
        review_dir / "review_summary.json",
        reviewer_name="Independent Chemist",
        independent_reviewer=True,
        reviewer_signoff=True,
    )

    assert summary["completed"]
    assert valid_completed_chemist_review(summary)
    assert summary["metrics"]["candidate_count"] == len(rows)
    assert summary["metrics"]["negative_control_rejection_rate"] == 1.0
    assert len(summary["provenance"]["review_form_sha256"]) == 64


def test_chemist_review_does_not_unblind_an_incomplete_form(
    tmp_path: Path,
) -> None:
    review_dir = tmp_path / "review"
    review_dir.mkdir()
    (review_dir / "review_form.csv").write_text(
        "case_id,candidate_id,chemist_decision\nCASE-1,CAND-1,\n",
        encoding="utf-8",
    )
    (review_dir / "answer_key.jsonl").write_text(
        "not valid JSON\n",
        encoding="utf-8",
    )
    (review_dir / "review_packet.jsonl").write_text("{}\n", encoding="utf-8")
    (review_dir / "review_packet.html").write_text("<html></html>", encoding="utf-8")
    (review_dir / "review_report.json").write_text(
        '{"candidate_count": 1}\n',
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="complete before unblinding"):
        adjudicate_chemist_review(
            review_dir,
            review_dir / "summary.json",
            reviewer_name="Independent Chemist",
            independent_reviewer=True,
            reviewer_signoff=True,
        )


def test_chemist_review_summary_requires_independent_signed_provenance() -> None:
    assert not valid_completed_chemist_review(
        {
            "schema_version": "1.0",
            "artifact_type": "generic_chemist_review_summary",
            "completed": True,
            "reviewer": {
                "name": "Reviewer",
                "independent": False,
                "signed_off": True,
            },
            "unresolved_systematic_defect_count": 0,
            "provenance": {},
        }
    )


def test_release_report_keeps_human_review_as_a_separate_gate(
    tmp_path: Path,
) -> None:
    conversion = {
        "failed_shard_count": 0,
        "input_row_count": 10,
        "output_row_count": 10,
        "integrity": {"valid": True},
        "index_eligibility_counts": {"ineligible": 2, "review_only": 3},
    }
    index = {
        "valid": True,
        "transformation_class_counts": {"generic": 5},
    }
    evaluation = {
        "split": {
            "reference_overlap_count": 0,
            "canonical_reaction_overlap_count": 0,
            "scaffold_overlap_count": 0,
        },
        "metrics": {
            "query_count": 2,
            "hard_incompatible_recommendation_count": 0,
        },
        "definition_versions": {"ranking": "1.0"},
    }
    calibration = {"promoted": True}
    values = {
        "conversion": conversion,
        "index": index,
        "grouped": evaluation,
        "scaffold": evaluation,
        "untouched": evaluation,
        "calibration": calibration,
    }
    paths = {}
    for name, value in values.items():
        path = tmp_path / f"{name}.json"
        path.write_text(json.dumps(value), encoding="utf-8")
        paths[name] = path

    report = build_release_readiness_report(
        paths["conversion"],
        paths["index"],
        paths["grouped"],
        paths["scaffold"],
        paths["untouched"],
        paths["calibration"],
        tmp_path / "release.json",
    )

    assert report["machine_release_ready"]
    assert not report["production_release_ready"]
    assert not report["human_review"]["summary_supplied"]
