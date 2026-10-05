"""Regressions for source-faithful consolidated dataset migration."""

import csv
import gzip
import json
from pathlib import Path

import pytest

from condition_recommender.conversion.input_schema import (
    discover_conversion_datasets,
    iter_conversion_records,
)
from condition_recommender.ingestion import (
    PreprocessingCancelled,
    detect_adapter,
    preprocess_file,
    preprocess_files,
)
from condition_recommender.ingestion.migration import (
    regenerate_intermediate_datasets,
    validate_intermediate_manifest,
)
from condition_recommender.ingestion.route_release import prepare_route_release


def _csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _routes(path: Path, reactions: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for index, reaction in enumerate(reactions):
            subtree = f"US1_{index}_0"
            handle.write(
                json.dumps(
                    {
                        "route_id": f"US1_{index}",
                        "original_tree": {
                            "depth": 1,
                            "num_reactions": 1,
                            "reaction_ids": ["42_0"],
                        },
                        "subtrees": [
                            {
                                "subtree_id": subtree,
                                "reactions": [
                                    {
                                        "_id": f"{subtree}_42_0",
                                        "reaction_smiles": reaction,
                                        "abstracted_reaction_smiles": (
                                            "[1CH3:1]>>[4CH3:1]"
                                        ),
                                    }
                                ],
                            }
                        ],
                    }
                )
                + "\n"
            )


def _read(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def test_route_extraction_keeps_agents_memberships_and_deduplicates_source_ids(
    tmp_path: Path,
) -> None:
    route = tmp_path / "uspto.higher-level.routes.jsonl.gz"
    original = tmp_path / "uspto_original.csv"
    reaction = "[CH3:1]Br.O>[Na+].[OH-]>[CH3:1]O"
    _routes(route, [reaction, reaction])
    _csv(original, [{"id": "uspto_42_0", "rxn_smiles": "[CH3:1]Br.O>>[CH3:1]O"}])
    output = tmp_path / "prepared"
    report = prepare_route_release(route, original, output)
    assert report["coverage_complete"]
    assert report["counts"]["route_count"] == 2
    assert report["counts"]["output_step_observations"] == 1
    assert len(_read(output / "routes.source.jsonl.gz")) == 2
    record = _read(output / "route_steps.observations.jsonl.gz")[0]
    assert record["reaction"]["reaction_smiles"] == "[CH3:1]Br.O>>[CH3:1]O"
    assert record["reaction"]["supplied_mapping_status"] == "supplied_unvalidated"
    assert record["source"]["source_groups"]["source_occurrence_count"] == "2"
    assert record["source"]["reference"] == "US1"
    assert {x["source_role_hint"] for x in record["conditions"]["components"]} == {None}
    raw = next(iter_conversion_records(output / "route_steps.observations.jsonl.gz"))
    assert {x.raw_identifier for x in raw.condition_component_inputs} == {
        "[Na+]",
        "[OH-]",
    }
    assert prepare_route_release(route, original, output)["reused"]


def test_route_conflicting_source_ids_keep_all_variants_for_review(
    tmp_path: Path,
) -> None:
    route = tmp_path / "routes.jsonl.gz"
    original = tmp_path / "uspto_original.csv"
    _routes(route, ["CBr>O>CO", "CBr>N>CN"])
    _csv(original, [{"id": "uspto_42_0", "rxn_smiles": "CBr>>CO"}])
    output = tmp_path / "prepared"
    report = prepare_route_release(route, original, output)
    assert report["counts"]["conflicting_source_ids"] == 1
    records = _read(output / "route_steps.observations.jsonl.gz")
    assert len(records) == 2
    assert len({r["observation_id"] for r in records}) == 2
    assert len({r["source"]["source_record_id"] for r in records}) == 2
    assert all(r["ingestion_status"] == "review" for r in records)
    assert all(
        "CONFLICTING_REACTION_STRINGS_FOR_SOURCE_ID" in r["warnings"] for r in records
    )


def test_route_csv_coverage_conflict_fails_gate(tmp_path: Path) -> None:
    route = tmp_path / "routes.jsonl.gz"
    original = tmp_path / "uspto_original.csv"
    _routes(route, ["CBr>O>CO"])
    _csv(original, [{"id": "uspto_42_0", "rxn_smiles": "CBr>>CN"}])
    with pytest.raises(ValueError, match="coverage gate failed"):
        prepare_route_release(route, original, tmp_path / "prepared")


def test_abstractions_are_archived_without_observed_reaction_structure(
    tmp_path: Path,
) -> None:
    source = tmp_path / "uspto_higher-level.csv"
    _csv(source, [{"id": "uspto_US1_0_0_42_0", "rxn_smiles": "[1CH3:1]>>[4CH3:1]"}])
    assert detect_adapter(source).adapter_id == "higher_level_abstraction_csv.v1"
    report = preprocess_file(source, tmp_path / "prepared")
    record = _read(Path(report["output_path"]))[0]
    assert record["observation_kind"] == "algorithmic_abstraction"
    assert record["reaction"]["reaction_smiles"] is None
    assert record["raw_fields"]["rxn_smiles"] == "[1CH3:1]>>[4CH3:1]"
    assert not discover_conversion_datasets(tmp_path / "prepared")
    with pytest.raises(ValueError, match="requires the abstraction adapter"):
        preprocess_file(
            source, tmp_path / "prepared", adapter_id="reaction_smiles_csv.v1"
        )


def test_ambiguous_reaction_export_requires_explicit_adapter(tmp_path: Path) -> None:
    source = tmp_path / "unknown.csv"
    _csv(source, [{"id": "1", "rxn_smiles": "CBr>>CO"}])
    with pytest.raises(ValueError, match="explicit reaction CSV adapter"):
        detect_adapter(source)
    report = preprocess_file(
        source, tmp_path / "prepared", adapter_id="reaction_smiles_csv.v1"
    )
    assert report["output_row_count"] == 1


def test_complete_migration_has_one_step_pool_and_reusable_manifest(
    tmp_path: Path,
) -> None:
    raw = tmp_path / "raw_datasets"
    release = raw / "routes" / "higher_level_retrosynthesis"
    original = release / "reactions" / "uspto_original.csv"
    abstraction = release / "reactions" / "uspto_higher-level.csv"
    route = release / "routes" / "uspto.higher-level.routes.jsonl.gz"
    _routes(route, ["CBr>O>CO", "CBr>O>CO"])
    _csv(original, [{"id": "uspto_42_0", "rxn_smiles": "CBr>>CO"}])
    _csv(abstraction, [{"id": "abstract_1", "rxn_smiles": "[1CH3:1]>>[4CH3:1]"}])
    output = tmp_path / "intermediate_datasets"
    report = regenerate_intermediate_datasets(raw, output)
    assert report["source_dataset_count"] == 3
    assert report["single_step_observation_count"] == 1
    assert report["abstraction_count"] == 1
    assert report["route_count"] == 2
    assert len(discover_conversion_datasets(output)) == 1
    assert validate_intermediate_manifest(output)["verified_artifact_count"] == 3
    regenerate_intermediate_datasets(raw, output)
    assert validate_intermediate_manifest(output)["coverage_complete"]
    path = discover_conversion_datasets(output)[0]
    path.write_bytes(b"broken")
    with pytest.raises(ValueError, match="Invalid intermediate artifact"):
        validate_intermediate_manifest(output)


def test_gui_batch_uses_same_layout_and_skips_duplicate_original_export(
    tmp_path: Path,
) -> None:
    raw = tmp_path / "raw_datasets"
    release = raw / "routes" / "higher_level_retrosynthesis"
    original = release / "reactions" / "uspto_original.csv"
    abstraction = release / "reactions" / "uspto_higher-level.csv"
    route = release / "routes" / "uspto.higher-level.routes.jsonl.gz"
    _routes(route, ["CBr>O>CO"])
    _csv(original, [{"id": "uspto_42_0", "rxn_smiles": "CBr>>CO"}])
    _csv(abstraction, [{"id": "abstract_1", "rxn_smiles": "[1CH3:1]>>[4CH3:1]"}])
    output = tmp_path / "intermediate_datasets"
    report = preprocess_files([original, abstraction, route], output, source_root=raw)
    assert report["file_count"] == 2
    assert len(discover_conversion_datasets(output)) == 1
    assert all(item["output_size_bytes"] > 0 for item in report["files"])
    assert (
        output / "routes/higher_level_retrosynthesis/routes.source.jsonl.gz"
    ).exists()


def test_cancelled_route_preparation_discards_temporary_files(tmp_path: Path) -> None:
    route = tmp_path / "routes.jsonl.gz"
    original = tmp_path / "uspto_original.csv"
    _routes(route, ["CBr>O>CO"])
    _csv(original, [{"id": "uspto_42_0", "rxn_smiles": "CBr>>CO"}])
    output = tmp_path / "prepared"
    with pytest.raises(PreprocessingCancelled):
        prepare_route_release(route, original, output, cancel_check=lambda: True)
    assert not list(output.glob("*.tmp"))
    assert not (output / "route_release.manifest.json").exists()
