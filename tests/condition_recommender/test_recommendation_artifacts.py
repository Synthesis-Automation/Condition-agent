import csv
import json
from pathlib import Path

import pytest

from condition_recommender.conversion import artifacts as artifact_module
from condition_recommender.conversion.artifacts import (
    build_recommendation_artifacts,
    combine_saved_recommendation_batches,
    discover_saved_conversion_batches,
    incomplete_saved_conversion_batches,
    recommendation_library_mode_dir,
    resume_saved_conversion_batch,
    save_recommendation_batch,
)
from condition_recommender.generic_indexing import load_generic_index
from condition_recommender.generic_api import GenericConditionRecommender


def test_recommendation_library_modes_use_isolated_directories(
    tmp_path: Path,
) -> None:
    assert recommendation_library_mode_dir(tmp_path, "Full") == tmp_path / "full"
    assert recommendation_library_mode_dir(tmp_path, "compact") == (
        tmp_path / "compact"
    )
    (tmp_path / "generic_index.sqlite").touch()
    assert recommendation_library_mode_dir(tmp_path, "full") == tmp_path
    assert recommendation_library_mode_dir(tmp_path, "compact") == (
        tmp_path / "compact"
    )
    with pytest.raises(ValueError, match="Unsupported recommendation library"):
        recommendation_library_mode_dir(tmp_path, "preview")


def _source_row(reaction_id: str) -> dict[str, str]:
    return {
        "reaction_id": reaction_id,
        "reaction_type": "Suzuki source",
        "reaction_smiles": (
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        ),
        "yield_pct": "80",
        "temperature_c": "",
        "time_h": "",
        "reference": f"reference-{reaction_id}",
        "reactant_cas": "",
        "product_cas": "",
        "reagent_cas": "584-08-7",
        "catalyst_cas": "",
        "solvent_cas": "108-88-3",
        "experimental_procedure": "",
        "stages": "1",
        "steps": "1",
        "notes": "",
    }


def test_artifact_workflow_builds_recommendation_data_without_review_csv(
    tmp_path: Path,
    monkeypatch,
) -> None:
    source = tmp_path / "source"
    source.mkdir()
    dataset = source / "reactions.csv"
    rows = [_source_row("one"), _source_row("two")]
    with dataset.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    output = tmp_path / "recommendation_data"
    output.mkdir()
    (output / "generic_index.json.gz").write_bytes(b"retired")
    (output / "generic_review_index.json.gz").write_bytes(b"retired")
    (output / "generic_review_index.sqlite").write_bytes(b"stale")
    progress = []
    index_builds = []
    real_index_builder = artifact_module.build_sqlite_generic_index

    def counting_index_builder(*args, **kwargs):
        index_builds.append(bool(kwargs.get("include_review")))
        return real_index_builder(*args, **kwargs)

    monkeypatch.setattr(
        artifact_module,
        "build_sqlite_generic_index",
        counting_index_builder,
    )

    report = build_recommendation_artifacts(
        source,
        output,
        shard_size=1,
        workers=1,
        build_fast_index=True,
        progress_callback=progress.append,
    )

    assert report["record_count"] == 2
    assert report["schema_version"] == "2.3"
    assert report["eligible_index_record_count"] == 2
    assert report["trusted_precedent_count"] == 2
    assert report["review_core_precedent_count"] == 0
    assert report["unrestricted_precedent_count"] is None
    assert report["query_core_eligible_count"] == 2
    assert report["shard_count"] == 2
    assert report["storage"]["shard_file_count"] == 2
    assert not (output / "records.jsonl.gz").exists()
    assert (output / "shard_manifest.json").is_file()
    assert not (output / "generic_index.json.gz").exists()
    assert not (output / "generic_review_index.json.gz").exists()
    assert (output / "generic_index.sqlite").is_file()
    assert not (output / "generic_review_index.sqlite").exists()
    assert (output / "recommendation_artifacts_report.json").is_file()
    assert not report["review_index_generated"]
    assert not report["review_index_reuses_trusted"]
    assert index_builds == [False]
    assert not any(
        name.startswith("legacy_json") for name in report["artifacts"]
    )
    trusted_recommender = GenericConditionRecommender.from_path(
        output / "generic_index.sqlite"
    )
    review_recommender = GenericConditionRecommender.from_path(
        output / "generic_index.sqlite",
        include_review=True,
    )
    assert not trusted_recommender.includes_review_precedents
    assert review_recommender.includes_review_precedents
    assert review_recommender.review_index_reuses_trusted
    assert review_recommender.index.precedent_scope.value == "trusted"
    review_result = review_recommender.recommend(rows[0]["reaction_smiles"])
    assert "UNRESTRICTED_MODE_REUSES_TRUSTED_INDEX" in review_result.warnings
    assert report["artifacts"]["fast_index"]["path"].endswith(
        "generic_index.sqlite"
    )
    assert "review_csv" not in report["artifacts"]
    assert not (output / "reaction_review.csv").exists()
    manifest_recommender = GenericConditionRecommender.from_path(
        output / "shard_manifest.json"
    )
    assert trusted_recommender.recommend(rows[0]["reaction_smiles"]) == (
        manifest_recommender.recommend(rows[0]["reaction_smiles"])
    )
    assert len(load_generic_index(output / "shard_manifest.json").rows) == 2
    assert progress[0].phase == "canonical_discovered"
    assert progress[-1].phase == "completed"


def test_artifact_workflow_reuses_completed_shards(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    dataset = source / "reactions.csv"
    row = _source_row("one")
    with dataset.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    output = tmp_path / "recommendation_data"

    build_recommendation_artifacts(
        source,
        output,
        shard_size=1,
        build_fast_index=False,
    )
    resumed = build_recommendation_artifacts(
        source,
        output,
        shard_size=1,
        build_fast_index=False,
    )

    assert resumed["reused_shard_count"] == 1


def test_artifact_workflow_converts_only_explicitly_selected_files(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source"
    source.mkdir()
    selected_paths = []
    for name in ("selected-one.csv", "selected-two.csv", "not-selected.csv"):
        path = source / name
        row = _source_row(name)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        if name.startswith("selected-"):
            selected_paths.append(path)

    output = tmp_path / "recommendation_data"
    report = build_recommendation_artifacts(
        selected_paths,
        output,
        shard_size=1,
        build_fast_index=False,
    )

    assert report["record_count"] == 2
    assert report["source_path"] is None
    assert report["source_paths"] == [
        str(path.resolve()) for path in selected_paths
    ]


def test_artifact_workflow_rejects_output_inside_source(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source"
    source.mkdir()

    with pytest.raises(ValueError, match="outside the source"):
        build_recommendation_artifacts(source, source / "generated")


def test_saved_batches_combine_into_one_active_recommender(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    first = source / "first.csv"
    second = source / "second.csv"
    for path, reaction_id in ((first, "one"), (second, "two")):
        row = _source_row(reaction_id)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
    library = tmp_path / "library"

    first_report = save_recommendation_batch(
        first,
        library,
        batch_name="first-batch",
        shard_size=1,
    )
    second_report = save_recommendation_batch(
        second,
        library,
        batch_name="second-batch",
        shard_size=1,
    )
    (library / "generic_review_index.sqlite").write_bytes(b"stale")
    report = combine_saved_recommendation_batches(library)

    assert first_report["batch_name"] == "first-batch"
    assert second_report["batch_name"] == "second-batch"
    assert len(discover_saved_conversion_batches(library)) == 2
    assert report["batch_count"] == 2
    assert report["input_record_count"] == 2
    assert report["record_count"] == 2
    assert report["duplicate_record_count"] == 0
    assert report["schema_version"] == "1.3"
    assert "review_csv" not in report["artifacts"]
    assert not (library / "reaction_review.csv").exists()
    assert (library / "combined_records.jsonl.gz").is_file()
    assert (library / "combined_batch_manifest.json").is_file()
    assert (library / "generic_index.sqlite").is_file()
    assert not (library / "generic_review_index.sqlite").exists()
    assert (library / "reference_catalog.jsonl.gz").is_file()
    recommender = GenericConditionRecommender.from_path(
        library / "generic_index.sqlite"
    )
    assert len(recommender.index.rows) == 2
    assert report["unrestricted_precedent_count"] is None
    assert not report["review_index_generated"]
    assert "review_core_index" not in report["artifacts"]
    with pytest.raises(FileNotFoundError, match="Review-core index is unavailable"):
        GenericConditionRecommender.from_path(
            library / "generic_index.sqlite",
            include_review=True,
        )


def test_combining_saved_batches_deduplicates_identical_observations(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.csv"
    row = _source_row("same")
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    library = tmp_path / "library"
    save_recommendation_batch(source, library, batch_name="one", shard_size=1)
    save_recommendation_batch(source, library, batch_name="two", shard_size=1)

    report = combine_saved_recommendation_batches(library)

    assert report["input_record_count"] == 2
    assert report["record_count"] == 1
    assert report["duplicate_record_count"] == 1
    assert len(load_generic_index(library / "generic_index.sqlite").rows) == 1


def test_overlapping_batches_share_converted_source_artifacts(
    tmp_path: Path,
) -> None:
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    paths = {}
    for name in ("a.csv", "b.csv", "c.csv"):
        path = source_dir / name
        row = _source_row(name)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        paths[name] = path
    library = tmp_path / "library"

    first = save_recommendation_batch(
        (paths["a.csv"], paths["b.csv"]),
        library,
        batch_name="ab",
        shard_size=1,
    )
    second = save_recommendation_batch(
        (paths["a.csv"], paths["c.csv"]),
        library,
        batch_name="ac",
        shard_size=1,
    )

    source_manifests = tuple(
        (library / "converted_sources").glob("*/*/shard_manifest.json")
    )
    assert len(source_manifests) == 3
    assert not (Path(first["batch_dir"]) / "shards").exists()
    assert not (Path(second["batch_dir"]) / "shards").exists()
    assert second["reused_source_file_count"] == 1
    first_manifest = json.loads(
        (Path(first["batch_dir"]) / "shard_manifest.json").read_text()
    )
    second_manifest = json.loads(
        (Path(second["batch_dir"]) / "shard_manifest.json").read_text()
    )
    first_a = next(
        item
        for item in first_manifest["source_manifests"]
        if item["source_path"] == str(paths["a.csv"].resolve())
    )
    second_a = next(
        item
        for item in second_manifest["source_manifests"]
        if item["source_path"] == str(paths["a.csv"].resolve())
    )
    assert first_a["path"] == second_a["path"]

    combined = combine_saved_recommendation_batches(library)
    assert combined["input_record_count"] == 4
    assert combined["record_count"] == 3
    assert combined["duplicate_record_count"] == 1


def test_automatic_batch_identity_ignores_selection_order(
    tmp_path: Path,
) -> None:
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    paths = []
    for name in ("a.csv", "b.csv"):
        path = source_dir / name
        row = _source_row(name)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        paths.append(path)
    library = tmp_path / "library"

    first = save_recommendation_batch(paths, library, shard_size=1)
    second = save_recommendation_batch(tuple(reversed(paths)), library, shard_size=1)

    assert first["batch_name"] == second["batch_name"]
    assert first["batch_dir"] == second["batch_dir"]
    assert second["reused_source_file_count"] == 2
    assert len(discover_saved_conversion_batches(library)) == 1
    assert len(
        load_generic_index(
            Path(second["batch_dir"]) / "shard_manifest.json"
        ).rows
    ) == 2


def test_changed_source_content_creates_a_new_source_artifact(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.csv"
    first_row = _source_row("first")
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(first_row))
        writer.writeheader()
        writer.writerow(first_row)
    library = tmp_path / "library"

    first = save_recommendation_batch(
        source, library, batch_name="before", shard_size=1
    )
    second_row = _source_row("second")
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(second_row))
        writer.writeheader()
        writer.writerow(second_row)
    second = save_recommendation_batch(
        source, library, batch_name="after", shard_size=1
    )

    first_artifact = first["artifacts"]["converted_sources"][0]
    second_artifact = second["artifacts"]["converted_sources"][0]
    assert first_artifact["path"] != second_artifact["path"]
    assert second["reused_source_file_count"] == 0
    assert len(
        tuple((library / "converted_sources").glob("*/*/shard_manifest.json"))
    ) == 2


def test_combining_rejects_an_incomplete_saved_batch(tmp_path: Path) -> None:
    source = tmp_path / "source.csv"
    row = _source_row("partial")
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    library = tmp_path / "library"
    saved = save_recommendation_batch(
        source,
        library,
        batch_name="partial",
        shard_size=1,
    )
    combine_saved_recommendation_batches(library)
    active_index = library / "generic_index.sqlite"
    original_index = active_index.read_bytes()
    manifest_path = Path(saved["batch_dir"]) / "shard_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source_files"][0]["coverage_complete"] = False
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="conversion is incomplete"):
        combine_saved_recommendation_batches(library)

    assert active_index.read_bytes() == original_index
    assert incomplete_saved_conversion_batches(library) == (manifest_path,)

    resumed = resume_saved_conversion_batch(manifest_path)
    rebuilt = combine_saved_recommendation_batches(library)

    assert resumed["reused_shard_count"] == 1
    assert incomplete_saved_conversion_batches(library) == ()
    assert rebuilt["record_count"] == 1
