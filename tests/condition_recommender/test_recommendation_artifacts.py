import csv
from pathlib import Path

import pytest

from condition_recommender.conversion.artifacts import (
    build_recommendation_artifacts,
)
from condition_recommender.conversion.concise_review import (
    CONCISE_REACTION_REVIEW_FIELDS,
)
from condition_recommender.generic_indexing import load_generic_index


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


def test_artifact_workflow_builds_recommendation_data_and_review_csv(
    tmp_path: Path,
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
    progress = []

    report = build_recommendation_artifacts(
        source,
        output,
        shard_size=1,
        workers=1,
        build_fast_index=True,
        progress_callback=progress.append,
    )

    assert report["record_count"] == 2
    assert report["eligible_index_record_count"] == 2
    assert report["shard_count"] == 2
    assert report["storage"]["shard_file_count"] == 2
    assert not (output / "records.jsonl.gz").exists()
    assert (output / "shard_manifest.json").is_file()
    assert (output / "generic_index.json.gz").is_file()
    assert (output / "recommendation_artifacts_report.json").is_file()
    assert len(load_generic_index(output / "generic_index.json.gz").rows) == 2
    assert len(load_generic_index(output / "shard_manifest.json").rows) == 2
    with (output / "reaction_review.csv").open(
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        review_rows = list(csv.DictReader(handle))
    assert len(review_rows) == 2
    assert tuple(review_rows[0]) == CONCISE_REACTION_REVIEW_FIELDS
    assert tuple(review_rows[0])[:10] == (
        "canonical_reaction_smiles",
        "reaction_label",
        "primary_reaction_pattern",
        "primary_reaction_pattern_count",
        "reaction_pattern_matches",
        "co_occurring_reaction_patterns",
        "identified_reaction_type",
        "compatible_reaction_types",
        "reaction_pattern_confidence",
        "reaction_pattern_requires_condition_evidence",
    )
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
