from pathlib import Path

from PyQt6 import QtWidgets

from app import reaction_converter_gui as gui
from condition_recommender.conversion.artifacts import (
    RecommendationArtifactProgress,
)


def test_review_window_accepts_folder_and_individual_file_inputs(
    qtbot, tmp_path: Path
) -> None:
    source = tmp_path / "datasets"
    nested = source / "nested"
    nested.mkdir(parents=True)
    (source / "first.csv").write_text("reaction_smiles\n", encoding="utf-8")
    (nested / "second.CSV").write_text("reaction_smiles\n", encoding="utf-8")
    selected = tmp_path / "selected.observations.jsonl.gz"
    selected.write_bytes(b"")
    window = gui.GenericReactionReviewWindow()
    qtbot.addWidget(window)

    window.add_source_inputs((str(source), str(selected), str(selected)))

    assert window.source_list.count() == 2
    assert window.source_inputs() == (
        str(source.resolve()),
        str(selected.resolve()),
    )
    assert "Found 3 conversion input file(s)" in window.source_summary.text()
    assert "1 selected file(s) and 1 folder(s)" in window.source_summary.text()
    assert window.status_box.isReadOnly()
    assert not window.cancel_button.isEnabled()
    assert Path(window.output_edit.text()) == gui.DEFAULT_OUTPUT_FOLDER
    assert window.shard_size_spin.value() == 1_000
    assert not window.build_index_check.isChecked()
    assert not window.build_retrosynthesis_check.isChecked()
    assert window.build_retrosynthesis_check.objectName() == "buildRetrosynthesis"
    assert window.combine_button.objectName() == "combineIndexButton"
    assert window.combine_button.isEnabled()
    assert window.batch_name_edit.objectName() == "batchName"
    assert "saved batch" in window.batch_summary.text()
    assert window.conversion_mode_combo.currentData() == "compact"
    assert window.conversion_mode_combo.objectName() == "conversionMode"
    assert not window.use_rxnmapper_check.isChecked()
    assert window.use_rxnmapper_check.objectName() == "useRxnMapper"
    options_layout = window.options_widget.layout()
    assert isinstance(options_layout, QtWidgets.QHBoxLayout)
    assert options_layout.indexOf(window.use_rxnmapper_check) >= 0
    assert options_layout.indexOf(window.build_index_check) >= 0
    assert options_layout.indexOf(window.build_retrosynthesis_check) >= 0
    assert not any(
        label.text().startswith("Outputs: shard_manifest.json")
        for label in window.findChildren(QtWidgets.QLabel)
    )
    assert window.worker_count_spin.isEnabled()
    window.use_rxnmapper_check.setChecked(True)
    assert window.worker_count_spin.value() == 1
    assert not window.worker_count_spin.isEnabled()

    window.build_retrosynthesis_check.setChecked(True)
    assert window.build_index_check.isChecked()
    assert not window.build_index_check.isEnabled()
    window.build_retrosynthesis_check.setChecked(False)
    assert not window.build_index_check.isChecked()
    assert window.build_index_check.isEnabled()


def test_review_worker_forwards_progress_and_result(monkeypatch) -> None:
    progress = RecommendationArtifactProgress(
        phase="canonical_completed",
        source_file_count=2,
        shard_count=2,
        row_count=10,
        message="Finished 10 reaction(s).",
    )

    def fake_conversion(
        source,
        output,
        *,
        shard_size,
        workers,
        build_fast_index,
        use_rxnmapper,
        progress_callback,
        cancel_check,
    ):
        assert source == "source"
        assert output == "output"
        assert shard_size == 500
        assert workers == 2
        assert build_fast_index
        assert use_rxnmapper
        assert not cancel_check()
        progress_callback(progress)
        return {"record_count": 10, "output_dir": output}

    monkeypatch.setattr(
        gui,
        "build_recommendation_artifacts",
        fake_conversion,
    )
    worker = gui.ReviewConversionWorker(
        "source",
        "output",
        shard_size=500,
        workers=2,
        build_fast_index=True,
        use_rxnmapper=True,
    )
    updates = []
    results = []
    worker.progress.connect(updates.append)
    worker.finished.connect(
        lambda success, report, error: results.append(
            (success, report, error)
        )
    )

    worker.run()

    assert updates == [progress]
    assert results == [
        (True, {"record_count": 10, "output_dir": "output"}, "")
    ]


def test_saved_batch_worker_can_save_then_combine(monkeypatch) -> None:
    calls = []

    def fake_save(source, output, **options):
        calls.append(
            (
                "save",
                source,
                output,
                options["batch_name"],
                options["conversion_mode"],
            )
        )
        return {
            "record_count": 3,
            "batch_name": "batch-a",
            "batch_dir": "output/batches/batch-a",
        }

    def fake_combine(output, **options):
        calls.append(("combine", output))
        return {"record_count": 5, "output_dir": output}

    monkeypatch.setattr(gui, "save_recommendation_batch", fake_save)
    monkeypatch.setattr(gui, "combine_saved_recommendation_batches", fake_combine)
    worker = gui.SavedBatchWorker(
        ("first.csv", "second.csv"),
        "output",
        batch_name="batch-a",
        combine_after_save=True,
        use_rxnmapper=False,
        conversion_mode="compact",
    )
    results = []
    worker.finished.connect(
        lambda success, report, error: results.append((success, report, error))
    )

    worker.run()

    assert calls == [
        (
            "save",
            ("first.csv", "second.csv"),
            "output",
            "batch-a",
            "compact",
        ),
        ("combine", "output"),
    ]
    assert results[0][0]
    assert results[0][1]["saved_batch"]["record_count"] == 3
    assert results[0][1]["combined"]["record_count"] == 5


def test_saved_batch_worker_builds_selected_retrosynthesis_mode_from_all_batches(
    monkeypatch,
) -> None:
    calls = []

    def fake_save(source, output, **options):
        calls.append(("save", output))
        return {
            "record_count": 3,
            "batch_name": "batch-a",
            "batch_dir": "output/batches/batch-a",
        }

    def fake_combine(output, **options):
        calls.append(("combine", output))
        return {"record_count": 5, "output_dir": output}

    def fake_retrosynthesis(source, output, **options):
        calls.append(
            (
                "retrosynthesis",
                source,
                output,
                options["library_mode"],
                options["workers"],
            )
        )
        return {"operator_count": 2, "library_mode": "compact"}

    monkeypatch.setattr(gui, "save_recommendation_batch", fake_save)
    monkeypatch.setattr(gui, "combine_saved_recommendation_batches", fake_combine)
    monkeypatch.setattr(
        gui,
        "build_retrosynthesis_operator_artifacts",
        fake_retrosynthesis,
    )
    worker = gui.SavedBatchWorker(
        ("first.csv",),
        "library/compact",
        workers=3,
        combine_after_save=False,
        conversion_mode="compact",
        build_retrosynthesis=True,
        retrosynthesis_output_folder="operators",
    )
    results = []
    worker.finished.connect(
        lambda success, report, error: results.append((success, report, error))
    )

    worker.run()

    assert calls == [
        ("save", "library/compact"),
        ("combine", "library/compact"),
        (
            "retrosynthesis",
            "library/compact",
            "operators",
            "compact",
            3,
        ),
    ]
    assert results[0][0]
    assert results[0][1]["retrosynthesis"]["operator_count"] == 2


def test_combine_worker_resumes_incomplete_batches_first(
    monkeypatch, tmp_path: Path
) -> None:
    manifest = tmp_path / "batches" / "partial" / "shard_manifest.json"
    calls = []

    monkeypatch.setattr(
        gui,
        "incomplete_saved_conversion_batches",
        lambda output: (manifest,),
    )

    def fake_resume(path, **options):
        calls.append(("resume", path, options["workers"]))
        return {"record_count": 3}

    def fake_combine(output, **options):
        calls.append(("combine", output))
        return {"record_count": 3, "output_dir": output}

    monkeypatch.setattr(gui, "resume_saved_conversion_batch", fake_resume)
    monkeypatch.setattr(gui, "combine_saved_recommendation_batches", fake_combine)
    worker = gui.CombineSavedBatchesWorker(
        "output",
        resume_incomplete=True,
        workers=4,
    )
    results = []
    worker.finished.connect(
        lambda success, report, error: results.append((success, report, error))
    )

    worker.run()

    assert calls == [("resume", manifest, 4), ("combine", "output")]
    assert results[0][0]
    assert results[0][1]["resumed_batch_count"] == 1


def test_retrosynthesis_builder_uses_unlimited_default_config(
    monkeypatch, tmp_path: Path
) -> None:
    calls = []
    updates = []

    def fake_build(source, output, **options):
        calls.append((Path(source), Path(output), options))
        options["progress_callback"](
            {
                "phase": "complete",
                "completed_shards": 2,
                "total_shards": 2,
                "source_rows": 12,
                "operator_count": 4,
            }
        )
        return object(), {
            "source_rows": 12,
            "operator_count": 4,
            "template_count": 6,
            "library_path": str(Path(output) / "operator_library_v3.json.gz"),
        }

    monkeypatch.setattr(gui, "build_full_scale_operator_library", fake_build)
    report = gui.build_retrosynthesis_operator_artifacts(
        tmp_path / "recommendations" / "full",
        tmp_path / "operators",
        library_mode="full",
        workers=4,
        progress_callback=updates.append,
    )

    source, output, options = calls[0]
    assert source == tmp_path / "recommendations" / "full"
    assert output == tmp_path / "operators" / "full"
    assert options["workers"] == 4
    assert options["config"].max_rows_per_shard is None
    assert report["library_mode"] == "full"
    assert report["output_dir"] == str(output.resolve())
    assert updates[0].phase == "retrosynthesis_complete"
    assert updates[0].row_count == 12
