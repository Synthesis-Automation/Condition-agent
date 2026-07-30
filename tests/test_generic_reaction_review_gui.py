from pathlib import Path

from app import reaction_converter_gui as gui
from condition_recommender.conversion.artifacts import (
    RecommendationArtifactProgress,
)


def test_review_window_reports_recursive_csv_count(qtbot, tmp_path: Path) -> None:
    source = tmp_path / "datasets"
    nested = source / "nested"
    nested.mkdir(parents=True)
    (source / "first.csv").write_text("reaction_smiles\n", encoding="utf-8")
    (nested / "second.CSV").write_text("reaction_smiles\n", encoding="utf-8")
    window = gui.GenericReactionReviewWindow()
    qtbot.addWidget(window)

    window.source_edit.setText(str(source))
    window.refresh_source_summary()

    assert "Found 2 CSV file(s)" in window.source_summary.text()
    assert window.status_box.isReadOnly()
    assert not window.cancel_button.isEnabled()
    assert Path(window.output_edit.text()) == gui.DEFAULT_OUTPUT_FOLDER
    assert window.shard_size_spin.value() == 1_000
    assert window.build_index_check.isChecked()
    assert window.use_rxnmapper_check.isChecked()
    assert window.use_rxnmapper_check.objectName() == "useRxnMapper"
    assert window.worker_count_spin.value() == 1
    assert not window.worker_count_spin.isEnabled()
    window.use_rxnmapper_check.setChecked(False)
    assert window.worker_count_spin.isEnabled()
    window.use_rxnmapper_check.setChecked(True)
    assert window.worker_count_spin.value() == 1
    assert not window.worker_count_spin.isEnabled()


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
