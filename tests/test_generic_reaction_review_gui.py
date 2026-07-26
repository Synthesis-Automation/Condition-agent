from pathlib import Path

from app import generic_reaction_review_gui as gui
from condition_recommender.conversion.concise_review import ConciseReviewProgress


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


def test_review_worker_forwards_progress_and_result(monkeypatch) -> None:
    progress = ConciseReviewProgress(
        phase="completed",
        file_index=2,
        file_count=2,
        row_count=10,
        current_file="",
        message="Finished 10 reaction(s).",
    )

    def fake_conversion(source, output, *, progress_callback, cancel_check):
        assert source == "source"
        assert output == "review.csv"
        assert not cancel_check()
        progress_callback(progress)
        return {"row_count": 10, "output_path": output}

    monkeypatch.setattr(
        gui,
        "convert_dataset_folder_to_concise_review_csv",
        fake_conversion,
    )
    worker = gui.ReviewConversionWorker("source", "review.csv")
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
        (True, {"row_count": 10, "output_path": "review.csv"}, "")
    ]
