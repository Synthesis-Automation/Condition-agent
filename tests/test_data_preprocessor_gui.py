from pathlib import Path

from app import data_preprocessor_gui as gui
from condition_recommender.ingestion import PreprocessingProgress


def test_preprocessor_window_adds_files_and_exposes_incremental_defaults(
    qtbot, tmp_path: Path
) -> None:
    source = tmp_path / "source.csv"
    source.write_text("unknown\nvalue\n", encoding="utf-8")
    window = gui.SourceDataPreprocessorWindow()
    qtbot.addWidget(window)

    window.add_source_files((str(source), str(source)))

    assert window.source_list.count() == 1
    assert window.source_files() == (str(source.resolve()),)
    assert "Selected 1 file(s)" in window.source_summary.text()
    assert "unrecognized: 1" in window.source_summary.text()
    assert Path(window.output_edit.text()) == gui.DEFAULT_OUTPUT_FOLDER
    assert window.adapter_combo.currentData() == "auto"
    assert not window.force_check.isChecked()
    assert window.status_box.isReadOnly()
    assert not window.cancel_button.isEnabled()


def test_preprocessor_worker_forwards_progress_and_result(monkeypatch) -> None:
    progress = PreprocessingProgress(
        phase="completed",
        source_file="source.csv",
        file_number=1,
        file_count=1,
        row_count=10,
        message="Completed source.csv: 10 observations.",
    )

    def fake_preprocess(
        sources,
        output,
        *,
        adapter_id,
        force,
        progress_callback,
        cancel_check,
    ):
        assert sources == ("source.csv",)
        assert output == "output"
        assert adapter_id == "auto"
        assert not force
        assert not cancel_check()
        progress_callback(progress)
        return {"file_count": 1, "row_count": 10, "output_dir": output}

    monkeypatch.setattr(gui, "preprocess_files", fake_preprocess)
    worker = gui.DataPreprocessingWorker(("source.csv",), "output")
    updates = []
    results = []
    worker.progress.connect(updates.append)
    worker.finished.connect(
        lambda success, report, error: results.append((success, report, error))
    )

    worker.run()

    assert updates == [progress]
    assert results == [
        (True, {"file_count": 1, "row_count": 10, "output_dir": "output"}, "")
    ]
