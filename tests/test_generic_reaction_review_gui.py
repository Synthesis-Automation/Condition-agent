"""The desktop converter delegates one complete resumable release build."""

from pathlib import Path
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from app import reaction_converter_gui as gui


def test_builder_defaults_to_all_data_with_no_mode_selector(qtbot):
    window = gui.GenericReactionReviewWindow()
    qtbot.addWidget(window)
    assert Path(window.source_edit.text()) == gui.DEFAULT_INPUT_FOLDER
    assert Path(window.output_edit.text()) == gui.DEFAULT_OUTPUT_FOLDER
    assert gui.DEFAULT_OUTPUT_FOLDER.name == "processed_datasets"
    assert not hasattr(window, "conversion_mode_combo")
    assert window.shard_size_spin.value() == 500
    assert window.log.isReadOnly()
    assert not window.cancel_button.isEnabled()


def test_worker_passes_complete_build_settings_and_reports_publication(monkeypatch, tmp_path):
    seen = {}
    result = {"release_id": "one", "coverage": {"observation_count": 15}}
    def build(source, output, **kwargs):
        seen.update(source=source, output=output, workers=kwargs["workers"], shard_size=kwargs["shard_size"])
        assert not kwargs["cancel_check"]()
        kwargs["progress"]({"phase": "conversion", "rows": 15})
        return result
    monkeypatch.setattr(gui, "build_processed_datasets", build)
    worker = gui.DatasetBuildWorker(str(tmp_path), str(tmp_path / "output"), workers=3, shard_size=500)
    finished, progress = [], []
    worker.finished.connect(lambda *args: finished.append(args))
    worker.progress.connect(progress.append)
    worker.run()
    assert seen["workers"] == 3 and seen["shard_size"] == 500
    assert progress[0]["rows"] == 15
    assert finished == [(True, result, "")]


def test_cancelled_worker_keeps_resume_information(monkeypatch, tmp_path):
    def build(*args, **kwargs):
        assert kwargs["cancel_check"]()
        raise InterruptedError("Completed shards remain reusable")
    monkeypatch.setattr(gui, "build_processed_datasets", build)
    worker = gui.DatasetBuildWorker(str(tmp_path), str(tmp_path / "out"), workers=1, shard_size=500)
    results = []
    worker.finished.connect(lambda *args: results.append(args))
    worker.request_cancel()
    worker.run()
    assert results[0][0] is False and results[0][1]["cancelled"]
    assert "reusable" in results[0][2]


def test_missing_source_does_not_start_worker(qtbot, tmp_path):
    window = gui.GenericReactionReviewWindow()
    qtbot.addWidget(window)
    window.source_edit.setText(str(tmp_path / "missing"))
    window.start_conversion()
    assert window.thread is None
    assert "does not exist" in window.status_label.text()
