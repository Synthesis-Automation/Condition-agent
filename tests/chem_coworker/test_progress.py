"""Tests for concise long-running CLI progress."""

from __future__ import annotations

import json
from io import StringIO

import pytest

from chem_coworker._cli import app as cli_app
from chem_coworker._cli.app import _parser
from chem_coworker._cli.progress import ProgressTask, concise_progress
from chem_coworker.contracts import (
    MultistepRetrosynthesisResponse,
    RetrosynthesisResponse,
)


class _TtyBuffer(StringIO):
    def isatty(self) -> bool:
        return True


def test_redirected_progress_emits_only_start_and_measured_result() -> None:
    stream = StringIO()

    with concise_progress("Searching bounded routes", stream=stream) as progress:
        progress.set_result("2 solved routes")

    lines = stream.getvalue().splitlines()
    assert lines[0] == "Searching bounded routes..."
    assert lines[1].startswith("Searching bounded routes - done in ")
    assert lines[1].endswith("; 2 solved routes")


def test_tty_progress_reuses_one_line_and_reports_failure() -> None:
    stream = _TtyBuffer()

    with pytest.raises(RuntimeError, match="search failed"):
        with ProgressTask(
            "Forward-validating disconnections",
            stream=stream,
            refresh_interval=60.0,
        ):
            raise RuntimeError("search failed")

    rendered = stream.getvalue()
    assert rendered.startswith("\rForward-validating disconnections - ")
    assert "\rForward-validating disconnections - failed after " in rendered
    assert rendered.endswith("\n")


def test_progress_can_be_disabled_and_cli_exposes_opt_out() -> None:
    stream = StringIO()

    with concise_progress("Hidden", enabled=False, stream=stream):
        pass

    assert stream.getvalue() == ""
    assert _parser().parse_args([]).progress is True
    assert _parser().parse_args(["--no-progress"]).progress is False


def _install_cli_fakes(monkeypatch: pytest.MonkeyPatch) -> None:
    class _ConditionCoworker:
        startup_warnings = ("using test index",)
        recommender = object()

    class _ConditionFactory:
        @staticmethod
        def from_default(**kwargs):
            del kwargs
            return _ConditionCoworker()

    class _RetrosynthesisCoworker:
        startup_warnings = ()

        def disconnect(self, request):
            return RetrosynthesisResponse(
                request=request,
                valid=True,
                answer="retro result",
            )

    class _RetrosynthesisFactory:
        @staticmethod
        def from_default(**kwargs):
            del kwargs
            return _RetrosynthesisCoworker()

    class _MultistepCoworker:
        def plan(self, request):
            return MultistepRetrosynthesisResponse(
                request=request,
                valid=False,
                error="test result",
                answer="test result",
            )

    class _MultistepFactory:
        @staticmethod
        def from_retrosynthesis_coworker(coworker, **kwargs):
            del coworker, kwargs
            return _MultistepCoworker()

    monkeypatch.setattr(cli_app, "ConditionCoworker", _ConditionFactory)
    monkeypatch.setattr(cli_app, "RetrosynthesisCoworker", _RetrosynthesisFactory)
    monkeypatch.setattr(
        cli_app,
        "MultistepRetrosynthesisCoworker",
        _MultistepFactory,
    )


@pytest.mark.parametrize(
    ("mode", "expected_exit", "progress_text", "system"),
    (
        (
            "retro",
            0,
            "Generating and forward-validating disconnections",
            "chem_coworker.retrosynthesis.v1",
        ),
        (
            "multistep",
            2,
            "Searching bounded routes",
            "chem_coworker.multistep_retrosynthesis.v1",
        ),
    ),
)
def test_one_shot_progress_stays_on_stderr_and_json_stays_clean(
    mode: str,
    expected_exit: int,
    progress_text: str,
    system: str,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _install_cli_fakes(monkeypatch)

    exit_code = cli_app.main(["--mode", mode, "--review", "off", "--json", "CCN"])

    captured = capsys.readouterr()
    assert exit_code == expected_exit
    assert json.loads(captured.out)["system"] == system
    assert "Warning: using test index" in captured.err
    assert progress_text in captured.err
