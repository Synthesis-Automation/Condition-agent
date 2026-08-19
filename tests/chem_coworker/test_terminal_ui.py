"""Tests for the enhanced prompt-toolkit and Rich terminal frontend."""

from __future__ import annotations

from io import StringIO

from prompt_toolkit.document import Document
from rich.console import Console

from condition_recommender import GenericRecommendationResult

from chem_coworker._cli.terminal_ui import (
    RichResponseRenderer,
    SlashCommandCompleter,
)
from chem_coworker.contracts import ConditionRequest, ConditionResponse


def test_slash_completion_covers_commands_and_profile_arguments() -> None:
    completer = SlashCommandCompleter(("default", "yield"))

    commands = list(completer.get_completions(Document("/pro"), None))
    profiles = list(completer.get_completions(Document("/profile yi"), None))

    assert [item.text for item in commands] == ["/profile"]
    assert [item.text for item in profiles] == ["yield"]


def test_rich_renderer_preserves_warning_text_as_literal_content() -> None:
    stream = StringIO()
    console = Console(file=stream, force_terminal=False, width=100)
    renderer = RichResponseRenderer(console)
    request = ConditionRequest(reaction_smiles="invalid")
    result = GenericRecommendationResult(
        query_reaction_smiles="invalid",
        valid=False,
        error="INVALID_REACTION",
        warnings=("Review [Pd] mapping",),
    )

    renderer(
        ConditionResponse(request=request, result=result, answer="unused"),
        False,
    )

    rendered = stream.getvalue()
    assert "Recommendation unavailable" in rendered
    assert "INVALID_REACTION" in rendered
    assert "Review [Pd] mapping" in rendered


def test_rich_renderer_supports_full_json_mode() -> None:
    stream = StringIO()
    console = Console(file=stream, force_terminal=False, width=100)
    renderer = RichResponseRenderer(console)
    request = ConditionRequest(reaction_smiles="C.N>>CN")
    result = GenericRecommendationResult(
        query_reaction_smiles="C.N>>CN",
        valid=True,
        named_family=None,
    )

    renderer(
        ConditionResponse(request=request, result=result, answer="unused"),
        True,
    )

    rendered = stream.getvalue()
    assert '"system": "condition_coworker.v1"' in rendered
    assert '"query_reaction_smiles": "C.N>>CN"' in rendered
