"""Tests for the reusable interactive condition coworker session."""

from __future__ import annotations

from dataclasses import dataclass, field
from condition_recommender import GenericRecommendationResult

from chem_coworker._cli.app import _parser
from chem_coworker._cli.interactive import InteractiveSession
from chem_coworker.contracts import ConditionRequest, ConditionResponse


@dataclass
class FakeCoworker:
    proposal: dict[str, object] = field(
        default_factory=lambda: {"status": "not_required", "requirements": ()}
    )
    requests: list[ConditionRequest] = field(default_factory=list)

    def prepare_reaction(self, reaction_smiles: str) -> dict[str, object]:
        return self.proposal

    def recommend(self, request: ConditionRequest) -> ConditionResponse:
        self.requests.append(request)
        result = GenericRecommendationResult(
            query_reaction_smiles=request.reaction_smiles,
            valid=True,
            named_family=None,
            retrieval_level="reaction_facet_exact",
        )
        return ConditionResponse(request=request, result=result, answer="test answer")


def _inputs(*values: str):
    remaining = iter(values)

    def read(_: str) -> str:
        return next(remaining)

    return read


def test_no_positional_reaction_selects_interactive_mode() -> None:
    args = _parser().parse_args([])

    assert args.reaction_smiles is None
    assert args.index is None


def test_session_updates_settings_and_recommends_repeatedly() -> None:
    coworker = FakeCoworker()
    output: list[str] = []
    session = InteractiveSession(
        coworker,
        input_fn=_inputs(
            "/top-k 2",
            "/profile yield",
            "C.N>>CN",
            "C.O>>CO",
            "/quit",
        ),
        output_fn=output.append,
    )

    exit_code = session.run()

    assert exit_code == 0
    assert [request.reaction_smiles for request in coworker.requests] == [
        "C.N>>CN",
        "C.O>>CO",
    ]
    assert all(request.top_k == 2 for request in coworker.requests)
    assert all(request.ranking_profile == "yield" for request in coworker.requests)
    assert output.count("test answer") == 2
    assert output[-1] == "Goodbye."


def test_session_guides_completion_option_selection() -> None:
    coworker = FakeCoworker(
        proposal={
            "status": "confirmation_recommended",
            "requirements": (
                {
                    "requirement_id": "source-1",
                    "rooted_fragment_smiles": "*N=[N+]=[N-]",
                    "options": (
                        {
                            "option_id": "sodium-azide",
                            "option_kind": "registered_substance",
                            "display_name": "Sodium azide",
                        },
                    ),
                },
            ),
        }
    )
    output: list[str] = []
    session = InteractiveSession(
        coworker,
        input_fn=_inputs("reaction-with-missing-source", "1", "/quit"),
        output_fn=output.append,
    )

    session.run()

    assert len(coworker.requests) == 1
    choices = coworker.requests[0].completion_choices
    assert len(choices) == 1
    assert choices[0].requirement_id == "source-1"
    assert choices[0].option_id == "sodium-azide"
    assert any("Sodium azide" in line for line in output)


def test_cancelled_completion_does_not_recommend() -> None:
    coworker = FakeCoworker(
        proposal={
            "requirements": (
                {
                    "requirement_id": "source-1",
                    "fragment_key": "azide",
                    "options": (),
                },
            )
        }
    )
    output: list[str] = []
    session = InteractiveSession(
        coworker,
        input_fn=_inputs("incomplete-reaction", "q", "/quit"),
        output_fn=output.append,
    )

    session.run()

    assert coworker.requests == []
    assert "Recommendation cancelled." in output


def test_json_toggle_and_save_previous_result(tmp_path) -> None:
    target = tmp_path / "recommendation.json"
    coworker = FakeCoworker()
    output: list[str] = []
    session = InteractiveSession(
        coworker,
        input_fn=_inputs(
            "/json on",
            "C.N>>CN",
            f"/save {target}",
            "/quit",
        ),
        output_fn=output.append,
    )

    session.run()

    assert target.is_file()
    saved = target.read_text(encoding="utf-8")
    assert '"query_reaction_smiles": "C.N>>CN"' in saved
    assert any('"system": "condition_coworker.v1"' in line for line in output)


def test_clear_command_uses_terminal_callback() -> None:
    cleared: list[bool] = []
    session = InteractiveSession(
        FakeCoworker(),
        input_fn=_inputs("/clear", "/quit"),
        output_fn=lambda _: None,
        clear_fn=lambda: cleared.append(True),
    )

    session.run()

    assert cleared == [True]


def test_session_model_and_review_selectors_flow_into_request(
    tmp_path,
    monkeypatch,
) -> None:
    from chem_coworker._cli import config

    monkeypatch.setattr(config, "CONFIG_PATH", tmp_path / "config.json")
    monkeypatch.setattr(
        "chem_coworker._cli.interactive.save_config",
        config.save_config,
    )
    coworker = FakeCoworker()
    output: list[str] = []
    session = InteractiveSession(
        coworker,
        input_fn=_inputs(
            "/model gpt-5.6-sol",
            "/review always",
            "/reasoning high",
            "/review-limit 3",
            "C.N>>CN",
            "/quit",
        ),
        output_fn=output.append,
    )

    session.run()

    settings = coworker.requests[0].review
    assert settings.model == "gpt-5.6-sol"
    assert settings.provider == "openai"
    assert settings.mode == "always"
    assert settings.reasoning_effort == "high"
    assert settings.max_candidates == 3
    assert (tmp_path / "config.json").is_file()
