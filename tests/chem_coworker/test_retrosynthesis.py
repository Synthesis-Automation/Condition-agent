"""Tests for one-step retrosynthesis in the coworker application layer."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from io import StringIO
import json

import pytest

from rich.console import Console

from condition_recommender.corpus_io import file_sha256

from core_retrosynthesis import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
    StrategyProposal,
    group_strategy_candidates,
    save_generic_library,
)

from chem_coworker._cli.interactive import InteractiveSession
from chem_coworker._cli.terminal_ui import RichResponseRenderer
from chem_coworker.contracts import (
    ConditionReviewSettings,
    RetrosynthesisRequest,
    RetrosynthesisResponse,
    RetrosynthesisReview,
)
from chem_coworker.retrosynthesis import RetrosynthesisCoworker
from chem_coworker.retrosynthesis_rendering import render_retrosynthesis
from chem_coworker.retrosynthesis_review import (
    LLMRetrosynthesisReviewer,
    RetrosynthesisCandidatePayload,
    RetrosynthesisReviewPayload,
    RetrosynthesisReviewTransportResult,
)
from tests.core_retrosynthesis_tests import test_focused_disconnection

library = test_focused_disconnection.library


def _published_retro_release(tmp_path):
    release = tmp_path / "releases" / "complete"
    release.mkdir(parents=True)
    artifact = release / "operator_library_v3.json.gz"
    save_generic_library(GenericTemplateLibrary(
        templates=(), source_row_count=7, accepted_observation_count=0,
        rejection_counts={}, definition={},
    ), artifact)
    manifest = {
        "schema_version": "processed_reaction_release.v1",
        "release_id": "complete", "build_complete": True,
        "artifacts": {"retro_library": {
            "relative_path": artifact.name,
            "size_bytes": artifact.stat().st_size,
            "sha256": file_sha256(artifact),
        }},
    }
    immutable = release / "manifest.json"
    immutable.write_text(json.dumps(manifest), encoding="utf-8")
    active = tmp_path / "manifest.json"
    active.write_text(json.dumps({
        "schema_version": "processed_reaction_active.v1",
        "release_manifest": "releases/complete/manifest.json",
        "sha256": file_sha256(immutable),
    }), encoding="utf-8")
    return active, immutable, artifact


def test_retrosynthesis_default_pins_declared_processed_library(tmp_path, monkeypatch):
    active, _, artifact = _published_retro_release(tmp_path)
    monkeypatch.setattr("chem_coworker.retrosynthesis.DEFAULT_RETROSYNTHESIS_LIBRARY_PATH", active)
    coworker = RetrosynthesisCoworker.from_default()
    assert coworker.library_path == artifact.resolve()
    assert coworker.library.source_row_count == 7
    active.write_text("{}", encoding="utf-8")
    assert coworker.library_path == artifact.resolve()
    assert coworker.library.source_row_count == 7


def test_retrosynthesis_default_rejects_unfinished_processed_release(tmp_path, monkeypatch):
    _, immutable, _ = _published_retro_release(tmp_path)
    manifest = json.loads(immutable.read_text(encoding="utf-8"))
    manifest["build_complete"] = False
    immutable.write_text(json.dumps(manifest), encoding="utf-8")
    monkeypatch.setattr("chem_coworker.retrosynthesis.DEFAULT_RETROSYNTHESIS_LIBRARY_PATH", immutable)
    with pytest.raises(ValueError, match="unfinished"):
        RetrosynthesisCoworker.from_default()


def _candidate(
    rank: int,
    *,
    precursor: str,
) -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles="CCN",
        precursor_smiles=precursor,
        proposed_reaction_smiles=f"{precursor}>>CCN",
        transformation_kind="substitution",
        abstraction_level="L2",
        compiler_engine="reaction_core",
        template_id=f"template-{rank}",
        score=0.9 - rank / 10,
        context_similarity=0.8,
        product_similarity=0.9,
        precursor_similarity=0.7,
        template_specificity=0.8,
        independent_reference_support=rank + 1,
        forward_validation_status="verified_signature",
        center_transition_key=f"center-{rank}",
        disconnection_site_key=f"SITE1:{rank}",
        precedent_reaction_ids=(f"precedent-{rank}",),
        operator_id=f"OP1:{rank}",
        realization_id=f"REAL2:{rank}",
        operator_signature=f"operator-{rank}",
        synthon_signature=f"SYN1:{rank}",
        diversity_rank=rank,
    )


def _strategy(rank: int, *, precursor: str) -> StrategyProposal:
    candidate = _candidate(rank, precursor=precursor)
    return group_strategy_candidates((candidate,), top_k_strategies=1)[0]


class FakeReviewer:
    def review(self, strategies, condition_evidence, settings):
        del condition_evidence, settings
        return RetrosynthesisReview(
            status="completed",
            provider="openai",
            model="gpt-5.6-terra",
            summary="reviewed",
            presentation_strategy_ids=tuple(
                item.strategy_id for item in reversed(strategies)
            ),
        )


def test_service_groups_validated_candidates_and_preserves_review_pool(
    library,
) -> None:
    coworker = RetrosynthesisCoworker(library=library, reviewer=FakeReviewer())
    response = coworker.disconnect(
        RetrosynthesisRequest(
            target_smiles="CCNC",
            top_k=1,
            include_conditions=False,
            review=ConditionReviewSettings(mode="always", max_candidates=2),
        )
    )

    assert response.valid is True
    assert len(response.strategies) == 2
    assert all(
        item.representative.forward_validation_status == "verified_signature"
        for item in response.strategies
    )
    assert response.review is not None
    assert response.review.presentation_strategy_ids == tuple(
        item.strategy_id for item in reversed(response.strategies)
    )
    assert [item.strategy_rank for item in response.strategies] == [1, 2]
    assert {item.representative.precursor_smiles for item in response.strategies} == {
        "CCBr.CN", "CBr.CCN",
    }
    diagnostics = response.to_dict()["search_diagnostics"]
    assert diagnostics["requested_strategy_count"] == diagnostics["returned_strategy_count"] == 2
    assert diagnostics["strategy_target_met"]
    assert diagnostics["scheduling"] == "similarity_first_l2_then_operator_site_round_robin"
    assert "One-step retrosynthesis strategies" in response.answer


@pytest.mark.parametrize("bond", [None, (1, 2)])
def test_service_reports_per_level_caps_and_strategy_shortfall(library, bond) -> None:
    response = RetrosynthesisCoworker(library).disconnect(
        RetrosynthesisRequest(
            target_smiles="CCNC", top_k=5,
            max_templates_to_apply=1, max_candidates_to_validate=1,
            max_realizations_per_strategy=1, include_l0=False, include_conditions=False,
            required_disconnection_bond=bond, focus_target_smiles="CCNC" if bond else None,
            review=ConditionReviewSettings(mode="off"),
        )
    )

    assert response.valid, response.error
    payload = response.to_dict()
    diagnostics = payload["search_diagnostics"]
    assert diagnostics["requested_strategy_count"] == 5
    assert diagnostics["returned_strategy_count"] == len(response.strategies) < 5
    assert not diagnostics["strategy_target_met"]
    assert diagnostics["budget_limited"]
    assert diagnostics["max_templates_per_level"] == diagnostics["max_validations_per_level"] == 1
    assert diagnostics["max_realizations_per_strategy"] == 1
    for level in diagnostics["level_diagnostics"].values():
        assert level["applied_template_count"] <= 1
        assert level["validation_attempt_count"] <= 1
    if bond:
        assert diagnostics["levels_attempted"] == ["L2"]
        assert payload["bond_focus"]["atom_ids"] == bond
        assert diagnostics["scheduling"] == "bond_focused_operator_ladder"
    else:
        assert diagnostics["levels_attempted"] == ["L2", "L1"]
        assert "bond_focus" not in payload


def test_service_returns_bounded_no_hit_with_search_evidence(library) -> None:
    response = RetrosynthesisCoworker(library).disconnect(
        RetrosynthesisRequest(
            target_smiles="C", top_k=1, include_conditions=False,
            review=ConditionReviewSettings(mode="off"),
        )
    )

    assert response.valid and not response.strategies
    assert "bounded search" in response.answer
    diagnostics = response.to_dict()["search_diagnostics"]
    assert diagnostics["levels_attempted"] == ["L2", "L1", "L0"]
    assert diagnostics["returned_strategy_count"] == diagnostics["validation_attempt_count"] == 0
    assert not diagnostics["strategy_target_met"]
    assert not diagnostics["budget_limited"]


@dataclass
class FakeTransport:
    payload: RetrosynthesisReviewPayload
    packet: dict[str, object] = field(default_factory=dict)

    def complete(self, evidence_packet, settings):
        del settings
        self.packet = dict(evidence_packet)
        return RetrosynthesisReviewTransportResult(
            payload=self.payload,
            response_id="retro-response",
        )


def test_llm_review_can_reorder_but_not_modify_validated_strategies() -> None:
    representative = _candidate(1, precursor="CCBr.N")
    alternate = replace(
        representative,
        precursor_smiles="CCCl.N",
        proposed_reaction_smiles="CCCl.N>>CCN",
        template_id="another-template",
        independent_reference_support=73,
        precedent_reaction_ids=("only-for-alternate",),
    )
    first, second = group_strategy_candidates(
        (
            representative,
            alternate,
            _candidate(2, precursor="CCI.N"),
        ),
        top_k_strategies=2,
    )
    payload = RetrosynthesisReviewPayload(
        summary="The second route has the cleaner selectivity case.",
        candidates=[
            RetrosynthesisCandidatePayload(
                review_id="strategy-1",
                suggested_rank=2,
                verdict="downrank",
                issue_codes=["chemoselectivity"],
                evidence_ids=["evidence.strategy-1"],
                rationale="Competing sites need evaluation.",
                confidence=0.8,
            ),
            RetrosynthesisCandidatePayload(
                review_id="strategy-2",
                suggested_rank=1,
                verdict="keep",
                issue_codes=[],
                evidence_ids=["evidence.strategy-2"],
                rationale="The supplied evidence supports this route.",
                confidence=0.75,
            ),
        ],
        questions=[],
    )
    transport = FakeTransport(payload)

    review = LLMRetrosynthesisReviewer(transport).review(
        (first, second),
        (),
        ConditionReviewSettings(mode="always"),
    )

    assert review.status == "completed"
    assert review.presentation_strategy_ids == (
        second.strategy_id,
        first.strategy_id,
    )
    assert transport.packet["review_strategy_ids"] == (
        "strategy-1",
        "strategy-2",
    )
    assert all(
        "STRAT1:" not in value
        for value in transport.packet["allowed_evidence_ids"]
    )
    assert first.representative.precursor_smiles == "CCBr.N"
    assert first.independent_reference_support == 73
    packet_candidate = transport.packet["strategies"][0]
    assert packet_candidate["independent_reference_support"] == 2
    assert packet_candidate["precedent_reaction_ids"] == ["precedent-1"]
    answer = render_retrosynthesis(RetrosynthesisResponse(
        request=RetrosynthesisRequest(target_smiles="CCN", include_conditions=False),
        valid=True, strategies=(first,),
    ))
    assert "2 template references" in answer and "precedent-1" in answer
    assert "only-for-alternate" not in answer


def test_llm_review_rejects_invented_evidence() -> None:
    strategy = _strategy(1, precursor="CCBr.N")
    payload = RetrosynthesisReviewPayload(
        summary="reviewed",
        candidates=[
            RetrosynthesisCandidatePayload(
                review_id="strategy-1",
                suggested_rank=1,
                verdict="flag",
                issue_codes=["other"],
                evidence_ids=["invented.evidence"],
                rationale="Unsupported claim.",
                confidence=0.9,
            )
        ],
        questions=[],
    )

    review = LLMRetrosynthesisReviewer(FakeTransport(payload)).review(
        (strategy,),
        (),
        ConditionReviewSettings(mode="always"),
    )

    assert review.status == "failed"
    assert review.presentation_strategy_ids == (strategy.strategy_id,)
    assert "unknown evidence IDs" in (review.warning or "")


class FakeConditionCoworker:
    @staticmethod
    def prepare_reaction(reaction_smiles):
        del reaction_smiles
        return {"requirements": ()}

    def recommend(self, request):
        raise AssertionError(f"condition mode was called unexpectedly: {request}")


@dataclass
class FakeRetrosynthesisCoworker:
    requests: list[RetrosynthesisRequest] = field(default_factory=list)

    def disconnect(self, request: RetrosynthesisRequest) -> RetrosynthesisResponse:
        self.requests.append(request)
        return RetrosynthesisResponse(
            request=request,
            valid=True,
            answer="retro answer",
        )


def test_interactive_mode_switch_dispatches_target_smiles() -> None:
    values = iter(("/mode retro", "CCN", "/quit"))
    retro = FakeRetrosynthesisCoworker()
    output: list[str] = []
    session = InteractiveSession(
        FakeConditionCoworker(),
        retrosynthesis_coworker=retro,
        input_fn=lambda _: next(values),
        output_fn=output.append,
    )

    assert session.run() == 0
    assert [item.target_smiles for item in retro.requests] == ["CCN"]
    assert "retro answer" in output


def test_rich_retro_details_identify_shown_and_original_ranks() -> None:
    first, second = group_strategy_candidates(
        (
            _candidate(1, precursor="CCBr.N"),
            _candidate(2, precursor="CCI.N"),
        ),
        top_k_strategies=2,
    )
    request = RetrosynthesisRequest(
        target_smiles="CCN",
        top_k=2,
        include_conditions=False,
    )
    response = RetrosynthesisResponse(
        request=request,
        valid=True,
        strategies=(first, second),
        review=RetrosynthesisReview(
            status="completed",
            provider="openai",
            model="gpt-5.6-terra",
            summary="Second strategy shown first.",
            presentation_strategy_ids=(second.strategy_id, first.strategy_id),
        ),
    )
    stream = StringIO()
    renderer = RichResponseRenderer(
        Console(file=stream, force_terminal=False, width=180)
    )

    renderer(response, False)

    rendered = stream.getvalue()
    assert "Shown" in rendered
    assert "Original" in rendered
    assert "Shown row 1 details · original rank 2" in rendered
    assert "Shown row 2 details · original rank 1" in rendered
