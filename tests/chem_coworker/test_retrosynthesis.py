"""Tests for one-step retrosynthesis in the coworker application layer."""

from __future__ import annotations

from dataclasses import dataclass, field
from io import StringIO

from rich.console import Console

from core_retrosynthesis import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
    StrategyProposal,
    group_strategy_candidates,
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
from chem_coworker.retrosynthesis_review import (
    LLMRetrosynthesisReviewer,
    RetrosynthesisCandidatePayload,
    RetrosynthesisReviewPayload,
    RetrosynthesisReviewTransportResult,
)


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
    monkeypatch,
) -> None:
    candidates = (
        _candidate(1, precursor="CCBr.N"),
        _candidate(2, precursor="CCI.N"),
    )
    monkeypatch.setattr(
        "chem_coworker.retrosynthesis.disconnect_operator_ladder",
        lambda *args, **kwargs: candidates,
    )
    coworker = RetrosynthesisCoworker(
        library=GenericTemplateLibrary(
            templates=(),
            source_row_count=0,
            accepted_observation_count=0,
            rejection_counts={},
            definition={},
        ),
        reviewer=FakeReviewer(),
    )

    response = coworker.disconnect(
        RetrosynthesisRequest(
            target_smiles="CCN",
            top_k=1,
            include_conditions=False,
            review=ConditionReviewSettings(mode="always"),
        )
    )

    assert response.valid is True
    assert len(response.strategies) == 2
    assert all(
        item.representative.forward_validation_status == "verified_signature"
        for item in response.strategies
    )
    assert response.review is not None
    assert "One-step retrosynthesis strategies" in response.answer


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
    first, second = group_strategy_candidates(
        (
            _candidate(1, precursor="CCBr.N"),
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
