"""Regressions for the bounded multistep coworker and route review."""

from __future__ import annotations

from contextlib import nullcontext
from types import SimpleNamespace

import pytest

from cas_tools import MoleculeIndexMatch, molecule_identity
from chem_coworker.contracts import (
    ConditionReviewSettings,
    MultistepRetrosynthesisRequest,
    MultistepRetrosynthesisResponse,
)
from chem_coworker._cli.interactive import InteractiveSession
from chem_coworker.multistep import MultistepRetrosynthesisCoworker
from chem_coworker.multistep_rendering import ordered_multistep_routes
from chem_coworker.multistep_review import (
    LLMMultistepReviewer,
    MultistepReviewPayload,
    MultistepReviewTransportResult,
    OpenAICompatibleMultistepReviewTransport,
)
from core_retrosynthesis import GenericDisconnectionCandidate, plan_multistep_routes


class _Stock:
    def __init__(self, known: tuple[str, ...] = ()) -> None:
        self.known = {
            molecule_identity(value).canonical_smiles  # type: ignore[union-attr]
            for value in known
        }

    def lookup(self, identity, *, provenance_limit: int = 5):
        if identity.canonical_smiles not in self.known:
            return None
        return MoleculeIndexMatch(
            canonical_smiles=identity.canonical_smiles,
            inchi_key=identity.inchi_key,
            occurrence_count=1,
            source_records=(
                {
                    "supplier": "test",
                    "source_role": "starting_material",
                    "stock_evidence": "supplier_in_stock",
                    "terminal_eligible": "true",
                },
            ),
        )


def _candidate(product: str, precursors: str, score: float):
    return GenericDisconnectionCandidate(
        target_smiles=product,
        precursor_smiles=precursors,
        proposed_reaction_smiles=f"{precursors}>>{product}",
        transformation_kind="test transformation",
        abstraction_level="L2",
        compiler_engine="reaction_core",
        template_id=f"template:{product}:{precursors}",
        score=score,
        context_similarity=score,
        product_similarity=score,
        precursor_similarity=score,
        template_specificity=1.0,
        independent_reference_support=2,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key=precursors,
        precedent_reaction_ids=("precedent",),
        operator_id="operator",
        realization_id=precursors,
        operator_signature="signature",
        synthon_signature=precursors,
    )


def _routes():
    expansions = {
        "CCCCCCCC": (
            _candidate("CCCCCCCC", "C.CCCCCC", 0.9),
            _candidate("CCCCCCCC", "CC.CCCCC", 0.8),
        ),
        "CCCCCC": (_candidate("CCCCCC", "CCN.CCO", 0.85),),
        "CCCCC": (_candidate("CCCCC", "CN.CCC", 0.75),),
    }

    def expand(product: str, top_k: int):
        return tuple(expansions.get(product, ()))[:top_k]

    return plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _Stock(("C", "CC", "CCN", "CCO", "CN", "CCC")),
        max_depth=2,
        molecular_weight_threshold=20.0,
        top_k_routes=2,
        expander=expand,
    )


def _protic_conflict_route():
    target = molecule_identity(
        "OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2"
    ).canonical_smiles  # type: ignore[union-attr]
    precursors = (
        "Brc1ccccc1.O=C(c1ccccc1)C12CC[N+](CCO)(CC1)CC2"
    )
    candidate = _candidate(target, precursors, 0.9)

    def expand(product: str, top_k: int):
        return (candidate,)[:top_k] if product == target else ()

    result = plan_multistep_routes(
        target,
        object(),
        _Stock(tuple(precursors.split("."))),
        max_depth=2,
        molecular_weight_threshold=20.0,
        top_k_routes=1,
        expander=expand,
    )
    assert len(result.routes) == 1
    return result.routes[0]


class _Transport:
    def __init__(self, payload: dict) -> None:
        self.payload = payload
        self.packet = None

    def complete(self, evidence_packet, settings):
        self.packet = evidence_packet
        return MultistepReviewTransportResult(
            payload=MultistepReviewPayload.model_validate(self.payload),
            response_id="response-1",
            input_tokens=10,
            output_tokens=20,
        )


def test_multistep_request_enforces_planner_depth_contract() -> None:
    with pytest.raises(ValueError, match="max_depth"):
        MultistepRetrosynthesisRequest(target_smiles="CC", max_depth=4)  # type: ignore[arg-type]


def test_llm_review_uses_short_aliases_and_only_changes_presentation_order() -> None:
    result = _routes()
    assert len(result.routes) == 2
    transport = _Transport(
        {
            "summary": "Second route better matches the preference.",
            "routes": [
                {
                    "review_id": "route-1",
                    "suggested_rank": 2,
                    "verdict": "downrank",
                    "issue_codes": ["user_preference_mismatch"],
                    "evidence_ids": ["evidence.route-1.step-1"],
                    "rationale": "The route is less aligned with the supplied preference.",
                    "confidence": 0.7,
                },
                {
                    "review_id": "route-2",
                    "suggested_rank": 1,
                    "verdict": "keep",
                    "issue_codes": [],
                    "evidence_ids": ["query.guidance", "evidence.route-2"],
                    "rationale": "The route better matches the supplied preference.",
                    "confidence": 0.8,
                },
            ],
            "questions": [],
        }
    )
    reviewer = LLMMultistepReviewer(transport)
    review = reviewer.review(
        result.routes,
        "solved",
        ConditionReviewSettings(mode="always", max_candidates=2),
        "Prefer a convergent route.",
    )

    assert review.status == "completed"
    assert review.presentation_route_ids == (
        result.routes[1].route_id,
        result.routes[0].route_id,
    )
    assert transport.packet["review_route_ids"] == ["route-1", "route-2"]
    assert transport.packet["query"]["strategic_guidance"]["text"] == (
        "Prefer a convergent route."
    )
    assert result.routes[0].route_id != "route-1"


def test_llm_review_rejects_unknown_evidence_ids_without_losing_routes() -> None:
    result = _routes()
    transport = _Transport(
        {
            "summary": "Review.",
            "routes": [
                {
                    "review_id": f"route-{index}",
                    "suggested_rank": index,
                    "verdict": "keep",
                    "issue_codes": [],
                    "evidence_ids": ["invented.evidence"],
                    "rationale": "No concern was found.",
                    "confidence": 0.5,
                }
                for index in range(1, len(result.routes) + 1)
            ],
            "questions": [],
        }
    )

    review = LLMMultistepReviewer(transport).review(
        result.routes,
        "solved",
        ConditionReviewSettings(mode="always", max_candidates=2),
    )

    assert review.status == "failed"
    assert "unknown evidence IDs" in (review.warning or "")
    assert review.presentation_route_ids == tuple(
        route.route_id for route in result.routes
    )


def test_graph_derived_route_issue_triggers_review_and_is_citable() -> None:
    route = _protic_conflict_route()
    transport = _Transport(
        {
            "summary": "The deterministic route issue needs resolution.",
            "routes": [
                {
                    "review_id": "route-1",
                    "suggested_rank": 1,
                    "verdict": "downrank",
                    "issue_codes": [
                        "cross_step_functional_group_compatibility"
                    ],
                    "evidence_ids": ["evidence.route-1.issue-1"],
                    "rationale": "A code-derived protic-quench warning is present.",
                    "confidence": 0.9,
                }
            ],
            "questions": [],
        }
    )

    review = LLMMultistepReviewer(transport).review(
        (route,),
        "solved",
        ConditionReviewSettings(mode="auto", max_candidates=1),
    )

    assert review.status == "completed"
    deterministic_issues = transport.packet["routes"][0][
        "deterministic_issues"
    ]
    assert deterministic_issues[0]["kind"] == "reaction_compatibility"
    assert deterministic_issues[0]["severity"] == "strong"
    assert deterministic_issues[0]["evidence_id"] in transport.packet[
        "allowed_evidence_ids"
    ]


def test_aliyun_multistep_message_explicitly_requests_json(monkeypatch) -> None:
    requests = []
    response = SimpleNamespace(
        choices=[
            SimpleNamespace(
                message=SimpleNamespace(
                    content=MultistepReviewPayload(
                        summary="No routes supplied.", routes=[], questions=[]
                    ).model_dump_json()
                ),
                finish_reason="stop",
            )
        ],
        usage=SimpleNamespace(prompt_tokens=10, completion_tokens=5),
        id="response-json",
    )

    def create(**kwargs):
        requests.append(kwargs)
        return response

    client = SimpleNamespace(
        chat=SimpleNamespace(completions=SimpleNamespace(create=create))
    )
    monkeypatch.setenv("ALIYUN_API_KEY", "test-key")
    transport = OpenAICompatibleMultistepReviewTransport(
        client_factory=lambda **_: client
    )

    result = transport.complete(
        {},
        ConditionReviewSettings(
            mode="always", provider="aliyun", model="glm-5.2"
        ),
    )

    assert result.response_id == "response-json"
    assert requests[0]["response_format"] == {"type": "json_object"}
    message_text = " ".join(
        str(message["content"]) for message in requests[0]["messages"]
    )
    assert "json" in message_text.casefold()


def test_coworker_preserves_search_result_and_applies_review_order(
    monkeypatch, tmp_path
) -> None:
    result = _routes()
    stock_path = tmp_path / "stock.sqlite"
    stock_path.touch()
    captured = {}

    def fake_plan(*args, **kwargs):
        captured.update(kwargs)
        return result

    class _Reviewer:
        def review(self, routes, route_kind, settings, guidance=""):
            transport = _Transport(
                {
                    "summary": "Reviewed.",
                    "routes": [
                        {
                            "review_id": f"route-{index}",
                            "suggested_rank": len(routes) - index + 1,
                            "verdict": "keep",
                            "issue_codes": [],
                            "evidence_ids": [f"evidence.route-{index}"],
                            "rationale": "Existing route retained.",
                            "confidence": 0.6,
                        }
                        for index in range(1, len(routes) + 1)
                    ],
                    "questions": [],
                }
            )
            return LLMMultistepReviewer(transport).review(
                routes, route_kind, settings, guidance
            )

    monkeypatch.setattr("chem_coworker.multistep.open_stock_lookup", lambda _: nullcontext(_Stock()))
    monkeypatch.setattr("chem_coworker.multistep.plan_multistep_routes", fake_plan)
    coworker = MultistepRetrosynthesisCoworker(
        library=object(),  # type: ignore[arg-type]
        stock_path=stock_path,
        reviewer=_Reviewer(),
    )
    response = coworker.plan(
        MultistepRetrosynthesisRequest(
            target_smiles="CCCCCCCC",
            top_k=2,
            max_depth=2,
            include_conditions=False,
            review=ConditionReviewSettings(mode="always", max_candidates=2),
        )
    )

    assert response.valid is True
    assert response.result is result
    assert captured["max_depth"] == 2
    assert ordered_multistep_routes(response)[0].route_id == result.routes[1].route_id
    assert "original rank 2" in response.answer


def test_interactive_multistep_controls_flow_into_request() -> None:
    class _ConditionCoworker:
        def prepare_reaction(self, reaction_smiles):
            return {"requirements": ()}

        def recommend(self, request):
            raise AssertionError("condition mode should not run")

    class _MultistepCoworker:
        request = None

        def plan(self, request):
            self.request = request
            return MultistepRetrosynthesisResponse(
                request=request,
                valid=False,
                error="test response",
                answer="test response",
            )

    values = iter(
        (
            "/mode multistep",
            "/depth 2",
            "/per-step 4",
            "/beam 9",
            "/expansions 12",
            "/guidance avoid palladium and minimize protecting groups",
            "CCCCCCCC",
            "/quit",
        )
    )
    planner = _MultistepCoworker()
    session = InteractiveSession(
        _ConditionCoworker(),
        multistep_coworker=planner,
        input_fn=lambda _: next(values),
        output_fn=lambda _: None,
    )

    assert session.run() == 0
    assert planner.request is not None
    assert planner.request.max_depth == 2
    assert planner.request.per_step_top_k == 4
    assert planner.request.beam_width == 9
    assert planner.request.max_expansions == 12
    assert planner.request.strategic_guidance.startswith("avoid palladium")
