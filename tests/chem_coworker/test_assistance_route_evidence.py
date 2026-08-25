"""Evidence projections for one-step and multistep deterministic outputs."""

from __future__ import annotations

from types import SimpleNamespace

from chem_coworker.assistance import (
    project_multistep_response,
    project_retrosynthesis_response,
)


class _Value:
    def __init__(self, payload):
        self.payload = payload

    def to_dict(self):
        return self.payload


def test_one_step_projection_keeps_graph_validation_and_conditions() -> None:
    candidate = SimpleNamespace(forward_validation_status="verified_signature")
    strategy = SimpleNamespace(
        strategy_id="strategy:internal",
        strategy_rank=1,
        representative_score=0.9,
        independent_reference_support=2,
        total_realization_count=1,
        representative=candidate,
        to_dict=lambda: {
            "schema_version": "1.0",
            "strategy_id": "strategy:internal",
            "operator_id": "operator:1",
            "representative": {
                "precursor_smiles": "CCBr.N",
                "forward_validation_status": "verified_signature",
                "selectivity_warnings": [],
            },
        },
    )
    response = SimpleNamespace(
        request=SimpleNamespace(target_smiles="CCN"),
        valid=True,
        error=None,
        warnings=(),
        library_path="operators.json.gz",
        strategies=(strategy,),
        condition_evidence=(
            SimpleNamespace(
                strategy_id="strategy:internal",
                evidence=_Value({"status": "supported", "recommendations": []}),
            ),
        ),
        to_dict=lambda: {"valid": True, "strategies": [strategy.to_dict()]},
    )

    projection = project_retrosynthesis_response(response)
    packet = projection.inspection_packet("strategy-1")
    payloads = {item["evidence_id"]: item["payload"] for item in packet["evidence"]}

    assert dict(projection.strategy_aliases) == {
        "strategy-1": "strategy:internal"
    }
    assert payloads["strategy-1.strategy"]["representative"][
        "forward_validation_status"
    ] == "verified_signature"
    assert payloads["strategy-1.conditions"]["condition_evidence"]["status"] == (
        "supported"
    )


def test_multistep_projection_never_presents_partial_route_as_solved() -> None:
    candidate = SimpleNamespace(
        precursor_compatibility_disposition="pass",
        precursor_compatibility_warning_strength=None,
        precursor_compatibility_assessments=(),
        precursor_compatibility_policy_definition_id="policy.v1",
        selectivity_warnings=(),
        strategic_class="unresolved",
        realization_id="realization:1",
        strategy_id="strategy:1",
    )
    condition = SimpleNamespace(
        status="insufficient_evidence",
        warnings=("no precedent",),
    )
    step = SimpleNamespace(
        step_id="step:1",
        candidate=candidate,
        condition_evidence=condition,
        to_dict=lambda: {
            "step_id": "step:1",
            "candidate": {"forward_validation_status": "verified_signature"},
            "condition_evidence": {"status": "insufficient_evidence"},
        },
    )
    leaf = SimpleNamespace(
        terminal=False,
        route_node_id="leaf:1",
        unresolved_reason="maximum_depth",
    )
    route = SimpleNamespace(
        route_id="route:internal",
        solved=False,
        route_cost=2.5,
        reaction_count=1,
        maximum_depth=2,
        evidence_summary=_Value({"condition_insufficient_step_count": 1}),
        warnings=("unresolved leaf",),
        steps=(step,),
        leaves=(leaf,),
        to_dict=lambda: {
            "route_id": "route:internal",
            "solved": False,
            "steps": [step.to_dict()],
            "leaves": [{"terminal": False, "unresolved_reason": "maximum_depth"}],
        },
    )
    result = SimpleNamespace(
        routes=(),
        partial_routes=(route,),
        diagnostics=_Value({"stopped_by_expansion_limit": True}),
        schema_version="2.0",
    )
    response = SimpleNamespace(
        request=SimpleNamespace(target_smiles="CCN"),
        valid=True,
        error=None,
        warnings=("No fully solved route.",),
        result=result,
        to_dict=lambda: {
            "valid": True,
            "result": {"routes": [], "partial_routes": [route.to_dict()]},
        },
    )

    projection = project_multistep_response(response)
    initial = projection.initial_packet()
    evidence = {item["evidence_id"]: item for item in initial["evidence"]}
    step_packet = projection.inspection_packet("route-1", step_index=1)

    assert evidence["query.multistep"]["payload"]["route_kind"] == "partial"
    assert evidence["route-1.summary"]["payload"]["solved"] is False
    assert evidence["route-1.summary"]["uncertainty"] == "This is a partial route."
    assert step_packet["evidence"][0]["payload"]["route_solved"] is False
    assert any(
        item["payload_type"] == "route_repair_proposal"
        for item in step_packet["evidence"]
    )
