"""Phase 4: Message thread cleanup tests — exhaustion guard and 7-tuple return."""

import pytest
from unittest.mock import MagicMock, patch
from chem_coworker.agent import ChemCoworker
from chem_coworker.tools import REGISTRY


def _make_agent():
    agent = object.__new__(ChemCoworker)
    agent.registry = REGISTRY
    agent.event_bus = MagicMock()
    agent.executor = MagicMock()
    agent.verbose = False
    return agent


class TestNativeToolLoopReturns7Tuple:
    def test_successful_loop_returns_7_elements(self):
        """Verify the return tuple has exactly 7 elements (adds messages)."""
        agent = _make_agent()

        # Build a minimal mock LLM that writes a final answer immediately
        mock_response = MagicMock()
        mock_response.tool_calls = []  # no tool calls → writes final answer
        mock_response.content = "The answer is 42."
        mock_llm = MagicMock()
        mock_llm.bind_tools.return_value.invoke.return_value = mock_response
        agent.llm = mock_llm

        from chem_coworker.workflow import WORKFLOW_REGISTRY
        workflow = WORKFLOW_REGISTRY.get_for_task("forward_chemistry")

        result = agent._run_native_tool_loop(
            query="test query",
            task_type="predict",
            smiles_list=[],
            workflow=workflow,
            primary_smiles="",
        )
        assert len(result) == 7, f"Expected 7-tuple, got {len(result)}-tuple"
        tool_results, hypothesis, confidence, warnings, llm_count, final_answer, messages = result
        assert final_answer == "The answer is 42."
        assert isinstance(messages, list)
        assert len(messages) >= 2  # at least SystemMessage + HumanMessage + response


class TestExhaustionGuard:
    def test_exhaustion_guard_fires_when_no_final_answer(self):
        """When max_iterations exhausted with no answer, closing call is made."""
        agent = _make_agent()

        # LLM that always returns tool calls → never writes an answer → loop exhausted
        always_tool_response = MagicMock()
        always_tool_response.tool_calls = [{"name": "nonexistent_tool", "args": {}, "id": "tc_1"}]
        always_tool_response.content = ""

        closing_response = MagicMock()
        closing_response.tool_calls = []
        closing_response.content = "Closing expert answer."

        call_count = {"n": 0}

        def mock_invoke(messages):
            call_count["n"] += 1
            if call_count["n"] <= 2:  # first 2 calls → always tool calls (exhaust loop)
                return always_tool_response
            return closing_response   # closing call

        from chem_coworker.workflow import WorkflowDefinition
        mini_workflow = WorkflowDefinition(
            name="test",
            system_prompt="test system",
            classifier_predicate=lambda t: True,
            max_iterations=2,   # exhaust after 2 iterations
        )

        mock_llm_bound = MagicMock()
        mock_llm_bound.invoke = mock_invoke
        mock_llm = MagicMock()
        mock_llm.bind_tools.return_value = mock_llm_bound
        agent.llm = mock_llm

        # Executor returns empty result for unknown tool
        agent.executor._run_parallel.return_value = {}

        _, _, _, warnings_out, llm_count, final_answer, messages = agent._run_native_tool_loop(
            query="test query",
            task_type="test",
            smiles_list=[],
            workflow=mini_workflow,
            primary_smiles="",
        )

        assert final_answer == "Closing expert answer."
        assert llm_count == 3   # 2 loop iterations + 1 closing call

    def test_no_exhaustion_guard_when_answer_present(self):
        """When the loop produces a final answer, the closing call is NOT made."""
        agent = _make_agent()

        normal_response = MagicMock()
        normal_response.tool_calls = []
        normal_response.content = "Normal answer."

        invoke_calls = {"n": 0}

        def mock_invoke(messages):
            invoke_calls["n"] += 1
            return normal_response

        from chem_coworker.workflow import WorkflowDefinition
        mini_workflow = WorkflowDefinition(
            name="test",
            system_prompt="test system",
            classifier_predicate=lambda t: True,
            max_iterations=5,
        )

        mock_llm_bound = MagicMock()
        mock_llm_bound.invoke = mock_invoke
        mock_llm = MagicMock()
        mock_llm.bind_tools.return_value = mock_llm_bound
        agent.llm = mock_llm

        _, _, _, _, llm_count, final_answer, _ = agent._run_native_tool_loop(
            query="test query",
            task_type="test",
            smiles_list=[],
            workflow=mini_workflow,
            primary_smiles="",
        )

        assert final_answer == "Normal answer."
        assert invoke_calls["n"] == 1   # only the one loop iteration call, no closing call


class TestNativeToolContractEnforcement:
    def test_partitions_tool_calls_and_defers_requires_until_provider_runs(self):
        """recommend_conditions should be deferred until detect_reaction_type has run."""
        agent = _make_agent()
        from chem_coworker.plan import ToolCall

        callables = agent.registry.get_callables()
        response_tool_calls = [
            {"name": "recommend_conditions", "args": {"reaction_smiles": "A>>B"}, "id": "tc1"},
            {"name": "detect_reaction_type", "args": {"reaction_smiles": "A>>B"}, "id": "tc2"},
        ]

        executable, blocked, warnings = agent._partition_native_tool_calls_by_contracts(
            response_tool_calls=response_tool_calls,
            callables=callables,
            tool_results={},
        )

        assert [c.name for c in executable] == ["detect_reaction_type"]
        assert "recommend_conditions" in blocked
        assert blocked["recommend_conditions"]["contract_violation"] is True
        assert blocked["recommend_conditions"]["deferred"] is True
        assert "reaction_type" in blocked["recommend_conditions"]["missing_requires"]
        assert any("recommend_conditions" in w for w in warnings)

    def test_allows_tool_when_required_context_key_exists(self):
        """recommend_conditions may run once reaction_type context is available."""
        agent = _make_agent()
        callables = agent.registry.get_callables()
        response_tool_calls = [
            {"name": "recommend_conditions", "args": {"reaction_smiles": "A>>B"}, "id": "tc1"},
        ]
        prior_results = {
            "detect_reaction_type": {
                "success": True,
                "reaction_type": "suzuki_miyaura",
                "reaction_type_id": "suzuki_miyaura",
            }
        }

        executable, blocked, warnings = agent._partition_native_tool_calls_by_contracts(
            response_tool_calls=response_tool_calls,
            callables=callables,
            tool_results=prior_results,
        )

        assert [c.name for c in executable] == ["recommend_conditions"]
        assert blocked == {}
        assert warnings == []


class TestStructuredExtractionCanonicalReactionType:
    def test_extract_structured_prefers_reaction_type_id_and_keeps_metadata(self):
        agent = _make_agent()
        structured = agent._extract_structured(
            {
                "detect_reaction_type": {
                    "success": True,
                    "reaction_type": "Suzuki",
                    "reaction_type_id": "suzuki_miyaura",
                    "family_label": "Suzuki-Miyaura Cross-Coupling",
                    "reaction_type_metadata": {
                        "id": "suzuki_miyaura",
                        "name": "Suzuki-Miyaura",
                        "category": "cross_coupling",
                        "aliases": ["Suzuki"],
                        "has_constraints": True,
                    },
                }
            }
        )
        assert structured["reaction_type"] == "suzuki_miyaura"
        assert structured["reaction_family"] == "Suzuki-Miyaura Cross-Coupling"
        assert structured["reaction_type_metadata"]["category"] == "cross_coupling"

    def test_extract_structured_uses_tool_declared_projection(self):
        agent = _make_agent()
        structured = agent._extract_structured(
            {
                "search_reaction_types": {
                    "success": True,
                    "matches": [{"id": "suzuki_miyaura", "score": 0.9}],
                }
            }
        )
        assert structured["taxonomy_matches"][0]["id"] == "suzuki_miyaura"

    def test_extract_structured_supports_custom_plugin_projection(self):
        from chem_coworker.tools._base import ToolPlugin

        agent = object.__new__(ChemCoworker)
        fake_plugin = ToolPlugin(
            name="custom_tool",
            category="test",
            description="d",
            fn=lambda: {},
            structured_projection=lambda result: {"custom_metric": result.get("value")},
        )
        agent.registry = type("R", (), {"_plugins": {"custom_tool": fake_plugin}})()
        agent.verbose = False

        structured = agent._extract_structured({"custom_tool": {"success": True, "value": 7}})
        assert structured == {"custom_metric": 7}


class TestAggregateConfidence:
    def test_aggregate_confidence_increases_with_strong_consistent_evidence(self):
        agent = _make_agent()
        score = agent._aggregate_confidence(
            tool_results={
                "detect_reaction_type": {
                    "success": True,
                    "reaction_type": "suzuki_miyaura",
                    "reaction_type_id": "suzuki_miyaura",
                    "confidence": 0.94,
                    "reaction_type_metadata": {"category": "cross_coupling"},
                    "reacted_motifs": ["Ar-Br"],
                    "formed_motifs": ["Biaryl"],
                },
                "analyze_bond_changes": {
                    "success": True,
                    "bonds_formed": [[1, 2]],
                    "key_bond_type": "C-C (Suzuki-type)",
                    "mapping_confidence": 0.88,
                },
                "recommend_conditions": {
                    "success": True,
                    "recommendations": [{"num_experiments": 12, "confidence": 0.83}],
                },
            },
            base_confidence=0.5,
            warnings=[],
            critic_findings=[],
        )
        assert score > 0.75

    def test_aggregate_confidence_decreases_with_unknown_and_warnings(self):
        agent = _make_agent()
        critic_finding = type("F", (), {"severity": type("S", (), {"value": "warning"})()})()
        score = agent._aggregate_confidence(
            tool_results={
                "detect_reaction_type": {
                    "success": True,
                    "reaction_type": "Unknown",
                    "confidence": 0.15,
                    "reaction_type_metadata": {},
                },
                "analyze_bond_changes": {
                    "success": True,
                    "bonds_formed": [],
                    "key_bond_type": "unknown",
                    "mapping_confidence": 0.18,
                },
                "recommend_conditions": {
                    "success": True,
                    "recommendations": [],
                },
                "recommend_conditions_blocked": {
                    "success": False,
                    "contract_violation": True,
                },
            },
            base_confidence=0.5,
            warnings=["Blocked 'recommend_conditions' due to unmet tool contracts.", "Tool 'x' failed: y"],
            critic_findings=[critic_finding],
        )
        assert score < 0.45
