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
