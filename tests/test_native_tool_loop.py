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

    def test_auto_inserts_missing_prerequisite_for_recommend_conditions(self):
        """If detect_reaction_type is missing and inferable, auto-insert it."""
        agent = _make_agent()
        callables = agent.registry.get_callables()
        response_tool_calls = [
            {"name": "recommend_conditions", "args": {"reaction_smiles": "A>>B"}, "id": "tc1"},
        ]

        executable, blocked, warnings = agent._partition_native_tool_calls_by_contracts(
            response_tool_calls=response_tool_calls,
            callables=callables,
            tool_results={},
        )

        assert [c.name for c in executable] == ["detect_reaction_type"]
        assert executable[0].args["reaction_smiles"] == "A>>B"
        assert "recommend_conditions" in blocked
        assert blocked["recommend_conditions"]["deferred"] is True
        assert any("Auto-inserted prerequisite(s): detect_reaction_type" in w for w in warnings)

    def test_dedupes_repeated_block_warnings_for_duplicate_calls(self):
        """Duplicate model tool calls should not spam identical contract warnings."""
        agent = _make_agent()
        callables = agent.registry.get_callables()
        response_tool_calls = [
            {"name": "recommend_conditions", "args": {"reaction_smiles": "A>>B"}, "id": "tc1"},
            {"name": "recommend_conditions", "args": {"reaction_smiles": "A>>B"}, "id": "tc2"},
        ]

        executable, blocked, warnings = agent._partition_native_tool_calls_by_contracts(
            response_tool_calls=response_tool_calls,
            callables=callables,
            tool_results={},
        )

        assert [c.name for c in executable] == ["detect_reaction_type"]
        assert "recommend_conditions" in blocked
        assert len(warnings) == 1


class TestDuplicateNativeToolCallsPreservePerCallResults:
    def test_duplicate_same_tool_calls_keep_results_by_call_id(self):
        agent = _make_agent()

        first_response = MagicMock()
        first_response.tool_calls = [
            {"name": "get_molecular_descriptors", "args": {"smiles": "CCO"}, "id": "call_1"},
            {"name": "get_molecular_descriptors", "args": {"smiles": "CCN"}, "id": "call_2"},
        ]
        first_response.content = ""

        final_response = MagicMock()
        final_response.tool_calls = []
        final_response.content = "done"

        invoke_state = {"n": 0}

        def _invoke(_messages):
            invoke_state["n"] += 1
            return first_response if invoke_state["n"] == 1 else final_response

        mock_bound = MagicMock()
        mock_bound.invoke = _invoke
        agent.llm = MagicMock()
        agent.llm.bind_tools.return_value = mock_bound

        agent.executor._run_parallel.return_value = (
            {"get_molecular_descriptors": {"success": True, "smiles": "CCN", "descriptors": {"MW": 45.0}}},
            {
                "call_1": {"success": True, "smiles": "CCO", "descriptors": {"MW": 46.0}},
                "call_2": {"success": True, "smiles": "CCN", "descriptors": {"MW": 45.0}},
            },
        )

        from chem_coworker.workflow import WORKFLOW_REGISTRY
        workflow = WORKFLOW_REGISTRY.get_for_task("forward_chemistry")

        tool_results, _, _, _, _, _, messages = agent._run_native_tool_loop(
            query="descriptor test",
            task_type="general",
            smiles_list=[],
            workflow=workflow,
            primary_smiles="",
            chemistry_state=agent._new_chemistry_run_state(),
        )

        assert "_tool_call_results_by_id" in tool_results
        by_id = tool_results["_tool_call_results_by_id"]
        assert by_id["call_1"]["smiles"] == "CCO"
        assert by_id["call_2"]["smiles"] == "CCN"

        tool_messages = [m for m in messages if m.__class__.__name__ == "ToolMessage"]
        assert len(tool_messages) >= 2
        payloads = {m.tool_call_id: m.content for m in tool_messages if getattr(m, "tool_call_id", "") in {"call_1", "call_2"}}
        assert "\"CCO\"" in payloads["call_1"]
        assert "\"CCN\"" in payloads["call_2"]


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


class TestSharedReactionContextCaching:
    def test_context_aware_wrappers_share_one_reaction_context(self):
        agent = _make_agent()
        state = agent._new_chemistry_run_state()

        agent._detect_reaction_type_from_context = lambda ctx: {"success": True, "ctx_id": id(ctx)}  # type: ignore[method-assign]
        agent._analyze_bond_changes_from_context = lambda ctx: {"success": True, "ctx_id": id(ctx)}  # type: ignore[method-assign]

        wrapped = agent._build_context_aware_callables(
            {
                "detect_reaction_type": lambda reaction_smiles: {"success": True},  # noqa: ARG005
                "analyze_bond_changes": lambda reaction_smiles: {"success": True},  # noqa: ARG005
            },
            state,
        )

        r1 = wrapped["detect_reaction_type"](reaction_smiles="A >> B")
        r2 = wrapped["analyze_bond_changes"](reaction_smiles="A>>B")

        assert len(state.reaction_contexts) == 1
        assert r1["ctx_id"] == r2["ctx_id"]

    def test_context_aware_wrappers_cache_conditions_by_top_k(self):
        agent = _make_agent()
        state = agent._new_chemistry_run_state()
        calls = {"n": 0}

        def _fake_recommend_conditions(reaction_smiles: str, top_k: int = 5):
            calls["n"] += 1
            return {"success": True, "reaction_smiles": reaction_smiles, "top_k": top_k}

        wrapped = agent._build_context_aware_callables(
            {"recommend_conditions": _fake_recommend_conditions},
            state,
        )

        a = wrapped["recommend_conditions"](reaction_smiles="A>>B", top_k=5)
        b = wrapped["recommend_conditions"](reaction_smiles="A >> B", top_k=5)
        c = wrapped["recommend_conditions"](reaction_smiles="A>>B", top_k=3)

        assert calls["n"] == 2
        assert a["top_k"] == 5 and b["top_k"] == 5 and c["top_k"] == 3

    def test_context_aware_wrappers_cache_molecule_tools(self):
        agent = _make_agent()
        state = agent._new_chemistry_run_state()
        counts = {"fg": 0, "desc": 0}

        def _fake_fg(smiles: str):
            counts["fg"] += 1
            return {"success": True, "smiles": smiles}

        def _fake_desc(smiles: str):
            counts["desc"] += 1
            return {"success": True, "smiles": smiles}

        wrapped = agent._build_context_aware_callables(
            {
                "inspect_functional_groups": _fake_fg,
                "get_molecular_descriptors": _fake_desc,
            },
            state,
        )

        wrapped["inspect_functional_groups"](smiles="c1ccccc1")
        wrapped["inspect_functional_groups"](smiles="c1ccccc1")
        wrapped["get_molecular_descriptors"](smiles="CCO")
        wrapped["get_molecular_descriptors"](smiles="CCO")

        assert counts == {"fg": 1, "desc": 1}
