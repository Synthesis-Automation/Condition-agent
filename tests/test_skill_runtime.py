from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import MagicMock

from chem_coworker.agent import ChemCoworker
from chem_coworker.skills import build_default_skill_registry
from chem_coworker.tools import REGISTRY
from chem_coworker.workflow import WORKFLOW_REGISTRY, WorkflowDefinition


def _make_runtime_agent() -> ChemCoworker:
    agent = object.__new__(ChemCoworker)
    agent.registry = REGISTRY
    agent.skill_registry = build_default_skill_registry()
    agent.event_bus = MagicMock()
    agent.executor = MagicMock()
    agent.verbose = False
    agent.provider = "openai"
    agent.model_name = "fake-model"
    return agent


def test_tool_registry_policy_resolves_expected_public_tools() -> None:
    names = REGISTRY.filtered_names_for_policy("general_chemistry", llm_exposed_only=True)
    assert "analyze_reaction" in names
    assert "retrosynthesis_step" in names
    assert "recommend_reaction_conditions" in names


def test_conditions_specialist_policy_exposes_atomic_condition_tools() -> None:
    names = REGISTRY.filtered_names_for_policy("conditions_specialist", llm_exposed_only=False)
    assert "build_reaction_context" in names
    assert "get_literature_condition_evidence" in names
    assert "get_motif_condition_evidence" in names
    assert "get_similarity_condition_evidence" in names
    assert "get_rule_condition_evidence" in names
    assert "compose_condition_candidates" in names
    assert "score_condition_candidates" in names


def test_workflow_defaults_expose_policy_and_default_skills() -> None:
    retro = WORKFLOW_REGISTRY.get_for_task("retrosynthesis")
    assert retro.tool_policy == "retro_specialist"
    assert "retrosynthesis_route_planning" in (retro.default_skill_ids or [])

    fwd = WORKFLOW_REGISTRY.get_for_task("forward_synthesis")
    assert fwd.tool_policy == "forward_specialist"
    assert "forward_prediction" in (fwd.default_skill_ids or [])
    assert "condition_recommendation" not in (fwd.default_skill_ids or [])


def test_agent_resolves_active_skills_from_defaults_and_query() -> None:
    agent = _make_runtime_agent()
    workflow = WORKFLOW_REGISTRY.get_for_task("forward_chemistry")

    active = agent._resolve_active_skill_records(
        query="Recommend conditions for A>>B",
        task_type="predict",
        workflow=workflow,
        smiles_present=True,
    )
    active_ids = [record.manifest.id for record in active]

    assert "reaction_analysis" in active_ids
    assert "condition_recommendation" in active_ids


def test_build_skill_system_messages_only_include_active_instructions() -> None:
    agent = _make_runtime_agent()
    workflow = WORKFLOW_REGISTRY.get_for_task("forward_chemistry")
    active = agent._resolve_active_skill_records(
        query="Recommend conditions for A>>B",
        task_type="predict",
        workflow=workflow,
        smiles_present=True,
    )

    messages = agent._build_skill_system_messages(workflow, active)
    contents = "\n".join(str(getattr(msg, "content", "")) for msg in messages)

    assert "Active skill instructions:" in contents
    assert "Condition Recommendation" in contents


def test_active_skill_allowlist_can_expand_workflow_tool_surface() -> None:
    agent = _make_runtime_agent()
    workflow = WORKFLOW_REGISTRY.get_for_task("forward_chemistry")
    fake_skill = SimpleNamespace(
        manifest=SimpleNamespace(
            id="retro_expander",
            tool_allowlist=["apply_hte_templates"],
        )
    )

    names = agent._resolve_native_tool_names(workflow=workflow, active_skill_records=[fake_skill])
    assert "apply_hte_templates" in names


def test_active_skill_tool_policy_can_expand_workflow_tool_surface() -> None:
    agent = _make_runtime_agent()
    workflow = WORKFLOW_REGISTRY.get_for_task("forward_chemistry")
    fake_skill = SimpleNamespace(
        manifest=SimpleNamespace(
            id="literature_expander",
            tool_policy="literature_curator",
            tool_allowlist=[],
        )
    )

    names = agent._resolve_native_tool_names(workflow=workflow, active_skill_records=[fake_skill])
    assert "search_notes" in names
    assert "read_notes" in names


def test_agent_only_adds_one_query_matched_skill() -> None:
    agent = _make_runtime_agent()
    active = agent._resolve_active_skill_records(
        query="Recommend conditions and prepare a literature summary for A>>B",
        task_type="predict",
        workflow=WORKFLOW_REGISTRY.get_for_task("forward_chemistry"),
        smiles_present=True,
    )
    active_ids = [record.manifest.id for record in active]

    assert active_ids == ["reaction_analysis", "condition_recommendation"]


def test_repetition_guard_stops_identical_tool_loop_without_progress() -> None:
    agent = _make_runtime_agent()
    agent._build_native_tools = lambda workflow=None, active_skill_records=None: [object()]
    agent.executor._run_parallel.return_value = ({}, {})
    agent._accumulate_token_usage = lambda *args, **kwargs: None

    repeated_response = MagicMock()
    repeated_response.tool_calls = [{"name": "nonexistent_tool", "args": {"x": 1}, "id": "tc1"}]
    repeated_response.content = ""

    closing_response = MagicMock()
    closing_response.tool_calls = []
    closing_response.content = "Final answer after loop bailout."

    bound_llm = MagicMock()
    bound_llm.invoke.side_effect = [
        repeated_response,
        repeated_response,
        repeated_response,
        closing_response,
    ]
    agent.llm = MagicMock()
    agent.llm.bind_tools.return_value = bound_llm

    workflow = WorkflowDefinition(
        name="test",
        system_prompt="system",
        classifier_predicate=lambda t: True,
        max_iterations=5,
        default_skill_ids=[],
        tool_policy="general_chemistry",
        llm_visible_tools=["analyze_reaction"],
    )

    result = agent._run_native_tool_loop(
        query="Loop test",
        task_type="general",
        smiles_list=[],
        workflow=workflow,
        primary_smiles="",
    )

    warnings = result[3]
    final_answer = result[5]
    llm_calls = result[4]

    assert any("repeated tool loop" in warning.lower() for warning in warnings)
    assert final_answer == "Final answer after loop bailout."
    assert llm_calls == 4


def test_native_tool_loop_binds_tools_once_when_skills_are_static() -> None:
    agent = _make_runtime_agent()
    agent._accumulate_token_usage = lambda *args, **kwargs: None

    calls_seen: list[list[str]] = []

    def _fake_build_native_tools(workflow=None, active_skill_records=None):
        names = agent._resolve_native_tool_names(workflow=workflow, active_skill_records=active_skill_records)
        calls_seen.append(names)
        return [SimpleNamespace(name=name) for name in names]

    agent._build_native_tools = _fake_build_native_tools

    first_response = MagicMock()
    first_response.tool_calls = [{"name": "analyze_reaction", "args": {"reaction_smiles": "A>>B"}, "id": "tc1"}]
    first_response.content = "Need literature summary and notes before concluding."

    second_response = MagicMock()
    second_response.tool_calls = []
    second_response.content = "done"

    bound_initial = MagicMock()
    bound_initial.invoke.side_effect = [first_response]
    agent.llm = MagicMock()
    agent.llm.bind_tools.return_value = bound_initial
    agent.executor._run_parallel.return_value = (
        {"analyze_reaction": {"success": True, "reaction_type": "suzuki_miyaura"}},
        {"tc1": {"success": True, "reaction_type": "suzuki_miyaura"}},
    )

    workflow = WorkflowDefinition(
        name="forward_chemistry",
        system_prompt="system",
        classifier_predicate=lambda t: True,
        max_iterations=3,
        default_skill_ids=["reaction_analysis"],
        tool_policy="general_chemistry",
        llm_visible_tools=["analyze_reaction"],
    )

    result = agent._run_native_tool_loop(
        query="Analyze route",
        task_type="general",
        smiles_list=[],
        workflow=workflow,
        primary_smiles="",
        chemistry_state=agent._new_chemistry_run_state(),
    )

    messages = result[6]
    system_contents = [
        str(getattr(msg, "content", ""))
        for msg in messages
        if msg.__class__.__name__ == "SystemMessage"
    ]

    assert agent.llm.bind_tools.call_count == 1
    assert all("search_notes" not in names for names in calls_seen)
    assert not any("Additional skills have become relevant." in content for content in system_contents)
