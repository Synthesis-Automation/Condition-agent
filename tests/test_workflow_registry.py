"""Phase 3: WorkflowRegistry unit tests."""

import pytest
from chem_coworker.workflow import WorkflowDefinition, WorkflowRegistry, WORKFLOW_REGISTRY


def _make_workflow(name: str, matches: str = None, is_fallback: bool = False) -> WorkflowDefinition:
    """Helper: build a minimal WorkflowDefinition."""
    predicate = (lambda t: t == matches) if matches else (lambda t: True)
    return WorkflowDefinition(name=name, system_prompt=f"System: {name}", classifier_predicate=predicate)


class TestWorkflowDefinitionDefaults:
    def test_max_iterations_default_is_8(self):
        w = _make_workflow("test")
        assert w.max_iterations == 8

    def test_critic_step_default_is_none(self):
        w = _make_workflow("test")
        assert w.critic_step is None

    def test_fields_stored(self):
        w = WorkflowDefinition(
            name="myflow",
            system_prompt="sys",
            classifier_predicate=lambda t: True,
            max_iterations=12,
        )
        assert w.name == "myflow"
        assert w.system_prompt == "sys"
        assert w.max_iterations == 12


class TestWorkflowRegistryRouting:
    def test_retrosynthesis_routes_to_retro_workflow(self):
        reg = WorkflowRegistry()
        reg.register(_make_workflow("retro", matches="retrosynthesis"))
        reg.register(_make_workflow("general"), is_fallback=True)
        assert reg.get_for_task("retrosynthesis").name == "retro"

    def test_non_retro_routes_to_fallback(self):
        reg = WorkflowRegistry()
        reg.register(_make_workflow("retro", matches="retrosynthesis"))
        reg.register(_make_workflow("general"), is_fallback=True)
        assert reg.get_for_task("forward_synthesis").name == "general"

    def test_unknown_task_routes_to_fallback(self):
        reg = WorkflowRegistry()
        reg.register(_make_workflow("retro", matches="retrosynthesis"))
        reg.register(_make_workflow("general"), is_fallback=True)
        assert reg.get_for_task("other_task").name == "general"

    def test_missing_fallback_raises_runtime_error(self):
        reg = WorkflowRegistry()
        reg.register(_make_workflow("retro", matches="retrosynthesis"))
        # no fallback registered
        with pytest.raises(RuntimeError, match="no fallback"):
            reg.get_for_task("something_else")

    def test_first_matching_predicate_wins(self):
        reg = WorkflowRegistry()
        reg.register(_make_workflow("first", matches="x"))
        reg.register(_make_workflow("second", matches="x"))
        reg.register(_make_workflow("fallback"), is_fallback=True)
        assert reg.get_for_task("x").name == "first"


class TestWorkflowRegistryNames:
    def test_names_includes_all_workflows(self):
        reg = WorkflowRegistry()
        reg.register(_make_workflow("retro", matches="retrosynthesis"))
        reg.register(_make_workflow("general"), is_fallback=True)
        names = reg.names()
        assert "retro" in names
        assert "general" in names

    def test_names_order_specific_before_fallback(self):
        reg = WorkflowRegistry()
        reg.register(_make_workflow("retro", matches="retrosynthesis"))
        reg.register(_make_workflow("general"), is_fallback=True)
        names = reg.names()
        assert names.index("retro") < names.index("general")

    def test_names_empty_registry(self):
        reg = WorkflowRegistry()
        assert reg.names() == []


class TestGlobalWorkflowRegistry:
    def test_retrosynthesis_workflow_exists(self):
        w = WORKFLOW_REGISTRY.get_for_task("retrosynthesis")
        assert w.name == "retrosynthesis"

    def test_retrosynthesis_max_iterations_is_10(self):
        w = WORKFLOW_REGISTRY.get_for_task("retrosynthesis")
        assert w.max_iterations == 10

    def test_retrosynthesis_has_system_prompt(self):
        w = WORKFLOW_REGISTRY.get_for_task("retrosynthesis")
        assert len(w.system_prompt) > 50   # non-trivial prompt

    def test_forward_chemistry_workflow_is_fallback(self):
        w = WORKFLOW_REGISTRY.get_for_task("predict_yield")
        assert w.name == "forward_chemistry"

    def test_forward_chemistry_max_iterations_is_8(self):
        w = WORKFLOW_REGISTRY.get_for_task("anything")
        assert w.max_iterations == 8

    def test_retrosynthesis_has_critic_step(self):
        """Phase 6: retrosynthesis workflow must have an enabled CriticStep."""
        from chem_coworker.workflow import CriticStep
        retro = WORKFLOW_REGISTRY.get_for_task("retrosynthesis")
        assert retro.critic_step is not None
        assert isinstance(retro.critic_step, CriticStep)
        assert retro.critic_step.enabled is True

    def test_forward_chemistry_has_no_critic_step(self):
        """Critic is retrosynthesis-specific; forward chemistry workflow has no critic."""
        fwd = WORKFLOW_REGISTRY.get_for_task("forward_chemistry")
        assert fwd.critic_step is None

    def test_names_lists_both_workflows(self):
        names = WORKFLOW_REGISTRY.names()
        assert "retrosynthesis" in names
        assert "forward_chemistry" in names

    def test_retrosynthesis_workflow_exposes_specialist_analysis_tools(self):
        w = WORKFLOW_REGISTRY.get_for_task("retrosynthesis")
        visible = set(w.llm_visible_tools or [])
        assert "featurize_molecule" in visible
        assert "assess_snar_feasibility" in visible
        assert "recommend_reaction_conditions" in visible
        assert "evaluate_synthesis_proposal" in visible
        assert "plan_route_candidates" in visible

    def test_forward_synthesis_workflow_exposes_specialist_analysis_tools(self):
        w = WORKFLOW_REGISTRY.get_for_task("forward_synthesis")
        visible = set(w.llm_visible_tools or [])
        assert "featurize_molecule" in visible
        assert "assess_snar_feasibility" in visible
        assert "recommend_reaction_conditions" in visible
        assert "evaluate_synthesis_proposal" in visible

    def test_fallback_workflow_exposes_forward_and_retro_facades(self):
        w = WORKFLOW_REGISTRY.get_for_task("anything_else")
        visible = set(w.llm_visible_tools or [])
        assert "retrosynthesis_step" in visible
        assert "forward_synthesis_step" in visible
