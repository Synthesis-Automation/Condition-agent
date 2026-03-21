"""
Phase 5 cleanup verification tests.

Each test asserts that a classic-path artifact has been removed or that
the new native-only API is in place. Tests are purely structural (import/attribute
checks) so they run instantly without a live LLM.
"""
import inspect
import pytest


# ---------------------------------------------------------------------------
# Classic prompt constants must be gone
# ---------------------------------------------------------------------------

class TestClassicPromptsRemoved:
    def test_reason_prompt_gone(self):
        import chem_coworker.prompts as m
        assert not hasattr(m, "REASON_PROMPT"), "REASON_PROMPT should have been deleted in Phase 5"

    def test_synthesize_prompt_gone(self):
        import chem_coworker.prompts as m
        assert not hasattr(m, "SYNTHESIZE_PROMPT"), "SYNTHESIZE_PROMPT should have been deleted in Phase 5"

    def test_observe_prompt_gone(self):
        import chem_coworker.prompts as m
        assert not hasattr(m, "OBSERVE_PROMPT"), "OBSERVE_PROMPT should have been deleted in Phase 5"

    def test_retro_reason_prompt_gone(self):
        import chem_coworker.retro_prompts as m
        assert not hasattr(m, "RETRO_REASON_PROMPT"), "RETRO_REASON_PROMPT should have been deleted in Phase 5"

    def test_retro_synthesize_prompt_gone(self):
        import chem_coworker.retro_prompts as m
        assert not hasattr(m, "RETRO_SYNTHESIZE_PROMPT"), "RETRO_SYNTHESIZE_PROMPT should have been deleted in Phase 5"

    def test_native_system_prompt_still_exists(self):
        """NATIVE_SYSTEM_PROMPT must survive — it drives the only remaining execution path."""
        from chem_coworker.prompts import NATIVE_SYSTEM_PROMPT
        assert isinstance(NATIVE_SYSTEM_PROMPT, str) and len(NATIVE_SYSTEM_PROMPT) > 50

    def test_native_retro_system_prompt_still_exists(self):
        from chem_coworker.retro_prompts import NATIVE_RETRO_SYSTEM_PROMPT
        assert isinstance(NATIVE_RETRO_SYSTEM_PROMPT, str) and len(NATIVE_RETRO_SYSTEM_PROMPT) > 50


# ---------------------------------------------------------------------------
# PlanParser must be gone from plan.py
# ---------------------------------------------------------------------------

class TestPlanParserRemoved:
    def test_plan_parser_class_gone(self):
        import chem_coworker.plan as m
        assert not hasattr(m, "PlanParser"), "PlanParser should have been deleted in Phase 5"

    def test_smiles_placeholder_re_gone(self):
        import chem_coworker.plan as m
        assert not hasattr(m, "_SMILES_PLACEHOLDER_RE"), "_SMILES_PLACEHOLDER_RE should have been deleted in Phase 5"

    def test_plan_data_classes_still_present(self):
        """Core plan data types must survive — ExecutionPlan is used by native path."""
        from chem_coworker.plan import ExecutionPlan, ToolCall, PlanRejected
        plan = ExecutionPlan(hypothesis="test", confidence=0.9, groups=[], rationale="r", raw_plan_text="")
        assert plan.hypothesis == "test"


# ---------------------------------------------------------------------------
# ToolExecutor must not have run_plan()
# ---------------------------------------------------------------------------

class TestExecutorRunPlanRemoved:
    def test_run_plan_gone(self):
        from chem_coworker.executor import ToolExecutor
        assert not hasattr(ToolExecutor, "run_plan"), "run_plan() should have been deleted from ToolExecutor in Phase 5"

    def test_run_parallel_still_present(self):
        from chem_coworker.executor import ToolExecutor
        assert hasattr(ToolExecutor, "_run_parallel"), "_run_parallel must survive for native tool loop"


# ---------------------------------------------------------------------------
# ChemCoworker API changes
# ---------------------------------------------------------------------------

class TestChemCoworkerNativeOnly:
    def test_no_native_tools_param(self):
        """`native_tools` init parameter was removed when classic path was deleted."""
        sig = inspect.signature(__import__("chem_coworker").ChemCoworker.__init__)
        assert "native_tools" not in sig.parameters, "native_tools param should be gone"

    def test_no_check_hypothesis(self):
        from chem_coworker.agent import ChemCoworker
        assert not hasattr(ChemCoworker, "_check_hypothesis"), "_check_hypothesis should be gone"

    def test_collect_caveats_present(self):
        from chem_coworker.agent import ChemCoworker
        assert hasattr(ChemCoworker, "_collect_caveats"), "_collect_caveats must replace _check_hypothesis"

    def test_no_run_observe_step(self):
        from chem_coworker.agent import ChemCoworker
        assert not hasattr(ChemCoworker, "_run_observe_step"), "_run_observe_step should be gone"

    def test_no_observe_threshold_constant(self):
        import chem_coworker.agent as m
        assert not hasattr(m, "_OBSERVE_THRESHOLD"), "_OBSERVE_THRESHOLD should be gone"

    def test_run_method_exists(self):
        from chem_coworker.agent import ChemCoworker
        assert callable(ChemCoworker.run)

    def test_native_tool_loop_returns_8_tuple(self):
        """Early-exit guard in _run_native_tool_loop must return 8-tuple including token usage."""
        from chem_coworker.agent import ChemCoworker
        src = inspect.getsource(ChemCoworker._run_native_tool_loop)
        assert '[], self._new_token_section("reason", "Reasoning")' in src, (
            "Early-exit must return messages plus a token-usage section."
        )


# ---------------------------------------------------------------------------
# ChemResponse field cleanup
# ---------------------------------------------------------------------------

class TestChemResponseFieldsRemoved:
    def test_plan_text_field_gone(self):
        from chem_coworker.response import ChemResponse
        assert not hasattr(ChemResponse, "plan_text") or \
            "plan_text" not in ChemResponse.__dataclass_fields__, \
            "plan_text should have been removed from ChemResponse in Phase 5"

    def test_plan_revised_field_gone(self):
        from chem_coworker.response import ChemResponse
        assert "plan_revised" not in ChemResponse.__dataclass_fields__, \
            "plan_revised should have been removed from ChemResponse in Phase 5"

    def test_observe_text_field_gone(self):
        from chem_coworker.response import ChemResponse
        assert "observe_text" not in ChemResponse.__dataclass_fields__, \
            "observe_text should have been removed from ChemResponse in Phase 5"

    def test_core_fields_still_present(self):
        from chem_coworker.response import ChemResponse
        for field in ("query", "answer", "hypothesis", "tools_called", "warnings"):
            assert field in ChemResponse.__dataclass_fields__
