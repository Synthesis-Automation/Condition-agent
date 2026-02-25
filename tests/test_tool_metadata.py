"""
Tests for Phase 1: ToolPlugin metadata fields (provides, requires, validators)
and data-contract SMILES resolution in ToolExecutor.
"""
from __future__ import annotations

import pytest
import sys
import types
from typing import Any, Dict, Optional
from unittest.mock import MagicMock


# ---------------------------------------------------------------------------
# ToolPlugin dataclass tests
# ---------------------------------------------------------------------------

class TestToolPluginMetadataDefaults:
    def test_provides_defaults_to_empty_list(self):
        from chem_coworker.tools._base import ToolPlugin
        p = ToolPlugin(name="x", category="c", description="d", fn=lambda: {})
        assert p.provides == []

    def test_requires_defaults_to_empty_list(self):
        from chem_coworker.tools._base import ToolPlugin
        p = ToolPlugin(name="x", category="c", description="d", fn=lambda: {})
        assert p.requires == []

    def test_validators_defaults_to_empty_list(self):
        from chem_coworker.tools._base import ToolPlugin
        p = ToolPlugin(name="x", category="c", description="d", fn=lambda: {})
        assert p.validators == []

    def test_structured_projection_defaults_to_none(self):
        from chem_coworker.tools._base import ToolPlugin
        p = ToolPlugin(name="x", category="c", description="d", fn=lambda: {})
        assert p.structured_projection is None

    def test_provides_stored_correctly(self):
        from chem_coworker.tools._base import ToolPlugin
        p = ToolPlugin(
            name="x", category="c", description="d", fn=lambda: {},
            provides=["resolved_smiles", "smiles"],
        )
        assert p.provides == ["resolved_smiles", "smiles"]

    def test_requires_stored_correctly(self):
        from chem_coworker.tools._base import ToolPlugin
        p = ToolPlugin(
            name="x", category="c", description="d", fn=lambda: {},
            requires=["reaction_type"],
        )
        assert p.requires == ["reaction_type"]

    def test_validators_stored_correctly(self):
        from chem_coworker.tools._base import ToolPlugin
        v = lambda r: None
        p = ToolPlugin(
            name="x", category="c", description="d", fn=lambda: {},
            validators=[v],
        )
        assert p.validators == [v]

    def test_structured_projection_stored_correctly(self):
        from chem_coworker.tools._base import ToolPlugin

        def proj(result):
            return {"x": result.get("y")}

        p = ToolPlugin(
            name="x", category="c", description="d", fn=lambda: {},
            structured_projection=proj,
        )
        assert p.structured_projection is proj


class TestToolPluginValidatorSemantics:
    def test_validator_receives_result_dict(self):
        from chem_coworker.tools._base import ToolPlugin
        seen = []
        def v(result):
            seen.append(result)
            return None
        p = ToolPlugin(name="x", category="c", description="d", fn=lambda: {}, validators=[v])
        result = {"success": True, "data": 42}
        for fn in p.validators:
            fn(result)
        assert seen == [result]

    def test_validator_returning_string_is_warning(self):
        from chem_coworker.tools._base import ToolPlugin
        def v(result):
            return "No precedents found"
        p = ToolPlugin(name="x", category="c", description="d", fn=lambda: {}, validators=[v])
        assert p.validators[0]({"success": True}) == "No precedents found"

    def test_validator_returning_none_means_pass(self):
        from chem_coworker.tools._base import ToolPlugin
        def v(result):
            return None
        p = ToolPlugin(name="x", category="c", description="d", fn=lambda: {}, validators=[v])
        assert p.validators[0]({"success": True}) is None


# ---------------------------------------------------------------------------
# conditions.py validator tests
# ---------------------------------------------------------------------------

class TestValidateRecommendConditions:
    def test_no_recs_returns_warning(self):
        from chem_coworker.tools.conditions import _validate_recommend_conditions
        result = {"success": True, "recommendations": []}
        msg = _validate_recommend_conditions(result)
        assert msg is not None
        assert "NO HTE" in msg

    def test_zero_exp_low_conf_returns_warning(self):
        from chem_coworker.tools.conditions import _validate_recommend_conditions
        result = {
            "success": True,
            "recommendations": [{"num_experiments": 0, "confidence": 0.1}],
        }
        msg = _validate_recommend_conditions(result)
        assert msg is not None
        assert "tentative" in msg

    def test_good_result_returns_none(self):
        from chem_coworker.tools.conditions import _validate_recommend_conditions
        result = {
            "success": True,
            "recommendations": [{"num_experiments": 15, "confidence": 0.85}],
        }
        assert _validate_recommend_conditions(result) is None

    def test_failed_tool_skipped(self):
        from chem_coworker.tools.conditions import _validate_recommend_conditions
        assert _validate_recommend_conditions({"success": False}) is None

    def test_zero_exp_high_conf_no_warning(self):
        """0 experiments but high confidence should NOT warn (threshold is conf < 0.3)."""
        from chem_coworker.tools.conditions import _validate_recommend_conditions
        result = {
            "success": True,
            "recommendations": [{"num_experiments": 0, "confidence": 0.7}],
        }
        assert _validate_recommend_conditions(result) is None

    def test_validator_registered_on_tool(self):
        """The validator must be registered on the ToolPlugin itself."""
        from chem_coworker.tools import REGISTRY
        plugin = REGISTRY._plugins.get("recommend_conditions")
        assert plugin is not None
        assert len(plugin.validators) > 0
        assert plugin.validators[0] is not None

    def test_recommend_conditions_surfaces_hte_timing(self, monkeypatch):
        from chem_coworker.tools.conditions import _recommend_conditions

        monkeypatch.setattr(
            "chemtools.recommend.hte_adapter.recommend_from_reaction",
            lambda *args, **kwargs: {
                "meta": {
                    "processing_time_ms": 1234.5,
                    "timing_ms": {
                        "input_parse_ms": 1.0,
                        "recommender_get_ms": 2.0,
                        "recommend_compute_ms": 1200.0,
                        "postprocess_ms": 31.5,
                        "total_ms": 1234.5,
                    },
                },
                "recommended_conditions": [],
                "recommendations": [],
            },
        )
        result = _recommend_conditions("CCO>>CC=O", top_k=3)
        assert result["success"] is True
        assert result["hte_processing_time_ms"] == 1234.5
        assert result["hte_timing_ms"]["recommend_compute_ms"] == 1200.0


# ---------------------------------------------------------------------------
# ToolPlugin provides/requires registration checks
# ---------------------------------------------------------------------------

class TestRegistryAnnotations:
    def test_detect_reaction_type_provides(self):
        from chem_coworker.tools import REGISTRY
        plugin = REGISTRY._plugins.get("detect_reaction_type")
        assert plugin is not None
        assert "reaction_type" in plugin.provides

    def test_recommend_conditions_provides(self):
        from chem_coworker.tools import REGISTRY
        plugin = REGISTRY._plugins.get("recommend_conditions")
        assert plugin is not None
        assert "recommendations" in plugin.provides

    def test_recommend_conditions_requires(self):
        from chem_coworker.tools import REGISTRY
        plugin = REGISTRY._plugins.get("recommend_conditions")
        assert plugin is not None
        assert "reaction_type" in plugin.requires

    def test_resolve_to_smiles_provides(self):
        from chem_coworker.tools import REGISTRY
        plugin = REGISTRY._plugins.get("resolve_to_smiles")
        assert plugin is not None
        assert "resolved_smiles" in plugin.provides or "smiles" in plugin.provides


class TestDetectReactionTypeCanonicalMetadata:
    def test_detect_reaction_type_returns_taxonomy_canonical_metadata(self, monkeypatch):
        from chem_coworker.tools.chemistry import _detect_reaction_type

        fake_unified = types.ModuleType("chemtools.featurizers.unified")
        fake_unified.featurize_reaction = lambda rxn: {
            "reaction_type": "Suzuki",
            "confidence": 0.91,
            "detection": {"validation": {}, "evidence": {}},
            "reaction_key": "rk1",
        }

        fake_catalog = types.ModuleType("chemtools.taxonomy.reaction_catalog")
        fake_catalog.resolve_reaction_type = lambda label: "suzuki_miyaura" if str(label).lower() == "suzuki" else None
        fake_catalog.get_reaction_type = lambda rid: types.SimpleNamespace(
            id="suzuki_miyaura",
            name="Suzuki-Miyaura",
            category="cross_coupling",
            aliases=["Suzuki"],
            constraints={"bond_formed": "C-C"},
        )

        monkeypatch.setitem(sys.modules, "chemtools.featurizers.unified", fake_unified)
        monkeypatch.setitem(sys.modules, "chemtools.taxonomy.reaction_catalog", fake_catalog)

        out = _detect_reaction_type("Brc1ccccc1.B(O)Oc1ccccc1>>c1ccc(-c2ccccc2)cc1")
        assert out["success"] is True
        assert out["reaction_type"] == "suzuki_miyaura"
        assert out["reaction_type_id"] == "suzuki_miyaura"
        assert out["reaction_type_metadata"]["id"] == "suzuki_miyaura"
        assert out["reaction_type_metadata"]["category"] == "cross_coupling"


class TestChemistryValidators:
    def test_validate_detect_reaction_type_warns_on_low_confidence(self):
        from chem_coworker.tools.chemistry import _validate_detect_reaction_type

        msg = _validate_detect_reaction_type(
            {
                "success": True,
                "reaction_type": "C_N_Coupling",
                "confidence": 0.21,
                "reaction_type_metadata": {"category": "cross_coupling"},
            }
        )
        assert msg is not None
        assert "Low reaction-type confidence" in msg

    def test_validate_detect_reaction_type_warns_on_unknown(self):
        from chem_coworker.tools.chemistry import _validate_detect_reaction_type

        msg = _validate_detect_reaction_type({"success": True, "reaction_type": "Unknown"})
        assert msg is not None
        assert "Unknown" in msg

    def test_validate_analyze_bond_changes_warns_on_missing_formed_bonds(self):
        from chem_coworker.tools.chemistry import _validate_analyze_bond_changes

        msg = _validate_analyze_bond_changes(
            {"success": True, "bonds_formed": [], "key_bond_type": "C-N", "mapping_confidence": 0.9}
        )
        assert msg is not None
        assert "no formed bonds" in msg.lower()

    def test_validate_analyze_bond_changes_warns_on_low_mapping_confidence(self):
        from chem_coworker.tools.chemistry import _validate_analyze_bond_changes

        msg = _validate_analyze_bond_changes(
            {"success": True, "bonds_formed": [[1, 2]], "key_bond_type": "C-C", "mapping_confidence": 0.2}
        )
        assert msg is not None
        assert "Low bond-mapping confidence" in msg


class TestRetrosynthesisHTETiming:
    def test_search_hte_precedent_surfaces_timing_breakdown(self, monkeypatch):
        from chem_coworker.tools.retrosynthesis import _search_hte_precedent

        monkeypatch.setattr(
            "chem_coworker.tools.retrosynthesis._map_reaction_to_family",
            lambda name: "suzuki_miyaura",
        )
        monkeypatch.setattr(
            "chem_coworker.tools.retrosynthesis._fast_load_hte_family_cached",
            lambda key: (
                {
                    "yield_value": 78.0,
                    "reaction_smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
                    "condition_core": "Pd/base",
                    "catalyst": {"name": "Pd"},
                    "base_uid": "K2CO3",
                    "solvent_uid": "dioxane",
                    "reagents": [],
                    "solvents": [],
                    "reference": "ref",
                    "rxn_type": "suzuki_miyaura",
                    "source_file": "demo.csv",
                },
            ),
        )

        result = _search_hte_precedent(
            target_smiles="c1ccc(-c2ccccc2)cc1",
            reaction_name="suzuki",
            top_k=1,
        )
        assert result["success"] is True
        timing = result["hte_search_timing_ms"]
        for key in ("input_parse_ms", "load_family_ms", "sort_ms", "drfp_rerank_ms", "format_ms", "total_ms"):
            assert key in timing
            assert isinstance(timing[key], (int, float))


# ---------------------------------------------------------------------------
# _extract_provides_smiles tests
# ---------------------------------------------------------------------------

class TestExtractProvidesSmiles:
    def _make_plugin(self, name: str, provides: list):
        from chem_coworker.tools._base import ToolPlugin
        return ToolPlugin(name=name, category="c", description="d", fn=lambda: {}, provides=provides)

    def test_finds_resolved_smiles(self):
        from chem_coworker.executor import _extract_provides_smiles
        plugin = self._make_plugin("resolve_to_smiles", ["resolved_smiles", "smiles"])
        plugins = {"resolve_to_smiles": plugin}
        group_results = {
            "resolve_to_smiles": {"success": True, "resolved_smiles": "c1ccccc1"}
        }
        assert _extract_provides_smiles(group_results, plugins) == "c1ccccc1"

    def test_falls_back_to_smiles_key(self):
        from chem_coworker.executor import _extract_provides_smiles
        plugin = self._make_plugin("resolver", ["smiles"])
        plugins = {"resolver": plugin}
        group_results = {
            "resolver": {"success": True, "smiles": "CCO"}
        }
        assert _extract_provides_smiles(group_results, plugins) == "CCO"

    def test_skips_failed_tool(self):
        from chem_coworker.executor import _extract_provides_smiles
        plugin = self._make_plugin("resolve_to_smiles", ["resolved_smiles"])
        plugins = {"resolve_to_smiles": plugin}
        group_results = {
            "resolve_to_smiles": {"success": False, "resolved_smiles": "c1ccccc1"}
        }
        assert _extract_provides_smiles(group_results, plugins) is None

    def test_skips_placeholder_value(self):
        from chem_coworker.executor import _extract_provides_smiles
        plugin = self._make_plugin("resolve_to_smiles", ["resolved_smiles"])
        plugins = {"resolve_to_smiles": plugin}
        group_results = {
            "resolve_to_smiles": {"success": True, "resolved_smiles": "(none)"}
        }
        assert _extract_provides_smiles(group_results, plugins) is None

    def test_no_matching_provides_key_returns_none(self):
        from chem_coworker.executor import _extract_provides_smiles
        plugin = self._make_plugin("detect_reaction_type", ["reaction_type"])
        plugins = {"detect_reaction_type": plugin}
        group_results = {
            "detect_reaction_type": {"success": True, "reaction_type": "Suzuki-Miyaura"}
        }
        # reaction_type is not in the SMILES keys — should return None
        assert _extract_provides_smiles(group_results, plugins) is None

    def test_empty_group_returns_none(self):
        from chem_coworker.executor import _extract_provides_smiles
        assert _extract_provides_smiles({}, {}) is None

    def test_plugin_not_in_registry_returns_none(self):
        from chem_coworker.executor import _extract_provides_smiles
        group_results = {
            "unknown_tool": {"success": True, "resolved_smiles": "c1ccccc1"}
        }
        assert _extract_provides_smiles(group_results, {}) is None


# ---------------------------------------------------------------------------
# _execute_one validator integration
# ---------------------------------------------------------------------------

class TestExecuteOneRunsValidators:
    def _make_executor(self, registry_plugins: dict):
        from chem_coworker.executor import ToolExecutor
        executor = ToolExecutor.__new__(ToolExecutor)
        executor.max_workers = 4
        executor.verbose = False
        executor.progress_cb = None
        executor.hooks = None
        executor.runtime_context = None
        mock_registry = MagicMock()
        mock_registry._plugins = registry_plugins
        executor.registry = mock_registry
        return executor

    def test_validator_warning_appended_to_warnings(self):
        from chem_coworker.executor import ToolExecutor
        from chem_coworker.tools._base import ToolPlugin
        from chem_coworker.plan import ToolCall

        warning_text = "No HTE precedents found"
        plugin = ToolPlugin(
            name="recommend_conditions", category="c", description="d",
            fn=lambda reaction_smiles: {"success": True, "recommendations": []},
            validators=[lambda r: warning_text if not r.get("recommendations") else None],
        )
        executor = self._make_executor({"recommend_conditions": plugin})

        call = ToolCall(name="recommend_conditions", args={"reaction_smiles": "CCO>>CC=O"})
        callables = {"recommend_conditions": plugin.fn}

        result = executor._execute_one(call, callables)
        assert "_warnings" in result
        assert warning_text in result["_warnings"]

    def test_passing_validator_adds_no_warnings(self):
        from chem_coworker.executor import ToolExecutor
        from chem_coworker.tools._base import ToolPlugin
        from chem_coworker.plan import ToolCall

        plugin = ToolPlugin(
            name="good_tool", category="c", description="d",
            fn=lambda: {"success": True, "data": "ok"},
            validators=[lambda r: None],
        )
        executor = self._make_executor({"good_tool": plugin})

        call = ToolCall(name="good_tool", args={})
        callables = {"good_tool": plugin.fn}

        result = executor._execute_one(call, callables)
        assert result.get("_warnings", []) == []

    def test_validator_exception_does_not_abort(self):
        from chem_coworker.executor import ToolExecutor
        from chem_coworker.tools._base import ToolPlugin
        from chem_coworker.plan import ToolCall

        def crashing_validator(r):
            raise RuntimeError("validator exploded")

        plugin = ToolPlugin(
            name="crashy_tool", category="c", description="d",
            fn=lambda: {"success": True},
            validators=[crashing_validator],
        )
        executor = self._make_executor({"crashy_tool": plugin})

        call = ToolCall(name="crashy_tool", args={})
        callables = {"crashy_tool": plugin.fn}

        # Should not raise
        result = executor._execute_one(call, callables)
        assert result.get("success") is True


class TestToolRuntimeContextInjection:
    def _make_executor(self):
        from chem_coworker.executor import ToolExecutor
        executor = ToolExecutor.__new__(ToolExecutor)
        executor.max_workers = 4
        executor.verbose = False
        executor.progress_cb = None
        executor.hooks = None
        executor.runtime_context = None
        executor.registry = MagicMock(_plugins={})
        return executor

    def test_executor_sets_contextvar_for_tool_call(self):
        from chem_coworker.plan import ToolCall
        from chem_coworker.tool_runtime import get_current_tool_runtime_context

        marker = object()

        def tool_fn():
            return {"success": True, "has_ctx": get_current_tool_runtime_context() is marker}

        executor = self._make_executor()
        call = ToolCall(name="ctx_tool", args={})
        result = executor._execute_one(call, {"ctx_tool": tool_fn}, runtime_context=marker)
        assert result["success"] is True
        assert result["has_ctx"] is True

    def test_chemistry_tool_prefers_runtime_context_short_circuit(self, monkeypatch):
        from chem_coworker.tools.chemistry import _detect_reaction_type

        class _FakeRuntimeContext:
            def detect_reaction_type(self, reaction_smiles: str):
                return {"success": True, "reaction_smiles": reaction_smiles, "cached": True}

        monkeypatch.setattr(
            "chem_coworker.tool_runtime.get_current_tool_runtime_context",
            lambda: _FakeRuntimeContext(),
        )
        result = _detect_reaction_type("A >> B")
        assert result["success"] is True
        assert result["cached"] is True

    def test_conditions_tool_prefers_runtime_context_cache(self, monkeypatch):
        from chem_coworker.tools.conditions import _recommend_conditions

        class _FakeRuntimeContext:
            def get_cached_conditions(self, reaction_smiles: str, top_k: int):
                return {"success": True, "reaction_smiles": reaction_smiles, "top_k": top_k, "cached": True}

        monkeypatch.setattr(
            "chem_coworker.tool_runtime.get_current_tool_runtime_context",
            lambda: _FakeRuntimeContext(),
        )
        result = _recommend_conditions("A>>B", top_k=7)
        assert result["success"] is True
        assert result["cached"] is True
        assert result["top_k"] == 7


class TestExecutorHTEHeavyConcurrencyGuard:
    def test_serializes_batch_when_multiple_hte_heavy_tools_present(self):
        from chem_coworker.executor import ToolExecutor
        from chem_coworker.plan import ToolCall

        executor = ToolExecutor.__new__(ToolExecutor)
        executor.max_workers = 4
        executor.verbose = False

        calls = [
            ToolCall(name="recommend_conditions", args={"reaction_smiles": "A>>B"}),
            ToolCall(name="search_hte_precedent", args={"target_smiles": "C"}),
            ToolCall(name="read_notes", args={"id": "suzuki"}),
        ]
        assert executor._compute_parallel_workers(calls) == 1

    def test_keeps_parallelism_for_non_heavy_batch(self):
        from chem_coworker.executor import ToolExecutor
        from chem_coworker.plan import ToolCall

        executor = ToolExecutor.__new__(ToolExecutor)
        executor.max_workers = 4
        executor.verbose = False

        calls = [
            ToolCall(name="read_notes", args={"id": "x"}),
            ToolCall(name="inspect_target", args={"smiles": "CCO"}),
            ToolCall(name="detect_reaction_type", args={"reaction_smiles": "A>>B"}),
        ]
        assert executor._compute_parallel_workers(calls) == 3
