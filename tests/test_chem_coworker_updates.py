"""
Unit tests for chem_coworker improvements:
  1. Diversity filter in _recommend_conditions
  2. _collect_caveats (replaced _check_hypothesis in Phase 5)
"""
import pytest

from chem_coworker.tools.conditions import _extract_metal, _diversify
from chem_coworker.agent import ChemCoworker
from chem_coworker.tools import REGISTRY


# ---------------------------------------------------------------------------
# Improvement 1: _extract_metal and _diversify
# ---------------------------------------------------------------------------

class TestExtractMetal:
    def test_pd_catalyst(self):
        assert _extract_metal("Pd(OAc)2") == "Pd"

    def test_pd2dba3(self):
        assert _extract_metal("Pd2dba3") == "Pd"

    def test_ni_catalyst(self):
        assert _extract_metal("NiCl2(dme)") == "Ni"

    def test_copper(self):
        assert _extract_metal("CuI") == "Cu"

    def test_plain_pd(self):
        assert _extract_metal("Pd") == "Pd"

    def test_unknown_returns_other(self):
        assert _extract_metal("TBAF") == "Other"

    def test_empty_string(self):
        assert _extract_metal("") == "Other"

    def test_iridium(self):
        assert _extract_metal("Ir(ppy)3") == "Ir"


class TestDiversify:
    def _make_rec(self, catalyst, rank):
        return {"rank": rank, "catalyst": catalyst, "avg_yield": 80.0}

    def test_single_family_preserves_top_k(self):
        recs = [self._make_rec("Pd(OAc)2", i) for i in range(1, 8)]
        result = _diversify(recs, top_k=5)
        assert len(result) == 5
        assert all("Pd" in r["catalyst"] for r in result)

    def test_two_families_interleaves(self):
        pd_recs = [self._make_rec("Pd(OAc)2", i) for i in range(1, 4)]
        cu_recs = [self._make_rec("CuI", i + 10) for i in range(1, 4)]
        result = _diversify(pd_recs + cu_recs, top_k=4)
        metals = [_extract_metal(r["catalyst"]) for r in result]
        assert "Pd" in metals
        assert "Cu" in metals

    def test_three_families_all_represented_in_top3(self):
        recs = (
            [self._make_rec("Pd(OAc)2", 1)]
            + [self._make_rec("CuI", 2)]
            + [self._make_rec("NiCl2", 3)]
            + [self._make_rec("Pd2dba3", 4)]
        )
        result = _diversify(recs, top_k=3)
        metals = {_extract_metal(r["catalyst"]) for r in result}
        assert len(metals) >= 3

    def test_top_k_never_exceeded(self):
        recs = [self._make_rec("Pd", i) for i in range(20)]
        result = _diversify(recs, top_k=5)
        assert len(result) == 5

    def test_empty_input_returns_empty(self):
        assert _diversify([], top_k=5) == []

    def test_fewer_recs_than_top_k(self):
        recs = [self._make_rec("Pd", 1), self._make_rec("Cu", 2)]
        result = _diversify(recs, top_k=5)
        assert len(result) == 2

    def test_no_catalyst_field_doesnt_crash(self):
        recs = [{"rank": 1, "avg_yield": 70.0}]  # no "catalyst" key
        result = _diversify(recs, top_k=3)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# Improvement 2: _collect_caveats (Phase 5 replacement for _check_hypothesis)
# ---------------------------------------------------------------------------

def _make_agent():
    """Create a ChemCoworker without an LLM (provider/key not needed for unit tests)."""
    agent = object.__new__(ChemCoworker)
    agent.registry = REGISTRY
    return agent


class TestCollectCaveats:
    def test_empty_results_returns_empty_string(self):
        agent = _make_agent()
        assert agent._collect_caveats({}, []) == ""

    def test_no_warnings_in_results_returns_empty_string(self):
        agent = _make_agent()
        results = {"recommend_conditions": {"success": True, "recommendations": []}}
        assert agent._collect_caveats(results, []) == ""

    def test_warnings_from_tool_results_included(self):
        agent = _make_agent()
        results = {
            "recommend_conditions": {
                "success": True,
                "recommendations": [],
                "_warnings": ["NO HTE precedents found — Do NOT invent conditions"],
            }
        }
        caveats = agent._collect_caveats(results, [])
        assert "NO HTE" in caveats
        assert "recommend_conditions" in caveats

    def test_existing_warnings_included(self):
        agent = _make_agent()
        caveats = agent._collect_caveats({}, ["Tool X failed: timeout"])
        assert "Tool X failed: timeout" in caveats

    def test_both_sources_combined(self):
        agent = _make_agent()
        results = {
            "detect_reaction_type": {
                "success": True,
                "_warnings": ["no match found"],
            }
        }
        caveats = agent._collect_caveats(results, ["LLM call failed"])
        assert "no match found" in caveats
        assert "LLM call failed" in caveats

    def test_exact_duplicate_lines_deduplicated(self):
        agent = _make_agent()
        # Same existing warning passed twice
        caveats = agent._collect_caveats({}, ["dup warning", "dup warning"])
        lines = caveats.splitlines()
        assert lines.count("• dup warning") == 1

    def test_non_dict_tool_result_skipped(self):
        agent = _make_agent()
        results = {"weird_tool": "just a string result"}
        assert agent._collect_caveats(results, []) == ""

    def test_multiple_tools_with_warnings(self):
        agent = _make_agent()
        results = {
            "tool_a": {"_warnings": ["w1"]},
            "tool_b": {"_warnings": ["w2", "w3"]},
        }
        caveats = agent._collect_caveats(results, [])
        assert "w1" in caveats
        assert "w2" in caveats
        assert "w3" in caveats

    def test_zero_experiment_warning_surfaced(self):
        """Validator-stored zero-experiment warning is surfaced via _collect_caveats."""
        agent = _make_agent()
        results = {
            "recommend_conditions": {
                "success": True,
                "recommendations": [{"rank": 1, "num_experiments": 0, "confidence": 0.1}],
                "_warnings": ["0 experiments backing this recommendation — treat as tentative"],
            }
        }
        caveats = agent._collect_caveats(results, [])
        assert "0 experiments" in caveats
        assert "tentative" in caveats
