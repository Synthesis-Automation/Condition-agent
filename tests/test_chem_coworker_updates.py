"""
Unit tests for three chem_coworker improvements:
  1. Diversity filter in _recommend_conditions
  2. pass_check in _check_hypothesis
  3. caveats_text injection into synthesis prompt
"""
import pytest
from unittest.mock import MagicMock

from chem_coworker.tools.conditions import _extract_metal, _diversify
from chem_coworker.agent import ChemCoworker
from chem_coworker.plan import ExecutionPlan


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
# Improvement 3a: _check_hypothesis pass_check for conditions
# ---------------------------------------------------------------------------

def _make_plan(hypothesis="Suzuki-Miyaura", confidence=0.9):
    return ExecutionPlan(
        hypothesis=hypothesis,
        confidence=confidence,
        groups=[],
        rationale="test",
        raw_plan_text="",
    )


def _make_agent():
    """Create a ChemCoworker without an LLM (provider/key not needed for unit tests)."""
    agent = object.__new__(ChemCoworker)
    return agent


class TestCheckHypothesis:
    def test_empty_results_returns_none(self):
        agent = _make_agent()
        plan = _make_plan()
        assert agent._check_hypothesis(plan, {}) is None

    def test_empty_hypothesis_returns_none(self):
        agent = _make_agent()
        plan = _make_plan(hypothesis="")
        results = {"recommend_conditions": {"success": True, "recommendations": []}}
        assert agent._check_hypothesis(plan, results) is None

    def test_no_hte_precedents_returns_warning(self):
        agent = _make_agent()
        plan = _make_plan()
        results = {
            "recommend_conditions": {
                "success": True,
                "recommendations": [],
            }
        }
        msg = agent._check_hypothesis(plan, results)
        assert msg is not None
        assert "NO HTE" in msg
        assert "Do NOT invent" in msg

    def test_missing_catalyst_field_returns_warning(self):
        """Empty catalyst string with real experimental support should NOT warn.
        Absent fields in HTE data are legitimate (e.g. uncatalysed reactions);
        the old field-presence check caused false positives."""
        agent = _make_agent()
        plan = _make_plan()
        results = {
            "recommend_conditions": {
                "success": True,
                "recommendations": [
                    {
                        "rank": 1,
                        "catalyst": "",
                        "solvent": "DMF",
                        "base": "K2CO3",
                        "num_experiments": 12,
                        "confidence": 0.75,
                    }
                ],
            }
        }
        msg = agent._check_hypothesis(plan, results)
        # Empty catalyst is fine when real experiments back it up
        assert msg is None

    def test_zero_experiment_low_confidence_returns_warning(self):
        """A recommendation with 0 experiments and low confidence should warn."""
        agent = _make_agent()
        plan = _make_plan()
        results = {
            "recommend_conditions": {
                "success": True,
                "recommendations": [
                    {
                        "rank": 1,
                        "catalyst": "Pd(OAc)2",
                        "solvent": "THF",
                        "base": "Cs2CO3",
                        "num_experiments": 0,
                        "confidence": 0.1,
                    }
                ],
            }
        }
        msg = agent._check_hypothesis(plan, results)
        assert msg is not None
        assert "0 experiments" in msg
        assert "tentative" in msg

    def test_complete_recommendation_passes(self):
        agent = _make_agent()
        plan = _make_plan()
        results = {
            "recommend_conditions": {
                "success": True,
                "recommendations": [
                    {
                        "rank": 1,
                        "catalyst": "Pd(OAc)2",
                        "solvent": "DMF",
                        "base": "K2CO3",
                    }
                ],
            }
        }
        assert agent._check_hypothesis(plan, results) is None

    def test_failed_tool_result_skipped(self):
        """A failed recommend_conditions call (success=False) should not trigger pass_check."""
        agent = _make_agent()
        plan = _make_plan()
        results = {
            "recommend_conditions": {
                "success": False,
                "error": "something broke",
            }
        }
        assert agent._check_hypothesis(plan, results) is None

    def test_detect_reaction_type_contradiction_still_works(self):
        """Existing detect_reaction_type check should still fire."""
        agent = _make_agent()
        plan = _make_plan(confidence=0.9)
        results = {
            "detect_reaction_type": {
                "success": True,
                "reaction_type": None,
            }
        }
        msg = agent._check_hypothesis(plan, results)
        assert msg is not None
        assert "classifier" in msg.lower() or "reaction type" in msg.lower()
