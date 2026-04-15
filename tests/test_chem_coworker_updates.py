"""
Unit tests for chem_coworker improvements:
  1. Diversity filter in _recommend_conditions
  2. _collect_caveats (replaced _check_hypothesis in Phase 5)
"""
import pandas as pd
import pytest

from chem_coworker.tools.conditions import _diversify, _extract_metal, _recommend_conditions
from chem_coworker.agent import ChemCoworker
from chem_coworker.tools import REGISTRY
from chemtools.recommend.recommender import (
    _apply_condition_plausibility_penalties,
    _extract_catalyst_metal_tokens,
    _row_contains_unexpected_catalyst_metals,
    _row_has_implausible_snar_annotation,
)


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


class TestRecommendConditionsSelectionMode:
    @staticmethod
    def _payload():
        return {
            "detection": {"detected_reaction_type": "Suzuki_miyaura", "confidence": 0.42},
            "extras": {
                "hte": {
                    "reaction_type_confidence": 0.42,
                    "is_fallback_match": True,
                    "matched_motifs": ["RXNEVT|sig=LGDisp+C-C|form=C-C"],
                    "total_matching_experiments": 3,
                    "database_coverage": 0.12,
                    "catalyst_requirement_enforced": True,
                    "required_catalyst_family": "Suzuki_miyaura",
                    "required_catalyst_classes": ["Pd"],
                    "filtered_missing_catalyst_rows": 1,
                    "retained_missing_catalyst_rows": 0,
                    "condition_quality_family": "Suzuki_miyaura",
                    "penalized_incomplete_condition_rows": 1,
                    "missing_required_condition_fields": {"base": 1},
                }
            },
            "recommended_conditions": [
                {
                    "rank": 1,
                    "conditions": {"catalyst": "Pd(OAc)2", "ligand": "SPhos", "base": "K3PO4", "solvent": "dioxane"},
                    "confidence": 0.91,
                    "component_scores": {"success_rate": 0.9, "avg_yield": 88.0, "median_yield": 89.0, "match_score": 0.95, "num_experiments": 6, "condition_quality_score": 1.0},
                    "reaction": "Suzuki_miyaura",
                },
                {
                    "rank": 2,
                    "conditions": {"catalyst": "Pd2(dba)3", "ligand": "XPhos", "base": "K3PO4", "solvent": "dioxane"},
                    "confidence": 0.84,
                    "component_scores": {"success_rate": 0.82, "avg_yield": 81.0, "median_yield": 82.0, "match_score": 0.92, "num_experiments": 5, "condition_quality_score": 1.0},
                    "reaction": "Suzuki_miyaura",
                },
                {
                    "rank": 3,
                    "conditions": {"catalyst": "CuI", "ligand": "phen", "base": "Cs2CO3", "solvent": "DMSO"},
                    "confidence": 0.66,
                    "component_scores": {"success_rate": 0.7, "avg_yield": 67.0, "median_yield": 66.0, "match_score": 0.75, "num_experiments": 1, "condition_quality_score": 1.0},
                    "reaction": "Suzuki_miyaura",
                },
            ],
            "recommendations": [],
        }

    def test_best_mode_preserves_rank_order(self, monkeypatch):
        monkeypatch.setattr(
            "chemtools.recommend.hte_adapter.recommend_from_reaction",
            lambda *args, **kwargs: self._payload(),
        )

        result = _recommend_conditions("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1", top_k=2, selection_mode="best")

        assert result["success"] is True
        assert [rec["catalyst"] for rec in result["recommendations"]] == ["Pd(OAc)2", "Pd2(dba)3"]
        assert result["selection_mode"] == "best"

    def test_diverse_mode_interleaves_metals_and_surfaces_evidence(self, monkeypatch):
        monkeypatch.setattr(
            "chemtools.recommend.hte_adapter.recommend_from_reaction",
            lambda *args, **kwargs: self._payload(),
        )

        result = _recommend_conditions("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1", top_k=2, selection_mode="diverse")

        assert result["success"] is True
        assert [rec["catalyst"] for rec in result["recommendations"]] == ["Pd(OAc)2", "CuI"]
        assert result["selection_mode"] == "diverse"
        assert result["reaction_type_confidence"] == pytest.approx(0.42)
        assert result["evidence"]["is_fallback_match"] is True
        assert result["evidence"]["filtered_missing_catalyst_rows"] == 1
        assert result["evidence"]["penalized_incomplete_condition_rows"] == 1
        assert "RXNEVT|" in result["evidence"]["matched_transformation"]
        assert result["_warnings"]

    def test_best_mode_warns_when_only_incomplete_catalyst_rows_exist(self, monkeypatch):
        payload = self._payload()
        payload["extras"]["hte"]["filtered_missing_catalyst_rows"] = 0
        payload["extras"]["hte"]["retained_missing_catalyst_rows"] = 2

        monkeypatch.setattr(
            "chemtools.recommend.hte_adapter.recommend_from_reaction",
            lambda *args, **kwargs: payload,
        )

        result = _recommend_conditions("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1", top_k=2, selection_mode="best")

        assert result["success"] is True
        assert result["evidence"]["retained_missing_catalyst_rows"] == 2
        assert any("missing required catalyst annotations" in warning for warning in result["_warnings"])

    def test_best_mode_warns_when_top_recommendation_is_incomplete(self, monkeypatch):
        payload = self._payload()
        payload["recommended_conditions"][0]["conditions"]["base"] = ""
        payload["recommended_conditions"][0]["missing_required_fields"] = ["base"]
        payload["recommended_conditions"][0]["component_scores"]["condition_quality_score"] = 0.55

        monkeypatch.setattr(
            "chemtools.recommend.hte_adapter.recommend_from_reaction",
            lambda *args, **kwargs: payload,
        )

        result = _recommend_conditions("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1", top_k=2, selection_mode="best")

        assert result["success"] is True
        assert result["recommendations"][0]["missing_required_fields"] == ["base"]
        assert any("Top condition recommendation is still missing required fields" in warning for warning in result["_warnings"])

    def test_best_mode_surfaces_plausibility_warnings(self, monkeypatch):
        payload = self._payload()
        payload["extras"]["hte"]["plausibility_family"] = "C_N_Coupling"
        payload["extras"]["hte"]["filtered_implausible_catalyst_rows"] = 2
        payload["extras"]["hte"]["penalized_implausible_mechanism_rows"] = 3
        payload["extras"]["hte"]["plausibility_issue_counts"] = {
            "unexpected_catalyst_metals": 2,
            "implausible_snar_annotation": 3,
        }

        monkeypatch.setattr(
            "chemtools.recommend.hte_adapter.recommend_from_reaction",
            lambda *args, **kwargs: payload,
        )

        result = _recommend_conditions("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1", top_k=2, selection_mode="best")

        assert result["success"] is True
        assert result["evidence"]["filtered_implausible_catalyst_rows"] == 2
        assert result["evidence"]["penalized_implausible_mechanism_rows"] == 3
        assert any("chemically implausible" in warning for warning in result["_warnings"])
        assert any("SNAr-style mechanism labels" in warning for warning in result["_warnings"])


class TestPlausibilityFilters:
    def test_extract_catalyst_metal_tokens_detects_mixed_system(self):
        assert _extract_catalyst_metal_tokens("Cu/[(S,S)-Teth-MsDPEN-RuCl]") == ("Cu", "Ru")

    def test_row_contains_unexpected_catalyst_metals(self):
        row = {"Catalyst": "Cu/[(S,S)-Teth-MsDPEN-RuCl]", "Catalyst_Type": ""}
        assert _row_contains_unexpected_catalyst_metals(row, allowed_classes=("Cu", "Pd", "Ni")) is True

    def test_row_has_implausible_snar_annotation_for_aryl_iodide(self):
        row = {
            "Reaction_Events": "context:LG=I, Nu=N, mech=SNAr | summary:fam=substitution",
            "Reaction_Key": "|Ar-I|AromN-H -> Ar-AromN | spectators: Ar-OR",
        }
        assert _row_has_implausible_snar_annotation(row) is True

    def test_apply_condition_plausibility_penalties_filters_and_penalizes(self):
        df = pd.DataFrame(
            [
                {
                    "Reaction_Type_Standardized": "C_N_Coupling",
                    "Source_Group": "literature",
                    "Catalyst": "Cu/[(S,S)-Teth-MsDPEN-RuCl]",
                    "Catalyst_Type": "",
                    "Reaction_Events": "context:LG=I, Nu=N, mech=SNAr | summary:fam=substitution",
                    "Reaction_Key": "|Ar-I|AromN-H -> Ar-AromN | spectators: Ar-OR",
                    "match_score": 1.0,
                },
                {
                    "Reaction_Type_Standardized": "C_N_Coupling",
                    "Source_Group": "literature",
                    "Catalyst": "CuI",
                    "Catalyst_Type": "",
                    "Reaction_Events": "context:LG=I, Nu=N, mech=SNAr | summary:fam=substitution",
                    "Reaction_Key": "|Ar-I|AromN-H -> Ar-AromN | spectators: Ar-OR",
                    "match_score": 1.0,
                },
                {
                    "Reaction_Type_Standardized": "C_N_Coupling",
                    "Source_Group": "literature",
                    "Catalyst": "CuI",
                    "Catalyst_Type": "",
                    "Reaction_Events": "context:LG=Br, Nu=N, mech=oa_based_coupling | summary:fam=substitution",
                    "Reaction_Key": "|Ar-Br|AromN-H -> Ar-AromN",
                    "match_score": 1.0,
                },
            ]
        )

        out_df, family, filtered_rows, penalized_rows, issue_counts = _apply_condition_plausibility_penalties(
            df,
            reaction_type_filter="C_N_Coupling",
            predicted_reaction_type="C_N_Coupling",
            reaction_type_confidence=0.95,
        )

        assert family == "C_N_Coupling"
        assert filtered_rows == 1
        assert penalized_rows == 1
        assert issue_counts["unexpected_catalyst_metals"] == 1
        assert issue_counts["implausible_snar_annotation"] == 1
        assert len(out_df) == 2
        assert out_df.iloc[0]["match_score"] == pytest.approx(0.35)
        assert out_df.iloc[1]["match_score"] == pytest.approx(1.0)


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
