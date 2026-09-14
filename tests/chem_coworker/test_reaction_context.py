"""Actual graph planning plus bounded condition enrichment of alternatives."""

from types import SimpleNamespace
from dataclasses import replace

import pytest

from chem_coworker.reaction_context import explore_reaction_context
from core_retrosynthesis import (
    build_generic_library,
    build_forward_library_from_generic,
)
from core_retrosynthesis import disconnect_strategies_detailed
import chem_coworker.reaction_context as context_module


@pytest.fixture(scope="module")
def libraries():
    retro = build_generic_library(
        (
            {
                "reaction_id": "bromide",
                "reference_id": "ref-1",
                "reaction_smiles": "CCBr.N>>CCN",
            },
        ),
        levels=("L1", "L2"),
        admission_mode="data_driven",
    )
    return retro, build_forward_library_from_generic(retro)


def test_actual_operators_propose_precursors_with_separate_conditions(libraries):
    retro, forward = libraries
    calls = []

    def recommend(reaction, **kwargs):
        calls.append((reaction, kwargs))
        return SimpleNamespace(to_dict=lambda: {"recommendations": []})

    result = explore_reaction_context(
        "CCI.N>>CCN",
        forward_library=lambda: forward,
        retro_library=lambda: retro,
        recommend=recommend,
    )
    assert result["advisory_only"]
    assert result["reactant_analysis"]["status"] == "complete"
    assert result["product_analysis"]["status"] == "complete"
    alternatives = result["product_analysis"]["alternatives"]
    assert alternatives and calls
    assert all(a["evidence_kind"] == "generated_hypothesis" for a in alternatives)
    assert any(a["relation"]["kind"] == "precursor_alternative" for a in alternatives)
    assert all(reaction != "CCI.N>>CCN" for reaction, _ in calls)
    assert all(kw["preferred_reaction_ids"] == ("bromide",) for _, kw in calls)
    assert len(calls) <= result["search_limits"]["alternatives"]


def test_missing_forward_preserves_retro_and_condition_failures(libraries):
    retro, _ = libraries

    def unavailable():
        raise FileNotFoundError("no forward artifact")

    def recommend(*args, **kwargs):
        raise RuntimeError("no condition artifact")

    result = explore_reaction_context(
        "CCI.N>>CCN",
        forward_library=unavailable,
        retro_library=lambda: retro,
        recommend=recommend,
    )
    assert result["reactant_analysis"]["status"] == "unavailable"
    assert result["product_analysis"]["status"] == "complete"
    assert any(
        a.get("condition_error") == "no condition artifact"
        for a in result["product_analysis"]["alternatives"]
    )


def test_multiple_products_do_not_silently_select_one():
    def must_not_run(*args, **kwargs):
        raise AssertionError("planner must not run")

    result = explore_reaction_context(
        "CCBr.N>>CCN.O",
        forward_library=must_not_run,
        retro_library=must_not_run,
        recommend=must_not_run,
    )
    assert result["reactant_analysis"]["status"] == "out_of_scope"
    assert result["product_analysis"]["error"] == "SINGLE_PRODUCT_REQUIRED"


def test_duplicate_proposals_merge_provenance_and_unverified_are_excluded(
    libraries, monkeypatch
):
    retro, _ = libraries
    candidate = (
        disconnect_strategies_detailed("CCN", retro).strategies[0].representative
    )
    duplicate = replace(candidate, precedent_reaction_ids=("second-source",))
    unverified = replace(candidate, forward_validation_status="unverified")
    search = SimpleNamespace(
        strategies=(
            SimpleNamespace(
                representative=candidate, alternate_realizations=(duplicate, unverified)
            ),
        ),
        to_dict=lambda: {"search_diagnostics": {}},
    )
    monkeypatch.setattr(
        context_module, "disconnect_strategies_detailed", lambda *a, **k: search
    )
    calls = []

    def unavailable():
        raise FileNotFoundError("absent")

    def recommend(reaction, **kwargs):
        calls.append(kwargs)
        return SimpleNamespace(to_dict=lambda: {"recommendations": []})

    result = explore_reaction_context(
        "CCI.N>>CCN",
        forward_library=unavailable,
        retro_library=lambda: retro,
        recommend=recommend,
    )["product_analysis"]
    assert len(result["alternatives"]) == len(calls) == 1
    assert calls[0]["preferred_reaction_ids"] == ("bromide", "second-source")
    assert result["excluded_unverified_count"] == 1


def test_unresolved_query_does_not_receive_alternative_recipe_calls(
    libraries, monkeypatch
):
    from condition_recommender.reaction_context import context_core

    retro, forward = libraries
    blocked = replace(
        context_core("CCI.N>>CCN"),
        levels=(),
        unavailable_reasons=("CONFLICTING_EDITS",),
    )
    monkeypatch.setattr(context_module, "context_core", lambda _: blocked)

    def must_not_recommend(*args, **kwargs):
        raise AssertionError("unresolved graph relation must not transfer conditions")

    result = explore_reaction_context(
        "CCI.N>>CCN",
        forward_library=lambda: forward,
        retro_library=lambda: retro,
        recommend=must_not_recommend,
    )
    assert result["product_analysis"]["alternatives"]
    assert all(
        a["relation"]["kind"] == "unresolved"
        for a in result["product_analysis"]["alternatives"]
    )
