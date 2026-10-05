"""Ordered starting-material decisions, exact indexed lookup and explicit uncertainty."""

import json
import sqlite3

import pytest

from condition_registry.models import Substance
from condition_registry.resolver import ConditionRegistry
from condition_recommender.fragment_index import build_fragment_index
from condition_recommender.fragment_search import lookup_exact_product
from condition_recommender.starting_materials import (
    POLICY_PATH, assess_starting_material, load_starting_material_policy,
)
from reactive_taxonomy.material_identity import identify_material


def registry(*entries):
    return ConditionRegistry(substances=[
        Substance(substance_id=identity, canonical_name=identity, cas=None, smiles=smiles)
        for identity, smiles in entries
    ])


@pytest.fixture
def index(tmp_path):
    products = ["CCO", "CCO", "CCCO", "C[C@H](O)F", "[13CH3]O", "C[NH3+]", "CC=O"]
    source = tmp_path / "source.jsonl"
    source.write_text("\n".join(json.dumps({
        "observation_id": f"obs-{i}", "reaction_id": f"rxn-{i}",
        "reference_id": f"ref-{i}", "reaction_smiles": f"C>>{product}",
    }) for i, product in enumerate(products)), encoding="utf-8")
    path = tmp_path / "fragments.sqlite"
    build_fragment_index(source, path)
    return path


def test_registry_has_priority_and_avoids_literature_lookup(index, monkeypatch):
    import condition_recommender.starting_materials as module

    def unexpected(*args, **kwargs):
        pytest.fail("Registry match must short-circuit literature lookup")

    monkeypatch.setattr(module, "lookup_exact_product", unexpected)
    result = assess_starting_material("[OH:3][CH2:2][CH3:1]", fragment_index=index,
                                      registry=registry(("ethanol", "CCO")))
    assert result.stop_expansion and result.stop_reason == "registry_match"
    assert result.registry["substance"]["substance_id"] == "ethanol"
    assert result.literature["status"] == "not_run"
    assert result.status == "assumed_terminal" and result.availability == "unknown"
    assert "unverified" in result.warnings[0]


def test_exact_literature_match_precedes_mass_and_preserves_ids(index):
    result = assess_starting_material("OCC", registry=registry(), fragment_index=index,
                                      mw_threshold=1.0)
    assert result.stop_reason == "exact_literature_match"
    assert result.literature["precedents"][0] == {
        "observation_id": "obs-0", "reaction_id": "rxn-0", "reference_id": "ref-0",
        "product_component_index": 0,
        "matched_side": "product", "match_extent": "whole_molecule",
    }
    assert result.literature["index_id"].startswith("FPI1:")
    assert result.availability == "unknown"


def test_exact_lookup_uses_sql_index_without_deserializing_library(index):
    with sqlite3.connect(index) as connection:
        connection.execute("UPDATE library SET payload=?", (b"not an RDKit library",))
        query_plan = connection.execute(
            "EXPLAIN QUERY PLAN SELECT id FROM products WHERE smiles=?", ("CCO",),
        ).fetchall()
    assert any("INDEX" in str(row) and "SEARCH" in str(row) for row in query_plan)
    result = lookup_exact_product(index, "[OH:1][CH2:2][CH3:3]", limit=1)
    assert result["status"] == "matched" and result["has_more"]
    assert len(result["precedents"]) == 1
    assert result == lookup_exact_product(index, "[OH:1][CH2:2][CH3:3]", limit=1)


@pytest.mark.parametrize("smiles", ["CO", "C[C@@H](O)F", "CC(O)F", "CN", "[12CH3]O",
                                    "C=CO", "CCO.[Na+]", "[OH-]"])
def test_exact_lookup_rejects_substructure_and_different_chemical_forms(index, smiles):
    assert lookup_exact_product(index, smiles)["status"] == "not_found"


def test_mass_cutoff_is_strict_and_configurable(index):
    mass = identify_material("CCCC").molecular_weight
    options = {"registry": registry(), "fragment_index": index}
    equal = assess_starting_material("CCCC", mw_threshold=mass, **options)
    assert not equal.stop_expansion and equal.stop_reason == "no_stopping_criterion"
    above = assess_starting_material("CCCC", mw_threshold=mass + 0.001, **options)
    assert above.stop_reason == "low_molecular_weight"
    assert above.registry["status"] == "unresolved"
    assert above.literature["status"] == "not_found"
    assert assess_starting_material("C" * 20, **options).stop_expansion is False


def test_disabled_stages_and_default_policy(index):
    result = assess_starting_material("CCO", registry=registry(("ethanol", "CCO")),
                                      fragment_index=index, allow_registry_stop=False,
                                      allow_literature_stop=False, allow_mw_stop=False)
    assert result.mw_threshold == 200.0 and not result.stop_expansion
    assert result.registry["reason"] == result.literature["reason"] == "disabled_by_policy"
    assert result.to_dict()["definition_version"] == "starting_material_policy.v1@1.0"


def test_explicit_exclusion_overrides_all_evidence(index):
    result = assess_starting_material("CCO", registry=registry(("ethanol", "CCO")),
                                      fragment_index=index, unavailable_starting_materials=["OCC"])
    assert not result.stop_expansion and result.stop_reason == "explicitly_unavailable"
    assert result.status == "excluded"
    assert result.registry["status"] == result.literature["status"] == "not_run"


def test_ambiguous_and_inactive_registry_records_do_not_prove_obtainability():
    ambiguous = assess_starting_material("CCO", registry=registry(("a", "CCO"), ("b", "OCC")),
                                         allow_literature_stop=False, allow_mw_stop=False)
    assert not ambiguous.stop_expansion and ambiguous.registry["status"] == "ambiguous"
    assert ambiguous.registry["candidates"] == ("a", "b")
    inactive_registry = ConditionRegistry(substances=[
        Substance("old", "old", None, "CCO", status="deprecated"),
    ])
    inactive = assess_starting_material("CCO", registry=inactive_registry,
                                        allow_literature_stop=False, allow_mw_stop=False)
    assert not inactive.stop_expansion and inactive.registry["status"] == "inactive"


@pytest.mark.parametrize("kind", ["unconfigured", "missing", "corrupt", "incompatible"])
def test_lookup_failures_remain_visible_during_mass_fallback(index, tmp_path, kind):
    path = index
    if kind == "unconfigured":
        path = None
    elif kind == "missing":
        path = tmp_path / "missing.sqlite"
    elif kind == "corrupt":
        index.write_bytes(b"invalid SQLite")
    else:
        with sqlite3.connect(index) as connection:
            manifest = json.loads(connection.execute("SELECT payload FROM metadata").fetchone()[0])
            manifest["rdkit_version"] = "incompatible"
            connection.execute("UPDATE metadata SET payload=?", (json.dumps(manifest),))
    result = assess_starting_material("CCCC", registry=registry(), fragment_index=path)
    assert result.stop_reason == "low_molecular_weight"
    assert result.literature["status"] == ("unavailable" if kind in {"unconfigured", "missing"} else "error")
    assert any("absence" in warning for warning in result.warnings)


def test_registry_failure_does_not_disappear_in_literature_fallback(index, monkeypatch):
    import condition_recommender.starting_materials as module

    def unavailable():
        raise OSError("Registry source is unavailable")

    monkeypatch.setattr(module, "get_registry", unavailable)
    result = assess_starting_material("CCO", fragment_index=index)
    assert result.stop_reason == "exact_literature_match"
    assert result.registry["status"] == "error"


@pytest.mark.parametrize("options", [
    {"mw_threshold": value} for value in (0, -1, True, float("nan"), float("inf"), "200")
] + [{"allow_mw_stop": 1}, {"unavailable_starting_materials": "CCO"}])
def test_bad_policy_overrides_are_rejected(options):
    with pytest.raises(ValueError):
        assess_starting_material("CCO", **options)


@pytest.mark.parametrize("field,value", [
    ("default_mw_threshold", -1), ("precedent_limit", 0), ("schema_version", "future"),
    ("check_order", ["molecular_weight", "registry", "literature"]),
    ("availability", "verified"), ("mw_comparison", "less_or_equal"),
    ("definition_version", ""), ("unknown_field", True),
])
def test_definition_loader_validates_policy(tmp_path, field, value):
    payload = json.loads(POLICY_PATH.read_text(encoding="utf-8"))
    payload[field] = value
    path = tmp_path / "policy.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError):
        load_starting_material_policy(path)
