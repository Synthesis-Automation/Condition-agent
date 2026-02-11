from __future__ import annotations

from chemtools.util.ingestion_routing import classify_reaction_for_taxonomy_benchmark


def test_classify_reaction_for_taxonomy_benchmark_eligible() -> None:
    result = classify_reaction_for_taxonomy_benchmark("CCO>>CCBr")
    assert result["excluded"] is False
    assert result["route"] == "eligible_taxonomy_benchmark"
    assert result["reason"] == "eligible"
    assert result["reactant_count"] == 1
    assert result["product_count"] == 1


def test_classify_reaction_for_taxonomy_benchmark_coordination_complex() -> None:
    result = classify_reaction_for_taxonomy_benchmark("CC[Ni]->CC>>CC")
    assert result["excluded"] is True
    assert result["route"] == "exclude_organometallic_or_coordination_complex"
    assert result["reason"] == "coordination_arrow_and_metal_token"
    assert "Ni" in result["metal_tokens"]


def test_classify_reaction_for_taxonomy_benchmark_unparseable_component() -> None:
    result = classify_reaction_for_taxonomy_benchmark("C1CC>>CC")
    assert result["excluded"] is True
    assert result["route"] == "exclude_unparseable_components"
    assert result["reason"] == "rdkit_component_parse_failures"
    assert result["unparseable_components"] >= 1
