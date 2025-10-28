"""Regression tests for the LangChain ChemTools wrappers."""

import lang_chain.chemtools_wrapper as wrapper
from lang_chain.chemtools_wrapper import (
    get_tool_descriptions,
    normalize_smiles_tool,
    normalize_reaction_tool,
)


def test_normalize_smiles_schema_contains_required_field():
    """Tool schema should expose the SMILES parameter for function-calling."""
    schema = normalize_smiles_tool.args_schema.model_json_schema()
    assert "smiles" in schema["properties"]
    assert "smiles" in schema.get("required", [])
    assert "normalize" in schema["properties"]["smiles"]["description"].lower()


def test_normalize_smiles_success_response():
    """Normalization should return a success flag and canonical SMILES."""
    result = normalize_smiles_tool.invoke({"smiles": "c1ccccc1"})
    assert result["success"], result
    assert result["smiles_norm"] == "c1ccccc1"
    assert "input" in result  # underlying metadata is preserved


def test_normalize_smiles_error_payload():
    """Invalid SMILES should produce a structured error payload."""
    result = normalize_smiles_tool.invoke({"smiles": "???"})
    assert not result["success"]
    assert "error" in result
    assert "details" in result


def test_normalize_reaction_returns_structured_payload():
    """Reaction normalization should keep success metadata."""
    payload = normalize_reaction_tool.invoke(
        {"reaction_smiles": "CCO.O>>CCO"}
    )
    assert payload["success"], payload
    assert "normalized" in payload
    assert payload["normalized"]


def test_get_tool_descriptions_includes_parameter_metadata():
    """Tool registry should surface parameter documentation for the CLI."""
    tool_meta = get_tool_descriptions()
    summary = next(item for item in tool_meta if item["name"] == "normalize_smiles_tool")
    smiles_param = next(param for param in summary["parameters"] if param["name"] == "smiles")
    assert smiles_param["required"] is True
    assert "normalize" in (smiles_param.get("description") or "").lower()


def test_recommendation_cache_reuse(monkeypatch):
    """recommend_conditions_tool and list_supported_cores_tool should share cache."""
    wrapper._recommendation_cache.clear()
    monkeypatch.setattr(wrapper, "_normalized_reaction_key", lambda _: "TEST_KEY")

    calls = []

    def fake_recommend(**kwargs):
        calls.append((kwargs["k"], kwargs["max_variants"]))
        return {
            "family": "Test",
            "detected_family": "Test",
            "recommendation": {"core": "Pd/XPhos"},
            "alternatives": {"cores": [("Pd/XPhos", 2)]},
            "precedent_pack": {"precedents": [1, 2]},
        }

    monkeypatch.setattr(wrapper, "recommend_from_reaction", fake_recommend)

    base_payload = {"reaction_smiles": "CCBr.CCO>>CCOCC"}

    first = wrapper.recommend_conditions_tool.invoke({**base_payload, "k": 5, "max_variants": 2})
    assert first["success"]

    second = wrapper.list_supported_cores_tool.invoke({**base_payload, "k": 5})
    assert second["success"]
    assert len(calls) == 1, "Expected cached recommendation to be reused"
    stats = wrapper.recommendation_cache_stats()
    assert stats["entries"] == 1
    wrapper.clear_recommendation_cache()
    assert wrapper.recommendation_cache_stats()["entries"] == 0
