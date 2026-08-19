"""Regression tests for the curated OpenAI model catalog."""

from chem_coworker._cli.config import DEFAULT_MODEL
from chem_coworker._cli.models import infer_provider, selectable_models
from llmtools.clients import AVAILABLE_MODELS, RECOMMENDED_MODELS


def test_current_openai_model_family_is_selectable() -> None:
    """Expose all GPT-5.6 tiers with their exact API identifiers."""
    expected = {"gpt-5.6-sol", "gpt-5.6-terra", "gpt-5.6-luna"}

    assert expected.issubset(AVAILABLE_MODELS["openai"])
    assert expected.issubset(
        {item["name"] for item in selectable_models() if item["provider"] == "openai"}
    )


def test_openai_recommendations_preserve_workload_tiers() -> None:
    """Map capability, balance, and efficiency workloads to distinct tiers."""
    recommendations = RECOMMENDED_MODELS["openai"]

    assert recommendations["reasoning"] == "gpt-5.6-sol"
    assert recommendations["advanced"] == "gpt-5.6-sol"
    assert recommendations["balanced"] == "gpt-5.6-terra"
    assert recommendations["fast"] == "gpt-5.6-luna"


def test_default_and_provider_inference_use_current_catalog() -> None:
    """Keep the CLI default current and infer every GPT-5.6 tier as OpenAI."""
    assert DEFAULT_MODEL == {"name": "gpt-5.6-terra", "provider": "openai"}
    for model in ("gpt-5.6-sol", "gpt-5.6-terra", "gpt-5.6-luna"):
        assert infer_provider(model) == "openai"
