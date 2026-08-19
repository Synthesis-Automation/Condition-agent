"""Curated model catalog for a future optional narration layer."""

from __future__ import annotations

from typing import Dict, List


_MODELS: List[Dict[str, str]] = [
    {"name": "gpt-5.6-terra", "provider": "openai"},
    {"name": "gpt-5.6-sol", "provider": "openai"},
    {"name": "gpt-5.6-luna", "provider": "openai"},
    {"name": "gpt-5.5", "provider": "openai"},
    {"name": "gpt-5.4", "provider": "openai"},
    {"name": "gpt-5.4-mini", "provider": "openai"},
    {"name": "gpt-5.4-nano", "provider": "openai"},
    {"name": "gpt-5.2", "provider": "openai"},
    {"name": "gpt-5-mini", "provider": "openai"},
    {"name": "gpt-5-nano", "provider": "openai"},
    {"name": "o3", "provider": "openai"},
    {"name": "gpt-4.1", "provider": "openai"},
    {"name": "glm-5.2", "provider": "aliyun"},
    {"name": "deepseek-v4-pro", "provider": "aliyun"},
    {"name": "deepseek-v4-flash-0731", "provider": "aliyun"},
    {"name": "qwen3.8-max", "provider": "aliyun"},
]


def selectable_models() -> List[Dict[str, str]]:
    """Return isolated copies of the curated model entries."""

    return [dict(item) for item in _MODELS]


def infer_provider(model: str) -> str:
    """Infer a provider for a known model, defaulting to OpenAI."""

    return next(
        (item["provider"] for item in _MODELS if item["name"] == model),
        "openai",
    )


def provider_model_set(provider: str) -> set[str]:
    """Return the known model names for one provider."""

    return {item["name"] for item in _MODELS if item["provider"] == provider}
