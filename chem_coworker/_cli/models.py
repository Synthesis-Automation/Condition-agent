"""Model catalog and selection helpers."""

from __future__ import annotations

from typing import Dict, List

from .style import C


_FALLBACK_MODELS: List[Dict[str, str]] = [
    {"name": "o4-mini", "provider": "openai"},
    {"name": "gpt-5.2", "provider": "openai"},
    {"name": "glm-5", "provider": "aliyun"},
    {"name": "glm-4.7", "provider": "aliyun"},
    {"name": "MiniMax-M2.1", "provider": "aliyun"},
    {"name": "deepseek-v3.2", "provider": "aliyun"},
    {"name": "qwen3-max", "provider": "aliyun"},
]

_PREFERRED_ORDER = [m["name"] for m in _FALLBACK_MODELS]


def selectable_models() -> List[Dict[str, str]]:
    """
    Build selectable models from shared catalog (`llmtools.clients`) when available.
    Falls back to the local curated list.
    """
    try:
        from llmtools.clients import AVAILABLE_MODELS  # local shared catalog
    except Exception:
        return list(_FALLBACK_MODELS)

    all_models: List[Dict[str, str]] = []
    seen = set()

    def _add(name: str, provider: str) -> None:
        key = (name, provider)
        if key not in seen:
            seen.add(key)
            all_models.append({"name": name, "provider": provider})

    for name in _PREFERRED_ORDER:
        for provider, names in AVAILABLE_MODELS.items():
            if name in names:
                _add(name, provider)
                break

    for provider in ("openai", "aliyun"):
        for name in AVAILABLE_MODELS.get(provider, []):
            _add(name, provider)

    return all_models or list(_FALLBACK_MODELS)


def infer_provider(model: str) -> str:
    """Infer provider from shared/fallback catalog; defaults to openai."""
    for item in selectable_models():
        if item["name"] == model:
            return item["provider"]
    return "openai"


def provider_model_set(provider: str) -> set[str]:
    """Return known models for a provider."""
    return {m["name"] for m in selectable_models() if m["provider"] == provider}


def select_model_interactive() -> Dict[str, str]:
    """Show numbered model menu and return selected {name, provider}."""
    models = selectable_models()
    print(f"\n  {C.LABEL}◆ Select model{C.R}")
    print(f"  {C.SECTION}{'─' * 44}{C.R}")
    for i, m in enumerate(models, 1):
        default_tag = f"  {C.DIM}← default{C.R}" if i == 1 else ""
        provider_color = C.CYAN if m["provider"] == "openai" else C.MAGENTA
        print(
            f"  {C.DIM}{i}.{C.R}  "
            f"{C.BOLD}{m['name']:<20}{C.R}"
            f"{provider_color}{m['provider']}{C.R}"
            f"{default_tag}"
        )
    print(f"  {C.SECTION}{'─' * 44}{C.R}")
    raw = input(f"  {C.PROMPT}>{C.R} Choice (1–{len(models)}, Enter=default): ").strip()
    if not raw:
        return models[0]
    try:
        idx = int(raw) - 1
        if 0 <= idx < len(models):
            return models[idx]
    except ValueError:
        pass
    print(f"  {C.DIM}Invalid — using default ({models[0]['name']}){C.R}")
    return models[0]
