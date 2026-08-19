"""Optional narration-model configuration kept separate from chemistry logic."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict


CONFIG_PATH = Path.home() / ".chemcoworker" / "config.json"
DEFAULT_MODEL: Dict[str, str] = {"name": "gpt-5.6-terra", "provider": "openai"}
DEFAULT_REVIEW_CONFIG: Dict[str, Any] = {
    "config_version": 3,
    **DEFAULT_MODEL,
    "review_mode": "auto",
    "reasoning_effort": "medium",
    "review_candidates": 10,
    "review_max_tokens": 8_000,
    "apply_review_order": True,
}


def load_config() -> Dict[str, Any]:
    """Load validated model/review settings, including older two-field files."""

    config = DEFAULT_REVIEW_CONFIG.copy()
    try:
        data = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        if isinstance(data, dict):
            config_version = data.get("config_version")
            if isinstance(data.get("name"), str) and data["name"].strip():
                config["name"] = data["name"]
            if data.get("provider") in {"openai", "aliyun"}:
                config["provider"] = data["provider"]
            if data.get("review_mode") in {"off", "auto", "always"}:
                config["review_mode"] = data["review_mode"]
            if data.get("reasoning_effort") in {
                "none",
                "low",
                "medium",
                "high",
                "xhigh",
                "max",
            }:
                config["reasoning_effort"] = data["reasoning_effort"]
            candidates = data.get("review_candidates")
            if isinstance(candidates, int) and 1 <= candidates <= 10:
                config["review_candidates"] = (
                    candidates if config_version in {2, 3} else max(10, candidates)
                )
            max_tokens = data.get("review_max_tokens")
            if isinstance(max_tokens, int) and 256 <= max_tokens <= 20_000:
                config["review_max_tokens"] = (
                    max_tokens if config_version == 3 else max(8_000, max_tokens)
                )
            if isinstance(data.get("apply_review_order"), bool):
                config["apply_review_order"] = data["apply_review_order"]
    except (OSError, ValueError, TypeError):
        pass
    return config


def save_config(
    model: str,
    provider: str,
    *,
    review_mode: str = "auto",
    reasoning_effort: str = "medium",
    review_candidates: int = 10,
    review_max_tokens: int = 8_000,
    apply_review_order: bool = True,
) -> None:
    """Persist interactive LLM review settings."""

    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CONFIG_PATH.write_text(
        json.dumps(
            {
                "config_version": 3,
                "name": model,
                "provider": provider,
                "review_mode": review_mode,
                "reasoning_effort": reasoning_effort,
                "review_candidates": review_candidates,
                "review_max_tokens": review_max_tokens,
                "apply_review_order": apply_review_order,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
