"""Explicit, model-independent execution profiles for scientific investigations."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any


PROFILE_NAMES = ("inherit", "quick", "research")
REASONING_EFFORTS = ("none", "minimal", "low", "medium", "high", "xhigh", "max", "ultra")
WEB_SEARCH_MODES = ("disabled", "cached", "live")


@dataclass(frozen=True)
class ResearchSettings:
    """Requested execution settings, without asserting provider acceptance."""

    profile: str
    reasoning_effort: str | None
    web_search: str | None
    timeout_seconds: float
    reasoning_source: str
    web_search_source: str
    timeout_source: str

    def to_dict(self) -> dict[str, Any]:
        """Serialize requests separately from where each value originated."""
        return {
            "profile": self.profile,
            "requested": {
                "reasoning_effort": self.reasoning_effort,
                "web_search": self.web_search,
                "timeout_seconds": self.timeout_seconds,
            },
            "sources": {
                "reasoning_effort": self.reasoning_source,
                "web_search": self.web_search_source,
                "timeout_seconds": self.timeout_source,
            },
        }


def resolve_research_profile(
    profile: str = "inherit", *, reasoning_effort: str | None = None,
    web_search: str | None = None, timeout_seconds: float | None = None,
) -> ResearchSettings:
    """Apply explicit overrides to a profile; never select or replace a model.

    ``inherit`` leaves provider reasoning and search configuration untouched while
    retaining the historical 900-second application deadline. A profile requests
    settings only: model/provider support and tool access must be observed at run
    time. Unsupported combinations fail through the runtime without downgrading.
    """
    if profile not in PROFILE_NAMES:
        raise ValueError(f"profile must be one of {PROFILE_NAMES}")
    if reasoning_effort is not None and reasoning_effort not in REASONING_EFFORTS:
        raise ValueError(f"reasoning_effort must be one of {REASONING_EFFORTS}")
    if web_search is not None and web_search not in WEB_SEARCH_MODES:
        raise ValueError(f"web_search must be one of {WEB_SEARCH_MODES}")
    defaults = {
        "inherit": (None, None, 900),
        "quick": ("medium", "cached", 300),
        "research": ("high", "live", 1800),
    }
    default_reasoning, default_search, default_timeout = defaults[profile]
    timeout = default_timeout if timeout_seconds is None else timeout_seconds
    if (
        isinstance(timeout, bool) or not isinstance(timeout, (int, float))
        or not math.isfinite(timeout) or not 1 <= timeout <= 7200
    ):
        raise ValueError("timeout_seconds must be a finite number between 1 and 7200")
    return ResearchSettings(
        profile=profile,
        reasoning_effort=default_reasoning if reasoning_effort is None else reasoning_effort,
        web_search=default_search if web_search is None else web_search,
        timeout_seconds=float(timeout),
        reasoning_source=("override" if reasoning_effort is not None else
                          "inherited" if default_reasoning is None else "profile"),
        web_search_source=("override" if web_search is not None else
                           "inherited" if default_search is None else "profile"),
        timeout_source="override" if timeout_seconds is not None else "profile",
    )
