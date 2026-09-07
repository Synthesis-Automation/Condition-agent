"""Versioned controls for bounded route exploration and guidance."""

from dataclasses import dataclass
from functools import lru_cache
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class SearchExplorationPolicy:
    """Scheduling limits; these values never authorize chemistry."""

    definition_id: str
    schema_version: str
    widening_interval: int
    guidance_cost_band: float
    ordering_definition_id: str


@lru_cache(maxsize=1)
def load_search_exploration_policy() -> SearchExplorationPolicy:
    """Load and validate the deterministic exploration definition."""

    path = Path(__file__).with_name("definitions") / "search_exploration.v1.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("definition_id") != "search_exploration.v1":
        raise ValueError("unexpected search exploration definition")
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported search exploration schema")
    interval = value.get("widening_interval")
    if type(interval) is not int or interval < 1:
        raise ValueError("widening interval must be a positive integer")
    band = value.get("guidance_cost_band")
    if type(band) not in (int, float) or not math.isfinite(band) or band <= 0:
        raise ValueError("guidance cost band must be finite and positive")
    if value.get("ordering_definition_id") != "dependency_route_ordering.v2":
        raise ValueError("unsupported dependency ordering definition")
    return SearchExplorationPolicy(**value)
