"""Lazy data/cache manager for recommendation backends."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from .models import SourceGroup


def _normalize_source_group_label(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text in {"", "any", "all"}:
        return "all"
    if text in {"literature", "dataset", "datasets", "lit"}:
        return "literature"
    if text in {"experiments", "experiment", "experiements"}:
        return "experiments"
    if text in {"protocols", "protocol"}:
        return "protocols"
    if text == "rules":
        return "rules"
    return text


def _resolve_db_path_for_source(db_path: str, source_group: str) -> str:
    normalized = _normalize_source_group_label(source_group)
    if normalized != "experiments":
        return db_path
    path = Path(db_path)
    if not path.is_dir():
        return db_path
    candidate = path / "experiments" / "HTE_canonical.csv"
    if candidate.exists():
        return str(candidate)
    return db_path


@dataclass
class LoadedResourceInfo:
    cache_key: str
    db_path: str
    source_group: str
    cache_hit: bool


class RecommendationDataManager:
    """Owns lazy recommender instances and explicit cache warming."""

    def __init__(self, base_db_path: str = "data/HTE_db") -> None:
        self.base_db_path = base_db_path
        self._recommender_cache: Dict[str, Any] = {}

    def _cache_key_for(self, db_path: str) -> str:
        try:
            return str(Path(db_path).resolve())
        except Exception:
            return str(db_path)

    def resolve_db_path(self, source_group: Optional[str | SourceGroup] = None) -> str:
        return _resolve_db_path_for_source(
            self.base_db_path,
            _normalize_source_group_label(source_group or "all"),
        )

    def get_recommender(
        self,
        *,
        source_group: Optional[str | SourceGroup] = None,
        db_path: Optional[str] = None,
    ) -> Tuple[Any, LoadedResourceInfo]:
        target_path = str(db_path or self.resolve_db_path(source_group))
        key = self._cache_key_for(target_path)
        recommender = self._recommender_cache.get(key)
        cache_hit = recommender is not None
        if recommender is None:
            from .recommender import HTERecommender

            recommender = HTERecommender(target_path)
            self._recommender_cache[key] = recommender
        return recommender, LoadedResourceInfo(
            cache_key=key,
            db_path=target_path,
            source_group=_normalize_source_group_label(source_group or "all"),
            cache_hit=cache_hit,
        )

    def warm(
        self,
        *,
        source_group: Optional[str | SourceGroup] = None,
        clear_memory_cache: bool = False,
    ) -> Dict[str, Any]:
        from .recommender import warm_hte_cache

        summary = warm_hte_cache(
            self.base_db_path,
            source_group=_normalize_source_group_label(source_group or "all"),
            clear_memory_cache=clear_memory_cache,
        )
        return summary

    def clear_local_cache(self) -> None:
        self._recommender_cache.clear()

    def cache_summary(self) -> Dict[str, Any]:
        return {
            "base_db_path": self.base_db_path,
            "cached_recommenders": len(self._recommender_cache),
            "cache_keys": list(self._recommender_cache.keys()),
        }

