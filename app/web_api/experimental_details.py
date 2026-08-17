"""Resolve observed experimental details from recommendation artifacts."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Dict, Mapping


EXPERIMENTAL_DETAIL_CATALOG_FILENAME = "experimental_detail_catalog.jsonl.gz"


def load_experimental_detail_catalog(
    index_path: str | Path,
) -> Dict[str, Dict[str, Any]]:
    """Load experimental records indexed by observation and reaction identity."""

    catalog_path = Path(index_path).parent / EXPERIMENTAL_DETAIL_CATALOG_FILENAME
    if not catalog_path.is_file():
        return {}

    records: Dict[str, Dict[str, Any]] = {}
    with gzip.open(catalog_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, Mapping):
                continue
            record = dict(value)
            observation_id = str(record.get("observation_id") or "")
            reaction_id = str(record.get("reaction_id") or "")
            if observation_id:
                records[f"observation:{observation_id}"] = record
            if reaction_id:
                records.setdefault(f"reaction:{reaction_id}", record)
    return records


def attach_recommendation_experimental_details(
    payload: Dict[str, Any],
    catalog: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Any]:
    """Attach available observed procedures in displayed precedent order."""

    for recommendation in payload.get("recommendations") or ():
        details = []
        seen = set()
        for reaction_id in recommendation.get("precedent_reaction_ids") or ():
            detail = catalog.get(f"reaction:{reaction_id}")
            if detail is None:
                continue
            observation_id = str(detail.get("observation_id") or reaction_id)
            if observation_id in seen:
                continue
            seen.add(observation_id)
            details.append(dict(detail))
        recommendation["precedent_experimental_details"] = details
    return payload


__all__ = [
    "EXPERIMENTAL_DETAIL_CATALOG_FILENAME",
    "attach_recommendation_experimental_details",
    "load_experimental_detail_catalog",
]
