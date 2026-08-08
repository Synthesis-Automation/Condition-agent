"""Resolve user-facing reference records from recommendation artifacts."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Dict, Mapping


REFERENCE_CATALOG_FILENAME = "reference_catalog.jsonl.gz"


def load_reference_catalog(index_path: str | Path) -> Dict[str, Dict[str, Any]]:
    """Load the reference catalog stored beside a recommendation index."""

    catalog_path = Path(index_path).parent / REFERENCE_CATALOG_FILENAME
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
            reference_id = str(value.get("reference_id") or "")
            if reference_id:
                records[reference_id] = dict(value)
    return records


def attach_recommendation_references(
    payload: Dict[str, Any],
    catalog: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Any]:
    """Attach readable reference records to each recommendation result."""

    for recommendation in payload.get("recommendations") or ():
        reference_ids = recommendation.get("precedent_reference_ids") or ()
        recommendation["precedent_references"] = [
            dict(catalog[reference_id])
            for reference_id in reference_ids
            if reference_id in catalog
        ]
    return payload


def attach_discovery_references(
    payload: Dict[str, Any],
    catalog: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Any]:
    """Attach the readable source reference to each discovery hit."""

    for hit in payload.get("hits") or ():
        reference = catalog.get(str(hit.get("reference_id") or ""))
        hit["reference_record"] = dict(reference) if reference else None
    return payload


__all__ = [
    "REFERENCE_CATALOG_FILENAME",
    "attach_discovery_references",
    "attach_recommendation_references",
    "load_reference_catalog",
]
