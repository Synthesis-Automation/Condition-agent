"""Read exact reference identities from the investigation's pinned catalogue."""

from __future__ import annotations

import gzip
import json
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .operations import ScientificOperations


def reference_records(
    operations: ScientificOperations, identities: set[str],
) -> tuple[dict[str, dict[str, Any]], str]:
    """Return matching catalogue records without a network lookup or inferred join."""
    try:
        path = operations._path("reference_catalog")
    except FileNotFoundError:
        return {}, "catalog_unavailable"
    found = {}
    if path.suffix == ".sqlite":
        from condition_recommender.processed_catalog import ProcessedCatalog
        return ProcessedCatalog(path).references(identities), "catalog_available"
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                record = json.loads(line)
                if record.get("reference_id") in identities:
                    found[record["reference_id"]] = record
    return found, "catalog_available"
