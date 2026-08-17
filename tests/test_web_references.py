"""Reference-catalog integration tests for the local web application."""

from __future__ import annotations

import gzip
import json

from app.web_api.references import (
    attach_recommendation_references,
    load_reference_catalog,
)


def test_reference_catalog_is_loaded_beside_index(tmp_path) -> None:
    index_path = tmp_path / "generic_index.sqlite"
    index_path.touch()
    record = {
        "reference_id": "REF1:paper",
        "raw_reference": "Example Journal (2024), 1, 1-5",
        "doi": "10.1000/example",
    }
    with gzip.open(
        tmp_path / "reference_catalog.jsonl.gz", "wt", encoding="utf-8"
    ) as handle:
        handle.write(json.dumps(record) + "\n")

    assert load_reference_catalog(index_path) == {"REF1:paper": record}


def test_reference_records_are_attached_to_recommendations() -> None:
    payload = {
        "recommendations": [
            {"precedent_reference_ids": ["REF1:paper", "REF1:missing"]}
        ]
    }
    catalog = {
        "REF1:paper": {
            "reference_id": "REF1:paper",
            "raw_reference": "Example Journal (2024), 1, 1-5",
        }
    }

    enriched = attach_recommendation_references(payload, catalog)

    assert enriched["recommendations"][0]["precedent_references"] == [
        catalog["REF1:paper"]
    ]
