"""Experimental-detail artifact integration tests for the local web app."""

from __future__ import annotations

import gzip
import json

from app.web_api.experimental_details import (
    attach_recommendation_experimental_details,
    load_experimental_detail_catalog,
)


def _detail() -> dict[str, str]:
    return {
        "observation_id": "observation:1",
        "reaction_id": "reaction:1",
        "procedure_text": "The mixture was stirred at 80 °C for 2 h.",
    }


def test_experimental_catalog_is_loaded_beside_index(tmp_path) -> None:
    index_path = tmp_path / "generic_index.sqlite"
    index_path.touch()
    with gzip.open(
        tmp_path / "experimental_detail_catalog.jsonl.gz",
        "wt",
        encoding="utf-8",
    ) as handle:
        handle.write(json.dumps(_detail()) + "\n")

    catalog = load_experimental_detail_catalog(index_path)

    assert catalog["observation:observation:1"] == _detail()
    assert catalog["reaction:reaction:1"] == _detail()


def test_experimental_details_are_attached_to_recommendations() -> None:
    payload = {"recommendations": [{"precedent_reaction_ids": ["reaction:1"]}]}
    catalog = {"reaction:reaction:1": _detail()}

    enriched = attach_recommendation_experimental_details(payload, catalog)

    assert enriched["recommendations"][0]["precedent_experimental_details"] == [
        _detail()
    ]
