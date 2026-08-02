"""Canonical notation and crosswalk consistency regressions."""

import json

from reactive_taxonomy import analyze_molecule, validate_taxonomy
from reactive_taxonomy.notation import (
    load_chemist_notation,
    render_fragment_notation,
)
from reactive_taxonomy.source_labels import validate_source_label_mappings


def test_canonical_fragment_symbols_have_one_meaning() -> None:
    payload = load_chemist_notation()
    symbols = {
        record["id"]: record["symbol"]
        for record in payload["fragment_notations"]
    }

    assert symbols["generic_organic"] == "R"
    assert symbols["alkyl"] == "Alk"
    assert symbols["aryl"] == "Ar"
    assert symbols["heteroaryl"] == "HetAr"
    assert symbols["halogen"] == "X"
    assert "HeteroAr" not in json.dumps(payload, sort_keys=True)


def test_notation_has_no_legacy_heteroaryl_alias() -> None:
    assert render_fragment_notation("heteroaryl") == "HetAr"
    assert validate_taxonomy() == []


def test_known_alkyl_sites_do_not_render_as_generic_r() -> None:
    labels = {
        site.chemist_label
        for site in analyze_molecule("CCCBr").reactive_site_hypotheses
    }

    assert "Alk–Br" in labels
    assert "R–Br" not in labels
    assert validate_source_label_mappings() == []
