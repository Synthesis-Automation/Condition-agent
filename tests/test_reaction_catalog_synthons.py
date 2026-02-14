from __future__ import annotations

import json
from pathlib import Path

from chemtools.taxonomy import reaction_catalog


def test_load_reaction_catalog_normalizes_synthon_slots(tmp_path: Path) -> None:
    payload = {
        "reaction_types": [
            {
                "id": "Toy_Coupling",
                "name": "Toy Coupling",
                "category": "test",
                "reactants": {
                    "electrophile": "@sp2_electrophiles",
                },
                "products": {},
                "constraints": {},
                "synthons": {
                    "electrophile": {"synthon_set": "electrophiles"},
                    "nucleophile": {
                        "include": ["amine_nucleophile"],
                        "min_hits": 1,
                        "min_reactants": 1,
                    },
                },
            }
        ],
        "motif_sets": {
            "sp2_electrophiles": {"members": ["Ar-Br"]},
        },
    }
    path = tmp_path / "reaction_types.test.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    reaction_catalog.load_reaction_catalog.cache_clear()
    definitions, _aliases = reaction_catalog.load_reaction_catalog(path)
    toy = definitions["Toy_Coupling"]

    assert "electrophile" in toy.synthons
    assert "nucleophile" in toy.synthons
    assert "sp2_electrophile" in toy.synthons["electrophile"].allowed
    assert toy.synthons["nucleophile"].allowed == ["amine_nucleophile"]
    assert toy.redox_neutral is None
