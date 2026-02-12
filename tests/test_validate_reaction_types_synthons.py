from __future__ import annotations

from chemtools.taxonomy.validate_reaction_types import _validate_constraints


def test_validate_reaction_types_flags_unknown_synthon_slot_keys() -> None:
    payload = {
        "reaction_types": [
            {
                "id": "Toy_Unknown_Key",
                "constraints": {},
                "synthons": {
                    "electrophile": {
                        "include": ["sp2_electrophile"],
                        "bad_key": True,
                    }
                },
            }
        ]
    }
    result = _validate_constraints(payload)
    joined = "\n".join(result["issues"])
    assert "synthons.electrophile unknown keys: bad_key" in joined


def test_validate_reaction_types_flags_unknown_synthon_ids() -> None:
    payload = {
        "reaction_types": [
            {
                "id": "Toy_Unknown_Id",
                "constraints": {},
                "synthons": {
                    "nucleophile": {
                        "include": ["not_a_real_synthon"],
                    }
                },
            }
        ]
    }
    result = _validate_constraints(payload)
    joined = "\n".join(result["issues"])
    assert "unknown synthon ids: not_a_real_synthon" in joined

