"""Unresolved molecular descriptors are absent from numeric template context."""

from core_retrosynthesis.context import _center_profiles
from reactive_taxonomy import analyze_molecule


def test_unsupported_nitrile_context_retains_observation_without_numeric_profile():
    smiles = "[CH3:1][C:2]#[N:3]"
    analysis = analyze_molecule(smiles)
    assert any(
        environment.center_atom_index == 1
        and environment.reactivity_profile.status == "unresolved"
        for environment in analysis.reactive_site_environments
    )
    observation = {
        "edits": [
            {
                "edit_type": "formed",
                "atom_1": {
                    "component_index": 0,
                    "atom_index": 1,
                    "atom_map_number": 2,
                    "element": "C",
                },
            }
        ],
    }
    assert _center_profiles(smiles + ">>", observation) == ()
