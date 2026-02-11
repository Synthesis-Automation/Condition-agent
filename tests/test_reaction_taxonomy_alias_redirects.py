from __future__ import annotations

from chemtools.taxonomy.reaction_catalog import resolve_reaction_type


def test_alcohol_halide_aliases_resolve_to_new_family() -> None:
    assert resolve_reaction_type("Appel_halogenation") == "Alcohol_to_Alkyl_Halide"
    assert resolve_reaction_type("Appel halogenation") == "Alcohol_to_Alkyl_Halide"
    assert (
        resolve_reaction_type("Chlorination_SOCl2_oxalyl_chloride")
        == "Alcohol_to_Alkyl_Halide"
    )


def test_aliphatic_halide_exchange_alias_maps_to_specific_family() -> None:
    assert (
        resolve_reaction_type("Aliphatic_Halide_Exchange")
        == "Aliphatic_Halide_Exchange"
    )
    assert (
        resolve_reaction_type("Aliphatic Halide Exchange")
        == "Aliphatic_Halide_Exchange"
    )


def test_balz_schiemann_alias_redirects_to_sandmeyer() -> None:
    assert resolve_reaction_type("BalzSchiemann") == "Sandmeyer"
