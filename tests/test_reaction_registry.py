"""
Tests for the unified reaction registry (chemtools/retro/reaction_registry.py).

Validates that all cross-references between the three retrosynthesis knowledge
systems (retron patterns, HTE templates, taxonomy) are consistent and that the
registry lookup API works correctly.
"""
from __future__ import annotations

import pytest


# ---------------------------------------------------------------------------
# Import guards
# ---------------------------------------------------------------------------

def test_registry_imports():
    """Registry module imports without error."""
    from chemtools.retro.reaction_registry import (  # noqa: F401
        REACTION_REGISTRY,
        get_by_taxonomy_id,
        get_taxonomy_id_for_retron,
        get_taxonomy_id_for_template,
        get_hte_families_for_retron,
        validate_registry,
    )


# ---------------------------------------------------------------------------
# Core validation: no broken cross-references
# ---------------------------------------------------------------------------

def test_registry_validates_clean():
    """validate_registry() must return zero warnings."""
    from chemtools.retro.reaction_registry import validate_registry
    warnings = validate_registry()
    assert warnings == [], (
        f"Registry has {len(warnings)} cross-reference issue(s):\n"
        + "\n".join(f"  - {w}" for w in warnings)
    )


# ---------------------------------------------------------------------------
# Coverage: all retrons and templates have taxonomy_id
# ---------------------------------------------------------------------------

def test_all_retrons_have_taxonomy_id():
    """Every retron in RETRONS must have a 'taxonomy_id' field."""
    from chemtools.retro.retron_patterns import RETRONS
    missing = [r["name"] for r in RETRONS if not r.get("taxonomy_id")]
    assert missing == [], f"Retrons missing taxonomy_id: {missing}"


def test_all_hte_templates_have_taxonomy_id():
    """Every HTE template must have a 'taxonomy_id' field."""
    from chemtools.retro.hte_templates import HTE_TEMPLATES
    missing = [t["name"] for t in HTE_TEMPLATES if not t.get("taxonomy_id")]
    assert missing == [], f"HTE templates missing taxonomy_id: {missing}"


# ---------------------------------------------------------------------------
# Lookup API: get_by_taxonomy_id
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("taxonomy_id,expected_retron", [
    ("Suzuki_miyaura", "biaryl_suzuki"),
    ("Heck", "aryl_alkene_heck"),
    ("Reductive_amination", "alpha_amino_reductive_amination"),
    ("Ring_Closing_Metathesis", "macrocycle_rce"),
])
def test_get_by_taxonomy_id_returns_retron(taxonomy_id, expected_retron):
    from chemtools.retro.reaction_registry import get_by_taxonomy_id
    entry = get_by_taxonomy_id(taxonomy_id)
    assert entry is not None, f"No registry entry for {taxonomy_id!r}"
    assert expected_retron in entry.retron_names, (
        f"{taxonomy_id}: expected retron '{expected_retron}' in {entry.retron_names}"
    )


@pytest.mark.parametrize("taxonomy_id,expected_template", [
    ("Suzuki_miyaura", "suzuki_miyaura"),
    ("C_N_Coupling", "buchwald_hartwig"),
    ("Amide_formation", "amide_coupling"),
    ("Click_azide_alkyne_cycloaddition", "cuaac_triazole"),
])
def test_get_by_taxonomy_id_returns_template(taxonomy_id, expected_template):
    from chemtools.retro.reaction_registry import get_by_taxonomy_id
    entry = get_by_taxonomy_id(taxonomy_id)
    assert entry is not None, f"No registry entry for {taxonomy_id!r}"
    assert expected_template in entry.template_names, (
        f"{taxonomy_id}: expected template '{expected_template}' in {entry.template_names}"
    )


def test_get_by_taxonomy_id_returns_none_for_unknown():
    from chemtools.retro.reaction_registry import get_by_taxonomy_id
    assert get_by_taxonomy_id("nonexistent_xyz_abc") is None


# ---------------------------------------------------------------------------
# Lookup API: get_taxonomy_id_for_retron
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("retron_name,expected_tid", [
    ("biaryl_suzuki", "Suzuki_miyaura"),
    ("aryl_amine_buchwald", "C_N_Coupling"),
    ("amide_direct", "Amide_formation"),
    ("aryl_ether_ullmann_o", "C_O_Coupling"),
    ("williamson_ether", "Alkyl_Nucleophilic_Substitution"),
    ("mitsunobu_inversion", "Mitsunobu_reaction"),
    ("cyclohexene_diels_alder", "Diels_Alder"),
    ("alkene_sharpless", "Dihydroxylation"),
    ("boc_protected_amine", "Boc_protection"),
])
def test_get_taxonomy_id_for_retron(retron_name, expected_tid):
    from chemtools.retro.reaction_registry import get_taxonomy_id_for_retron
    result = get_taxonomy_id_for_retron(retron_name)
    assert result == expected_tid, (
        f"get_taxonomy_id_for_retron({retron_name!r}): "
        f"expected {expected_tid!r}, got {result!r}"
    )


def test_get_taxonomy_id_for_retron_returns_none_for_unknown():
    from chemtools.retro.reaction_registry import get_taxonomy_id_for_retron
    assert get_taxonomy_id_for_retron("no_such_retron") is None


# ---------------------------------------------------------------------------
# Lookup API: get_taxonomy_id_for_template
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("template_name,expected_tid", [
    ("suzuki_miyaura", "Suzuki_miyaura"),
    ("buchwald_hartwig", "C_N_Coupling"),
    ("snar_amination", "C_N_Coupling"),
    ("snar_co", "C_O_Coupling"),
    ("cuaac_triazole", "Click_azide_alkyne_cycloaddition"),
    ("wacker_oxidation", "Wacker_oxidation"),
    ("liebeskind_srogl", "Liebeskind_Srogl"),
    ("giese_radical", "Giese_radical_addition"),
])
def test_get_taxonomy_id_for_template(template_name, expected_tid):
    from chemtools.retro.reaction_registry import get_taxonomy_id_for_template
    result = get_taxonomy_id_for_template(template_name)
    assert result == expected_tid, (
        f"get_taxonomy_id_for_template({template_name!r}): "
        f"expected {expected_tid!r}, got {result!r}"
    )


# ---------------------------------------------------------------------------
# Lookup API: get_hte_families_for_retron
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("retron_name,expected_family", [
    ("aryl_amine_buchwald", "C_N_Coupling"),
    ("amide_direct", "Amide_formation"),
    ("alpha_amino_reductive_amination", "Reductive_amination"),
])
def test_get_hte_families_for_retron(retron_name, expected_family):
    from chemtools.retro.reaction_registry import get_hte_families_for_retron
    families = get_hte_families_for_retron(retron_name)
    assert expected_family in families, (
        f"get_hte_families_for_retron({retron_name!r}): "
        f"expected '{expected_family}' in {families}"
    )


# ---------------------------------------------------------------------------
# Difficulty: unified difficulty is ≤ both contributing values
# ---------------------------------------------------------------------------

def test_registry_difficulty_is_minimum():
    from chemtools.retro.reaction_registry import REACTION_REGISTRY
    from chemtools.retro.retron_patterns import RETRONS
    from chemtools.retro.hte_templates import HTE_TEMPLATES

    retron_diff = {r["name"]: r.get("difficulty", 0.5) for r in RETRONS}
    template_diff = {t["name"]: t.get("difficulty", 0.5) for t in HTE_TEMPLATES}

    for tid, entry in REACTION_REGISTRY.items():
        for rn in entry.retron_names:
            if rn in retron_diff:
                assert entry.difficulty <= retron_diff[rn] + 1e-6, (
                    f"Registry entry '{tid}' difficulty {entry.difficulty} "
                    f"> retron '{rn}' difficulty {retron_diff[rn]}"
                )
        for tn in entry.template_names:
            if tn in template_diff:
                assert entry.difficulty <= template_diff[tn] + 1e-6, (
                    f"Registry entry '{tid}' difficulty {entry.difficulty} "
                    f"> template '{tn}' difficulty {template_diff[tn]}"
                )


# ---------------------------------------------------------------------------
# Registry size sanity
# ---------------------------------------------------------------------------

def test_registry_has_expected_minimum_entries():
    """Registry should have at least 60 entries (original 60 + 29 new)."""
    from chemtools.retro.reaction_registry import REACTION_REGISTRY
    assert len(REACTION_REGISTRY) >= 60, (
        f"Registry has only {len(REACTION_REGISTRY)} entries, expected ≥ 60"
    )


def test_new_taxonomy_entries_in_registry():
    """Spot-check that the newly added taxonomy entries appear in the registry."""
    from chemtools.retro.reaction_registry import REACTION_REGISTRY
    new_ids = [
        "Aldol_addition",
        "Grignard_addition",
        "Diels_Alder",
        "Mitsunobu_reaction",
        "Boc_protection",
        "Sulfonamide_synthesis",
        "Wacker_oxidation",
        "Giese_radical_addition",
        "Thiol_ene",
    ]
    for tid in new_ids:
        assert tid in REACTION_REGISTRY, (
            f"New taxonomy entry '{tid}' is missing from the registry"
        )
