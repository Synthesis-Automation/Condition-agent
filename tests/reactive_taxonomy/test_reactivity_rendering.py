import json

from reactive_taxonomy import analyze_molecule
from reactive_taxonomy.descriptors import (
    reactivity_profile_tokens,
    render_reactivity_profile,
    render_reactivity_profile_expanded,
)


def _profile(smiles: str, site_type: str):
    result = analyze_molecule(smiles)
    assert result.valid
    site = next(site for site in result.reactive_site_hypotheses if site.site_type == site_type)
    environment = next(
        item
        for item in result.reactive_site_environments
        if item.hypothesis_id == site.hypothesis_id
    )
    assert environment.reactivity_profile is not None
    return environment.reactivity_profile


def test_compact_renderer_prioritizes_chemist_facing_factors() -> None:
    aromatic = _profile(
        "Cc1c(Br)c(C)ccc1",
        "leaving_group",
    )
    text = render_reactivity_profile(aromatic)

    assert "benzene (6-membered)" in text
    assert "ortho burden high (2/2)" in text
    assert "electron demand" in text
    assert "score=" not in text
    assert "atom" not in text
    assert "heavy" not in text
    assert "method" not in text


def test_heteroatom_renderer_exposes_attached_group_class() -> None:
    profile = _profile("CC(C)(C)N", "pronucleophile_XH")
    text = render_reactivity_profile(profile)
    assert "primary N" in text
    assert "attached Alkyl: tertiary" in text
    assert "lone-pair availability high" in text


def test_expanded_renderer_retains_scores_methods_and_contributors() -> None:
    profile = _profile("Brc1ccccn1", "leaving_group")
    text = render_reactivity_profile_expanded(profile)
    assert "steric score=" in text
    assert "electronic score=" in text
    assert "aromatic_steric_graph_v1" in text
    assert "contributors=" in text


def test_signature_tokens_exclude_raw_scores_indices_and_display_labels() -> None:
    profile = _profile("Brc1ccccn1", "leaving_group")
    tokens = reactivity_profile_tokens(profile)
    serialized = json.dumps(tokens)

    assert "context:aromatic" in tokens
    assert "aromatic_family:pyridine" in tokens
    assert "accessibility_score" not in serialized
    assert "activation_score" not in serialized
    assert "atom_index" not in serialized
    assert "Brc1ccccn1" not in serialized


def test_profile_tokens_are_invariant_to_equivalent_smiles_ordering() -> None:
    first = _profile("Brc1ccccn1", "leaving_group")
    second = _profile("n1ccccc1Br", "leaving_group")
    assert reactivity_profile_tokens(first) == reactivity_profile_tokens(second)


def test_profile_serialization_is_deterministic() -> None:
    profile = _profile("CC(C)Br", "leaving_group")
    first = json.dumps(profile.to_dict(), sort_keys=True)
    second = json.dumps(profile.to_dict(), sort_keys=True)
    assert first == second
