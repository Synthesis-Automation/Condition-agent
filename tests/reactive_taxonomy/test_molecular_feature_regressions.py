"""Chemical counterexamples and graph-serialization regressions from review."""

from collections import Counter

import pytest
from rdkit import Chem

from reactive_taxonomy import analyze_molecule, detect_reactive_site_hypotheses
from reactive_taxonomy.descriptors import (
    reactivity_profile_tokens,
    render_reactivity_profile,
)


def _profiles(smiles: str) -> dict:
    analysis = analyze_molecule(smiles)
    assert analysis.valid
    sites = {site.hypothesis_id: site for site in analysis.reactive_site_hypotheses}
    return {
        sites[env.hypothesis_id].canonical_signature: env.reactivity_profile
        for env in analysis.reactive_site_environments
    }


def _fingerprint(smiles: str) -> Counter:
    analysis = analyze_molecule(smiles)
    environments = {
        env.hypothesis_id: env.reactivity_profile
        for env in analysis.reactive_site_environments
    }
    return Counter(
        (
            site.canonical_signature,
            reactivity_profile_tokens(environments[site.hypothesis_id]),
            environments[site.hypothesis_id].electronic.activation_score,
            environments[site.hypothesis_id].steric.accessibility_score,
        )
        for site in analysis.reactive_site_hypotheses
    )


@pytest.mark.parametrize(
    "smiles",
    [
        "CC=CCO",
        "CC#CCO",
        "CC=C(C)O",
        "COc1ccc(Br)cc1",
        "Clc1cccc[n+]1C",
    ],
)
def test_features_are_invariant_to_randomized_smiles(smiles: str) -> None:
    expected = _fingerprint(smiles)
    molecule = Chem.MolFromSmiles(smiles)
    for variant in Chem.MolToRandomSmilesVect(molecule, 20, randomSeed=1234):
        assert _fingerprint(variant) == expected, variant


@pytest.mark.parametrize(
    "smiles,allylic,propargylic",
    [
        ("BrCC=O", False, False),
        ("BrCC#N", False, False),
        ("BrCC=C", True, False),
        ("BrCC#C", False, True),
        ("BrCOC=C", False, False),
    ],
)
def test_alkyl_activation_requires_a_carbon_carbon_pi_system(
    smiles: str,
    allylic: bool,
    propargylic: bool,
) -> None:
    context = _profiles(smiles)["LG|Alkyl|Br"].context
    assert context.allylic is allylic
    assert context.propargylic is propargylic


@pytest.mark.parametrize("smiles", ["Clc1cccc[nH+]1", "Clc1cccc[n+]1C"])
def test_cationic_aromatic_nitrogen_is_not_a_pyrrole_donor(smiles: str) -> None:
    profile = _profiles(smiles)["LG|HetAr|Cl"]
    assert profile.context.heteroatoms[0].aromatic_role == "cationic_aromatic"
    assert all(item.contribution >= 0 for item in profile.electronic.contributions)


@pytest.mark.parametrize(
    "smiles,signature",
    [
        ("CC#N", "PI|Nitrile"),
        ("CC=NC", "PI|PolarizedC=N"),
    ],
)
def test_unsupported_profiles_do_not_claim_neutral_or_open(
    smiles: str,
    signature: str,
) -> None:
    profile = _profiles(smiles)[signature]
    assert profile.status == "unresolved"
    assert profile.steric.accessibility_score is None
    assert profile.electronic.activation_score is None
    assert profile.steric.evidence.confidence == 0
    assert reactivity_profile_tokens(profile) == ()
    text = render_reactivity_profile(profile)
    assert "unknown" in text and "balanced" not in text and "access open" not in text


def test_rdkit_helper_keeps_original_indices_including_disconnected_components() -> (
    None
):
    molecule = Chem.MolFromSmiles("NCC.BrCCO")
    molecule = Chem.RenumberAtoms(
        molecule, list(reversed(range(molecule.GetNumAtoms())))
    )
    sites = detect_reactive_site_hypotheses(molecule)
    for site in sites:
        roles = site.details["atom_roles"]
        if site.canonical_signature.startswith("XH|N"):
            assert molecule.GetAtomWithIdx(roles["center"][0]).GetSymbol() == "N"
        if site.site_type == "leaving_group":
            assert molecule.GetAtomWithIdx(roles["handle"][0]).GetSymbol() == "Br"
            assert (
                molecule.GetBondBetweenAtoms(roles["anchor"][0], roles["handle"][0])
                is not None
            )
        component = Chem.GetMolFrags(molecule)[site.component_index]
        assert set(site.atom_indices) <= set(component)
        assert "environment" not in site.context_features


def test_context_opt_out_skips_descriptor_construction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def forbidden(*args, **kwargs):
        raise AssertionError("descriptor construction must be skipped")

    monkeypatch.setattr("reactive_taxonomy.api.build_site_environment", forbidden)
    result = analyze_molecule("CCN", include_context_features=False)
    assert result.valid and result.reactive_site_hypotheses
    assert result.reactive_site_environments == ()
    assert all(site.context_features == {} for site in result.reactive_site_hypotheses)
    assert analyze_molecule("CCN", site_types=[]).reactive_site_hypotheses == ()


@pytest.mark.parametrize("smiles", ["[2H]OC", "[2H]N(C)C"])
def test_isotopic_hydrogen_counts_agree_with_the_detected_handle(smiles: str) -> None:
    profile = next(iter(_profiles(smiles).values()))
    assert profile.reactive_center.hydrogen_count == 1


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Brc1ccccc1", "balanced"),
        ("COc1ccc(Br)cc1", "slightly_rich"),
        ("FC(F)(F)c1ccc(Br)cc1", "slightly_poor"),
        ("O=[N+]([O-])c1ccc(Br)cc1", "electron_poor"),
    ],
)
def test_aromatic_substituent_electronic_matrix(smiles: str, expected: str) -> None:
    electronic = _profiles(smiles)["LG|Ar|Br"].electronic
    assert electronic.activation_class == expected
    assert electronic.evidence.confidence < 1


def test_aromatic_resonance_requires_direct_attachment_and_correct_position() -> None:
    para = _profiles("COc1ccc(Br)cc1")["LG|Ar|Br"].electronic
    meta = _profiles("COc1cccc(Br)c1")["LG|Ar|Br"].electronic
    remote = _profiles("COCc1ccc(Br)cc1")["LG|Ar|Br"].electronic
    assert any(
        c.pathway == "resonance" and c.positional_relation == "para"
        for c in para.contributions
    )
    assert not any(c.pathway == "resonance" for c in meta.contributions)
    assert not any(
        c.source_id == "aromatic_substituent:alkoxy" for c in remote.contributions
    )
    assert para.activation_score < meta.activation_score


def test_aromatic_rule_validation_rejects_missing_roles_and_invalid_weights() -> None:
    from copy import deepcopy
    from reactive_taxonomy.descriptors.registry import (
        descriptor_rules,
        validate_aromatic_substituent_rules,
    )

    rules = deepcopy(descriptor_rules()["electronic"]["aromatic_substituents"])
    rules["rules"][0]["smarts"] = "CO"
    with pytest.raises(ValueError, match="roles"):
        validate_aromatic_substituent_rules(rules)
    rules = deepcopy(descriptor_rules()["electronic"]["aromatic_substituents"])
    rules["rules"][0]["resonance"]["para"] = 2
    with pytest.raises(ValueError, match="weights"):
        validate_aromatic_substituent_rules(rules)
