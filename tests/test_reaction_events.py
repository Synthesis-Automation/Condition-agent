from __future__ import annotations

import pytest

from chemtools.featurizers.formatters.reaction_events import (
    format_multi_event_signature,
    summarize_reaction_events,
)
from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


def _event_kinds(payload: dict) -> set[str]:
    kinds: set[str] = set()
    for ev in payload.get("events") or []:
        if isinstance(ev, dict):
            kind = ev.get("kind")
            if kind:
                kinds.add(str(kind))
    return kinds


def test_reaction_events_detect_displacement_and_cn_formation() -> None:
    summary = summarize_reaction_events(
        reaction_smiles="Ic1ccccc1.N>>Nc1ccccc1",
        bond_key="break: C-I | form: C-N",
        fallback_bond_key=None,
        reacted_motifs=["Ar-I", "Ar-NH2"],
        formed_motifs=["Ar-NR2"],
        mapping_warning=None,
    )
    kinds = _event_kinds(summary)
    assert "c_n_bond_formation" in kinds
    assert "leaving_group_displacement" in kinds
    quality = summary.get("reaction_key_quality") or {}
    assert quality.get("level") in {"high", "medium", "low"}
    assert isinstance(quality.get("score_0_1"), float)
    redox = summary.get("redox_assessment") or {}
    assert redox.get("classification") in {"redox_neutral", "net_oxidation", "net_reduction", "uncertain"}
    assert isinstance(redox.get("confidence"), float)
    profile = summary.get("transformation_profile") or {}
    assert profile.get("molecularity") in {"intramolecular", "intermolecular_or_multi_component", "unknown"}
    assert "C-N" in (profile.get("formed_bond_classes") or [])
    assert "I" in (profile.get("leaving_groups") or [])
    assert "N" in (profile.get("nucleophile_elements") or [])


def test_reaction_events_flag_amidation_without_activation() -> None:
    summary = summarize_reaction_events(
        reaction_smiles="O=C(O)c1ccccc1.NCC>>O=C(NCC)c1ccccc1",
        bond_key="form: C-N",
        fallback_bond_key=None,
        reacted_motifs=["Ar-CO2H", "Bn-NHR"],
        formed_motifs=["Ar-CONR2"],
        mapping_warning=None,
    )
    assert "amidation_without_explicit_activation_marker" in (summary.get("anomalies") or [])
    assert "amidation_like" in _event_kinds(summary)


def test_reaction_events_redox_assessment_marks_balanced_substitution_as_neutral() -> None:
    summary = summarize_reaction_events(
        reaction_smiles="Ic1ccccc1.N>>Nc1ccccc1",
        bond_key="form: C-N | break: C-I",
        fallback_bond_key=None,
        reacted_motifs=["Ar-I", "Ar-NH2"],
        formed_motifs=["Ar-NR2"],
        mapping_warning=None,
    )
    redox = summary.get("redox_assessment") or {}
    assert redox.get("classification") == "redox_neutral"


def test_reaction_events_redox_assessment_marks_c_h_formation_as_reduction() -> None:
    summary = summarize_reaction_events(
        reaction_smiles="C=C>>CC",
        bond_key="form: C-H",
        fallback_bond_key=None,
        reacted_motifs=[],
        formed_motifs=[],
        mapping_warning=None,
    )
    redox = summary.get("redox_assessment") or {}
    assert redox.get("classification") == "net_reduction"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reaction_events_for_intramolecular_annulation_example() -> None:
    rxn = (
        "CCCCc1nnn(-c2ccccc2S(=O)(=O)Nc2ccc(C)cc2)c1I"
        ">>"
        "CCCCc1nnn2c1N(c1ccc(C)cc1)S(=O)(=O)c1ccccc1-2"
    )
    result = featurize_reaction(rxn, options={"detailed": True, "confirm_coupling_products": True})
    summary = result.get("reaction_events") or {}
    kinds = _event_kinds(summary)
    assert "c_n_bond_formation" in kinds
    assert "leaving_group_displacement" in kinds
    assert "ring_closure_or_annulation" in kinds
    quality = summary.get("reaction_key_quality") or {}
    assert quality.get("level") in {"high", "medium", "low"}


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reaction_events_for_multi_transform_example_flags_anomalies() -> None:
    rxn = (
        "c1ccc2c3c([nH]c2c1)CNCC3.ClCc1ccc(OCc2ccccc2)cc1.O=C(O)c1ccc(O)cc1"
        ">>"
        "O=C(c1ccc(O)cc1)N1CCc2c(n(Cc3ccc(O)cc3)c3ccccc23)C1"
    )
    result = featurize_reaction(rxn, options={"detailed": True, "confirm_coupling_products": True})
    summary = result.get("reaction_events") or {}
    kinds = _event_kinds(summary)
    assert "c_n_bond_formation" in kinds
    assert "amidation_like" in kinds
    anomalies = set(summary.get("anomalies") or [])
    assert "amidation_without_explicit_activation_marker" in anomalies
    assert "multi_transform_record_possible" in anomalies


def test_event_signature_formatter_for_multi_event_payload() -> None:
    signature = format_multi_event_signature(
        [
            {"kind": "c_o_bond_formation", "confidence": 0.9},
            {"kind": "benzyl_o_alkylation_like", "confidence": 0.84},
            {"kind": "ester_hydrolysis_like", "confidence": 0.82},
        ]
    )
    assert signature.startswith("BzOAlk+EsterHyd")


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reaction_key_displays_multi_event_signature_for_benzyl_alkylation_hydrolysis() -> None:
    rxn = (
        "Oc1c(Br)cc(Br)c2cccnc12.CCOC(=O)c1cc2cc3c(cc2nc1CBr)OCO3"
        ">>"
        "O=C(O)c1cc2cc3c(cc2nc1COc1c(Br)cc(Br)c2cccnc12)OCO3"
    )
    result = featurize_reaction(rxn, options={"detailed": True, "confirm_coupling_products": True})
    reaction_key = str(result.get("reaction_key") or "")
    reaction_type = result.get("reaction_type")
    summary = result.get("reaction_events") or {}
    kinds = _event_kinds(summary)
    assert "benzyl_o_alkylation_like" in kinds
    assert "ester_hydrolysis_like" in kinds
    assert "events:" in reaction_key
    assert "BzOAlk" in reaction_key
    assert "EsterHyd" in reaction_key
    if isinstance(reaction_type, dict):
        rt_text = str(reaction_type.get("reaction_type") or "")
    else:
        rt_text = str(reaction_type or "")
    assert rt_text.startswith("Multi-Event:")
