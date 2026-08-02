"""Single-renderer regressions for graph-derived reaction labels."""

from reactive_taxonomy import featurize_reaction


CLICK = (
    "C#Cc1ccc(Cl)cc1.[N-]=[N+]=NC1COc2c(Br)cc([N+](=O)[O-])c3cccc1c23"
    ">>O=[N+]([O-])c1cc(Br)c2c3c(cccc13)C(n1cc(-c3ccc(Cl)cc3)nn1)CO2"
)


def test_reaction_label_is_derived_from_normalized_edits() -> None:
    result = featurize_reaction("Brc1ccccc1.CN>>CNc1ccccc1")
    assert result.reaction_label is not None
    assert result.reaction_label.text == "Ar–Br + Alk–NH₂ → Ar–NH–Alk"
    assert result.reaction_label.basis == "reaction_sites"
    assert set(result.to_dict()["reaction_label"]) == {
        "text",
        "status",
        "basis",
        "evidence",
        "confidence",
        "warnings",
        "style",
        "definition_version",
        "event_count",
        "schema_version",
    }


def test_cycloaddition_receives_chemist_friendly_topology_rendering() -> None:
    result = featurize_reaction(CLICK)
    assert result.reaction_label is not None
    assert "ring" in result.reaction_label.text.lower()
    assert "C₂N₃" in result.reaction_label.text
    assert result.reaction_label.basis == "ring_topology"
    assert result.interpretation is not None
    assert result.interpretation.primary_pattern_id == "cycloaddition_like"


def test_ascii_style_changes_display_not_structural_identity() -> None:
    reaction = "CCBr>>C=C"
    unicode = featurize_reaction(reaction)
    ascii_result = featurize_reaction(reaction, label_style="ascii")
    assert unicode.reaction_label is not None
    assert ascii_result.reaction_label is not None
    assert unicode.reaction_label.text != ascii_result.reaction_label.text
    assert unicode.reaction_signature is not None
    assert ascii_result.reaction_signature is not None
    assert (
        unicode.reaction_signature.signature_id
        == ascii_result.reaction_signature.signature_id
    )


def test_unresolved_reaction_is_explicitly_unavailable() -> None:
    result = featurize_reaction("CC.O>>N#N")
    assert result.reaction_label is not None
    assert result.reaction_label.status == "unavailable"
    assert result.reaction_label.text == "Unavailable"
