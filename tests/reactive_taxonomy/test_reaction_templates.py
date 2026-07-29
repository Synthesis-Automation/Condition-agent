"""Contracts for the expandable single-event reaction-template registry."""

from __future__ import annotations

import json

import pytest

from reactive_taxonomy import (
    ReactionTemplateError,
    derive_reaction_template,
    load_reaction_template_registry,
    match_reaction_templates,
    upsert_reaction_template,
    validate_reaction_template_registry,
)
from reactive_taxonomy.reaction_template_cli import main as template_cli_main


ACETAL_REFERENCE = (
    "[CH3:7][CH:1]=[O:2].[CH3:4][OH:3].[CH3:6][OH:5]"
    ">>"
    "[CH3:7][CH:1]([O:3][CH3:4])[O:5][CH3:6]"
)

RENORMALIZED_ACETAL_REFERENCE = (
    "[OH:31][CH3:41].[CH3:71][CH:11]=[O:21].[OH:51][CH3:61]"
    ">>"
    "[CH3:71][CH:11]([O:51][CH3:61])[O:31][CH3:41]"
)

INCOMPLETE_REPORTED_ACETAL = (
    "CCO.COc1cccc(C=O)c1>>CCOC(OCC)c1cccc(OC)c1"
)


def _acetal_template(*, status: str = "active"):
    return derive_reaction_template(
        ACETAL_REFERENCE,
        template_id="carbonyl_to_dialkoxy",
        display_name="Carbonyl to dialkoxy",
        family_id="acetalization",
        aliases=("acetal formation",),
        transformation_class="carbonyl_diheteroatom_condensation",
        status=status,  # type: ignore[arg-type]
    )


def test_mapped_reference_compiles_to_one_deterministic_event() -> None:
    template = _acetal_template()
    renormalized = derive_reaction_template(
        RENORMALIZED_ACETAL_REFERENCE,
        template_id="renormalized_acetal",
        display_name="Renormalized acetal",
    )

    assert template.edit_component_count == 1
    assert template.edit_archetype == "substitution"
    assert len(template.edits) == 5
    assert [participant.explicit_count for participant in template.participants] == [
        1,
        1,
        2,
    ]
    assert template.edit_fingerprint == renormalized.edit_fingerprint
    assert template.edit_fingerprint.startswith("RTE1:")
    assert template.definition_hash.startswith("RTD1:")


def test_reference_requires_complete_heavy_atom_mapping() -> None:
    with pytest.raises(
        ReactionTemplateError,
        match="Every reactant heavy atom",
    ):
        derive_reaction_template(
            "CC=O.[CH3:4][OH:3]>>CC([O:3][CH3:4])O",
            template_id="invalid_partial_map",
            display_name="Invalid partial map",
        )


def test_importer_rejects_disconnected_multi_event_reference() -> None:
    reaction = (
        "[CH3:1][Br:2].[NH2:3].[CH3:4][Cl:5].[OH:6]"
        ">>"
        "[CH3:1][NH2:3].[CH3:4][OH:6]"
    )
    with pytest.raises(
        ReactionTemplateError,
        match="exactly one connected edit event",
    ):
        derive_reaction_template(
            reaction,
            template_id="two_independent_substitutions",
            display_name="Two independent substitutions",
        )


def test_registry_round_trip_and_explicit_replace(tmp_path) -> None:
    path = tmp_path / "reaction_templates.v1.json"
    template = _acetal_template()
    upsert_reaction_template(template, path)

    loaded = load_reaction_template_registry(path)
    assert loaded == (template,)
    assert validate_reaction_template_registry(path) == ()
    with pytest.raises(ReactionTemplateError, match="already exists"):
        upsert_reaction_template(template, path)

    replacement = derive_reaction_template(
        ACETAL_REFERENCE,
        template_id=template.template_id,
        display_name="Acetal and ketal formation",
        family_id="acetalization",
        status="draft",
    )
    upsert_reaction_template(replacement, path, replace_existing=True)
    assert load_reaction_template_registry(path) == (replacement,)


def test_query_signature_is_derived_and_template_is_interpretive(tmp_path) -> None:
    path = tmp_path / "reaction_templates.v1.json"
    template = _acetal_template()
    upsert_reaction_template(template, path)

    result = match_reaction_templates(RENORMALIZED_ACETAL_REFERENCE, path=path)

    assert result.valid
    assert result.signature_id is not None
    assert result.signature_id.startswith("RS3:")
    assert result.edit_fingerprint == template.edit_fingerprint
    assert [match.template_id for match in result.matches] == [
        "carbonyl_to_dialkoxy"
    ]
    assert result.matches[0].family_id == "acetalization"


def test_minimum_acetal_template_matches_requested_incomplete_example(
    tmp_path,
) -> None:
    path = tmp_path / "reaction_templates.v1.json"
    template = _acetal_template()
    upsert_reaction_template(template, path)

    result = match_reaction_templates(
        INCOMPLETE_REPORTED_ACETAL,
        path=path,
    )

    assert result.valid
    assert result.signature_id is None
    assert result.evidence == "template_center_transition_hypothesis"
    assert result.edit_fingerprint == template.edit_fingerprint
    assert [match.template_id for match in result.matches] == [
        "carbonyl_to_dialkoxy"
    ]
    assert result.matches[0].provisional is True
    assert result.matches[0].confidence == 0.7
    assert result.matches[0].evidence == (
        "template_center_transition_hypothesis"
    )
    assert (
        "PROVISIONAL_TEMPLATE_MATCH_WITHOUT_ATOM_PROVENANCE"
        in result.warnings
    )
    assert "UNACCOUNTED_PRODUCT_HEAVY_ATOMS" in result.warnings


def test_acetal_template_does_not_match_hemiacetal_or_reduction(
    tmp_path,
) -> None:
    path = tmp_path / "reaction_templates.v1.json"
    upsert_reaction_template(_acetal_template(), path)

    hemiacetal = match_reaction_templates(
        "CC=O.CO>>CC(O)OC",
        path=path,
    )
    reduction = match_reaction_templates("CC=O>>CCO", path=path)

    assert hemiacetal.matches == ()
    assert reduction.matches == ()


def test_drafts_require_explicit_query_opt_in(tmp_path) -> None:
    path = tmp_path / "reaction_templates.v1.json"
    upsert_reaction_template(_acetal_template(status="draft"), path)

    assert match_reaction_templates(ACETAL_REFERENCE, path=path).matches == ()
    assert len(
        match_reaction_templates(
            ACETAL_REFERENCE,
            path=path,
            include_drafts=True,
        ).matches
    ) == 1


def test_reaction_template_cli_import_list_validate_show_and_match(
    tmp_path, capsys
) -> None:
    path = tmp_path / "reaction_templates.v1.json"
    common = ["--registry", str(path)]
    assert template_cli_main(
        [
            *common,
            "import",
            "--mapped-reaction",
            ACETAL_REFERENCE,
            "--id",
            "carbonyl_to_dialkoxy",
            "--name",
            "Carbonyl to dialkoxy",
            "--family",
            "acetalization",
            "--status",
            "active",
            "--format",
            "json",
        ]
    ) == 0
    imported = json.loads(capsys.readouterr().out)
    assert imported["saved"] is True
    assert imported["template"]["edit_component_count"] == 1

    assert template_cli_main([*common, "validate", "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)["valid"] is True

    assert template_cli_main([*common, "list", "--format", "json"]) == 0
    listed = json.loads(capsys.readouterr().out)
    assert listed["count"] == 1

    assert template_cli_main(
        [
            *common,
            "show",
            "carbonyl_to_dialkoxy",
            "--format",
            "json",
        ]
    ) == 0
    shown = json.loads(capsys.readouterr().out)
    assert shown["family_id"] == "acetalization"

    assert template_cli_main(
        [
            *common,
            "match",
            RENORMALIZED_ACETAL_REFERENCE,
            "--format",
            "json",
        ]
    ) == 0
    matched = json.loads(capsys.readouterr().out)
    assert matched["signature_id"].startswith("RS3:")
    assert matched["matches"][0]["template_id"] == "carbonyl_to_dialkoxy"

    assert template_cli_main(
        [
            *common,
            "match",
            INCOMPLETE_REPORTED_ACETAL,
            "--format",
            "json",
        ]
    ) == 0
    incomplete = json.loads(capsys.readouterr().out)
    assert incomplete["signature_id"] is None
    assert incomplete["evidence"] == "template_center_transition_hypothesis"
    assert incomplete["matches"][0]["provisional"] is True
