"""Contracts for the expandable single-event reaction-template registry."""

from __future__ import annotations

import json
from dataclasses import replace

import pytest

from reactive_taxonomy import (
    ReactionTemplateError,
    derive_reaction_template,
    load_reaction_template_registry,
    match_reaction_templates,
    upsert_reaction_template,
    validate_reaction_template,
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

DARZENS_REFERENCE = (
    "[Cl:1][CH2:2][C:3](=[O:4])[O:5][CH2:6][CH3:7]."
    "[CH3:8][CH:9]=[O:10]"
    ">>"
    "[CH:2]1([C:3](=[O:4])[O:5][CH2:6][CH3:7])"
    "[CH:9]([CH3:8])[O:10]1"
)

DARZENS_QUERY = "CCOC(=O)CCl.CC=O>>CCOC(=O)C1OC1C"


def _acetal_template(*, status: str = "active"):
    return derive_reaction_template(
        ACETAL_REFERENCE,
        template_id="carbonyl_to_dialkoxy",
        display_name="Carbonyl to dialkoxy",
        family_id="acetalization",
        aliases=("acetal formation",),
        template_label="Acetal formation",
        product_label="acetal",
        transformation_class="carbonyl_diheteroatom_condensation",
        status=status,  # type: ignore[arg-type]
    )


def _darzens_template(*, status: str = "active"):
    return derive_reaction_template(
        DARZENS_REFERENCE,
        template_id="darzens_epoxide_formation",
        display_name="Darzens epoxide formation",
        family_id="darzens_reaction",
        aliases=("Darzens condensation",),
        template_label="Darzens reaction",
        product_label="glycidic ester",
        role_labels={"activated_sp3_carbon": "α-halo ester"},
        role_required_tokens={
            "activated_sp3_carbon": ("alpha_to:ester",)
        },
        atom_element_alternatives={1: ("Cl", "Br", "I")},
        transformation_class="carbonyl_epoxide_condensation",
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
    assert [
        (
            role.role_id,
            role.site_type,
            role.atom_map_numbers,
        )
        for role in template.roles
    ] == [
        ("carbonyl", "electrophilic_center", (1,)),
        ("alcohol", "pronucleophile_XH", (3, 5)),
    ]


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


def test_template_requires_compact_taxonomy_role_annotations() -> None:
    template = _acetal_template()

    errors = validate_reaction_template(replace(template, roles=()))

    assert "carbonyl_to_dialkoxy:missing_roles" in errors
    assert (
        "carbonyl_to_dialkoxy:reference_contract_mismatch:roles"
        in errors
    )


def test_darzens_reference_derives_connected_roles_and_curated_halides() -> None:
    template = _darzens_template()

    assert template.edit_component_count == 1
    assert template.edit_archetype == "addition"
    assert len(template.edits) == 5
    assert [
        (
            role.role_id,
            role.site_type,
            role.atom_map_numbers,
            role.display_label,
            role.required_context_tokens,
        )
        for role in template.roles
    ] == [
        (
            "activated_sp3_carbon",
            "pronucleophile_XH",
            (2,),
            "α-halo ester",
            ("alpha_to:ester",),
        ),
        ("carbonyl", "electrophilic_center", (9,), None, ()),
    ]
    assert [
        (item.atom_map_number, item.elements)
        for item in template.atom_element_alternatives
    ] == [(1, ("Br", "Cl", "I"))]


def test_element_alternatives_must_include_the_reference_element() -> None:
    with pytest.raises(
        ReactionTemplateError,
        match="must include reference element Cl",
    ):
        derive_reaction_template(
            DARZENS_REFERENCE,
            template_id="invalid_darzens_halides",
            display_name="Invalid Darzens halides",
            atom_element_alternatives={1: ("Br", "I")},
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
    path = tmp_path / "reaction_templates.v2.json"
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
    path = tmp_path / "reaction_templates.v2.json"
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
    path = tmp_path / "reaction_templates.v2.json"
    template = _acetal_template()
    upsert_reaction_template(template, path)

    result = match_reaction_templates(
        INCOMPLETE_REPORTED_ACETAL,
        path=path,
    )

    assert result.valid
    assert result.signature_id is None
    assert result.evidence == (
        "exact_template_reconstruction_with_inferred_multiplicity"
    )
    assert result.edit_fingerprint == template.edit_fingerprint
    assert [match.template_id for match in result.matches] == [
        "carbonyl_to_dialkoxy"
    ]
    assert result.matches[0].provisional is True
    assert result.matches[0].confidence == 0.85
    assert result.matches[0].evidence == (
        "exact_template_reconstruction_with_inferred_multiplicity"
    )
    assert result.matches[0].predicted_product_smiles == (
        "CCOC(OCC)c1cccc(OC)c1"
    )
    assert result.matches[0].inferred_multiplicity is True
    interpretation = result.matches[0].interpretation
    assert interpretation is not None
    assert interpretation.template_label == "Acetal formation"
    assert interpretation.structural_label == (
        "R–CH=O + 2 × R–OH → acetal"
    )
    assert interpretation.predicted_product_smiles == (
        "CCOC(OCC)c1cccc(OC)c1"
    )
    assert [
        (
            binding.role_id,
            binding.site_type,
            binding.multiplicity,
            binding.inferred_multiplicity,
        )
        for binding in interpretation.roles
    ] == [
        ("carbonyl", "electrophilic_center", 1, False),
        ("alcohol", "pronucleophile_XH", 2, True),
    ]
    carbonyl = interpretation.roles[0]
    alcohol = interpretation.roles[1]
    assert carbonyl.steric_class == "moderate"
    assert carbonyl.steric_score == pytest.approx(0.278)
    assert carbonyl.electronic_class == "slightly_poor"
    assert alcohol.steric_class == "open"
    assert alcohol.steric_score == pytest.approx(0.1)
    assert alcohol.electronic_class == "medium"
    assert (
        "EXACT_MAIN_PRODUCT_RECONSTRUCTION_FROM_TEMPLATE"
        in result.warnings
    )
    assert "INFERRED_REACTANT_MULTIPLICITY" in result.warnings
    assert "UNACCOUNTED_PRODUCT_HEAVY_ATOMS" in result.warnings


def test_explicit_second_alcohol_upgrades_reconstruction_evidence(
    tmp_path,
) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    template = _acetal_template()
    upsert_reaction_template(template, path)

    result = match_reaction_templates(
        "CCO.CCO.COc1cccc(C=O)c1>>CCOC(OCC)c1cccc(OC)c1",
        path=path,
    )

    assert result.signature_id is None
    assert result.evidence == "exact_template_reconstruction"
    assert result.edit_fingerprint == template.edit_fingerprint
    assert result.matches[0].confidence == 0.95
    assert result.matches[0].provisional is False
    assert result.matches[0].inferred_multiplicity is False
    interpretation = result.matches[0].interpretation
    assert interpretation is not None
    assert [
        binding.component_index
        for binding in interpretation.roles
        if binding.role_id == "alcohol"
    ] == [0, 1]
    assert all(
        not binding.inferred_multiplicity
        for binding in interpretation.roles
    )
    assert "INFERRED_REACTANT_MULTIPLICITY" not in result.warnings


def test_same_template_reconstructs_a_ketal_with_repeated_methanol(
    tmp_path,
) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    upsert_reaction_template(_acetal_template(), path)

    result = match_reaction_templates(
        "CO.CC(=O)C>>COC(C)(C)OC",
        path=path,
    )

    assert result.evidence == (
        "exact_template_reconstruction_with_inferred_multiplicity"
    )
    assert result.matches[0].family_id == "acetalization"
    assert result.matches[0].confidence == 0.85


def test_center_match_does_not_claim_exact_reconstruction_for_wrong_alcohol(
    tmp_path,
) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    upsert_reaction_template(_acetal_template(), path)

    result = match_reaction_templates(
        "CO.COc1cccc(C=O)c1>>CCOC(OCC)c1cccc(OC)c1",
        path=path,
    )

    assert result.evidence == "template_center_transition_hypothesis"
    assert result.matches[0].confidence == 0.7
    assert result.matches[0].interpretation is None
    assert (
        "EXACT_MAIN_PRODUCT_RECONSTRUCTION_FROM_TEMPLATE"
        not in result.warnings
    )


def test_equivalent_template_families_preserve_interpretation_ambiguity(
    tmp_path,
) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    first = _acetal_template()
    second = derive_reaction_template(
        ACETAL_REFERENCE,
        template_id="alternative_dialkoxy_interpretation",
        display_name="Alternative dialkoxy interpretation",
        family_id="alternative_acetal_family",
        status="active",
    )
    upsert_reaction_template(first, path)
    upsert_reaction_template(second, path)

    result = match_reaction_templates(
        INCOMPLETE_REPORTED_ACETAL,
        path=path,
    )

    assert result.edit_fingerprint == first.edit_fingerprint
    assert result.signature_id is None
    assert {match.family_id for match in result.matches} == {
        "acetalization",
        "alternative_acetal_family",
    }
    assert all(
        match.evidence
        == "exact_template_reconstruction_with_inferred_multiplicity"
        for match in result.matches
    )


def test_acetal_template_does_not_match_hemiacetal_or_reduction(
    tmp_path,
) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    upsert_reaction_template(_acetal_template(), path)

    hemiacetal = match_reaction_templates(
        "CC=O.CO>>CC(O)OC",
        path=path,
    )
    reduction = match_reaction_templates("CC=O>>CCO", path=path)

    assert hemiacetal.matches == ()
    assert reduction.matches == ()


def test_darzens_template_labels_and_profiles_exact_chloro_query(
    tmp_path,
) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    template = _darzens_template()
    upsert_reaction_template(template, path)

    result = match_reaction_templates(DARZENS_QUERY, path=path)

    assert result.signature_id is not None
    assert result.edit_fingerprint == template.edit_fingerprint
    assert [match.template_id for match in result.matches] == [
        "darzens_epoxide_formation"
    ]
    match = result.matches[0]
    assert match.evidence == "query_derived_edit_fingerprint"
    assert match.confidence == 1.0
    assert match.interpretation is not None
    assert match.interpretation.template_label == "Darzens reaction"
    assert match.interpretation.structural_label == (
        "α-halo ester + R–CH=O → glycidic ester"
    )
    assert match.interpretation.predicted_product_smiles == (
        "CCOC(=O)C1OC1C"
    )
    assert [
        (
            binding.role_id,
            binding.steric_class,
            binding.electronic_class,
        )
        for binding in match.interpretation.roles
    ] == [
        ("activated_sp3_carbon", "open", "slightly_poor"),
        ("carbonyl", "open", "slightly_poor"),
    ]


@pytest.mark.parametrize(
    ("reaction", "carbonyl_label"),
    (
        (
            "COC(=O)CBr.O=Cc1ccccc1"
            ">>COC(=O)C1OC1c1ccccc1",
            "R–CH=O",
        ),
        (
            "CCOC(=O)CI.CC(C)=O"
            ">>CCOC(=O)C1OC1(C)C",
            "R2C=O",
        ),
    ),
)
def test_darzens_curated_halide_and_carbonyl_generalization(
    tmp_path,
    reaction: str,
    carbonyl_label: str,
) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    upsert_reaction_template(_darzens_template(), path)

    result = match_reaction_templates(reaction, path=path)

    assert [match.template_id for match in result.matches] == [
        "darzens_epoxide_formation"
    ]
    match = result.matches[0]
    assert match.evidence == "exact_template_reconstruction"
    assert match.confidence == 0.95
    assert match.interpretation is not None
    assert match.interpretation.structural_label == (
        f"α-halo ester + {carbonyl_label} → glycidic ester"
    )
    assert "EXACT_MAIN_PRODUCT_RECONSTRUCTION_FROM_TEMPLATE" in result.warnings


def test_darzens_match_is_reactant_order_invariant(tmp_path) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    upsert_reaction_template(_darzens_template(), path)

    forward = match_reaction_templates(DARZENS_QUERY, path=path)
    reversed_order = match_reaction_templates(
        "CC=O.CCOC(=O)CCl>>CCOC(=O)C1OC1C",
        path=path,
    )

    assert forward.matches[0].interpretation is not None
    assert reversed_order.matches[0].interpretation is not None
    assert (
        forward.matches[0].interpretation.structural_label
        == reversed_order.matches[0].interpretation.structural_label
    )


def test_darzens_rejects_missing_activation_and_wrong_product(tmp_path) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    upsert_reaction_template(_darzens_template(), path)

    missing_activation = match_reaction_templates(
        "CCl.CC=O>>CC1OC1C",
        path=path,
    )
    wrong_product = match_reaction_templates(
        "CCOC(=O)CCl.CC=O>>CCOC(=O)CC(O)C",
        path=path,
    )
    wrong_activation_family = match_reaction_templates(
        "N#CCCl.CC=O>>N#CC1OC1C",
        path=path,
    )

    assert missing_activation.matches == ()
    assert wrong_product.matches == ()
    assert wrong_activation_family.matches == ()


def test_drafts_require_explicit_query_opt_in(tmp_path) -> None:
    path = tmp_path / "reaction_templates.v2.json"
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
    path = tmp_path / "reaction_templates.v2.json"
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
            "--template-label",
            "Acetal formation",
            "--product-label",
            "acetal",
            "--status",
            "active",
            "--format",
            "json",
        ]
    ) == 0
    imported = json.loads(capsys.readouterr().out)
    assert imported["saved"] is True
    assert imported["template"]["edit_component_count"] == 1
    assert imported["template"]["roles"][0]["role_id"] == "carbonyl"
    assert imported["template"]["template_label"] == "Acetal formation"

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
    assert incomplete["evidence"] == (
        "exact_template_reconstruction_with_inferred_multiplicity"
    )
    assert incomplete["matches"][0]["provisional"] is True
    assert incomplete["matches"][0]["interpretation"]["structural_label"] == (
        "R–CH=O + 2 × R–OH → acetal"
    )


def test_reaction_template_cli_imports_generalized_darzens_template(
    tmp_path,
    capsys,
) -> None:
    path = tmp_path / "reaction_templates.v2.json"
    common = ["--registry", str(path)]

    assert template_cli_main(
        [
            *common,
            "import",
            "--mapped-reaction",
            DARZENS_REFERENCE,
            "--id",
            "darzens_epoxide_formation",
            "--name",
            "Darzens epoxide formation",
            "--family",
            "darzens_reaction",
            "--template-label",
            "Darzens reaction",
            "--product-label",
            "glycidic ester",
            "--role-label",
            "activated_sp3_carbon=α-halo ester",
            "--role-tokens",
            "activated_sp3_carbon=alpha_to:ester",
            "--atom-elements",
            "1=Cl,Br,I",
            "--status",
            "active",
            "--format",
            "json",
        ]
    ) == 0
    imported = json.loads(capsys.readouterr().out)
    assert imported["template"]["roles"][0]["display_label"] == (
        "α-halo ester"
    )
    assert imported["template"]["roles"][0][
        "required_context_tokens"
    ] == ["alpha_to:ester"]
    assert imported["template"]["atom_element_alternatives"] == [
        {
            "atom_map_number": 1,
            "elements": ["Br", "Cl", "I"],
        }
    ]

    bromo_query = (
        "COC(=O)CBr.O=Cc1ccccc1"
        ">>COC(=O)C1OC1c1ccccc1"
    )
    assert template_cli_main(
        [
            *common,
            "match",
            bromo_query,
            "--format",
            "json",
        ]
    ) == 0
    matched = json.loads(capsys.readouterr().out)
    assert matched["matches"][0]["template_id"] == (
        "darzens_epoxide_formation"
    )
    assert matched["matches"][0]["evidence"] == (
        "exact_template_reconstruction"
    )
