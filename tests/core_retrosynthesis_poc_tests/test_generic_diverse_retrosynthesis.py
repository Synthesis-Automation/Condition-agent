"""Generic graph-operator compilation and diverse-search regressions."""

from __future__ import annotations

import pytest

from reactive_taxonomy import featurize_reaction
from core_retrosynthesis_poc import (
    GenericTemplateLibrary,
    build_generic_library,
    compile_generic_templates,
    disconnect_generic_target,
)
from core_retrosynthesis_poc.generic_compiler import classify_reaction_smiles


REACTIONS = {
    "acyl_substitution": (
        "O=C(O)c1ccccc1.NCc1ccccc1>>"
        "O=C(NCc1ccccc1)c1ccccc1"
    ),
    "c_c_coupling": (
        "Brc1ccccc1.OB(O)c1ccccc1>>"
        "c1ccc(-c2ccccc2)cc1"
    ),
    "carbonyl_oxidation": "OCCc1ccccc1>>O=CCc1ccccc1",
    "carbonyl_reduction": "O=Cc1ccccc1>>OCc1ccccc1",
    "conjugate_addition": (
        "CC1C(=O)CCCC1=O.O=[N+]([O-])/C=C/c1ccccc1Cl>>"
        "CC1(C(C[N+](=O)[O-])c2ccccc2Cl)C(=O)CCCC1=O"
    ),
    "carbonyl_condensation": (
        "CCN.O=Cc1ccc(Oc2ccccc2)cc1>>"
        "CCNCc1ccc(Oc2ccccc2)cc1"
    ),
    "ring_formation": "C#CC.[N-]=[N+]=NCC>>CCn1cc(C)nn1",
}


def _row(kind: str, *, ordinal: int = 1) -> dict:
    reaction_smiles = REACTIONS[kind]
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.valid and analysis.reaction_core is not None
    value = analysis.to_dict()
    return {
        "reaction_id": f"{kind}-{ordinal}",
        "reference_id": f"reference-{kind}-{ordinal}",
        "reaction_smiles": reaction_smiles,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
        "reaction_completeness": value["reaction_completeness"],
        "_sampling_cohort": kind,
    }


@pytest.mark.parametrize("kind", tuple(REACTIONS))
@pytest.mark.parametrize("engine", ("reaction_core", "rdchiral"))
def test_compiles_source_round_tripped_diverse_archetypes(
    kind: str,
    engine: str,
) -> None:
    result = compile_generic_templates(_row(kind), engine=engine)

    assert result.rejection_reason is None
    assert result.templates
    assert {template.transformation_kind for template in result.templates} == {kind}
    assert all(template.reaction_smarts.count(">>") == 1 for template in result.templates)


def test_generic_reduction_generalizes_and_hard_filters_archetype() -> None:
    library = build_generic_library(
        (_row("carbonyl_reduction"),),
        engine="reaction_core",
    )

    candidates = disconnect_generic_target(
        "OCc1ccc(F)cc1",
        library,
        transformations=("carbonyl_reduction",),
        top_k=5,
    )
    incompatible = disconnect_generic_target(
        "COc1ccc(F)cc1",
        library,
        transformations=("carbonyl_reduction",),
        top_k=5,
    )

    assert candidates
    assert candidates[0].precursor_smiles == "O=Cc1ccc(F)cc1"
    assert candidates[0].forward_validation_status == "verified_signature"
    assert incompatible == ()


def test_generic_library_serialization_round_trip() -> None:
    library = build_generic_library(
        (_row("carbonyl_reduction"), _row("ring_formation")),
        engine="reaction_core",
    )

    loaded = GenericTemplateLibrary.from_dict(library.to_dict())

    assert loaded == library


def test_structural_classifier_does_not_use_family_names() -> None:
    row = _row("acyl_substitution")
    row["source_declared_family"] = "incorrect_family_label"

    result = compile_generic_templates(row, engine="reaction_core")

    assert result.templates[0].transformation_kind == "acyl_substitution"
    assert classify_reaction_smiles(REACTIONS["acyl_substitution"]) == (
        "acyl_substitution"
    )
    assert classify_reaction_smiles("CCO>>CCO") is None
