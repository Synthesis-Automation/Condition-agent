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
from core_retrosynthesis_poc.diverse_benchmark import run_diverse_benchmark
from core_retrosynthesis_poc.generic_compiler import (
    classify_reaction_smiles,
    classify_reaction_with_site,
)


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


def test_disconnection_site_ignores_suzuki_handle_form() -> None:
    product = "Cc1ccc(-c2ccccc2)cc1"
    boronic_acid = f"Brc1ccccc1.Cc1ccc(B(O)O)cc1>>{product}"
    pinacol_chloride = (
        "CC1(C)OB(c2ccc(C)cc2)OC1(C)C.Clc1ccccc1>>"
        f"{product}"
    )

    acid_kind, acid_site = classify_reaction_with_site(boronic_acid)
    pinacol_kind, pinacol_site = classify_reaction_with_site(pinacol_chloride)

    assert acid_kind == pinacol_kind == "c_c_coupling"
    assert acid_site.startswith("SITE1:")
    assert acid_site == pinacol_site


def test_disconnection_site_is_invariant_to_mapped_source_serialization() -> None:
    template = compile_generic_templates(
        _row("conjugate_addition"),
        engine="reaction_core",
        levels=("L1",),
    ).templates[0]
    precedent = template.precedents[0]

    mapped = classify_reaction_with_site(precedent.mapped_reaction_smiles)
    canonical = classify_reaction_with_site(
        f"{precedent.precursor_smiles}>>{precedent.product_smiles}"
    )

    assert mapped == canonical


def test_site_diverse_ranking_exposes_distinct_product_bonds() -> None:
    library = build_generic_library((_row("c_c_coupling"),), engine="reaction_core")

    candidates = disconnect_generic_target(
        "Cc1ccc(-c2ccc(-c3ccc(F)cc3)cc2)cc1",
        library,
        transformations=("c_c_coupling",),
        top_k=4,
        max_candidates_to_validate=20,
        diversify_sites=True,
    )

    unique_sites = {candidate.disconnection_site_key for candidate in candidates}
    assert len(unique_sites) >= 2
    assert len({candidate.disconnection_site_key for candidate in candidates[:2]}) == 2


def test_diverse_benchmark_writes_leakage_safe_artifacts(tmp_path) -> None:
    rows = tuple(
        _row(kind, ordinal=ordinal)
        for kind in ("carbonyl_reduction", "ring_formation")
        for ordinal in range(1, 7)
    )

    report = run_diverse_benchmark(
        rows,
        tmp_path,
        test_fraction=0.3,
        max_targets_per_transformation=2,
        top_k=3,
        max_candidates_to_validate=5,
    )

    assert report["split"]["evaluated_targets"] >= 2
    assert set(report["metrics"]) == {
        "baseline",
        "core_context",
        "core_l1_context",
        "core_l2_context",
        "core_site_diverse",
        "ensemble_baseline_core_context",
        "ensemble_baseline_core_site_diverse",
    }
    assert report["metrics"]["core_context"]["valid_candidate_fraction"] == 1.0
    assert report["metrics"]["core_context"]["top1_site_recall"] > 0.0
    assert (
        report["metrics"]["core_context"]["top10_site_recall"]
        >= report["metrics"]["core_context"]["top1_site_recall"]
    )
    assert report["target_results"][0]["expected_disconnection_site_key"]
    assert report["target_results"][0]["difficulty"] in {
        "single_site",
        "multi_site",
    }
    assert report["metrics_by_difficulty"]
    assert (tmp_path / "baseline_templates.json.gz").is_file()
    assert (tmp_path / "core_templates.json.gz").is_file()
    assert (tmp_path / "comparison.json").is_file()
