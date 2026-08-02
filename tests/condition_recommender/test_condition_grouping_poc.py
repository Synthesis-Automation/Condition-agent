"""Tests for the isolated weak-label condition grouping proof of concept."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from condition_recommender.condition_grouping_poc import (
    build_condition_cores,
    load_weak_label_materials,
    parse_condition_display,
    run_condition_grouping_poc,
)


_FIELDS = (
    "source_reaction_type",
    "condition_recipe_id",
    "condition_display",
    "yield_pct",
    "temperature_c",
    "time_h",
    "condition_identity_uncertainty",
)


def _write_fixture(path: Path) -> None:
    rows = []
    recipes = (
        ("Pd(OAc)2", "XPhos", "K2CO3", "Dioxane"),
        ("Pd(OAc)2", "XPhos", "Cs2CO3", "Dioxane"),
        ("Pd2(dba)3", "SPhos", "K2CO3", "THF"),
        ("CuI", "No ligand", "K3PO4", "DMSO"),
        ("CuBr", "No ligand", "K2CO3", "DMF"),
        ("NiCl2", "BINAP", "NaOtBu", "PhMe"),
    )
    for index, (catalyst, ligand, base, solvent) in enumerate(recipes):
        ligand_value = (
            ligand
            if ligand == "No ligand"
            else f"{ligand} [ligand]"
        )
        display = (
            f"{catalyst} [catalyst]; {ligand_value}; "
            f"{base} [base]; {solvent} [solvent]"
        )
        rows.append(
            {
                "source_reaction_type": "coupling_a" if index < 3 else "coupling_b",
                "condition_recipe_id": f"RCR1:{index}",
                "condition_display": display,
                "yield_pct": str(50 + index),
                "temperature_c": "80" if index % 2 else "",
                "time_h": "12" if index % 2 else "",
                "condition_identity_uncertainty": "false",
            }
        )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def test_parse_condition_display_is_order_invariant_and_retains_absence() -> None:
    left = parse_condition_display(
        "Pd(OAc)2 [catalyst]; K2CO3 [base]; No ligand"
    )
    right = parse_condition_display(
        "No ligand; K2CO3 [base]; Pd(OAc)2 [catalyst]"
    )

    assert left.material_id == right.material_id
    assert left.declared_absences == ("ligand",)
    assert "core_absence:ligand" in left.core_feature_tokens
    assert not left.warnings


def test_solvent_and_additive_are_cross_references_not_core_identity() -> None:
    left = parse_condition_display(
        "Pd(OAc)2 [catalyst]; XPhos [ligand]; K2CO3 [base]; "
        "Dioxane [solvent]"
    )
    right = parse_condition_display(
        "Pd(OAc)2 [catalyst]; XPhos [ligand]; K2CO3 [base]; "
        "THF [solvent]; LiCl [additive]"
    )

    assert left.material_id != right.material_id
    assert left.core_id == right.core_id
    assert left.core_display == right.core_display
    assert all("solvent" not in token for token in left.core_feature_tokens)
    assert any("core_pair:catalyst_base=" in token for token in left.core_feature_tokens)


def test_no_ligand_only_changes_a_catalyst_core() -> None:
    amidation = parse_condition_display(
        "HATU [coupling_reagent]; DIPEA [base]; No ligand; DMF [solvent]"
    )
    catalyst = parse_condition_display(
        "CuI [catalyst]; K3PO4 [base]; No ligand; DMSO [solvent]"
    )

    assert "No ligand" not in amidation.core_display
    assert "core_absence:ligand" not in amidation.core_feature_tokens
    assert any(
        "core_pair:coupling_reagent_base=" in token
        for token in amidation.core_feature_tokens
    )
    assert "No ligand" in catalyst.core_display
    assert "core_absence:ligand" in catalyst.core_feature_tokens


def test_material_aggregation_excludes_operating_variants(tmp_path: Path) -> None:
    source = tmp_path / "weak.csv"
    _write_fixture(source)
    materials, audit = load_weak_label_materials(source)

    assert len(materials) == 6
    assert audit["row_count"] == 6
    assert audit["unique_recipe_id_count"] == 6
    assert audit["parse_warning_counts"] == {}
    cores, unresolved = build_condition_cores(materials)
    assert len(cores) == 6
    assert unresolved == ()


def test_grouping_poc_writes_reproducible_review_artifacts(tmp_path: Path) -> None:
    source = tmp_path / "weak.csv"
    _write_fixture(source)
    first = run_condition_grouping_poc(
        source,
        tmp_path / "first",
        cluster_count=2,
        latent_dimensions=3,
        seed=7,
        silhouette_sample_size=6,
    )
    second = run_condition_grouping_poc(
        source,
        tmp_path / "second",
        cluster_count=2,
        latent_dimensions=3,
        seed=7,
        silhouette_sample_size=6,
    )

    assert first.report["feature_policy"]["cross_referenced_not_clustered"] == (
        "solvents",
        "additives",
        "other_context_components",
        "temperature_c",
        "time_h",
    )
    assert first.report["feature_policy"]["excluded"] == (
        "source_reaction_type",
        "reactive_site_labels_and_signatures",
        "yield_pct",
        "z_score",
        "procedure_text",
    )
    assert first.report["model"]["populated_cluster_count"] == 2
    first_groups = [
        json.loads(line)
        for line in first.groups_path.read_text(encoding="utf-8").splitlines()
    ]
    second_groups = [
        json.loads(line)
        for line in second.groups_path.read_text(encoding="utf-8").splitlines()
    ]
    assert [group["condition_group_id"] for group in first_groups] == [
        group["condition_group_id"] for group in second_groups
    ]
    assert first.assignments_path.exists()
    assert first.model_path.exists()
    assert first.report_markdown_path.exists()
    assert (
        sum(first.report["quality"]["core_assignment_status_counts"].values())
        == 6
    )
    assert first_groups[0]["context_cross_references"]["solvent_systems"]
