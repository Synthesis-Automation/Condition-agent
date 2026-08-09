"""Core-to-SMARTS compilation, application, and comparison regressions."""

from __future__ import annotations

from dataclasses import asdict

import pytest

from reactive_taxonomy import featurize_reaction
from retrosynthesis_poc import build_library as build_baseline_library
from retrosynthesis_poc import save_library as save_baseline_library
from core_retrosynthesis_poc import (
    CoreTemplateLibrary,
    build_library,
    compile_core_templates,
    disconnect_target,
    disconnect_ensemble,
    load_library,
    run_comparison,
    save_library,
    split_by_reference,
)
from core_retrosynthesis_poc.cli import main


MAPPED_REACTIONS = {
    "C-N": (
        "[c:1]1(Br)[cH:3][cH:4][cH:5][cH:6][cH:7]1."
        "[NH2:8][CH3:9]>>"
        "[c:1]1([NH:8][CH3:9])[cH:3][cH:4][cH:5][cH:6][cH:7]1"
    ),
    "C-O": (
        "[c:1]1(Br)[cH:3][cH:4][cH:5][cH:6][cH:7]1."
        "[OH:8][CH3:9]>>"
        "[c:1]1([O:8][CH3:9])[cH:3][cH:4][cH:5][cH:6][cH:7]1"
    ),
    "C-S": (
        "[c:1]1(Br)[cH:3][cH:4][cH:5][cH:6][cH:7]1."
        "[SH:8][CH3:9]>>"
        "[c:1]1([S:8][CH3:9])[cH:3][cH:4][cH:5][cH:6][cH:7]1"
    ),
}


def _row(bond_kind: str, *, ordinal: int = 1) -> dict:
    reaction_smiles = MAPPED_REACTIONS[bond_kind]
    analysis = featurize_reaction(reaction_smiles)
    value = analysis.to_dict()
    assert analysis.reaction_core is not None
    return {
        "reaction_id": f"reaction-{bond_kind}-{ordinal}",
        "observation_id": f"observation-{bond_kind}-{ordinal}",
        "reference_id": f"reference-{bond_kind}-{ordinal}",
        "reaction_smiles": reaction_smiles,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
        "reaction_completeness": value["reaction_completeness"],
    }


def _c_n_chloride_row() -> dict:
    reaction_smiles = MAPPED_REACTIONS["C-N"].replace("(Br)", "(Cl)")
    analysis = featurize_reaction(reaction_smiles)
    value = analysis.to_dict()
    return {
        "reaction_id": "reaction-C-N-chloride",
        "observation_id": "observation-C-N-chloride",
        "reference_id": "reference-C-N-chloride",
        "reaction_smiles": reaction_smiles,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
        "reaction_completeness": value["reaction_completeness"],
    }


@pytest.mark.parametrize("bond_kind", ("C-N", "C-O", "C-S"))
def test_compiles_l1_l2_and_round_trips_source(bond_kind: str) -> None:
    result = compile_core_templates(_row(bond_kind))

    assert result.rejection_reason is None
    assert {template.abstraction_level for template in result.templates} == {
        "L1",
        "L2",
    }
    assert {template.bond_kind for template in result.templates} == {bond_kind}
    assert all(template.template_id.startswith("CRT1:") for template in result.templates)
    assert all(template.operator_id.startswith("CRO1:") for template in result.templates)
    assert all(">>" in template.reaction_smarts for template in result.templates)


def test_l1_core_is_map_and_partner_order_invariant() -> None:
    first = _row("C-N")
    reaction = MAPPED_REACTIONS["C-N"]
    reactants, product = reaction.split(">>")
    reversed_reaction = ".".join(reversed(reactants.split("."))) + ">>" + product
    analysis = featurize_reaction(reversed_reaction)
    value = analysis.to_dict()
    second = {
        **first,
        "reaction_smiles": reversed_reaction,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
    }

    first_l1 = compile_core_templates(first, levels=("L1",)).templates[0]
    second_l1 = compile_core_templates(second, levels=("L1",)).templates[0]

    assert first_l1.reaction_smarts == second_l1.reaction_smarts
    assert first_l1.operator_id == second_l1.operator_id


def test_compiled_smarts_is_invariant_to_source_atom_map_numbers() -> None:
    first = _row("C-N")
    replacements = {1: 101, 3: 17, 4: 88, 5: 42, 6: 73, 7: 9, 8: 55, 9: 2}
    remapped_reaction = MAPPED_REACTIONS["C-N"]
    for old, new in replacements.items():
        remapped_reaction = remapped_reaction.replace(f":{old}]", f":x{new}]")
    remapped_reaction = remapped_reaction.replace(":x", ":")
    analysis = featurize_reaction(remapped_reaction)
    value = analysis.to_dict()
    second = {
        **first,
        "reaction_smiles": remapped_reaction,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
    }

    first_templates = compile_core_templates(first).templates
    second_templates = compile_core_templates(second).templates

    assert [template.reaction_smarts for template in first_templates] == [
        template.reaction_smarts for template in second_templates
    ]
    assert [template.operator_id for template in first_templates] == [
        template.operator_id for template in second_templates
    ]


def test_context_is_separate_from_smarts_and_contains_site_profiles() -> None:
    template = compile_core_templates(_row("C-N"), levels=("L1",)).templates[0]
    context = template.precedents[0].context

    assert context.center_profiles
    assert {profile.role for profile in context.center_profiles} == {
        "carbon",
        "heteroatom",
    }
    assert "accessibility" not in template.reaction_smarts
    assert "activation" not in template.reaction_smarts


def test_core_library_and_search_recover_cx_precursors() -> None:
    library, report = build_library(
        _row(kind) for kind in ("C-N", "C-O", "C-S")
    )

    assert report.accepted_observation_count == 3
    assert report.unique_template_count == 6
    examples = {
        "CNc1ccccc1": "Brc1ccccc1.CN",
        "COc1ccc(F)cc1": "CO.Fc1ccc(Br)cc1",
        "CSc1ccccc1": "Brc1ccccc1.CS",
    }
    for target, expected in examples.items():
        candidates = disconnect_target(target, library, top_k=5)
        assert candidates
        assert expected in {candidate.precursor_smiles for candidate in candidates}
        assert all(0.0 <= candidate.context_similarity <= 1.0 for candidate in candidates)


def test_library_gzip_round_trip(tmp_path) -> None:
    library, _ = build_library((_row("C-N"), _row("C-O")))
    destination = tmp_path / "core.json.gz"

    save_library(library, destination)
    loaded = load_library(destination)

    assert loaded == library
    assert CoreTemplateLibrary.from_dict(asdict(library)) == library


def test_reference_split_has_no_leakage() -> None:
    rows = tuple(_row("C-N", ordinal=index) for index in range(1, 30))

    train, test = split_by_reference(rows, test_fraction=0.3)

    train_references = {row["reference_id"] for row in train}
    test_references = {row["reference_id"] for row in test}
    assert train and test
    assert train_references.isdisjoint(test_references)


def test_comparison_builds_both_methods_without_reference_leakage(tmp_path) -> None:
    rows = tuple(_row("C-N", ordinal=index) for index in range(1, 8))

    report = run_comparison(rows, tmp_path, test_fraction=0.3, top_k=5)

    assert report["split"]["reference_leakage_count"] == 0
    assert report["split"]["evaluated_targets"] >= 1
    assert report["builds"]["baseline"]["accepted_observations"] >= 1
    assert report["builds"]["core"]["accepted_observations"] >= 1
    assert (tmp_path / "baseline_templates.json.gz").is_file()
    assert (tmp_path / "core_templates.json.gz").is_file()
    assert (tmp_path / "comparison.json").is_file()


def test_concise_cli_honors_top_k(tmp_path, capsys) -> None:
    library, _ = build_library((_row("C-N"), _c_n_chloride_row()))
    destination = tmp_path / "core.json.gz"
    save_library(library, destination)

    exit_code = main(
        [
            "disconnect",
            str(destination),
            "CNc1ccccc1",
            "--concise",
            "--top-k",
            "2",
        ]
    )

    assert exit_code == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2
    assert all(">>" in line for line in lines)


def test_baseline_first_ensemble_preserves_order_and_adds_fallback(
    tmp_path,
    capsys,
) -> None:
    row = _row("C-N")
    baseline_library, _ = build_baseline_library((row,))
    core_library, _ = build_library((row, _c_n_chloride_row()))

    candidates = disconnect_ensemble(
        "CNc1ccccc1",
        baseline_library,
        core_library,
        top_k=2,
    )

    assert len(candidates) == 2
    assert candidates[0].source_method == "rdchiral_baseline"
    assert candidates[0].candidate.precursor_smiles == "Brc1ccccc1.CN"
    assert candidates[1].source_method == "core_l1_context"
    assert len({candidate.candidate.precursor_smiles for candidate in candidates}) == 2

    baseline_path = tmp_path / "baseline.json.gz"
    core_path = tmp_path / "core.json.gz"
    save_baseline_library(baseline_library, baseline_path)
    save_library(core_library, core_path)
    exit_code = main(
        [
            "disconnect-ensemble",
            str(baseline_path),
            str(core_path),
            "CNc1ccccc1",
            "--concise",
            "--top-k",
            "2",
        ]
    )
    assert exit_code == 0
    assert len(capsys.readouterr().out.splitlines()) == 2


def test_rejects_multi_event_core() -> None:
    row = _row("C-S")
    row["reaction_core"] = {**row["reaction_core"], "event_count": 2}

    result = compile_core_templates(row)

    assert result.templates == ()
    assert result.rejection_reason == "not_single_event"


def test_connected_precursor_only_handle_round_trips() -> None:
    reaction = (
        "O=Cc1cccc(O)c1."
        "Brc1ccc([I+]c2ccc(Br)cc2Br)c(Br)c1.[Cl-]"
        ">>O=Cc1cccc(Oc2ccc(Br)cc2Br)c1"
    )
    analysis = featurize_reaction(reaction)
    value = analysis.to_dict()
    row = {
        "reaction_id": "diaryliodonium-c-o",
        "observation_id": "diaryliodonium-c-o",
        "reference_id": "diaryliodonium-c-o",
        "reaction_smiles": reaction,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
        "reaction_completeness": value["reaction_completeness"],
    }

    result = compile_core_templates(row, levels=("L1",))

    assert result.rejection_reason is None
    assert len(result.templates) == 1
    assert "[#53+]" in result.templates[0].precursor_smarts
