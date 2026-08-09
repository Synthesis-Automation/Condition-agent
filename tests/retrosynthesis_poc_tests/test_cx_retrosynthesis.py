"""C-X template extraction and one-step disconnection regressions."""

from __future__ import annotations

from dataclasses import asdict
import json

import pytest

from reactive_taxonomy import featurize_reaction
from retrosynthesis_poc import (
    RetrosynthesisLibrary,
    build_library,
    disconnect_target,
    load_library,
    save_library,
)
from retrosynthesis_poc.extraction import extract_cx_template
from retrosynthesis_poc.cli import main
from retrosynthesis_poc.library import iter_rows


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
    assert analysis.valid and analysis.observation is not None
    assert analysis.reaction_core is not None
    value = analysis.to_dict()
    return {
        "reaction_id": f"reaction-{bond_kind}-{ordinal}",
        "observation_id": f"observation-{bond_kind}-{ordinal}",
        "reference_id": f"reference-{ordinal}",
        "reaction_smiles": reaction_smiles,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
        "reaction_completeness": value["reaction_completeness"],
    }


def _hydrogen_coupled_pyrazole_row() -> dict:
    reaction_smiles = (
        "COc1ccccc1I.c1cn[nH]c1>>"
        "COc1ccccc1-n1cccn1"
    )
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.reaction_core is not None
    assert analysis.reaction_core.event_count == 2
    value = analysis.to_dict()
    return {
        "reaction_id": "hydrogen-coupled-pyrazole",
        "observation_id": "hydrogen-coupled-pyrazole",
        "reference_id": "hydrogen-coupled-pyrazole",
        "reaction_smiles": reaction_smiles,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
        "reaction_completeness": value["reaction_completeness"],
    }


@pytest.mark.parametrize("bond_kind", ("C-N", "C-O", "C-S"))
def test_extracts_and_round_trips_each_supported_cx_bond(bond_kind: str) -> None:
    result = extract_cx_template(_row(bond_kind))

    assert result.rejection_reason is None
    assert result.template is not None
    assert result.template.bond_kind == bond_kind
    assert result.template.template_id.startswith("CXT1:")
    assert ">>" in result.template.reaction_smarts


def test_library_deduplicates_templates_and_counts_references() -> None:
    first = _row("C-O", ordinal=1)
    second = _row("C-O", ordinal=2)
    library, report = build_library((first, second))

    assert report.accepted_observation_count == 2
    assert report.unique_template_count == 1
    assert len(library.templates) == 1
    template = library.templates[0]
    assert template.observation_support == 2
    assert template.independent_reference_support == 2
    assert len(template.precedents) == 2


def test_target_disconnection_generalizes_to_a_substituted_product() -> None:
    library, _ = build_library((_row("C-O"),))

    candidates = disconnect_target(
        "COc1ccc(F)cc1",
        library,
        allowed_bonds=("C-O",),
        top_k=5,
    )

    assert candidates
    assert candidates[0].bond_kind == "C-O"
    assert candidates[0].precursor_smiles == "CO.Fc1ccc(Br)cc1"
    assert candidates[0].forward_validation_status == "verified_signature"
    assert candidates[0].proposed_reaction_smiles.endswith(">>COc1ccc(F)cc1")


def test_search_respects_bond_filter_and_is_deterministic() -> None:
    library, _ = build_library(_row(kind) for kind in ("C-N", "C-O", "C-S"))

    first = disconnect_target("CNc1ccccc1", library, allowed_bonds=("C-N",))
    second = disconnect_target("CNc1ccccc1", library, allowed_bonds=("C-N",))

    assert first == second
    assert first
    assert {candidate.bond_kind for candidate in first} == {"C-N"}
    assert first[0].precursor_smiles == "Brc1ccccc1.CN"


def test_rejects_unsigned_or_multi_event_observations() -> None:
    row = _row("C-S")
    row["reaction_core"] = {
        **row["reaction_core"],
        "event_count": 2,
    }

    result = extract_cx_template(row)

    assert result.template is None
    assert result.rejection_reason == "not_single_event"


def test_accepts_local_hydrogen_satellite_and_generalizes_pyrazole() -> None:
    row = _hydrogen_coupled_pyrazole_row()

    library, report = build_library((row,))
    candidates = disconnect_target(
        "N#Cc1ccccc1-n1cccn1",
        library,
        allowed_bonds=("C-N",),
        top_k=5,
    )

    assert report.accepted_observation_count == 1
    assert candidates
    assert candidates[0].precursor_smiles == "N#Cc1ccccc1I.c1cn[nH]c1"


def test_rejects_nonlocal_hydrogen_satellite() -> None:
    row = _hydrogen_coupled_pyrazole_row()
    relations = row["reaction_core"]["event_relations"]
    row["reaction_core"] = {
        **row["reaction_core"],
        "event_relations": [
            {
                **relation,
                "shortest_paths": [
                    {**path, "bond_count": 2}
                    for path in relation["shortest_paths"]
                ],
            }
            for relation in relations
        ],
    }

    result = extract_cx_template(row)

    assert result.template is None
    assert result.rejection_reason == "not_single_event"


def test_unmapped_inferred_core_is_materialized_and_revalidated() -> None:
    reaction_smiles = (
        "Cc1ccc(O)cc1.COc1ccc(N)cc1>>"
        "COc1ccc(Nc2ccc(C)cc2)cc1"
    )
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.reaction_core is not None
    assert analysis.reaction_core.quality.status == "review"
    value = analysis.to_dict()
    row = {
        "reaction_id": "inferred-reaction",
        "observation_id": "inferred-observation",
        "reaction_smiles": reaction_smiles,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
    }

    result = extract_cx_template(row)

    assert result.rejection_reason is None
    assert result.template is not None
    assert result.template.bond_kind == "C-N"
    assert result.template.precedent.mapping_evidence == (
        "global_atom_correspondence"
    )
    assert result.template.precedent.mapping_confidence == 0.8
    assert ":" in result.template.precedent.mapped_reaction_smiles


def test_library_json_gzip_round_trip(tmp_path) -> None:
    library, _ = build_library((_row("C-N"), _row("C-S")))
    destination = tmp_path / "templates.json.gz"

    save_library(library, destination)
    loaded = load_library(destination)

    assert loaded == library
    assert RetrosynthesisLibrary.from_dict(asdict(library)) == library


def test_directory_row_globs_select_requested_cohorts(tmp_path) -> None:
    included = tmp_path / "C_N_Coupling.part-00000.jsonl"
    excluded = tmp_path / "Suzuki.part-00000.jsonl"
    included.write_text(json.dumps({"row": "included"}) + "\n", encoding="utf-8")
    excluded.write_text(json.dumps({"row": "excluded"}) + "\n", encoding="utf-8")

    rows = tuple(iter_rows(tmp_path, include=("C_N_Coupling*.jsonl",)))

    assert rows == ({"row": "included"},)


def test_disconnect_cli_concise_prints_requested_reaction_smiles(
    tmp_path,
    capsys,
) -> None:
    library, _ = build_library((_row("C-N"), _row("C-O")))
    destination = tmp_path / "templates.json.gz"
    save_library(library, destination)

    exit_code = main(
        [
            "disconnect",
            str(destination),
            "COc1ccc(NC)cc1",
            "--concise",
            "--top-k",
            "2",
        ]
    )

    assert exit_code == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2
    assert all(">>" in line for line in lines)
    assert all(not line.startswith("[") for line in lines)
