"""Reusable canonical molecule-index regressions."""

from __future__ import annotations

import csv

from cas_tools import (
    CanonicalMoleculeIndex,
    build_canonical_molecule_index,
    molecule_identity,
)
from cas_tools.molecule_index_cli import main


def _catalog(path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("cas_no", "compound_smiles", "reaction_id"),
        )
        writer.writeheader()
        writer.writerows(
            (
                {
                    "cas_no": "64-17-5",
                    "compound_smiles": "OCC",
                    "reaction_id": "reaction-1",
                },
                {
                    "cas_no": "64-17-5",
                    "compound_smiles": "CCO",
                    "reaction_id": "reaction-2",
                },
                {
                    "cas_no": "67-56-1",
                    "compound_smiles": "[CH3:7][OH:8]",
                    "reaction_id": "reaction-3",
                },
                {
                    "cas_no": "invalid",
                    "compound_smiles": "not-smiles",
                    "reaction_id": "reaction-4",
                },
            )
        )


def test_build_and_lookup_normalizes_serialization_and_atom_maps(tmp_path) -> None:
    source = tmp_path / "catalog.csv"
    output = tmp_path / "molecules.sqlite"
    _catalog(source)

    report = build_canonical_molecule_index(
        source,
        output,
        provenance_columns=("cas_no", "reaction_id"),
    )

    assert report.source_rows == 4
    assert report.accepted_rows == 3
    assert report.invalid_smiles_rows == 1
    assert report.unique_molecules == 2
    with CanonicalMoleculeIndex(output) as index:
        ethanol = index.lookup("C(C)O")
        methanol = index.lookup(molecule_identity("CO"))

    assert ethanol is not None
    assert ethanol.canonical_smiles == "CCO"
    assert ethanol.occurrence_count == 2
    assert ethanol.source_records == (
        {"cas_no": "64-17-5", "reaction_id": "reaction-1"},
        {"cas_no": "64-17-5", "reaction_id": "reaction-2"},
    )
    assert methanol is not None
    assert methanol.canonical_smiles == "CO"


def test_module_cli_builds_configurable_index(tmp_path, capsys) -> None:
    source = tmp_path / "catalog.csv"
    output = tmp_path / "molecules.sqlite"
    _catalog(source)

    status = main(
        (
            str(source),
            str(output),
            "--provenance-column",
            "cas_no",
        )
    )

    assert status == 0
    assert output.is_file()
    assert '"unique_molecules": 2' in capsys.readouterr().out

