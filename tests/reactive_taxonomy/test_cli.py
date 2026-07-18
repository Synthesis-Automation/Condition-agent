import csv
import json

from reactive_taxonomy.cli import main


def test_validate_and_self_test(capsys) -> None:
    assert main(["validate", "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)["valid"] is True
    assert main(["self-test"]) == 0
    assert "overall: PASS" in capsys.readouterr().out


def test_molecule_and_reaction_json(capsys) -> None:
    assert main(["molecule", "Brc1ccccc1", "--format", "json"]) == 0
    molecule = json.loads(capsys.readouterr().out)
    assert molecule["sites"][0]["chemist_label"] == "Ar–Br"

    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    assert main(["reaction", reaction, "--format", "json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["evidence_quality"] == "exact_product_reconstruction"
    assert payload["reaction_label"] == "Ar–Br + Ar–B(OH)2 → Ar–Ar"


def test_concise_molecule_and_reaction_output(capsys) -> None:
    assert main(["molecule", "Brc1ccccc1", "--concise"]) == 0
    molecule_output = capsys.readouterr().out
    assert "Reactive sites:" in molecule_output
    assert "Ar–Br — leaving_group, available" in molecule_output
    assert "atom0" not in molecule_output

    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    assert main(["reaction", reaction, "--concise"]) == 0
    reaction_output = capsys.readouterr().out
    assert "Reaction: Ar–Br + Ar–B(OH)2 → Ar–Ar" in reaction_output
    assert "Evidence: exact_product_reconstruction" in reaction_output
    assert "Product connection: Ar–Ar (C_C)" in reaction_output
    assert "selected grammar" not in reaction_output

    assert main(["molecule", "CC(C)(C)N", "--concise"]) == 0
    tert_butylamine_output = capsys.readouterr().out
    assert "attached Alkyl: tertiary" in tert_butylamine_output


def test_batch_autodetects_column_and_writes_jsonl(tmp_path, capsys) -> None:
    source = tmp_path / "molecules.csv"
    output = tmp_path / "results.jsonl"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "smiles"])
        writer.writeheader()
        writer.writerow({"name": "bromobenzene", "smiles": "Brc1ccccc1"})
        writer.writerow({"name": "aniline", "smiles": "Nc1ccccc1"})

    assert main(["batch", str(source), "--mode", "molecule", "--output", str(output)]) == 0
    assert "valid: 2" in capsys.readouterr().out
    records = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]
    assert len(records) == 2
    assert records[0]["source"]["name"] == "bromobenzene"


def test_batch_writes_molecule_summary_csv(tmp_path, capsys) -> None:
    source = tmp_path / "molecules.csv"
    output = tmp_path / "results.csv"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "smiles", "source_id"])
        writer.writeheader()
        writer.writerow({"name": "acetaldehyde", "smiles": "CC=O", "source_id": "mol-1"})
        writer.writerow({"name": "styrene", "smiles": "C=Cc1ccccc1", "source_id": "mol-2"})

    assert main(["batch", str(source), "--mode", "molecule", "--output", str(output)]) == 0
    assert "CSV output:" in capsys.readouterr().out
    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        smiles_index = reader.fieldnames.index("smiles")
        assert reader.fieldnames[smiles_index + 1:smiles_index + 3] == [
            "reactive_site_labels",
            "functional_group_labels",
        ]
        rows = list(reader)
    assert len(rows) == 2
    assert rows[0]["electrophilic_center_count"] == "1"
    assert rows[1]["unsaturated_bond_count"] == "1"
    assert "analysis_json" not in rows[1]


def test_batch_autodetects_clean_reaction_column(tmp_path, capsys) -> None:
    source = tmp_path / "reactions.csv"
    output = tmp_path / "results.csv"
    reaction = (
        "Brc1ccc(C#N)cc1.OB(O)c1ccccc1"
        ">>N#Cc1ccc(-c2ccccc2)cc1"
    )
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rxn_smiles_clean", "source_id"])
        writer.writeheader()
        writer.writerow({"rxn_smiles_clean": reaction, "source_id": "rxn-1"})

    assert main([
        "batch", str(source), "--mode", "reaction", "--format", "json",
        "--output", str(output),
    ]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["valid"] == 1
    assert summary["evidence_counts"] == {"exact_product_reconstruction": 1}
    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        smiles_index = reader.fieldnames.index("rxn_smiles_clean")
        assert reader.fieldnames[smiles_index + 1:smiles_index + 3] == [
            "reaction_label",
            "spectator_groups",
        ]
        rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["spectator_groups"] == "nitrile"
