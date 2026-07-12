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


def test_batch_autodetects_clean_reaction_column(tmp_path, capsys) -> None:
    source = tmp_path / "reactions.csv"
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rxn_smiles_clean"])
        writer.writeheader()
        writer.writerow({"rxn_smiles_clean": reaction})

    assert main(["batch", str(source), "--mode", "reaction", "--format", "json"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["valid"] == 1
    assert summary["evidence_counts"] == {"exact_product_reconstruction": 1}
