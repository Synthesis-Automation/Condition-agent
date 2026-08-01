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
    assert payload["reaction_label"]["concise"] == (
        "Ar1–Br + Ar2–B(OH)2 → Ar1–Ar2"
    )


def test_concise_molecule_and_reaction_output(capsys) -> None:
    assert main(["molecule", "Brc1ccccc1", "--concise"]) == 0
    molecule_output = capsys.readouterr().out
    assert "Reactive sites:" in molecule_output
    assert "Ar–Br — leaving_group, available" in molecule_output
    assert "atom0" not in molecule_output

    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    assert main(["reaction", reaction, "--concise"]) == 0
    reaction_output = capsys.readouterr().out
    assert "Reaction: Ar1–Br + Ar2–B(OH)2 → Ar1–Ar2" in reaction_output
    assert "Evidence: exact_product_reconstruction" in reaction_output
    assert "Product connection: Ar1–Ar2 (C_C)" in reaction_output
    assert (
        "Partner analysis: electrophile=benzene (6-membered); "
        "ortho burden none (0/2); electron demand balanced"
        in reaction_output
    )
    assert (
        "transfer_partner=benzene (6-membered); "
        "ortho burden none (0/2); electron demand balanced"
        in reaction_output
    )
    assert "selected grammar" not in reaction_output

    cn_reaction = "Brc1ccccn1.Nc1ccccc1>>c1ccc(Nc2ccccn2)cc1"
    assert main(["reaction", cn_reaction, "--concise"]) == 0
    cn_output = capsys.readouterr().out
    assert (
        "Partner analysis: electrophile=pyridine (6-membered), "
        "heteroatoms N@ortho" in cn_output
    )
    assert "nucleophile=primary N, aryl delocalized" in cn_output

    alkyl_reaction = "CC(C)Br.N>>CC(C)N"
    assert main(["reaction", alkyl_reaction, "--concise"]) == 0
    alkyl_output = capsys.readouterr().out
    assert (
        "Partner analysis: electrophile=secondary alkyl; "
        "access hindered (alpha-branched); center polarization balanced; "
        "beta-H present"
        in alkyl_output
    )

    assert main(["molecule", "CC(C)(C)N", "--concise"]) == 0
    tert_butylamine_output = capsys.readouterr().out
    assert "attached Alkyl: tertiary" in tert_butylamine_output


def test_concise_reaction_output_exposes_minimized_core(capsys) -> None:
    reaction = (
        "[CH3:1][OH:2].O[CH3:5]."
        "[CH:3](=[O:4])[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
        ">>[CH3:1][O:2][CH:3]([O:4][CH3:5])"
        "[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
    )

    assert main(["reaction", reaction, "--concise"]) == 0
    output = capsys.readouterr().out

    assert "Reaction minimization:" in output
    assert (
        "Minimized reaction: C(H)(Ar)(=O) + O(H)(R) "
        "→ C(H)(Ar)(O-R)2"
    ) in output
    assert "Core evidence: verified (validated_atom_mapping" in output
    assert "Core shape (retrieval): RSH2:" in output
    assert "Center transition (diagnostic only): RCS2:" in output
    assert "1 event(s); 1 primary center(s); 4 active atom(s)" in output
    assert "retained aryl [Fc1ccccc1] (1 port)" in output
    assert "departing heteroatom [O] (1 port)" in output


def test_concise_ambiguous_reaction_explains_structural_evidence(capsys) -> None:
    reaction = (
        "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
        ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
    )

    assert main(["reaction", reaction, "--concise"]) == 0
    output = capsys.readouterr().out

    assert "Evidence: ambiguous_atom_correspondence" in output
    assert "Atom accounting: 17 reactant → 14 product heavy atoms" in output
    assert "Unaccounted product atoms: none" in output
    assert "Atoms not in the main product: Cl × 1, N × 1, O × 1" in output
    assert (
        "Correspondence ambiguity: 2 distinct edit hypotheses; "
        "verified edits withheld"
        in output
    )
    assert output.count("Edit hypothesis ") == 2
    assert output.count("REH1:") == 2
    assert "global_correspondence; 4 correspondences; unverified" in output
    assert "formed: C(" in output
    assert (
        "Detected functional groups: reactants "
        "[aryl halide, ketone, primary amine, secondary amine] → "
        "product [aryl halide, pyrrolic nitrogen]"
        in output
    )
    assert "Net bond inventory (unmapped, not verified edits):" in output
    assert "Grammar checks: 9 candidates (9 product mismatch)" in output
    assert "Rejected interpretations:" in output
    assert "reductive carbonyl–amine coupling" in output
    assert (
        "Retrieval: not eligible — contradicted or incomplete structure "
        "(candidate hypotheses, confidence 0.40)"
        in output
    )
    assert (
        "registered grammars did not reconstruct the supplied product; "
        "this does not invalidate the product structure"
        in output
    )


def test_batch_autodetects_column_and_writes_jsonl(tmp_path, capsys) -> None:
    source = tmp_path / "molecules.csv"
    output = tmp_path / "results.jsonl"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "smiles"])
        writer.writeheader()
        writer.writerow({"name": "bromobenzene", "smiles": "Brc1ccccc1"})
        writer.writerow({"name": "aniline", "smiles": "Nc1ccccc1"})

    assert (
        main(["batch", str(source), "--mode", "molecule", "--output", str(output)]) == 0
    )
    assert "valid: 2" in capsys.readouterr().out
    records = [
        json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()
    ]
    assert len(records) == 2
    assert records[0]["source"]["name"] == "bromobenzene"


def test_batch_writes_molecule_summary_csv(tmp_path, capsys) -> None:
    source = tmp_path / "molecules.csv"
    output = tmp_path / "results.csv"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "smiles", "source_id"])
        writer.writeheader()
        writer.writerow(
            {"name": "acetaldehyde", "smiles": "CC=O", "source_id": "mol-1"}
        )
        writer.writerow(
            {"name": "styrene", "smiles": "C=Cc1ccccc1", "source_id": "mol-2"}
        )

    assert (
        main(["batch", str(source), "--mode", "molecule", "--output", str(output)]) == 0
    )
    assert "CSV output:" in capsys.readouterr().out
    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        smiles_index = reader.fieldnames.index("smiles")
        assert reader.fieldnames[smiles_index + 1 : smiles_index + 3] == [
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
    reaction = "Brc1ccc(C#N)cc1.OB(O)c1ccccc1>>N#Cc1ccc(-c2ccccc2)cc1"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rxn_smiles_clean", "source_id"])
        writer.writeheader()
        writer.writerow({"rxn_smiles_clean": reaction, "source_id": "rxn-1"})

    assert (
        main(
            [
                "batch",
                str(source),
                "--mode",
                "reaction",
                "--format",
                "json",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    summary = json.loads(capsys.readouterr().out)
    assert summary["valid"] == 1
    assert summary["evidence_counts"] == {"exact_product_reconstruction": 1}
    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        smiles_index = reader.fieldnames.index("rxn_smiles_clean")
        assert reader.fieldnames[smiles_index + 1 : smiles_index + 4] == [
            "reaction_core_label",
            "reaction_display_label",
            "spectator_groups",
        ]
        rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["spectator_groups"] == "nitrile"


def test_batch_writes_concise_reaction_csv(tmp_path, capsys) -> None:
    source = tmp_path / "reactions.csv"
    output = tmp_path / "results.csv"
    reaction = "Brc1ccc(C#N)cc1.OB(O)c1ccccc1>>N#Cc1ccc(-c2ccccc2)cc1"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["reaction_smiles", "source_id", "yield_pct"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "reaction_smiles": reaction,
                "source_id": "rxn-1",
                "yield_pct": "91",
            }
        )

    assert (
        main(
            [
                "batch",
                str(source),
                "--mode",
                "reaction",
                "--concise",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    assert "CSV output:" in capsys.readouterr().out
    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames is not None
        assert reader.fieldnames[:5] == [
            "reaction_smiles",
            "reaction_core_label",
            "reaction_display_label",
            "partner_analysis",
            "spectator_groups",
        ]
        assert "reaction_core_shape_key" in reader.fieldnames
        assert "reaction_core_quality_status" in reader.fieldnames
        assert "reaction_core_state_changes" in reader.fieldnames
        assert "reaction_core_remote_subgraphs" in reader.fieldnames
        assert "reaction_display_label_detailed" in reader.fieldnames
        assert "formed_ring_sizes" in reader.fieldnames
        assert "ring_count_delta" in reader.fieldnames
        assert "ring_change_count" in reader.fieldnames
        assert "ring_changes_json" in reader.fieldnames
        rows = list(reader)

    assert len(rows) == 1
    assert rows[0]["reaction_smiles"] == reaction
    assert rows[0]["partner_analysis"] == (
        "electrophile=benzene (6-membered); ortho burden none (0/2); "
        "electron demand balanced; transfer_partner=benzene (6-membered); "
        "ortho burden none (0/2); electron demand balanced"
    )
    assert rows[0]["spectator_groups"] == "nitrile"
    assert rows[0]["reaction_core_available"] == "False"


def test_batch_reaction_csv_exposes_minimized_core(tmp_path, capsys) -> None:
    source = tmp_path / "mapped_reactions.csv"
    output = tmp_path / "mapped_results.csv"
    reaction = (
        "[CH3:1][OH:2].O[CH3:5]."
        "[CH:3](=[O:4])[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
        ">>[CH3:1][O:2][CH:3]([O:4][CH3:5])"
        "[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
    )
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reaction_smiles"])
        writer.writeheader()
        writer.writerow({"reaction_smiles": reaction})

    assert (
        main(
            [
                "batch",
                str(source),
                "--mode",
                "reaction",
                "--concise",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    assert "CSV output:" in capsys.readouterr().out
    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        row = next(csv.DictReader(handle))

    assert row["reaction_core_available"] == "True"
    assert list(row)[:3] == [
        "reaction_smiles",
        "reaction_core_label",
        "reaction_display_label",
    ]
    assert row["reaction_core_label"] == (
        "C(H)(Ar)(=O) + O(H)(R) → C(H)(Ar)(O-R)2"
    )
    assert row["reaction_core_evidence_status"] == "verified"
    assert row["reaction_core_exact_key"].startswith("RCX2:")
    assert row["reaction_core_typed_key"].startswith("RCT2:")
    assert row["reaction_core_shape_key"].startswith("RSH2:")
    assert row["reaction_core_mapping_equivalence_key"].startswith("RME1:")
    assert row["reaction_core_quality_status"] in {"pass", "review"}
    assert row["reaction_core_center_transition_key"].startswith("RCS2:")
    assert row["reaction_core_event_count"] == "1"
    assert row["reaction_core_primary_center_count"] == "1"
    assert row["reaction_core_remote_classes"] == "alkyl; aryl; heteroatom"
    assert "reactant:retained:aryl:Fc1ccccc1:1" in row[
        "reaction_core_remote_subgraphs"
    ]
