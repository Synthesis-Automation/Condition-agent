import csv
import json
from pathlib import Path

from condition_recommender import recommend_conditions_from_labels
from condition_recommender.recommend_cli import main as recommend_cli_main
from reactive_taxonomy import resolve_source_label


def _row(
    reaction_type: str,
    fg_a: str,
    fg_b: str,
    *,
    yield_pct: float,
    catalyst: str,
    z_score: float = 0.0,
) -> dict[str, str]:
    row = {
        "yield_pct": str(yield_pct),
        "source_reaction_type": reaction_type,
        "z_score": str(z_score),
        "base": "K2CO3",
        "catalyst": catalyst,
        "primary_solvent": "dioxane",
        "ligand": "",
        "additive": "",
        "coupling_reagent": "",
        "secondary_solvent": "water",
        "tertiary_solvent": "",
        "procedure_text": "2 h at 80 °C",
    }
    for label, prefix in (
        (fg_a, "reactive_site_1"),
        (fg_b, "reactive_site_2"),
    ):
        columns = resolve_source_label(label).to_columns(prefix)
        row.update(
            {
                f"{prefix}_{suffix}": columns[f"{prefix}_{suffix}"]
                for suffix in (
                    "signature",
                    "center_class",
                    "attachment_class",
                    "alpha_branched",
                )
            }
        )
    return row


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_label_recommender_uses_family_and_unordered_signature_pair(
    tmp_path: Path,
) -> None:
    path = tmp_path / "labels.csv"
    _write(
        path,
        [
            _row(
                "Suzuki-Miyaura", "ArBr", "ArB(OH)2", yield_pct=80, catalyst="exact-a"
            ),
            _row(
                "Suzuki-Miyaura",
                "ArB(OH)2",
                "ArBr",
                yield_pct=75,
                catalyst="exact-reversed",
            ),
            _row(
                "Suzuki-Miyaura", "ArCl", "ArB(OH)2", yield_pct=100, catalyst="partial"
            ),
            _row(
                "Buchwald-Hartwig",
                "ArBr",
                "ArB(OH)2",
                yield_pct=100,
                catalyst="wrong-family",
            ),
        ],
    )

    result = recommend_conditions_from_labels(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        records_path=path,
        top_k=3,
    )

    assert result.valid
    assert result.grammar_id == "boron_transfer_coupling"
    assert result.candidate_count == 3
    assert len(result.recommendations) == 3
    assert result.recommendations[0].conditions["catalyst"] == "exact-a"
    assert result.recommendations[0].signature_similarity == 1.0
    assert result.recommendations[1].conditions["catalyst"] == "exact-reversed"
    assert result.recommendations[1].signature_similarity == 1.0
    assert all(
        "wrong-family" not in item.conditions.values()
        for item in result.recommendations
    )


def test_label_recommender_rewards_matching_alpha_branch_qualifier(
    tmp_path: Path,
) -> None:
    path = tmp_path / "labels.csv"
    _write(
        path,
        [
            _row(
                "Buchwald-Hartwig",
                "ArBr",
                "RNH2 a-branch",
                yield_pct=80,
                catalyst="qualified",
            ),
            _row("Buchwald-Hartwig", "ArBr", "RNH2", yield_pct=80, catalyst="generic"),
        ],
    )

    result = recommend_conditions_from_labels(
        "Brc1ccccc1.CC(C)N>>CC(C)Nc1ccccc1",
        records_path=path,
        top_k=2,
    )

    assert result.valid
    assert result.recommendations[0].conditions["catalyst"] == "qualified"
    assert (
        result.recommendations[0].qualifier_similarity
        > result.recommendations[1].qualifier_similarity
    )


def test_label_recommender_supports_alkylation_as_sp3_c_n_substitution(
    tmp_path: Path,
) -> None:
    path = tmp_path / "labels.csv"
    _write(
        path,
        [
            _row(
                "Alkylation",
                "Alkyl-Br",
                "RNH2",
                yield_pct=82,
                catalyst="sp3-cn",
            ),
            _row(
                "Buchwald-Hartwig",
                "Alkyl-Br",
                "RNH2",
                yield_pct=99,
                catalyst="wrong-source-type",
            ),
        ],
    )

    result = recommend_conditions_from_labels(
        "CCBr.CN>>CCNC",
        records_path=path,
        top_k=2,
    )

    assert result.valid
    assert result.grammar_id == "sp3_c_n_substitution"
    assert result.candidate_count == 1
    assert result.recommendations[0].signature_similarity == 1.0
    assert result.recommendations[0].conditions["catalyst"] == "sp3-cn"


def test_label_recommender_supports_heck_with_alkene_label(
    tmp_path: Path,
) -> None:
    path = tmp_path / "labels.csv"
    _write(
        path,
        [
            _row("Heck", "ArBr", "alkene", yield_pct=78, catalyst="heck"),
            _row(
                "Suzuki-Miyaura",
                "ArBr",
                "alkene",
                yield_pct=99,
                catalyst="wrong-source-type",
            ),
        ],
    )

    result = recommend_conditions_from_labels(
        "Brc1ccccc1.C=C>>C=Cc1ccccc1",
        records_path=path,
        top_k=2,
    )

    assert result.valid
    assert result.grammar_id == "terminal_alkene_heck_coupling"
    assert result.candidate_count == 1
    assert result.recommendations[0].signature_similarity == 1.0
    assert result.recommendations[0].conditions["catalyst"] == "heck"


def test_label_recommender_supports_acid_or_carboxylate_amide_label(
    tmp_path: Path,
) -> None:
    path = tmp_path / "labels.csv"
    _write(
        path,
        [
            _row(
                "Amide coupling",
                "RCO2H or M",
                "RNH2",
                yield_pct=84,
                catalyst="amide",
            ),
        ],
    )

    result = recommend_conditions_from_labels(
        "CC(=O)O.CN>>CNC(C)=O",
        records_path=path,
        top_k=1,
    )

    assert result.valid
    assert result.grammar_id == "amide_formation"
    assert result.candidate_count == 1
    assert result.recommendations[0].signature_similarity == 1.0
    assert result.recommendations[0].conditions["catalyst"] == "amide"


def test_label_recommender_supports_activated_carbon_arylation(
    tmp_path: Path,
) -> None:
    path = tmp_path / "labels.csv"
    _write(
        path,
        [
            _row(
                "Arylation, acidic C-H",
                "ArBr",
                "Alkyl-H acidic",
                yield_pct=84,
                catalyst="activated-carbon",
            ),
            _row(
                "Arylation, acidic C-H",
                "ArBr",
                "RNH2",
                yield_pct=99,
                catalyst="incompatible-nitrogen",
            ),
        ],
    )

    result = recommend_conditions_from_labels(
        "Brc1ccccc1.CC#N>>N#CCc1ccccc1",
        records_path=path,
        top_k=2,
    )

    assert result.valid
    assert result.grammar_id == "sp2_c_activated_c_substitution"
    assert result.candidate_count == 2
    assert len(result.recommendations) == 1
    assert result.recommendations[0].conditions["catalyst"] == (
        "activated-carbon"
    )
    assert 0.0 < result.recommendations[0].signature_similarity < 1.0


def test_label_recommender_supports_aromatic_ch_arylation_without_alkyne_leakage(
    tmp_path: Path,
) -> None:
    path = tmp_path / "labels.csv"
    _write(
        path,
        [
            _row(
                "CH-Activation",
                "ArBr",
                "ArH",
                yield_pct=81,
                catalyst="direct-arylation",
            ),
            _row(
                "CH-Activation",
                "ArH",
                "alkyne",
                yield_pct=99,
                catalyst="incompatible-alkyne",
            ),
        ],
    )

    result = recommend_conditions_from_labels(
        "Brc1ccccc1.c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        records_path=path,
        top_k=2,
    )

    assert result.valid
    assert result.grammar_id == "sp2_c_aromatic_ch_substitution"
    assert result.candidate_count == 2
    assert len(result.recommendations) == 1
    assert result.recommendations[0].signature_similarity == 1.0
    assert result.recommendations[0].conditions["catalyst"] == (
        "direct-arylation"
    )


def test_label_recommender_requires_positive_top_k() -> None:
    result = recommend_conditions_from_labels("not-used", top_k=0)
    assert not result.valid
    assert result.error == "TOP_K_MUST_BE_POSITIVE"


def test_label_recommender_rejects_unrepresented_intramolecular_topology() -> None:
    result = recommend_conditions_from_labels("NCCc1ccccc1Br>>c1ccc2c(c1)CCN2")

    assert not result.valid
    assert result.grammar_id == "sp2_c_n_substitution"
    assert result.error == "QUERY_TOPOLOGY_NOT_SUPPORTED_BY_LABEL_DATASET"
    assert result.warnings == ("LABEL_DATASET_HAS_NO_REACTION_TOPOLOGY",)


def test_recommend_cli_returns_requested_top_k(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    path = tmp_path / "labels.csv"
    _write(
        path,
        [
            _row("Suzuki-Miyaura", "ArBr", "ArB(OH)2", yield_pct=80, catalyst="one"),
            _row("Suzuki-Miyaura", "ArBr", "ArB(OH)2", yield_pct=70, catalyst="two"),
        ],
    )
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    monkeypatch.setattr(
        "sys.argv",
        ["recommend_cli", reaction, "--records", str(path), "--top-k", "2"],
    )

    recommend_cli_main()

    payload = json.loads(capsys.readouterr().out)
    assert payload["valid"] is True
    assert len(payload["recommendations"]) == 2
