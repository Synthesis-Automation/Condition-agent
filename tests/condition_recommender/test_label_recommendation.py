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
    row.update(resolve_source_label(fg_a).to_columns("reactive_site_1"))
    row.update(resolve_source_label(fg_b).to_columns("reactive_site_2"))
    return row


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_label_recommender_uses_family_and_unordered_signature_pair(tmp_path: Path) -> None:
    path = tmp_path / "labels.csv"
    _write(path, [
        _row("Suzuki-Miyaura", "ArBr", "ArB(OH)2", yield_pct=80, catalyst="exact-a"),
        _row("Suzuki-Miyaura", "ArB(OH)2", "ArBr", yield_pct=75, catalyst="exact-reversed"),
        _row("Suzuki-Miyaura", "ArCl", "ArB(OH)2", yield_pct=100, catalyst="partial"),
        _row("Buchwald-Hartwig", "ArBr", "ArB(OH)2", yield_pct=100, catalyst="wrong-family"),
    ])

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
    assert all("wrong-family" not in item.conditions.values() for item in result.recommendations)


def test_label_recommender_rewards_matching_alpha_branch_qualifier(tmp_path: Path) -> None:
    path = tmp_path / "labels.csv"
    _write(path, [
        _row("Buchwald-Hartwig", "ArBr", "RNH2 a-branch", yield_pct=80, catalyst="qualified"),
        _row("Buchwald-Hartwig", "ArBr", "RNH2", yield_pct=80, catalyst="generic"),
    ])

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


def test_label_recommender_requires_positive_top_k() -> None:
    result = recommend_conditions_from_labels("not-used", top_k=0)
    assert not result.valid
    assert result.error == "TOP_K_MUST_BE_POSITIVE"


def test_recommend_cli_returns_requested_top_k(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    path = tmp_path / "labels.csv"
    _write(path, [
        _row("Suzuki-Miyaura", "ArBr", "ArB(OH)2", yield_pct=80, catalyst="one"),
        _row("Suzuki-Miyaura", "ArBr", "ArB(OH)2", yield_pct=70, catalyst="two"),
    ])
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    monkeypatch.setattr(
        "sys.argv",
        ["recommend_cli", reaction, "--records", str(path), "--top-k", "2"],
    )

    recommend_cli_main()

    payload = json.loads(capsys.readouterr().out)
    assert payload["valid"] is True
    assert len(payload["recommendations"]) == 2
