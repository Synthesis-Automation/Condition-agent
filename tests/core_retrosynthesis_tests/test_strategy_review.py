"""Review generation preserves source failures and existing human judgments."""

import json

from core_retrosynthesis import strategy_review


def test_review_is_unscored_and_does_not_replace_human_judgments(
    tmp_path, monkeypatch,
) -> None:
    monkeypatch.setattr(strategy_review, "molecule_svg", lambda *args: "<svg/>")
    monkeypatch.setattr(strategy_review, "reaction_svg", lambda *args: "<svg/>")
    report = tmp_path / "report.json"
    report.write_text(json.dumps({"cases": [{
        "observation_id": "unusable", "reaction_smiles": "invalid",
        "source_rejection": "unresolved_correspondence", "engines": {},
    }, {
        "observation_id": "usable", "reaction_smiles": "CC=O>>CCO",
        "target_smiles": "CCO", "engines": {"strategy": {"result": {
            "strategies": [{"representative": {
                "proposed_reaction_smiles": "CC=O>>CCO",
                "strategy_id": "STRAT1:example",
            }}],
        }}},
    }]}), encoding="utf-8")
    page = strategy_review.render_strategy_review(report)
    text = page.read_text(encoding="utf-8")
    assert "unblinded" in text
    assert "unresolved_correspondence" in text
    sheet = page.parent / "review.csv"
    assert "STRAT1:example,,,," in sheet.read_text(encoding="utf-8")
    sheet.write_text("preserve human judgments", encoding="utf-8")
    strategy_review.render_strategy_review(report)
    assert sheet.read_text(encoding="utf-8") == "preserve human judgments"
