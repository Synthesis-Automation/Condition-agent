"""Human review artifacts for frozen single-step search comparisons."""

from __future__ import annotations

import csv
import html
import json
from pathlib import Path
from typing import Any

from .html_report import molecule_svg, reaction_svg


def render_strategy_review(report_path: str | Path) -> Path:
    """Render an explicitly unblinded packet and an unscored review sheet.

    Automated reconstruction checks do not establish synthetic feasibility.
    Existing human scores are never overwritten when refreshing the HTML.
    """

    report_path = Path(report_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    output = report_path.parent / "chemist_review"
    output.mkdir(parents=True, exist_ok=True)
    sections = []
    review_rows = []

    def candidate_card(candidate: dict[str, Any], rank: int) -> str:
        reaction = candidate.get("proposed_reaction_smiles", "")
        return (
            f'<h4>Choice {rank}</h4>{reaction_svg(reaction)}'
            f'<code>{html.escape(reaction)}</code>'
            '<p>Supporting observations: '
            + html.escape(", ".join(candidate.get("precedent_reaction_ids", [])))
            + "</p>"
        )

    for number, case in enumerate(report["cases"], start=1):
        case_id = f"R{number:03d}"
        target = case.get("target_smiles", "")
        columns = []
        for engine in ("flat", "strategy"):
            evidence = case["engines"].get(engine, {})
            result = evidence.get("result", {})
            groups = (
                result.get("strategies", [])
                if engine == "strategy"
                else [
                    {"representative": candidate, "alternate_realizations": []}
                    for candidate in result.get("candidates", [])
                ]
            )
            cards = []
            for rank, group in enumerate(groups, start=1):
                representative = group["representative"]
                cards.append(
                    '<article>'
                    + candidate_card(representative, rank)
                    + "".join(
                        '<details><summary>Alternate precursor choice</summary>'
                        + candidate_card(alternate, index)
                        + "</details>"
                        for index, alternate in enumerate(
                            group.get("alternate_realizations", []), start=2
                        )
                    )
                    + "</article>"
                )
                review_rows.append({
                    "case_id": case_id,
                    "observation_id": case["observation_id"],
                    "engine": engine,
                    "rank": rank,
                    "strategy_id": representative.get("strategy_id", ""),
                    "chemically_plausible": "",
                    "useful_disconnection": "",
                    "distinct_useful_alternative": "",
                    "comment": "",
                })
            columns.append(
                f'<section><h3>{engine.title()}</h3>'
                + ("".join(cards) or "<p>No proposals returned.</p>")
                + "</section>"
            )
        rejection = case.get("source_rejection")
        sections.append(
            f'<section class="case" id="{case_id}"><h2>{case_id}</h2>'
            + (molecule_svg(target) if target else "")
            + f'<p><code>{html.escape(target)}</code></p>'
            + (
                f'<p>Source compilation failure: {html.escape(rejection)}</p>'
                if rejection else ""
            )
            + '<details><summary>Observed reaction (reference answer)</summary>'
            + reaction_svg(case["reaction_smiles"])
            + "</details>"
            + '<div class="columns">' + "".join(columns) + "</div></section>"
        )
    sheet = output / "review.csv"
    if not sheet.exists():
        with sheet.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=[
                "case_id", "observation_id", "engine", "rank", "strategy_id",
                "chemically_plausible", "useful_disconnection",
                "distinct_useful_alternative", "comment",
            ])
            writer.writeheader()
            writer.writerows(review_rows)
    page = output / "review_packet.html"
    page.write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<title>Single-step strategy review</title><style>'
        'body{font:16px system-ui;max-width:1500px;margin:auto;padding:24px;'
        'color:#172c38;background:#f4f7f9}h1,h2,h3{color:#123f55}'
        '.columns{display:grid;grid-template-columns:1fr 1fr;gap:20px}'
        'article,.case{background:white;padding:16px;margin:16px 0;'
        'border:1px solid #d4dce2;border-radius:8px;min-width:0}'
        'svg{max-width:100%;height:auto}code{overflow-wrap:anywhere}'
        'summary{cursor:pointer;padding:8px}nav a{display:inline-block;margin:5px}'
        '@media(max-width:800px){.columns{grid-template-columns:1fr}}'
        '</style><h1>Single-step strategy review</h1>'
        '<p>This is an <strong>unblinded diagnostic review</strong> of frozen '
        'held-out queries. It is not a completed chemist assessment. '
        'Signature verification checks graph reconstruction; it does not prove '
        'reaction feasibility or useful synthesis.</p>'
        '<p>Use <a href="review.csv">the unscored review sheet</a>. Mark each '
        'judgment yes, no, or uncertain, and explain disagreements. Assess '
        'the proposed choice before opening the observed reference reaction. '
        'Missing source compilations remain visible. Search settings were '
        'frozen before outcomes; do not tune against these cases.</p><nav>'
        + " ".join(
            f'<a href="#R{i:03d}">R{i:03d}</a>'
            for i in range(1, len(sections) + 1)
        )
        + "</nav>" + "".join(sections) + "</html>",
        encoding="utf-8",
    )
    return page
