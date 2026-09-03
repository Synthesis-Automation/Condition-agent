"""Self-contained structure-rich review artifacts for external route admission."""

from __future__ import annotations

from html import escape
import json
from pathlib import Path
from typing import Any

from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

from .external_route_admission import ExternalRouteAssessment


EXTERNAL_PROPOSAL_REVIEW_VERSION = "external_proposal_review.v1"


def _reaction_svg(reaction_smiles: str) -> str:
    try:
        reaction = rdChemReactions.ReactionFromSmarts(
            reaction_smiles, useSmiles=True
        )
    except Exception:
        reaction = None
    if reaction is None:
        return "<div class='missing'>Reaction depiction unavailable</div>"
    drawer = rdMolDraw2D.MolDraw2DSVG(980, 230)
    drawer.drawOptions().clearBackground = False
    drawer.DrawReaction(reaction)
    drawer.FinishDrawing()
    value = drawer.GetDrawingText()
    start = value.find("<svg")
    return value[start:] if start >= 0 else value


def _badge(label: str, value: str, tone: str) -> str:
    return (
        f"<span class='badge {escape(tone)}'><small>{escape(label)}</small>"
        f"<strong>{escape(value)}</strong></span>"
    )


def _gate_rows(gates: tuple) -> str:
    return "".join(
        "<tr>"
        f"<td><code>{escape(gate.gate_id)}</code></td>"
        f"<td><span class='status {escape(gate.status)}'>{escape(gate.status)}</span></td>"
        f"<td>{escape(gate.summary)}</td>"
        "</tr>"
        for gate in gates
    )


def render_external_route_assessment_html(
    assessment: ExternalRouteAssessment,
    *,
    title: str = "External route admission review",
) -> str:
    """Render route topology, structures, gates, and precedents for review."""

    topology = _gate_rows(assessment.topology_gates)
    step_cards = []
    for step in assessment.step_assessments:
        item = step.assessment
        reaction = (
            item.normalized_mapped_reaction_smiles
            or (
                f"{item.canonical_precursor_smiles}>>{item.canonical_target_smiles}"
                if item.canonical_precursor_smiles and item.canonical_target_smiles
                else ""
            )
        )
        precedent_cards = "".join(
            "<article class='precedent'>"
            f"{_reaction_svg(match.mapped_reaction_smiles)}"
            f"<p><strong>{escape(match.reference_id or 'unidentified reference')}</strong>"
            f" · product {match.product_similarity:.2f}"
            f" · precursors {match.precursor_similarity:.2f}</p>"
            "</article>"
            for match in item.precedent_matches[:3]
        ) or "<p class='missing'>No admitted precedent attached.</p>"
        source_claims = "".join(
            f"<li>{escape(source.source_kind)} · {escape(source.source_id)}"
            + (
                f" — {escape(source.quoted_claim)}"
                if source.quoted_claim
                else ""
            )
            + "</li>"
            for source in item.sources
        ) or "<li>no source metadata</li>"
        tone = "good" if item.admission_eligible else "warn"
        step_cards.append(
            "<section class='step-card'>"
            f"<h2>{escape(step.external_step_id)}</h2>"
            "<div class='badges'>"
            f"{_badge('status', item.status, tone)}"
            f"{_badge('evidence', item.strongest_evidence_tier, tone)}"
            f"{_badge('review eligible', str(item.admission_eligible).lower(), tone)}"
            "</div>"
            f"<div class='reaction'>{_reaction_svg(reaction) if reaction else '<div class=missing>Missing structures</div>'}</div>"
            f"<details><summary>Canonical reaction</summary><code>{escape(reaction)}</code></details>"
            "<h3>Evidence gates</h3>"
            "<table><thead><tr><th>Gate</th><th>Result</th><th>Evidence</th></tr></thead>"
            f"<tbody>{_gate_rows(item.gates)}</tbody></table>"
            "<h3>Top admitted precedents</h3>"
            f"<div class='precedents'>{precedent_cards}</div>"
            "<h3>External provenance (not validation)</h3>"
            f"<ul>{source_claims}</ul>"
            "</section>"
        )
    leaf_items = "".join(
        f"<li><code>{escape(smiles)}</code></li>" for smiles in assessment.leaf_smiles
    ) or "<li>none</li>"
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<style>
:root {{ color-scheme: light; --ink:#16202a; --muted:#667085; --line:#dbe2ea; --paper:#fff; --wash:#f4f7fa; --good:#087a55; --warn:#ad5b00; --bad:#b42318; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--wash); color:var(--ink); font:15px/1.5 Inter,Segoe UI,sans-serif; }}
main {{ width:min(1180px,calc(100% - 32px)); margin:32px auto 80px; }}
header,.step-card {{ background:var(--paper); border:1px solid var(--line); border-radius:16px; padding:24px; box-shadow:0 8px 24px #1d293908; }}
.step-card {{ margin-top:20px; }}
h1,h2,h3 {{ margin:0 0 12px; }} h3 {{ margin-top:22px; font-size:1rem; }}
.lede {{ color:var(--muted); max-width:78ch; }}
.badges {{ display:flex; flex-wrap:wrap; gap:8px; margin:12px 0 18px; }}
.badge {{ display:grid; gap:1px; border-radius:9px; padding:7px 11px; background:#eef2f6; }}
.badge.good {{ color:var(--good); background:#e8f7f0; }} .badge.warn {{ color:var(--warn); background:#fff3df; }}
.badge small {{ text-transform:uppercase; letter-spacing:.06em; font-size:.65rem; }}
.reaction {{ overflow-x:auto; min-height:170px; display:grid; place-items:center; border:1px solid var(--line); border-radius:12px; background:#fff; }}
.reaction svg,.precedent svg {{ max-width:100%; height:auto; }}
table {{ width:100%; border-collapse:collapse; }} th,td {{ padding:9px 10px; text-align:left; vertical-align:top; border-bottom:1px solid var(--line); }}
th {{ color:var(--muted); font-size:.75rem; text-transform:uppercase; letter-spacing:.05em; }}
.status {{ display:inline-block; border-radius:999px; padding:2px 8px; font-weight:700; }}
.status.pass {{ color:var(--good); background:#e8f7f0; }} .status.warn,.status.unresolved {{ color:var(--warn); background:#fff3df; }} .status.fail {{ color:var(--bad); background:#feeceb; }} .status.not_run,.status.out_of_scope {{ color:var(--muted); background:#eef2f6; }}
.precedents {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(280px,1fr)); gap:12px; }}
.precedent {{ border:1px solid var(--line); border-radius:12px; padding:10px; overflow:hidden; }}
.precedent svg {{ max-height:155px; }} code {{ overflow-wrap:anywhere; }}
.missing,.muted {{ color:var(--muted); }} details {{ margin-top:8px; }}
</style>
</head>
<body><main>
<header>
<h1>{escape(title)}</h1>
<p class="lede">Every physical step is reconstructed from molecular graphs. External names and claims are shown only as provenance. This artifact is review-only and cannot change route ranking.</p>
<div class="badges">
{_badge('route status', assessment.status, 'good' if assessment.admission_eligible else 'warn')}
{_badge('step count', str(len(assessment.step_assessments)), 'neutral')}
{_badge('route eligible', str(assessment.admission_eligible).lower(), 'good' if assessment.admission_eligible else 'warn')}
{_badge('actionable', 'false', 'warn')}
</div>
<h3>Route topology gates</h3>
<table><thead><tr><th>Gate</th><th>Result</th><th>Evidence</th></tr></thead><tbody>{topology}</tbody></table>
<h3>Unverified route leaves</h3><ul>{leaf_items}</ul>
</header>
{''.join(step_cards)}
</main></body></html>"""


def write_external_route_assessment_review(
    assessment: ExternalRouteAssessment,
    output_json: str | Path,
    output_html: str | Path,
    *,
    title: str = "External route admission review",
) -> dict[str, Any]:
    """Write deterministic JSON evidence and a self-contained HTML review."""

    json_path = Path(output_json)
    html_path = Path(output_html)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(assessment.to_dict(), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    html_path.write_text(
        render_external_route_assessment_html(assessment, title=title),
        encoding="utf-8",
    )
    return {
        "assessment_id": assessment.assessment_id,
        "route_id": assessment.route_id,
        "status": assessment.status,
        "admission_eligible": assessment.admission_eligible,
        "ranking_influence": assessment.ranking_influence,
        "output_json": str(json_path.resolve()),
        "output_html": str(html_path.resolve()),
        "review_version": EXTERNAL_PROPOSAL_REVIEW_VERSION,
    }


__all__ = [
    "EXTERNAL_PROPOSAL_REVIEW_VERSION",
    "render_external_route_assessment_html",
    "write_external_route_assessment_review",
]
