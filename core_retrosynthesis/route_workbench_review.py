"""Self-contained SVG review for provider-backed route workbench reports."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Mapping

from .html_report import molecule_svg, reaction_svg


def _text(value: Any, fallback: str = "unavailable") -> str:
    rendered = str(value or "").strip()
    return rendered or fallback


def _metric(value: Any) -> str:
    try:
        return f"{float(value):.3f}"
    except (TypeError, ValueError):
        return "—"


def _badge(label: str, state: str) -> str:
    normalized = state.lower().replace("_", "-")
    return (
        f'<span class="badge state-{html.escape(normalized)}">'
        f"{html.escape(label)}</span>"
    )


def _target_panel(case: Mapping[str, Any]) -> str:
    target = _text(case.get("target_smiles"))
    return (
        '<div class="target-panel"><div class="target-drawing">'
        f"{molecule_svg(target, width=480, height=260)}</div>"
        '<div class="target-meta">'
        f'<h2>{html.escape(_text(case.get("target_name"), "Unnamed target"))}</h2>'
        f'<code>{html.escape(target)}</code>'
        f'<p>{html.escape(_text(case.get("challenge_focus"), "No challenge note"))}</p>'
        "</div></div>"
    )


def _verification_table(verification: Mapping[str, Any]) -> str:
    rows = []
    for gate in verification.get("gates") or ():
        status = _text(gate.get("status"))
        rows.append(
            "<tr>"
            f'<td>{html.escape(_text(gate.get("gate")))}</td>'
            f'<td>{_badge(status, status)}</td>'
            f'<td>{html.escape(_text(gate.get("message") or gate.get("summary"), ""))}</td>'
            "</tr>"
        )
    return (
        '<table class="verification"><thead><tr><th>Verification gate</th>'
        "<th>Status</th><th>Evidence</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
    )


def _issue_list(route_report: Mapping[str, Any]) -> str:
    weakest = route_report.get("weakest_issue_id")
    issues = route_report.get("issues") or ()
    if not issues:
        return '<p class="empty">No deterministic route issues.</p>'
    values = []
    for issue in issues:
        is_weakest = issue.get("issue_id") == weakest
        values.append(
            f'<li class="{"weakest" if is_weakest else ""}">'
            f'{_badge(_text(issue.get("severity")), _text(issue.get("severity")))}'
            f'<strong>{html.escape(_text(issue.get("kind")).replace("_", " ").title())}</strong> '
            f'{html.escape(_text(issue.get("message"), ""))}'
            f'{" <em>Weakest issue</em>" if is_weakest else ""}</li>'
        )
    return '<ul class="issues">' + "".join(values) + "</ul>"


def _step_panel(
    step: Mapping[str, Any],
    evidence: Mapping[str, Any],
    index: int,
) -> str:
    precursors = tuple(str(item) for item in step.get("precursor_smiles") or ())
    product = _text(step.get("product_smiles"))
    reaction = f'{".".join(precursors)}>>{product}'
    precedent = evidence.get("precedent_evidence") or {}
    matches = precedent.get("matches") or ()
    best = matches[0] if matches else None
    precedent_line = "No admitted precedent match"
    if best is not None:
        precedent_line = (
            f'{_text(best.get("reaction_id"))}; product similarity '
            f'{_metric(best.get("product_similarity"))}; precursor similarity '
            f'{_metric(best.get("precursor_similarity"))}'
        )
    compatibility = _text(
        evidence.get("reaction_compatibility_disposition"), "unknown"
    )
    selectivity_count = int(evidence.get("selectivity_warning_count") or 0)
    return (
        '<section class="step">'
        f'<div class="step-heading"><h3>Step {index}</h3>'
        f'{_badge(compatibility, compatibility)}'
        f'{_badge(f"{selectivity_count} selectivity warning(s)", "advisory" if selectivity_count else "pass")}'
        "</div>"
        f'<div class="reaction-drawing">{reaction_svg(reaction)}</div>'
        '<dl class="step-evidence">'
        f'<div><dt>Provider</dt><dd>{html.escape(_text(evidence.get("provider_id")))} '
        f'rank {html.escape(_text(evidence.get("provider_rank")))}</dd></div>'
        f'<div><dt>Strategy</dt><dd>{html.escape(_text(evidence.get("strategic_class")).replace("_", " "))}</dd></div>'
        f'<div><dt>Operator</dt><dd><code>{html.escape(_text(evidence.get("operator_id")))}</code></dd></div>'
        f'<div><dt>Precedent</dt><dd>{html.escape(precedent_line)} '
        f'({len(matches)} shown matches)</dd></div>'
        "</dl></section>"
    )


def _leaf_panel(leaf: Mapping[str, Any]) -> str:
    smiles = _text(leaf.get("smiles") or leaf.get("canonical_smiles"))
    terminal = bool(leaf.get("terminal"))
    reason = _text(leaf.get("unresolved_reason"), "terminal")
    return (
        f'<article class="leaf {"terminal" if terminal else "unresolved"}">'
        f'<div class="leaf-drawing">{molecule_svg(smiles, width=300, height=170)}</div>'
        f'{_badge("terminal" if terminal else "unresolved", "pass" if terminal else "fail")}'
        f'<code>{html.escape(smiles)}</code><small>{html.escape(reason)}</small></article>'
    )


def _case_card(case: Mapping[str, Any], index: int) -> str:
    if case.get("execution_status") != "completed":
        return (
            '<article class="case failed">'
            f'<h2>{html.escape(_text(case.get("target_name")))}</h2>'
            f'<p>{html.escape(_text(case.get("error")))}</p></article>'
        )
    report = case.get("report") or {}
    routes = report.get("routes") or ()
    if not routes:
        return (
            '<article class="case failed">'
            f"{_target_panel(case)}<p>No retained route was available.</p></article>"
        )
    best = routes[0]
    route = best.get("route") or {}
    verification = best.get("verification") or {}
    steps = route.get("steps") or ()
    evidence = best.get("step_evidence") or ()
    step_panels = "".join(
        _step_panel(step, evidence[position] if position < len(evidence) else {}, position + 1)
        for position, step in enumerate(steps)
    )
    leaf_panels = "".join(_leaf_panel(item) for item in route.get("leaves") or ())
    route_kind = _text(report.get("route_kind"))
    verification_status = _text(verification.get("status"))
    return (
        f'<article class="case" id="case-{index}">'
        f"{_target_panel(case)}"
        '<div class="route-summary">'
        f'{_badge(route_kind, route_kind)}{_badge(verification_status, verification_status)}'
        f'<span>{len(steps)} retained step(s)</span>'
        f'<span>{sum(not bool(item.get("terminal")) for item in route.get("leaves") or ())} unresolved leaf/leaves</span>'
        "</div>"
        f'<section><h3>Retrosynthetic route</h3>{step_panels}</section>'
        f'<section><h3>Route leaves</h3><div class="leaves">{leaf_panels}</div></section>'
        f'<section><h3>Deterministic issues</h3>{_issue_list(best)}</section>'
        f'<details><summary>Verification gates</summary>{_verification_table(verification)}</details>'
        "</article>"
    )


def render_provider_route_workbench_review_html(
    report: Mapping[str, Any],
    *,
    title: str = "Provider-backed retrosynthesis review",
) -> str:
    """Render best retained routes from a workbench panel as inline SVG HTML."""

    cases = tuple(report.get("cases") or ())
    if not cases:
        raise ValueError("workbench report contains no cases")
    cards = "".join(_case_card(case, index) for index, case in enumerate(cases, 1))
    summary = report.get("summary") or {}
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title><style>
:root{{--paper:#f4f6f5;--card:#fff;--ink:#17201d;--muted:#66716c;--line:#d5ddd9;--green:#1d6a4b;--amber:#976519;--red:#9b4038;--blue:#315f78}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}}
header,main{{max-width:1500px;margin:auto;padding:20px}}header{{padding-bottom:6px}}header h1{{margin:0 0 4px}}header p{{margin:0;color:var(--muted)}}
.case{{background:var(--card);border:1px solid var(--line);border-radius:10px;margin:18px 0 30px;padding:16px}}.case.failed{{border-left:5px solid var(--red)}}
.target-panel{{display:grid;grid-template-columns:minmax(300px,480px) 1fr;gap:20px;align-items:center;border-bottom:1px solid var(--line);padding-bottom:12px}}
.target-drawing svg{{display:block;width:100%;height:250px}}.target-meta h2{{margin:0 0 8px}}code{{overflow-wrap:anywhere}}.target-meta code{{display:block}}.target-meta p{{color:var(--muted)}}
.route-summary,.step-heading{{display:flex;align-items:center;gap:8px;flex-wrap:wrap}}.route-summary{{padding:12px 0}}.route-summary span{{color:var(--muted)}}
.badge{{display:inline-block;border-radius:999px;padding:3px 8px;background:#e9eeeb;color:var(--muted);font-size:12px}}.state-pass,.state-solved,.state-verified-with-cautions{{background:#dceee5;color:var(--green)}}
.state-fail,.state-failed,.state-partial{{background:#f2dfdc;color:var(--red)}}.state-warn,.state-advisory,.state-strong{{background:#f5ead5;color:var(--amber)}}
.step{{border-top:1px solid var(--line);padding:10px 0 16px}}.step-heading h3{{margin:0 auto 0 0}}.reaction-drawing{{overflow-x:auto;background:#fafbfa}}.reaction-drawing svg{{display:block;width:100%;min-width:700px;height:240px}}
.step-evidence{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:4px 18px;margin:8px 0}}.step-evidence div{{border-top:1px solid var(--line);padding:5px 0;min-width:0}}dt{{color:var(--muted);font-size:12px}}dd{{margin:2px 0;overflow-wrap:anywhere}}
.leaves{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:10px}}.leaf{{border-left:4px solid var(--green);background:#f8faf9;padding:8px;min-width:0}}.leaf.unresolved{{border-left-color:var(--red)}}.leaf-drawing svg{{display:block;width:100%;height:165px}}.leaf code,.leaf small{{display:block;margin-top:5px}}.leaf small{{color:var(--muted)}}
.issues{{padding-left:20px}}.issues li{{margin:7px 0}}.issues .badge{{margin-right:7px}}.issues li.weakest{{border-left:4px solid var(--red);padding-left:8px}}.issues em{{color:var(--red)}}.empty{{color:var(--muted)}}
details{{margin-top:14px}}summary{{cursor:pointer;font-weight:600}}.verification{{border-collapse:collapse;width:100%;margin-top:8px}}.verification th,.verification td{{border-top:1px solid var(--line);padding:7px;text-align:left;vertical-align:top}}
@media(max-width:780px){{.target-panel{{grid-template-columns:1fr}}.step-evidence{{grid-template-columns:1fr}}.reaction-drawing svg{{min-width:580px}}}}@media print{{body{{background:#fff}}.case{{break-inside:avoid}}}}
</style></head><body><header><h1>{html.escape(title)}</h1>
<p>{html.escape(_text(summary.get("completed_count"), "0"))} completed case(s); best retained route shown for each target. Solved status reflects the configured terminal policy, not experimental proof.</p>
</header><main>{cards}</main></body></html>"""


def write_provider_route_workbench_review_html(
    source_report: str | Path,
    output_html: str | Path,
    *,
    title: str = "Provider-backed retrosynthesis review",
) -> dict[str, Any]:
    """Load a workbench JSON report and write a self-contained SVG review."""

    source = Path(source_report).resolve()
    report = json.loads(source.read_text(encoding="utf-8"))
    document = render_provider_route_workbench_review_html(report, title=title)
    output = Path(output_html).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8", newline="\n")
    return {
        "source_report": str(source),
        "output_html": str(output),
        "case_count": len(report.get("cases") or ()),
        "html_bytes": output.stat().st_size,
    }


__all__ = [
    "render_provider_route_workbench_review_html",
    "write_provider_route_workbench_review_html",
]
