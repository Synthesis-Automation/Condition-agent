"""Self-contained review HTML for held-out v1 coupled-strategy evaluation."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Mapping

from .html_report import molecule_svg, reaction_svg


def load_v1_coupled_strategy_evaluation(
    source: str | Path,
) -> dict[str, Any]:
    """Load and validate one held-out v1 evaluation report."""

    value = json.loads(Path(source).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("evaluation report must be an object")
    if value.get("artifact_type") != "v1_coupled_strategy_heldout_evaluation":
        raise ValueError("unexpected evaluation artifact type")
    if not isinstance(value.get("results"), list):
        raise ValueError("evaluation report requires results")
    return value


def _rank(value: Any) -> str:
    return "not recovered" if value is None else f"rank {int(value)}"


def _candidate_rows(result: Mapping[str, Any]) -> str:
    rows = []
    expected = (result.get("case") or {}).get("expected_intermediate_smiles")
    for rank, candidate in enumerate(
        result.get("baseline_top_level_candidates") or (), 1
    ):
        precursor = str(candidate.get("precursor_smiles") or "")
        marker = " expected" if expected and expected in precursor.split(".") else ""
        rows.append(
            f'<tr class="{marker.strip()}"><td>{rank}</td>'
            f'<td>{molecule_svg(precursor, width=260, height=120)}</td>'
            f'<td><code>{html.escape(str(candidate.get("operator_id") or ""))}</code></td>'
            f'<td>{float(candidate.get("score") or 0):.3f}</td></tr>'
        )
    return "".join(rows) or '<tr><td colspan="4">No validated candidate</td></tr>'


def _action_rows(result: Mapping[str, Any]) -> str:
    rows = []
    expected = (result.get("case") or {}).get("expected_intermediate_smiles")
    for rank, action in enumerate(result.get("promoted_actions") or (), 1):
        intermediate = str(action.get("intermediate_smiles") or "")
        marker = " expected" if intermediate == expected else ""
        physical = (
            reaction_svg(
                str(action.get("first_reaction_smiles") or ""),
                sub_image_size=(140, 120),
            )
            + reaction_svg(
                str(action.get("second_reaction_smiles") or ""),
                sub_image_size=(140, 120),
            )
        )
        rows.append(
            f'<tr class="{marker.strip()}"><td>{rank}</td><td>{physical}</td>'
            f'<td>{html.escape(intermediate)}</td>'
            f'<td>{float(action.get("score") or 0):.3f}</td></tr>'
        )
    return "".join(rows) or '<tr><td colspan="4">No validated macro action</td></tr>'


def _case_card(result: Mapping[str, Any], ordinal: int) -> str:
    case = result.get("case") or {}
    target = str(case.get("target_smiles") or "")
    expected_route = reaction_svg(
        str(case.get("observed_first_reaction_smiles") or ""),
        sub_image_size=(143, 130),
    ) + reaction_svg(
        str(case.get("observed_second_reaction_smiles") or ""),
        sub_image_size=(143, 130),
    )
    return f"""
    <section class="card">
      <h2>{ordinal}. {html.escape(str(case.get('relationship_class') or '').replace('_', ' '))}</h2>
      <div class="meta">held-out {html.escape(str(case.get('split') or ''))} patent
      <code>{html.escape(str(case.get('patent_id') or ''))}</code> · v2 annotation:
      {html.escape(str(case.get('v2_dependency_class') or '').replace('_', ' '))}</div>
      <div class="target">{molecule_svg(target, width=360, height=180)}<code>{html.escape(target)}</code></div>
      <h3>Observed two-step ground truth</h3><div class="route">{expected_route}</div>
      <div class="scoreline">
        Ordinary intermediate: <b>{_rank(result.get('baseline_intermediate_rank'))}</b> ·
        ordinary depth-two pair: <b>{_rank(result.get('baseline_operator_pair_rank'))}</b> ·
        promoted v1 pair: <b>{_rank(result.get('promoted_operator_pair_rank'))}</b>
      </div>
      <details><summary>Ordinary one-step candidates</summary>
        <table><thead><tr><th>Rank</th><th>Precursors</th><th>Operator</th><th>Score</th></tr></thead>
        <tbody>{_candidate_rows(result)}</tbody></table>
      </details>
      <details open><summary>Promoted logical actions (two physical steps)</summary>
        <table><thead><tr><th>Rank</th><th>Validated physical reactions</th><th>Intermediate</th><th>Score</th></tr></thead>
        <tbody>{_action_rows(result)}</tbody></table>
      </details>
      <div class="audit">Validation attempts — ordinary: {int(result.get('baseline_validation_attempt_count') or 0)},
      promoted: {int(result.get('promoted_validation_attempt_count') or 0)}. One-step fallback preserved:
      {str(bool(result.get('one_step_fallback_preserved'))).lower()}.</div>
    </section>"""


def render_v1_coupled_strategy_evaluation_html(
    report: Mapping[str, Any],
    *,
    title: str = "V1 two-step strategy held-out evaluation",
) -> str:
    """Render a compact, self-contained chemistry comparison."""

    metrics = report.get("metrics") or {}
    results = report.get("results") or ()
    panel_id = str(report.get("frozen_panel_id") or "legacy selection")
    library_source = str(report.get("operator_library_source") or "in memory")
    cards = "".join(
        _case_card(result, ordinal)
        for ordinal, result in enumerate(results, 1)
    )
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>
body{{font:14px/1.45 system-ui,sans-serif;margin:0;background:#f4f1ea;color:#17201d}}
main{{max-width:1180px;margin:auto;padding:28px}} h1,h2,h3{{margin:.35em 0}}
.summary,.card{{background:white;border:1px solid #d8d4ca;border-radius:12px;padding:18px;margin:0 0 18px}}
.metrics{{display:flex;gap:12px;flex-wrap:wrap}} .metric{{background:#e8f2ec;border-radius:9px;padding:12px;min-width:170px}}
.metric b{{font-size:24px;display:block}} .meta,.audit{{color:#56625d;margin:6px 0}}
.target,.route{{display:flex;align-items:center;gap:12px;flex-wrap:wrap}} .scoreline{{background:#edf4ff;padding:10px;border-radius:8px;margin:10px 0}}
table{{border-collapse:collapse;width:100%}} th,td{{border-top:1px solid #ddd;padding:7px;text-align:left;vertical-align:top}}
th{{background:#f7f6f2}} tr.expected{{background:#e5f6e9}} code{{font-size:11px;word-break:break-all}}
details{{margin:10px 0}} summary{{cursor:pointer;font-weight:650}} svg{{max-width:100%;height:auto}}
</style></head><body><main>
<section class="summary"><h1>{html.escape(title)}</h1>
<p>Patent-held-out real route segments. V1 structural relationships determine promotion; v2 labels are displayed only as annotations. Every returned macro retains two independently forward-validated physical reactions.</p>
<p class="meta">Panel: <code>{html.escape(panel_id)}</code><br>Library: <code>{html.escape(library_source)}</code></p>
<div class="metrics">
<div class="metric"><b>{int(report.get('panel_case_count') or 0)}</b>held-out cases</div>
<div class="metric"><b>{float(metrics.get('baseline_pair_recall') or 0):.1%}</b>ordinary depth-two recall</div>
<div class="metric"><b>{float(metrics.get('promoted_pair_recall') or 0):.1%}</b>promoted v1 recall</div>
<div class="metric"><b>{int(metrics.get('promoted_validation_attempt_count') or 0)}</b>promoted validations</div>
</div></section>{cards}
</main></body></html>"""


def write_v1_coupled_strategy_evaluation_html(
    report: Mapping[str, Any],
    output_path: str | Path,
    *,
    title: str = "V1 two-step strategy held-out evaluation",
) -> dict[str, Any]:
    """Write the self-contained review page."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        render_v1_coupled_strategy_evaluation_html(report, title=title),
        encoding="utf-8",
        newline="\n",
    )
    return {
        "output_html": str(destination.resolve()),
        "html_bytes": destination.stat().st_size,
        "panel_case_count": int(report.get("panel_case_count") or 0),
    }


__all__ = [
    "load_v1_coupled_strategy_evaluation",
    "render_v1_coupled_strategy_evaluation_html",
    "write_v1_coupled_strategy_evaluation_html",
]
