"""Human-readable and blind-review artifacts for Phase 4 assessments."""

from __future__ import annotations

from html import escape
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

from .chemistry import molecule_without_maps
from .partition_assessment import PartitionAssessmentResult
from .synthetic_partition import SyntheticPartition, analyze_partition_target


PARTITION_ASSESSMENT_REVIEW_VERSION = "partition_assessment_review.v1"
_MODULE_COLORS = (
    (0.24, 0.49, 0.91),
    (0.91, 0.35, 0.28),
    (0.20, 0.68, 0.42),
    (0.62, 0.37, 0.82),
    (0.95, 0.62, 0.16),
    (0.12, 0.68, 0.72),
)
_MODULE_COLOR_HEX = ("#3d7de8", "#e85947", "#33ad6b", "#9e5ed1", "#f29e29", "#1fadb8")
_MAX_DEPICTED_PRECEDENTS = 3


def _items(values: Iterable[str]) -> str:
    items = tuple(values)
    if not items:
        return "<span class='muted'>none</span>"
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def _embedded_svg(value: str) -> str:
    start = value.find("<svg")
    return value[start:] if start >= 0 else value


def _molecule_svg(smiles: str, *, width: int = 260, height: int = 170) -> str:
    molecule = molecule_without_maps(smiles)
    if molecule is None:
        return "<div class='structure-error'>Structure unavailable</div>"
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().clearBackground = False
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    return _embedded_svg(drawer.GetDrawingText())


def _reaction_svg(
    reaction_smiles: str,
    *,
    width: int = 940,
    height: int = 230,
) -> str:
    try:
        reaction = rdChemReactions.ReactionFromSmarts(
            reaction_smiles,
            useSmiles=True,
        )
    except Exception:
        reaction = None
    if reaction is None:
        return "<div class='structure-error'>Reaction depiction unavailable</div>"
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().clearBackground = False
    drawer.DrawReaction(reaction)
    drawer.FinishDrawing()
    return _embedded_svg(drawer.GetDrawingText())


def _target_partition_svg(partition: SyntheticPartition) -> str:
    molecule = Chem.MolFromSmiles(partition.target_smiles)
    if molecule is None:
        return "<div class='structure-error'>Target depiction unavailable</div>"
    _, _, target_atoms, _ = analyze_partition_target(partition.target_smiles)
    index_by_map = {item.atom_map: item.atom_index for item in target_atoms}
    atom_colors: dict[int, tuple[float, float, float]] = {}
    atom_radii: dict[int, float] = {}
    highlighted: list[int] = []
    for module_index, module in enumerate(partition.modules):
        color = _MODULE_COLORS[module_index % len(_MODULE_COLORS)]
        for atom_map in module.target_atom_maps:
            atom_index = index_by_map.get(atom_map)
            if atom_index is None:
                continue
            highlighted.append(atom_index)
            atom_colors[atom_index] = color
            atom_radii[atom_index] = 0.34
    drawer = rdMolDraw2D.MolDraw2DSVG(560, 330)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.addAtomIndices = False
    drawer.DrawMolecule(
        molecule,
        highlightAtoms=highlighted,
        highlightAtomColors=atom_colors,
        highlightAtomRadii=atom_radii,
    )
    drawer.FinishDrawing()
    return _embedded_svg(drawer.GetDrawingText())


def _badge(label: str, value: str, tone: str = "neutral") -> str:
    return (
        f"<span class='badge {escape(tone)}'><span>{escape(label)}</span>"
        f"<strong>{escape(value)}</strong></span>"
    )


def _recipe_summary(recommendation: Mapping[str, Any]) -> str:
    recipe = recommendation.get("resolved_recipe")
    if not isinstance(recipe, Mapping):
        return str(recommendation.get("recipe_id") or "unresolved recipe")
    components = recipe.get("components") or ()
    labels = []
    for component in components:
        if not isinstance(component, Mapping):
            continue
        labels.append(
            str(
                component.get("canonical_name")
                or component.get("substance_id")
                or component.get("raw_identifier")
                or "component"
            )
        )
    return ", ".join(labels) or str(
        recommendation.get("recipe_core_id")
        or recommendation.get("recipe_id")
        or "resolved recipe"
    )


def _module_legend(partition: SyntheticPartition) -> str:
    return "".join(
        "<div class='module-key'>"
        f"<span style='background:{_MODULE_COLOR_HEX[index % len(_MODULE_COLOR_HEX)]}'></span>"
        f"<strong>M{index + 1}</strong>"
        f"<small>{len(module.target_atom_maps)} atoms · maps "
        f"{escape(', '.join(map(str, module.target_atom_maps)))}</small>"
        "</div>"
        for index, module in enumerate(partition.modules)
    )


def _frontier_gallery(
    realization: Any,
    partition: SyntheticPartition,
) -> str:
    if not realization.frontier_states:
        return "<p class='empty'>No frontier structures retained.</p>"
    module_labels = {
        module.module_id: f"M{index}"
        for index, module in enumerate(partition.modules, start=1)
    }
    return "".join(
        "<article class='molecule-tile'>"
        f"<div class='molecule-graphic'>{_molecule_svg(state.mapped_smiles)}</div>"
        f"<div class='molecule-label'><strong>{escape(' + '.join(module_labels.get(item, item) for item in state.module_ids) or 'unassigned')}</strong>"
        f"<span>{escape(state.evidence_status)}</span></div>"
        f"<code>{escape(state.mapped_smiles)}</code>"
        "</article>"
        for state in realization.frontier_states
    )


def _precedent_gallery(step: Any) -> str:
    matches = step.precedent_lookup.matches if step.precedent_lookup else ()
    if not matches:
        return "<div class='empty'>No admitted structure-backed precedent.</div>"
    depicted = matches[:_MAX_DEPICTED_PRECEDENTS]
    cards = "".join(
        "<article class='precedent-card'>"
        "<div class='precedent-heading'>"
        f"<strong>Precedent {index}</strong>"
        f"{_badge('product', f'{item.product_similarity:.2f}', 'evidence')}"
        f"{_badge('precursor', f'{item.precursor_similarity:.2f}', 'evidence')}"
        "</div>"
        f"<div class='reaction-graphic precedent-graphic'>{_reaction_svg(item.mapped_reaction_smiles, width=860, height=210)}</div>"
        "<details><summary>Source identifiers and reaction SMILES</summary>"
        f"<p><strong>Reference:</strong> <code>{escape(item.reference_id or '—')}</code></p>"
        f"<p><strong>Reaction:</strong> <code>{escape(item.reaction_id or '—')}</code></p>"
        f"<p><code>{escape(item.mapped_reaction_smiles)}</code></p></details>"
        "</article>"
        for index, item in enumerate(depicted, start=1)
    )
    remainder = ""
    if len(matches) > len(depicted):
        remainder = (
            f"<p class='muted'>+ {len(matches) - len(depicted)} additional admitted "
            "matches are retained in the embedded JSON.</p>"
        )
    return cards + remainder


def _condition_and_caution_panel(step: Any) -> str:
    recommendations = step.condition_evidence.recommendations
    recipes = "".join(
        "<li><strong>Recipe</strong> " + escape(_recipe_summary(item)) + "</li>"
        for item in recommendations
    )
    caution_values = (*step.cautions, *step.warnings, *step.structural_issues)
    return (
        "<div class='evidence-grid'>"
        "<section class='evidence-box'><h4>Conditions</h4>"
        f"{_badge('support', step.condition_evidence.status, 'condition')}"
        f"<ul>{recipes or '<li class="muted">No canonical recipe retrieved</li>'}</ul>"
        "</section>"
        "<section class='evidence-box'><h4>Compatibility and selectivity</h4>"
        f"{_items(caution_values)}"
        "</section></div>"
    )


def _route_panel(
    route: Any,
    realization: Any,
    partition: SyntheticPartition,
    route_index: int,
) -> str:
    step_cards = []
    for step_index, step in enumerate(route.step_assessments, start=1):
        weakest = step.step_id == route.weakest_step_id
        step_cards.append(
            f"<article class='step-card{' weakest' if weakest else ''}'>"
            "<div class='step-heading'>"
            f"<div><span class='step-number'>{step_index}</span><h3>Route step {step_index}</h3>"
            f"<small>Retrosynthetic depth {step.depth}</small></div>"
            "<div class='badges'>"
            f"{_badge('structure', step.structural_status, 'good' if step.structural_status == 'validated' else 'danger')}"
            f"{_badge('precedent', step.precedent_evidence_level, 'evidence')}"
            f"{_badge('conditions', step.condition_evidence.status, 'condition')}"
            f"{'<span class="weakest-flag">weakest step</span>' if weakest else ''}"
            "</div></div>"
            "<p class='scheme-label'>Forward view · precursors → product</p>"
            f"<div class='reaction-graphic'>{_reaction_svg(step.reaction_smiles)}</div>"
            "<details class='smiles-details'><summary>Reaction SMILES and identifiers</summary>"
            f"<p><code>{escape(step.reaction_smiles)}</code></p>"
            f"<p>Operator <code>{escape(step.operator_id)}</code><br>"
            f"Template <code>{escape(step.template_id)}</code></p></details>"
            "<details class='evidence-details' open>"
            f"<summary>Depicted precedents ({step.independent_precedent_count} independent references)</summary>"
            f"{_precedent_gallery(step)}</details>"
            f"{_condition_and_caution_panel(step)}"
            "</article>"
        )
    interface_chips = "".join(
        "<div class='interface-chip'>"
        f"<strong>{escape(item.realization_status)}</strong>"
        f"<span>{escape(item.weakest_precedent_evidence)} precedent</span>"
        f"<small>{escape(', '.join(item.condition_statuses) or 'conditions unavailable')}</small>"
        "</div>"
        for item in route.interface_assessments
    )
    route_tone = "danger" if route.status == "hard_incompatible" else "condition"
    return (
        f"<section class='route-panel{' active' if route_index == 1 else ''}' data-route-panel='{route_index}'>"
        "<div class='route-summary card'>"
        "<div><p class='eyebrow'>ROUTE OVERVIEW</p>"
        f"<h2>Route {route_index}</h2><div class='badges'>"
        f"{_badge('assessment', route.status, route_tone)}"
        f"{_badge('realization', route.source_realization_status, 'neutral')}"
        f"{_badge('steps', str(len(route.step_assessments)), 'neutral')}"
        "</div></div>"
        "<div class='route-metrics'>"
        f"<div><strong>{route.precedent_supported_step_count}/{len(route.step_assessments)}</strong><span>precedent-backed steps</span></div>"
        f"<div><strong>{route.condition_supported_step_count}/{len(route.step_assessments)}</strong><span>condition-backed steps</span></div>"
        f"<div><strong>{route.protection_burden_count}</strong><span>protection / auxiliary atoms</span></div>"
        f"<div><strong>{route.unresolved_latent_atom_count}</strong><span>unclassified latent atoms</span></div>"
        "</div></div>"
        "<section class='card'><div class='section-heading'><div><p class='eyebrow'>RETROSYNTHETIC FRONTIER</p>"
        "<h2>Target-derived modules in their latent forms</h2></div>"
        "<p>Structures are shown independently; tactical atoms do not change module ownership.</p></div>"
        f"<div class='frontier-gallery'>{_frontier_gallery(realization, partition)}</div></section>"
        "<section class='card'><div class='section-heading'><div><p class='eyebrow'>STRATEGIC INTERFACES</p>"
        "<h2>Coverage at a glance</h2></div></div>"
        f"<div class='interface-strip'>{interface_chips or '<p class="empty">No interfaces.</p>'}</div></section>"
        f"{''.join(step_cards) or '<section class="card empty">No executable route steps retained.</section>'}"
        "<section class='card review'>"
        "<div class='section-heading'><div><p class='eyebrow'>YOUR ASSESSMENT</p><h2>Chemist review</h2></div>"
        "<p>Your entries are saved locally in this browser.</p></div>"
        "<div class='review-grid'>"
        f"<label>Overall route judgment<select data-field='judgment' data-case='{escape(route.assessment_id)}'>"
        "<option value=''>Select…</option><option>credible</option>"
        "<option>credible_with_changes</option><option>not_credible</option>"
        "<option>insufficient_information</option></select></label>"
        f"<label>Weakest step<select data-field='weakest-step' data-case='{escape(route.assessment_id)}'>"
        "<option value=''>Select…</option>"
        + "".join(
            f"<option value='{index}'>{index} — depth {step.depth}</option>"
            for index, step in enumerate(route.step_assessments, start=1)
        )
        + "</select></label></div>"
        f"<label>Conditions, selectivity, handles, and protection notes<textarea data-field='notes' data-case='{escape(route.assessment_id)}'></textarea></label>"
        "<details><summary>System warnings</summary>"
        f"{_items(route.warnings)}</details></section></section>"
    )


def render_partition_assessment_html(result: PartitionAssessmentResult) -> str:
    """Render a structure-first static evidence and chemist-review report."""

    realization_by_id = {
        item.realization_id: item for item in result.source.realizations
    }
    route_panels = "".join(
        _route_panel(
            route,
            realization_by_id[route.realization_id],
            result.source.partition,
            route_index,
        )
        for route_index, route in enumerate(result.route_assessments, start=1)
    )
    route_tabs = "".join(
        f"<button class='route-tab{' active' if index == 1 else ''}' data-route-tab='{index}'>"
        f"Route {index}<small>{escape(route.status)}</small></button>"
        for index, route in enumerate(result.route_assessments, start=1)
    )
    partition = result.source.partition
    payload = json.dumps(
        result.to_dict(), sort_keys=True, separators=(",", ":")
    ).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Synthetic partition structure review</title>
<style>
:root {{ --ink:#17212b; --muted:#687583; --line:#dce3e8; --paper:#fff; --canvas:#eef2f4; --navy:#173b56; --cyan:#26a6a1; --amber:#b86a13; --red:#b83c3c; }}
* {{ box-sizing:border-box; }} body {{ margin:0; color:var(--ink); background:var(--canvas); font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif; }}
header,.card,.step-card {{ max-width:1240px; margin:0 auto 1.25rem; background:var(--paper); border:1px solid var(--line); border-radius:16px; box-shadow:0 5px 20px #1c34420b; }}
header {{ margin-top:1.5rem; padding:1.5rem; }} .hero {{ display:grid; grid-template-columns:minmax(420px,1.1fr) minmax(320px,.9fr); gap:1.5rem; align-items:center; }}
h1 {{ margin:.2rem 0 .5rem; font-size:clamp(1.8rem,3vw,2.8rem); letter-spacing:-.035em; }} h2 {{ margin:.15rem 0 .65rem; }} h3,h4 {{ margin:.2rem 0 .65rem; }} p {{ line-height:1.5; }}
.eyebrow {{ color:var(--cyan); font-size:.74rem; font-weight:800; letter-spacing:.13em; margin:0; }} .muted,.empty,small {{ color:var(--muted); }} code {{ overflow-wrap:anywhere; font-size:.78rem; }}
.target-graphic,.reaction-graphic,.molecule-graphic {{ display:flex; align-items:center; justify-content:center; overflow:auto; background:#fff; }} .target-graphic svg,.reaction-graphic svg,.molecule-graphic svg {{ max-width:100%; height:auto; }}
.module-legend {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:.55rem; margin:1rem 0; }} .module-key {{ display:grid; grid-template-columns:12px auto; column-gap:.55rem; align-items:center; padding:.6rem; border:1px solid var(--line); border-radius:10px; }} .module-key>span {{ width:12px; height:32px; border-radius:8px; grid-row:1/3; }} .module-key small {{ display:block; overflow-wrap:anywhere; }}
.notice {{ border-left:4px solid var(--amber); background:#fff8ed; padding:.8rem 1rem; border-radius:7px; }} .route-nav {{ position:sticky; top:0; z-index:20; display:flex; gap:.6rem; max-width:1240px; margin:0 auto 1.25rem; padding:.65rem; background:#eef2f4eF; backdrop-filter:blur(10px); }}
.route-tab {{ border:1px solid var(--line); background:white; color:var(--ink); border-radius:10px; padding:.7rem 1rem; font:inherit; font-weight:750; cursor:pointer; }} .route-tab small {{ display:block; font-weight:500; }} .route-tab.active {{ color:white; border-color:var(--navy); background:var(--navy); }} .route-tab.active small {{ color:#d8e9f3; }}
.route-panel {{ display:none; }} .route-panel.active {{ display:block; }} .card,.step-card {{ padding:1.35rem; }} .route-summary {{ display:flex; align-items:center; justify-content:space-between; gap:1rem; }} .route-metrics {{ display:grid; grid-template-columns:repeat(4,minmax(120px,1fr)); gap:.65rem; }} .route-metrics div {{ background:#f5f8fa; border-radius:10px; padding:.7rem; }} .route-metrics strong,.route-metrics span {{ display:block; }} .route-metrics strong {{ font-size:1.35rem; }}
.badges {{ display:flex; gap:.45rem; flex-wrap:wrap; align-items:center; }} .badge {{ display:inline-flex; gap:.35rem; align-items:center; border:1px solid var(--line); border-radius:999px; padding:.28rem .55rem; font-size:.73rem; }} .badge span {{ color:var(--muted); }} .badge.good {{ border-color:#9ad7b4; background:#effbf4; }} .badge.evidence {{ border-color:#b8caea; background:#f2f6fd; }} .badge.condition {{ border-color:#efd09f; background:#fff8ec; }} .badge.danger {{ border-color:#efaaaa; background:#fff1f1; }}
.section-heading {{ display:flex; justify-content:space-between; gap:1rem; align-items:end; }} .section-heading>p {{ max-width:460px; color:var(--muted); margin:.2rem 0; }} .frontier-gallery {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:.85rem; }} .molecule-tile {{ border:1px solid var(--line); border-radius:12px; overflow:hidden; min-width:0; }} .molecule-label {{ display:flex; justify-content:space-between; padding:.5rem .7rem 0; }} .molecule-tile code {{ display:block; padding:.5rem .7rem .8rem; }}
.interface-strip {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(210px,1fr)); gap:.7rem; }} .interface-chip {{ border-top:4px solid var(--cyan); background:#f5fafb; border-radius:9px; padding:.75rem; }} .interface-chip>* {{ display:block; }}
.step-card {{ border-left:6px solid var(--navy); }} .step-card.weakest {{ border-left-color:var(--amber); }} .step-heading,.precedent-heading {{ display:flex; justify-content:space-between; align-items:center; gap:.8rem; }} .step-heading>div:first-child {{ display:grid; grid-template-columns:42px auto; column-gap:.7rem; align-items:center; }} .step-heading small {{ grid-column:2; }} .step-number {{ grid-row:1/3; display:grid; place-items:center; width:42px; height:42px; color:white; background:var(--navy); border-radius:50%; font-weight:850; }} .weakest-flag {{ color:#7c4306; background:#ffe8c4; border-radius:999px; padding:.35rem .65rem; font-size:.76rem; font-weight:800; }}
.scheme-label {{ margin:1rem 0 -.55rem; color:var(--muted); font-size:.75rem; font-weight:750; letter-spacing:.04em; text-transform:uppercase; }} .reaction-graphic {{ min-height:190px; margin:1rem 0; border:1px solid var(--line); border-radius:12px; }} details {{ margin:.75rem 0; }} summary {{ cursor:pointer; font-weight:750; }} .evidence-details {{ background:#f7f9fa; border:1px solid var(--line); border-radius:12px; padding:.8rem; }} .precedent-card {{ background:white; border:1px solid var(--line); border-radius:10px; padding:.75rem; margin:.7rem 0; }} .precedent-graphic {{ min-height:170px; margin:.5rem 0; border:0; }}
.evidence-grid,.review-grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:.8rem; }} .evidence-box {{ border:1px solid var(--line); border-radius:10px; padding:.85rem; }} .review {{ border-top:5px solid var(--cyan); }} label {{ display:block; font-weight:700; margin:.7rem 0; }} select,textarea {{ display:block; width:100%; margin-top:.4rem; border:1px solid #aebbc5; border-radius:8px; padding:.7rem; font:inherit; background:white; }} textarea {{ min-height:7rem; }} .structure-error {{ padding:2rem; color:var(--red); }}
@media (max-width:800px) {{ .hero,.evidence-grid,.review-grid {{ grid-template-columns:1fr; }} .route-summary,.section-heading {{ display:block; }} .route-metrics {{ grid-template-columns:repeat(2,1fr); margin-top:1rem; }} .module-legend {{ grid-template-columns:1fr; }} header,.card,.step-card {{ border-radius:0; border-left:0; border-right:0; }} body {{ padding:0; }} }}
</style></head><body>
<header><div class="hero"><div class="target-graphic">{_target_partition_svg(partition)}</div>
<div><p class="eyebrow">K-WAY SYNTHETIC PARTITION · PHASE 4</p><h1>Structure-first route review</h1>
<p>Review the molecular structures and transformations first. Evidence details remain available directly beneath each scheme.</p>
<div class="badges">{_badge("modules", str(partition.k), "neutral")}{_badge("routes", str(len(result.route_assessments)), "neutral")}{_badge("assessment", result.status, "condition")}</div>
<div class="module-legend">{_module_legend(partition)}</div>
<div class="notice">Review-only: the evidence shown here does not change Phase 3 route ordering or realization status.</div></div></div>
<details><summary>Global system warnings</summary>{_items(result.warnings)}</details></header>
<nav class="route-nav" aria-label="Route selection">{route_tabs or '<span class="empty">No routes available</span>'}</nav>
{route_panels or '<section class="card"><p>No route realization was available to assess.</p></section>'}
<script type="application/json" id="partition-assessment-data">{payload}</script>
<script>
const storageKey='partition-assessment-review:' + {json.dumps(partition.partition_id)};
const saved=JSON.parse(localStorage.getItem(storageKey)||'{{}}');
document.querySelectorAll('[data-case]').forEach(el=>{{
  const fieldKey=el.dataset.case+':'+el.dataset.field;
  if(saved[fieldKey]!==undefined) el.value=saved[fieldKey];
  el.addEventListener('input',()=>{{saved[fieldKey]=el.value;localStorage.setItem(storageKey,JSON.stringify(saved));}});
}});
document.querySelectorAll('[data-route-tab]').forEach(button=>button.addEventListener('click',()=>{{
  document.querySelectorAll('[data-route-tab]').forEach(item=>item.classList.toggle('active',item===button));
  document.querySelectorAll('[data-route-panel]').forEach(panel=>panel.classList.toggle('active',panel.dataset.routePanel===button.dataset.routeTab));
  window.scrollTo({{top:document.querySelector('.route-nav').offsetTop,behavior:'smooth'}});
}}));
</script></body></html>"""


def write_partition_assessment_review(
    result: PartitionAssessmentResult,
    json_path: str | Path,
    html_path: str | Path,
) -> None:
    """Write the self-contained Phase 4 JSON and HTML artifacts."""

    output_json = Path(json_path)
    output_html = Path(html_path)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    output_html.write_text(render_partition_assessment_html(result), "utf-8")


def _review_source_kind(route_status: str, source_status: str) -> str:
    if route_status == "hard_incompatible":
        return "negative_control"
    if source_status != "fully_realized":
        return "abstention"
    return "realization"


def build_partition_blind_review_packet(
    results: Iterable[PartitionAssessmentResult],
    *,
    seed: int = 47,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build separated chemistry cases and an answer key with stable blinding."""

    cases: list[dict[str, Any]] = []
    answers: list[dict[str, Any]] = []
    for result in results:
        by_realization = {
            item.realization_id: item for item in result.source.realizations
        }
        for assessment in result.route_assessments:
            realization = by_realization[assessment.realization_id]
            case_id = (
                "SPREV1:"
                + hashlib.sha256(
                    f"{seed}\0{assessment.assessment_id}".encode("utf-8")
                ).hexdigest()[:20]
            )
            steps = []
            for step in assessment.step_assessments:
                matches = step.precedent_lookup.matches if step.precedent_lookup else ()
                steps.append(
                    {
                        "reaction_smiles": step.reaction_smiles,
                        "precedents": [
                            {
                                "reference_id": item.reference_id,
                                "reaction_smiles": item.mapped_reaction_smiles,
                            }
                            for item in matches
                        ],
                        "condition_recipes": [
                            {
                                "recipe_id": item.get("recipe_id"),
                                "resolved_recipe": item.get("resolved_recipe"),
                                "precedent_reference_ids": item.get(
                                    "precedent_reference_ids"
                                ),
                            }
                            for item in step.condition_evidence.recommendations
                        ],
                        "compatibility_cautions": list(step.cautions),
                    }
                )
            cases.append(
                {
                    "case_id": case_id,
                    "target_smiles": result.source.target_smiles,
                    "module_target_atom_maps": [
                        list(module.target_atom_maps)
                        for module in result.source.partition.modules
                    ],
                    "frontier_smiles": [
                        state.mapped_smiles for state in realization.frontier_states
                    ],
                    "steps": steps,
                    "review_questions": {
                        "route_credible": None,
                        "weakest_step_index": None,
                        "conditions_plausible": None,
                        "handle_or_protection_conflict": None,
                        "notes": "",
                    },
                    "schema_version": PARTITION_ASSESSMENT_REVIEW_VERSION,
                }
            )
            answers.append(
                {
                    "case_id": case_id,
                    "assessment_id": assessment.assessment_id,
                    "realization_id": assessment.realization_id,
                    "source_kind": _review_source_kind(
                        assessment.status,
                        assessment.source_realization_status,
                    ),
                    "system_assessment_status": assessment.status,
                    "source_realization_status": assessment.source_realization_status,
                    "weakest_step_id": assessment.weakest_step_id,
                }
            )

    def order(value: Mapping[str, Any]) -> str:
        return hashlib.sha256(f"{seed}\0{value['case_id']}".encode("utf-8")).hexdigest()

    cases.sort(key=order)
    answers.sort(key=lambda value: value["case_id"])
    kinds = sorted({item["source_kind"] for item in answers})
    packet = {
        "artifact_type": "synthetic_partition_blind_review_packet",
        "seed": seed,
        "case_count": len(cases),
        "cases": cases,
        "warnings": [
            code
            for kind, code in (
                ("negative_control", "NO_NEGATIVE_CONTROL_PRESENT"),
                ("abstention", "NO_ABSTENTION_CASE_PRESENT"),
            )
            if kind not in kinds
        ],
        "schema_version": PARTITION_ASSESSMENT_REVIEW_VERSION,
    }
    answer_key = {
        "artifact_type": "synthetic_partition_blind_review_answer_key",
        "case_count": len(answers),
        "source_kinds": kinds,
        "answers": answers,
        "schema_version": PARTITION_ASSESSMENT_REVIEW_VERSION,
    }
    return packet, answer_key


def write_partition_blind_review_packet(
    results: Iterable[PartitionAssessmentResult],
    packet_path: str | Path,
    answer_key_path: str | Path,
    *,
    seed: int = 47,
) -> None:
    """Write blind chemistry cases separately from system labels and ranks."""

    packet, answer_key = build_partition_blind_review_packet(results, seed=seed)
    review_output = Path(packet_path)
    key_output = Path(answer_key_path)
    review_output.parent.mkdir(parents=True, exist_ok=True)
    key_output.parent.mkdir(parents=True, exist_ok=True)
    review_output.write_text(
        json.dumps(packet, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    key_output.write_text(
        json.dumps(answer_key, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


__all__ = [
    "PARTITION_ASSESSMENT_REVIEW_VERSION",
    "build_partition_blind_review_packet",
    "render_partition_assessment_html",
    "write_partition_assessment_review",
    "write_partition_blind_review_packet",
]
