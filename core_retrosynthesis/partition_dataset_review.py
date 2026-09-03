"""Deterministic K-way partition evaluation packets from observed route trees."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from html import escape
import hashlib
import json
from pathlib import Path
import random
from typing import Any, Iterable, Sequence

from .chemistry import digest
from .html_report import molecule_svg, reaction_svg
from .partition_projection import RoutePartitionProjection, project_route_partitions
from .partition_review import render_partition_svg
from .route_contract import MoleculeOccurrenceNode, ReactionRouteTree, RouteReactionNode
from .route_conversion import iter_route_trees


PARTITION_DATASET_REVIEW_SCHEMA_VERSION = "1.0"
PARTITION_DATASET_REVIEW_ALGORITHM_VERSION = "observed_route_k_way_review.v2"
DEFAULT_PARTITION_DATASET_REVIEW_SEED = 20_260_903
_MODULE_COLORS = (
    "#4f8fe8",
    "#e8735f",
    "#61b87a",
    "#b482db",
    "#efb24a",
    "#4bbec0",
)


@dataclass(frozen=True)
class PartitionDatasetReviewCase:
    """One observed route and all of its resolved K-way frontiers."""

    case_id: str
    selection_rank: int
    route_tree: ReactionRouteTree
    projection: RoutePartitionProjection
    maximum_k: int
    fully_projected: bool

    def to_dict(self) -> dict[str, Any]:
        """Return a self-contained JSON-compatible case."""

        return {
            "case_id": self.case_id,
            "selection_rank": self.selection_rank,
            "route_tree": self.route_tree.to_dict(),
            "projection": self.projection.to_dict(),
            "maximum_k": self.maximum_k,
            "fully_projected": self.fully_projected,
        }


@dataclass(frozen=True)
class PartitionDatasetReview:
    """An enriched, reproducible K-way functionality test packet."""

    source_path: str
    source_sha256: str
    source_route_count: int
    projection_count_by_maximum_k: tuple[tuple[int, int], ...]
    minimum_selected_k: int
    eligible_route_count: int
    sample_size: int
    seed: int
    selection_method: str
    cases: tuple[PartitionDatasetReviewCase, ...]
    warnings: tuple[str, ...]
    algorithm_version: str = PARTITION_DATASET_REVIEW_ALGORITHM_VERSION
    schema_version: str = PARTITION_DATASET_REVIEW_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != PARTITION_DATASET_REVIEW_SCHEMA_VERSION:
            raise ValueError("unsupported partition dataset-review schema")
        if self.algorithm_version != PARTITION_DATASET_REVIEW_ALGORITHM_VERSION:
            raise ValueError("unsupported partition dataset-review algorithm")
        if self.sample_size != len(self.cases):
            raise ValueError("partition dataset-review sample count disagrees")
        if len({case.case_id for case in self.cases}) != len(self.cases):
            raise ValueError("partition dataset-review cases must be unique")

    def to_dict(self) -> dict[str, Any]:
        """Return a self-contained JSON-compatible review packet."""

        return {
            "source_path": self.source_path,
            "source_sha256": self.source_sha256,
            "source_route_count": self.source_route_count,
            "projection_count_by_maximum_k": [
                [maximum_k, count]
                for maximum_k, count in self.projection_count_by_maximum_k
            ],
            "minimum_selected_k": self.minimum_selected_k,
            "eligible_route_count": self.eligible_route_count,
            "sample_size": self.sample_size,
            "seed": self.seed,
            "selection_method": self.selection_method,
            "cases": [case.to_dict() for case in self.cases],
            "warnings": list(self.warnings),
            "algorithm_version": self.algorithm_version,
            "schema_version": self.schema_version,
        }


def _file_sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _stable_key(tree_id: str, seed: int) -> str:
    return hashlib.sha256(f"{seed}\0{tree_id}".encode("utf-8")).hexdigest()


def build_partition_dataset_review(
    source: str | Path,
    *,
    sample_size: int = 10,
    seed: int = DEFAULT_PARTITION_DATASET_REVIEW_SEED,
    minimum_k: int = 3,
) -> PartitionDatasetReview:
    """Project a corpus and select a reproducible K-way-enriched sample."""

    if sample_size < 1:
        raise ValueError("partition review sample size must be positive")
    if minimum_k < 1:
        raise ValueError("minimum selected k must be positive")
    source_path = Path(source).resolve()
    projected: list[tuple[ReactionRouteTree, RoutePartitionProjection, int]] = []
    counts: Counter[int] = Counter()
    source_count = 0
    for tree in iter_route_trees(source_path):
        source_count += 1
        projection = project_route_partitions(tree)
        maximum_k = max(
            (frontier.partition.k for frontier in projection.frontiers),
            default=0,
        )
        counts[maximum_k] += 1
        if maximum_k >= minimum_k:
            projected.append((tree, projection, maximum_k))
    if sample_size > len(projected):
        raise ValueError(
            f"sample size {sample_size} exceeds {len(projected)} eligible routes"
        )
    fully_projected = [
        item for item in projected if not item[1].unresolved_occurrence_ids
    ]
    anchor_pool = fully_projected or projected
    anchor = min(anchor_pool, key=lambda item: _stable_key(item[0].tree_id, seed))
    remaining = sorted(
        (item for item in projected if item[0].tree_id != anchor[0].tree_id),
        key=lambda item: item[0].tree_id,
    )
    selected = [anchor]
    if sample_size > 1:
        selected.extend(random.Random(seed).sample(remaining, sample_size - 1))
    cases = tuple(
        PartitionDatasetReviewCase(
            case_id=digest(
                "SPDATA1",
                tree.tree_id,
                projection.algorithm_version,
                str(seed),
            ),
            selection_rank=index,
            route_tree=tree,
            projection=projection,
            maximum_k=maximum_k,
            fully_projected=not projection.unresolved_occurrence_ids,
        )
        for index, (tree, projection, maximum_k) in enumerate(selected, start=1)
    )
    warnings = {"K_WAY_ENRICHED_SAMPLE_NOT_CORPUS_REPRESENTATIVE"}
    if any(not case.fully_projected for case in cases):
        warnings.add("SAMPLE_CONTAINS_PARTIALLY_PROJECTED_ROUTES")
    return PartitionDatasetReview(
        source_path=str(source_path),
        source_sha256=_file_sha256(source_path),
        source_route_count=source_count,
        projection_count_by_maximum_k=tuple(sorted(counts.items())),
        minimum_selected_k=minimum_k,
        eligible_route_count=len(projected),
        sample_size=sample_size,
        seed=seed,
        selection_method=(
            "fully_projected_anchor_then_seeded_random_without_replacement"
        ),
        cases=cases,
        warnings=tuple(sorted(warnings)),
    )


def _reaction_nodes(
    root: MoleculeOccurrenceNode,
) -> tuple[RouteReactionNode, ...]:
    values: list[RouteReactionNode] = []
    if root.reaction is not None:
        values.append(root.reaction)
        for child in root.reaction.children:
            values.extend(_reaction_nodes(child))
    return tuple(sorted(values, key=lambda item: (item.depth, item.step_id)))


def _warnings(values: Iterable[str]) -> str:
    items = tuple(values)
    if not items:
        return "<span class='muted'>none</span>"
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def _module_legend(frontier: Any) -> str:
    return "".join(
        "<span class='module-chip'>"
        f"<i style='background:{_MODULE_COLORS[index % len(_MODULE_COLORS)]}'></i>"
        f"M{index + 1} · {len(module.target_atom_maps)} target atoms"
        "</span>"
        for index, module in enumerate(frontier.partition.modules)
    )


def _latent_state_cards(frontier: Any) -> str:
    module_labels = {
        module.module_id: f"M{index}"
        for index, module in enumerate(frontier.partition.modules, start=1)
    }
    return "".join(
        "<article class='latent-card'>"
        f"<div class='latent-svg'>{molecule_svg(state.mapped_smiles, width=300, height=190)}</div>"
        f"<h4>{escape(' + '.join(module_labels.get(value, value) for value in state.module_ids) or 'partial module')}</h4>"
        f"<code>{escape(state.mapped_smiles)}</code>"
        f"<small>{len(state.target_atom_maps)} target atoms · "
        f"{len(state.non_target_atoms)} tactical/non-target atoms</small>"
        "</article>"
        for state in frontier.latent_states
    )


def _frontier_panel(frontier: Any) -> str:
    interfaces = "".join(
        "<tr>"
        f"<td>{escape(interface.interface_kind)}</td>"
        f"<td>{escape(interface.evidence_level)}</td>"
        f"<td>{escape('; '.join(f'{left}-{right}:{kind}' for left, right, kind in interface.target_bonds) or 'unary')}</td>"
        "</tr>"
        for interface in frontier.partition.interfaces
    )
    return (
        "<section class='frontier-panel'>"
        "<div class='frontier-heading'>"
        f"<div><span>Observed frontier depth {frontier.depth}</span>"
        f"<h3>k = {frontier.partition.k}</h3></div>"
        f"<code>{escape(frontier.partition.partition_id)}</code></div>"
        f"<div class='module-legend'>{_module_legend(frontier)}</div>"
        "<div class='partition-layout'>"
        f"<div class='partition-svg'>{render_partition_svg(frontier.partition)}</div>"
        f"<div class='latent-grid'>{_latent_state_cards(frontier)}</div>"
        "</div>"
        "<details><summary>Strategic interfaces and warnings</summary>"
        "<table><thead><tr><th>Interface</th><th>Evidence</th><th>Target bonds</th></tr></thead>"
        f"<tbody>{interfaces or '<tr><td colspan=3>none</td></tr>'}</tbody></table>"
        f"{_warnings((*frontier.warnings, *frontier.partition.warnings))}</details>"
        "</section>"
    )


def _case_panel(case: PartitionDatasetReviewCase) -> str:
    tree = case.route_tree
    unresolved = len(case.projection.unresolved_occurrence_ids)
    reactions = "".join(
        "<article class='reaction-card'>"
        f"<h4>Observed step · retrosynthetic depth {reaction.depth}</h4>"
        f"<div>{reaction_svg(reaction.reaction_smiles)}</div>"
        f"<code>{escape(reaction.reaction_smiles)}</code>"
        "</article>"
        for reaction in _reaction_nodes(tree.root)
    )
    return (
        f"<article class='case' data-case='{escape(case.case_id)}' "
        f"data-k='{case.maximum_k}' data-resolved='{str(case.fully_projected).lower()}'>"
        "<div class='case-heading'>"
        f"<div><span class='case-number'>Case {case.selection_rank}</span>"
        f"<h2>Observed route · maximum successfully projected k={case.maximum_k}</h2>"
        f"<p>{tree.reaction_count} steps · {len(case.projection.frontiers)} resolved frontiers · "
        f"{unresolved} unresolved occurrences</p></div>"
        f"<span class='status {'good' if case.fully_projected else 'caution'}'>"
        f"{'fully projected' if case.fully_projected else 'partial projection'}</span>"
        "</div>"
        "<section class='target-panel'><div>"
        f"{molecule_svg(tree.target_smiles, width=430, height=260)}</div><div>"
        f"<h3>Target</h3><code>{escape(tree.target_smiles)}</code>"
        f"<p>Tree <code>{escape(tree.tree_id)}</code><br>"
        f"Patent <code>{escape(tree.patent_id or 'unavailable')}</code> · "
        f"split <strong>{escape(tree.split or 'unknown')}</strong></p>"
        f"<div class='warning'>{_warnings(case.projection.warnings)}</div></div></section>"
        f"<div class='frontier-stack'>{''.join(_frontier_panel(item) for item in case.projection.frontiers)}</div>"
        "<details class='observed-route'><summary>Observed reaction sequence used for projection</summary>"
        f"{reactions}</details>"
        "<section class='review'><label>Is the deepest partition useful?"
        f"<select data-field='useful' data-case='{escape(case.case_id)}'>"
        "<option value=''>Select…</option><option>yes</option><option>partly</option>"
        "<option>no</option><option>cannot_assess</option></select></label>"
        "<label>Projection correctness"
        f"<select data-field='correctness' data-case='{escape(case.case_id)}'>"
        "<option value=''>Select…</option><option>correct</option><option>minor_issue</option>"
        "<option>incorrect</option><option>cannot_assess</option></select></label>"
        f"<label class='notes'>Notes<textarea data-field='notes' data-case='{escape(case.case_id)}'></textarea></label>"
        "</section></article>"
    )


def render_partition_dataset_review_html(review: PartitionDatasetReview) -> str:
    """Render a self-contained graphical K-way corpus review."""

    coverage = "".join(
        f"<div><strong>{count:,}</strong><span>routes with maximum successfully projected k={maximum_k}</span></div>"
        for maximum_k, count in review.projection_count_by_maximum_k
    )
    k_options = "".join(
        f'<option value="{maximum_k}">k={maximum_k}</option>'
        for maximum_k in sorted({case.maximum_k for case in review.cases})
    )
    payload = json.dumps(review.to_dict(), sort_keys=True, separators=(",", ":")).replace(
        "</", "<\\/"
    )
    storage_key = f"partition-dataset-review:{review.source_sha256[:16]}:{review.seed}"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>K-way observed-route fragmentation review</title><style>
:root{{--ink:#1f2d35;--muted:#687983;--paper:#eef3f5;--card:#fff;--line:#d5e0e5;--blue:#3178b9;--green:#267a52;--amber:#9b650d}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}}header{{background:linear-gradient(120deg,#163b52,#21665e);color:#fff;padding:28px max(20px,calc((100% - 1500px)/2))}}header h1{{margin:.2rem 0;font-size:30px}}header p{{max-width:1050px;color:#dcebef}}main{{max-width:1500px;margin:auto;padding:18px}}code{{overflow-wrap:anywhere;font-size:12px}}.metrics{{display:grid;grid-template-columns:repeat(4,minmax(130px,1fr));gap:8px;margin-top:18px}}.metrics div{{background:#ffffff17;border:1px solid #ffffff2b;border-radius:10px;padding:10px}}.metrics strong,.metrics span{{display:block}}.metrics strong{{font-size:22px}}.toolbar{{position:sticky;top:0;z-index:10;display:flex;gap:12px;align-items:end;background:#fff;border:1px solid var(--line);border-radius:10px;padding:10px;margin-bottom:18px;box-shadow:0 4px 18px #1737481c}}.toolbar label{{display:grid;gap:3px;color:var(--muted)}}select,textarea,button{{font:inherit;border:1px solid #aebdc5;border-radius:7px;background:#fff;padding:7px}}button{{background:#245c78;color:#fff;border:0;cursor:pointer}}#visible{{margin-left:auto}}.case{{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:18px;margin:0 0 22px;box-shadow:0 4px 16px #17374812}}.case-heading,.frontier-heading{{display:flex;justify-content:space-between;align-items:start;gap:12px}}.case-heading h2,.frontier-heading h3{{margin:.15rem 0}}.case-number,.frontier-heading span{{color:var(--blue);font-weight:750;text-transform:uppercase;letter-spacing:.06em}}.status{{border-radius:999px;padding:6px 10px;font-weight:700}}.status.good{{background:#e8f7ef;color:var(--green)}}.status.caution{{background:#fff3dc;color:var(--amber)}}.target-panel{{display:grid;grid-template-columns:minmax(360px,500px) 1fr;align-items:center;gap:18px;border:1px solid var(--line);border-radius:12px;padding:10px;margin:14px 0}}.target-panel svg{{width:100%;height:auto;max-height:270px}}.frontier-stack{{display:grid;gap:14px}}.frontier-panel{{border-left:5px solid var(--blue);background:#f8fbfc;border-radius:10px;padding:14px}}.module-legend{{display:flex;flex-wrap:wrap;gap:8px;margin:8px 0}}.module-chip{{background:#fff;border:1px solid var(--line);border-radius:999px;padding:5px 9px}}.module-chip i{{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px}}.partition-layout{{display:grid;grid-template-columns:minmax(440px,560px) 1fr;gap:14px;align-items:start}}.partition-svg{{background:#fff;border:1px solid var(--line);border-radius:9px;overflow:auto;text-align:center}}.partition-svg svg{{max-width:100%;height:auto}}.latent-grid{{display:grid;grid-template-columns:repeat(3,minmax(180px,1fr));gap:8px}}.latent-card{{background:#fff;border:1px solid var(--line);border-radius:9px;padding:8px;min-width:0}}.latent-card h4{{margin:3px 0}}.latent-card small{{display:block;color:var(--muted);margin-top:5px}}.latent-svg svg{{width:100%;height:auto;max-height:190px}}details{{margin-top:10px}}summary{{cursor:pointer;font-weight:700}}table{{width:100%;border-collapse:collapse;margin-top:8px}}th,td{{text-align:left;border:1px solid var(--line);padding:6px}}th{{background:#edf3f6}}.observed-route{{border-top:1px solid var(--line);padding-top:12px}}.reaction-card{{border:1px solid var(--line);border-radius:8px;padding:8px;margin:8px 0;overflow:auto}}.reaction-card svg{{max-width:100%;height:auto;max-height:210px}}.review{{display:grid;grid-template-columns:1fr 1fr;gap:10px;border-top:1px solid var(--line);margin-top:14px;padding-top:12px}}.review label{{display:grid;gap:4px;font-weight:700}}.review .notes{{grid-column:1/-1}}textarea{{min-height:80px;resize:vertical}}.muted{{color:var(--muted)}}.warning{{color:var(--amber)}}
@media(max-width:950px){{.partition-layout,.target-panel{{grid-template-columns:1fr}}.latent-grid{{grid-template-columns:repeat(2,1fr)}}.metrics{{grid-template-columns:repeat(2,1fr)}}.toolbar{{position:static;flex-wrap:wrap}}#visible{{width:100%;margin:0}}}}
</style></head><body><header><p>K-WAY SYNTHETIC PARTITION · OBSERVED ROUTE DATASET</p><h1>{review.sample_size}-route graphical fragmentation review</h1><p>This is an enriched functionality sample from observed routes with maximum successfully projected k≥{review.minimum_selected_k}; it is not an unbiased estimate of corpus accuracy. Module colors describe target-atom ownership, while separate structures show the corresponding latent precursor states.</p><div class="metrics"><div><strong>{review.source_route_count:,}</strong><span>source routes</span></div><div><strong>{review.eligible_route_count}</strong><span>eligible projected-k routes</span></div><div><strong>{review.sample_size}</strong><span>selected cases</span></div>{coverage}</div></header><main><div class="toolbar"><label>Projection status<select id="resolved"><option value="all">All</option><option value="true">Fully projected</option><option value="false">Partial projection</option></select></label><label>Maximum projected k<select id="max-k"><option value="all">All</option>{k_options}</select></label><button id="export">Export review JSON</button><span id="visible"></span></div><div class="warning">{_warnings(review.warnings)}</div>{''.join(_case_panel(case) for case in review.cases)}</main><script type="application/json" id="partition-dataset-review-data">{payload}</script><script>
const key={json.dumps(storage_key)},cards=[...document.querySelectorAll('.case')];let saved={{}};try{{saved=JSON.parse(localStorage.getItem(key)||'{{}}')}}catch(_){{saved={{}}}};
for(const el of document.querySelectorAll('[data-field]')){{const id=el.dataset.case,field=el.dataset.field;el.value=(saved[id]||{{}})[field]||'';el.addEventListener(el.tagName==='TEXTAREA'?'input':'change',()=>{{saved[id]=saved[id]||{{}};saved[id][field]=el.value;try{{localStorage.setItem(key,JSON.stringify(saved))}}catch(_){{}}}})}}
function filter(){{const resolved=document.querySelector('#resolved').value,k=document.querySelector('#max-k').value;let n=0;for(const card of cards){{const show=(resolved==='all'||card.dataset.resolved===resolved)&&(k==='all'||card.dataset.k===k);card.hidden=!show;if(show)n++}}document.querySelector('#visible').textContent=`${{n}} / ${{cards.length}} cases`}}document.querySelector('#resolved').addEventListener('change',filter);document.querySelector('#max-k').addEventListener('change',filter);document.querySelector('#export').addEventListener('click',()=>{{const out={{schema_version:'1.0',source_sha256:{json.dumps(review.source_sha256)},seed:{review.seed},reviews:cards.map(card=>({{case_id:card.dataset.case,...(saved[card.dataset.case]||{{}})}}))}};const blob=new Blob([JSON.stringify(out,null,2)],{{type:'application/json'}}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='k-way-partition-dataset-review.json';a.click();URL.revokeObjectURL(a.href)}});filter();
</script></body></html>"""


def write_partition_dataset_review(
    review: PartitionDatasetReview,
    json_path: str | Path,
    html_path: str | Path,
) -> None:
    """Write the deterministic evidence packet and graphical review."""

    output_json = Path(json_path)
    output_html = Path(html_path)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(review.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    output_html.write_text(
        render_partition_dataset_review_html(review),
        encoding="utf-8",
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a graphical K-way review from observed route trees"
    )
    parser.add_argument("source_route_trees")
    parser.add_argument("output_json")
    parser.add_argument("output_html")
    parser.add_argument("--sample-size", type=int, default=10)
    parser.add_argument("--seed", type=int, default=DEFAULT_PARTITION_DATASET_REVIEW_SEED)
    parser.add_argument("--minimum-k", type=int, default=3)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the observed-route K-way review builder."""

    arguments = _parser().parse_args(argv)
    review = build_partition_dataset_review(
        arguments.source_route_trees,
        sample_size=arguments.sample_size,
        seed=arguments.seed,
        minimum_k=arguments.minimum_k,
    )
    write_partition_dataset_review(
        review,
        arguments.output_json,
        arguments.output_html,
    )
    print(
        json.dumps(
            {
                "source_route_count": review.source_route_count,
                "eligible_route_count": review.eligible_route_count,
                "sample_size": review.sample_size,
                "output_json": str(Path(arguments.output_json).resolve()),
                "output_html": str(Path(arguments.output_html).resolve()),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "DEFAULT_PARTITION_DATASET_REVIEW_SEED",
    "PARTITION_DATASET_REVIEW_ALGORITHM_VERSION",
    "PARTITION_DATASET_REVIEW_SCHEMA_VERSION",
    "PartitionDatasetReview",
    "PartitionDatasetReviewCase",
    "build_partition_dataset_review",
    "render_partition_dataset_review_html",
    "write_partition_dataset_review",
]
