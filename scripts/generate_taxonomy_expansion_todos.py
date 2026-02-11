#!/usr/bin/env python3
"""
Generate patch-ready taxonomy expansion TODOs from discovery clusters.

Input:
  JSON produced by scripts/discover_reaction_coverage_gaps.py

Outputs:
  - JSON with ranked TODO entries and patch stubs
  - Markdown triage document for human curation

Example:
  python scripts/generate_taxonomy_expansion_todos.py \
    --input results/reaction_coverage_discovery.json \
    --output-json results/taxonomy_expansion_todos.json \
    --output-markdown results/taxonomy_expansion_todos.md
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.taxonomy.reaction_catalog import load_reaction_catalog


DEFAULT_DISCOVERY_INPUT = "results/reaction_coverage_discovery.json"
DEFAULT_TODOS_JSON = "results/taxonomy_expansion_todos.json"
DEFAULT_TODOS_MD = "results/taxonomy_expansion_todos.md"
DEFAULT_ORGANIC_COMPOUNDS = (
    REPO_ROOT / "chemtools" / "taxonomy" / "data" / "organic_compounds.v1.3.json"
)
DEFAULT_SCAFFOLD_MOTIFS = (
    REPO_ROOT / "chemtools" / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"
)
DEFAULT_ORGANIC_GROUPS = (
    REPO_ROOT / "chemtools" / "taxonomy" / "data" / "organic_groups.v1.3.json"
)


def _to_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except Exception:
        return default


def _to_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _pair_list_to_dict(value: Any) -> Dict[str, int]:
    out: Dict[str, int] = {}
    if not isinstance(value, list):
        return out
    for item in value:
        if (
            isinstance(item, (list, tuple))
            and len(item) == 2
            and str(item[0]).strip()
        ):
            out[str(item[0]).strip()] = _to_int(item[1], 0)
    return out


def _split_signature_side(text: str) -> Set[str]:
    value = str(text or "").strip()
    if not value or value == "none":
        return set()
    return {
        token.strip()
        for token in value.split("|")
        if token and token.strip() and token.strip() != "none"
    }


def _parse_motif_signature(signature: str) -> Tuple[Set[str], Set[str]]:
    value = str(signature or "").strip()
    if "->" not in value:
        return set(), set()
    left, right = value.split("->", 1)
    return _split_signature_side(left), _split_signature_side(right)


def _is_unclassified(value: Any) -> bool:
    text = str(value or "").strip()
    if not text:
        return True
    low = text.lower()
    return low == "unknown" or text.startswith("Unclassified-")


def _generated_compound_id(a_group: str, b_group: str) -> str:
    if b_group.startswith("-"):
        return f"{a_group}{b_group}"
    return f"{a_group}-{b_group}"


def _load_json(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload if isinstance(payload, dict) else {}


def _load_compound_ids(paths: Sequence[Path]) -> Set[str]:
    motif_ids: Set[str] = set()
    for path in paths:
        if not path.exists():
            continue
        payload = _load_json(path)
        compounds = payload.get("compounds")
        if not isinstance(compounds, list):
            continue
        for entry in compounds:
            if not isinstance(entry, dict):
                continue
            comp_id = str(entry.get("id") or "").strip()
            if comp_id:
                motif_ids.add(comp_id)
                continue
            a_group = str(entry.get("A") or "").strip()
            b_group = str(entry.get("B") or "").strip()
            if a_group and b_group:
                motif_ids.add(_generated_compound_id(a_group, b_group))
    return motif_ids


def _load_group_ids(path: Path) -> Set[str]:
    if not path.exists():
        return set()
    payload = _load_json(path)
    groups = payload.get("groups")
    if not isinstance(groups, list):
        return set()
    out: Set[str] = set()
    for item in groups:
        if not isinstance(item, dict):
            continue
        group_id = str(item.get("id") or "").strip()
        if group_id:
            out.add(group_id)
    return out


def load_taxonomy_context(
    *,
    compounds_path: Path = DEFAULT_ORGANIC_COMPOUNDS,
    scaffold_path: Path = DEFAULT_SCAFFOLD_MOTIFS,
    groups_path: Path = DEFAULT_ORGANIC_GROUPS,
) -> Dict[str, Any]:
    definitions, _ = load_reaction_catalog()
    reacted_universe: Set[str] = set()
    formed_universe: Set[str] = set()
    for defn in definitions.values():
        for slot in defn.reactants.values():
            reacted_universe.update(str(v) for v in (slot.allowed or []) if str(v).strip())
        for slot in defn.products.values():
            formed_universe.update(str(v) for v in (slot.allowed or []) if str(v).strip())

    return {
        "reaction_definitions": definitions,
        "reacted_slot_motifs": reacted_universe,
        "formed_slot_motifs": formed_universe,
        "compound_ids": _load_compound_ids([compounds_path, scaffold_path]),
        "group_ids": _load_group_ids(groups_path),
    }


def _slugify_token(text: str, *, default: str = "cluster") -> str:
    value = str(text or "").strip().lower()
    if not value:
        return default
    value = re.sub(r"[^a-z0-9]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or default


def _make_family_stub_id(
    *,
    cluster_index: int,
    event_signature: str,
    reacted: Iterable[str],
    formed: Iterable[str],
) -> str:
    event_part = _slugify_token(event_signature, default="unknown_event")
    reacted_part = _slugify_token("_".join(sorted(reacted)[:2]), default="reacted")
    formed_part = _slugify_token("_".join(sorted(formed)[:2]), default="formed")
    return f"todo_family_{cluster_index:03d}_{event_part}_{reacted_part}_{formed_part}"


def _infer_compound_stub(
    motif_id: str,
    *,
    group_ids: Set[str],
) -> Dict[str, Any]:
    motif = str(motif_id or "").strip()
    if not motif:
        return {
            "id": "",
            "todo": "empty_motif_id",
            "notes": "motif id is empty; inspect source rows",
        }

    if "-" not in motif:
        return {
            "id": motif,
            "todo": "define_compound_id_explicitly",
            "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json",
        }

    left, right = motif.split("-", 1)
    a_group = left.strip()
    b_payload = right.strip()
    b_group = b_payload if b_payload.startswith("-") else f"-{b_payload}"

    notes: List[str] = []
    if a_group and a_group not in group_ids:
        notes.append(f"missing_group:{a_group}")
    if b_group and b_group not in group_ids:
        notes.append(f"missing_group:{b_group}")

    return {
        "template": "single_bond",
        "A": a_group,
        "B": b_group,
        "description": f"TODO: inferred from unresolved cluster for motif '{motif}'.",
        "notes": notes,
    }


def _reaction_slot_additions(
    *,
    reaction_id: str,
    reacted_missing: List[str],
    formed_missing: List[str],
) -> Dict[str, Any]:
    return {
        "reaction_id": reaction_id,
        "reactant_slot_additions_todo": {
            "electrophile_or_substrate": reacted_missing,
        },
        "product_slot_additions_todo": {
            "product": formed_missing,
        },
        "constraints_review_todo": [
            "review include/exclude reacted constraints",
            "review include/exclude formed constraints",
            "confirm min_reactant_slot_matches / min_product_slot_matches",
        ],
    }


def _rank_priority(
    *,
    cluster_count: int,
    reason_counts: Dict[str, int],
) -> int:
    score = int(cluster_count)
    score += 3 * int(reason_counts.get("unknown_reaction_type", 0) > 0)
    score += 2 * int(reason_counts.get("motif_outside_reaction_taxonomy", 0) > 0)
    score += 2 * int(reason_counts.get("unclassified_motif", 0) > 0)
    score += 1 * int(reason_counts.get("low_reaction_key_quality", 0) > 0)
    score += 1 * int(reason_counts.get("empty_reaction_key", 0) > 0)
    return score


def build_todo_payload(
    discovery_payload: Dict[str, Any],
    *,
    taxonomy_context: Dict[str, Any],
    top_clusters: int = 40,
    top_motifs_per_cluster: int = 8,
) -> Dict[str, Any]:
    clusters = discovery_payload.get("clusters")
    if not isinstance(clusters, list):
        clusters = []
    clusters = list(clusters)[: max(1, int(top_clusters))]

    reacted_slot_motifs = set(taxonomy_context.get("reacted_slot_motifs") or set())
    formed_slot_motifs = set(taxonomy_context.get("formed_slot_motifs") or set())
    compound_ids = set(taxonomy_context.get("compound_ids") or set())
    group_ids = set(taxonomy_context.get("group_ids") or set())

    todos: List[Dict[str, Any]] = []
    for index, cluster in enumerate(clusters, start=1):
        if not isinstance(cluster, dict):
            continue
        motif_signature = str(cluster.get("motif_signature") or "")
        event_signature = str(cluster.get("event_signature") or "")
        reacted, formed = _parse_motif_signature(motif_signature)
        reason_counts = _pair_list_to_dict(cluster.get("reasons"))
        count = _to_int(cluster.get("count"), 0)
        top_candidates = cluster.get("top_taxonomy_candidates")
        if not isinstance(top_candidates, list):
            top_candidates = []
        top_candidate = top_candidates[0] if top_candidates else {}
        top_reaction_id = str((top_candidate or {}).get("reaction_id") or "").strip()
        top_score = _to_float((top_candidate or {}).get("score"), 0.0)

        reacted_missing = sorted(
            [
                motif
                for motif in reacted
                if not _is_unclassified(motif) and motif not in reacted_slot_motifs
            ]
        )
        formed_missing = sorted(
            [
                motif
                for motif in formed
                if not _is_unclassified(motif) and motif not in formed_slot_motifs
            ]
        )
        unresolved_motifs = sorted(set(reacted_missing + formed_missing))

        existing_compounds_missing_slots = sorted(
            [motif for motif in unresolved_motifs if motif in compound_ids]
        )
        unknown_compounds = sorted(
            [motif for motif in unresolved_motifs if motif not in compound_ids]
        )

        if top_reaction_id and top_score >= 0.70:
            action_track = "expand_existing_reaction_type"
        elif top_reaction_id and top_score >= 0.45:
            action_track = "review_constraints_or_product_patterns"
        else:
            action_track = "propose_new_reaction_family"

        proposed_reaction_patch = None
        proposed_family_stub = None
        if action_track in {"expand_existing_reaction_type", "review_constraints_or_product_patterns"} and top_reaction_id:
            proposed_reaction_patch = _reaction_slot_additions(
                reaction_id=top_reaction_id,
                reacted_missing=reacted_missing[: max(1, int(top_motifs_per_cluster))],
                formed_missing=formed_missing[: max(1, int(top_motifs_per_cluster))],
            )
        else:
            proposed_family_stub = {
                "id": _make_family_stub_id(
                    cluster_index=index,
                    event_signature=event_signature,
                    reacted=reacted,
                    formed=formed,
                ),
                "aliases": [],
                "description": "TODO: new reaction family candidate from unresolved cluster.",
                "reactants": {
                    "electrophile_or_substrate": sorted(list(reacted))[: max(1, int(top_motifs_per_cluster))],
                    "nucleophile_or_partner": [],
                },
                "products": {
                    "product": sorted(list(formed))[: max(1, int(top_motifs_per_cluster))],
                },
                "constraints": {
                    "include_reacted": [],
                    "exclude_reacted": [],
                    "include_formed": [],
                    "exclude_formed": [],
                    "min_reactant_slot_matches": 0,
                    "min_product_slot_matches": 0,
                },
            }

        compound_stubs = [
            _infer_compound_stub(motif, group_ids=group_ids)
            for motif in unknown_compounds[: max(1, int(top_motifs_per_cluster))]
        ]

        todo_entry = {
            "rank": index,
            "priority_score": _rank_priority(cluster_count=count, reason_counts=reason_counts),
            "cluster_id": cluster.get("cluster_id"),
            "cluster_count": count,
            "motif_signature": motif_signature,
            "event_signature": event_signature,
            "recommended_action_from_discovery": cluster.get("recommended_action"),
            "action_track": action_track,
            "reason_counts": reason_counts,
            "top_candidate": top_candidate,
            "gap_analysis": {
                "missing_reacted_slot_motifs": reacted_missing,
                "missing_formed_slot_motifs": formed_missing,
                "existing_compounds_missing_slots": existing_compounds_missing_slots,
                "unknown_compound_ids": unknown_compounds,
            },
            "patch_stubs": {
                "reaction_type_update": proposed_reaction_patch,
                "new_reaction_family": proposed_family_stub,
                "new_compound_entries": compound_stubs,
            },
            "samples": list(cluster.get("samples") or []),
        }
        todos.append(todo_entry)

    todos.sort(
        key=lambda item: (
            int(item.get("priority_score", 0)),
            int(item.get("cluster_count", 0)),
            -int(item.get("rank", 0)),
        ),
        reverse=True,
    )
    for idx, item in enumerate(todos, start=1):
        item["priority_rank"] = idx

    summary = {
        "input_clusters": len(clusters),
        "todo_entries": len(todos),
        "tracks": {
            "expand_existing_reaction_type": sum(1 for row in todos if row.get("action_track") == "expand_existing_reaction_type"),
            "review_constraints_or_product_patterns": sum(
                1 for row in todos if row.get("action_track") == "review_constraints_or_product_patterns"
            ),
            "propose_new_reaction_family": sum(1 for row in todos if row.get("action_track") == "propose_new_reaction_family"),
        },
    }
    return {"summary": summary, "todos": todos}


def write_markdown(path: Path, payload: Dict[str, Any]) -> None:
    summary = payload.get("summary") if isinstance(payload, dict) else {}
    todos = payload.get("todos") if isinstance(payload, dict) else []
    if not isinstance(summary, dict):
        summary = {}
    if not isinstance(todos, list):
        todos = []

    lines: List[str] = []
    lines.append("# Taxonomy Expansion TODOs")
    lines.append("")
    lines.append("Generated from unresolved discovery clusters.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Input clusters: {summary.get('input_clusters', 0)}")
    lines.append(f"- TODO entries: {summary.get('todo_entries', 0)}")
    tracks = summary.get("tracks")
    if isinstance(tracks, dict):
        lines.append(f"- Expand existing reaction type: {tracks.get('expand_existing_reaction_type', 0)}")
        lines.append(
            f"- Review constraints/product patterns: {tracks.get('review_constraints_or_product_patterns', 0)}"
        )
        lines.append(f"- Propose new reaction family: {tracks.get('propose_new_reaction_family', 0)}")
    lines.append("")

    for todo in todos:
        if not isinstance(todo, dict):
            continue
        lines.append(f"## Priority {todo.get('priority_rank', '?')}: Cluster {todo.get('cluster_id')}")
        lines.append("")
        lines.append(f"- Cluster count: {todo.get('cluster_count', 0)}")
        lines.append(f"- Priority score: {todo.get('priority_score', 0)}")
        lines.append(f"- Action track: `{todo.get('action_track', '')}`")
        lines.append(f"- Motif signature: `{todo.get('motif_signature', '')}`")
        lines.append(f"- Event signature: `{todo.get('event_signature', '')}`")
        top_candidate = todo.get("top_candidate")
        if isinstance(top_candidate, dict) and top_candidate:
            lines.append(
                f"- Top taxonomy candidate: `{top_candidate.get('reaction_id', '')}` (score={top_candidate.get('score', 0)})"
            )
        reasons = todo.get("reason_counts")
        if isinstance(reasons, dict) and reasons:
            reason_text = ", ".join(f"{k}:{v}" for k, v in sorted(reasons.items()))
            lines.append(f"- Reasons: {reason_text}")

        gap = todo.get("gap_analysis")
        if isinstance(gap, dict):
            lines.append("")
            lines.append("### Gap Analysis")
            lines.append("")
            lines.append(
                f"- Missing reacted-slot motifs: {', '.join(gap.get('missing_reacted_slot_motifs') or []) or '(none)'}"
            )
            lines.append(
                f"- Missing formed-slot motifs: {', '.join(gap.get('missing_formed_slot_motifs') or []) or '(none)'}"
            )
            lines.append(
                f"- Existing compound IDs missing in reaction slots: {', '.join(gap.get('existing_compounds_missing_slots') or []) or '(none)'}"
            )
            lines.append(
                f"- Unknown compound IDs: {', '.join(gap.get('unknown_compound_ids') or []) or '(none)'}"
            )

        stubs = todo.get("patch_stubs")
        if isinstance(stubs, dict):
            lines.append("")
            lines.append("### Patch Stubs")
            lines.append("")
            for key in ("reaction_type_update", "new_reaction_family", "new_compound_entries"):
                value = stubs.get(key)
                if not value:
                    continue
                lines.append(f"- `{key}`")
                lines.append("```json")
                lines.append(json.dumps(value, indent=2, ensure_ascii=True))
                lines.append("```")

        lines.append("")
        lines.append("### Validation Checklist")
        lines.append("")
        lines.append("- Confirm representative reactions are true single-step transforms.")
        lines.append("- Verify motifs against product/reactant structures with RDKit inspection.")
        lines.append("- If updating existing family, ensure constraints still disambiguate close families.")
        lines.append("- Run taxonomy validators and reaction reliability report after edits.")
        lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate taxonomy expansion TODOs from discovery output")
    parser.add_argument(
        "--input",
        default=DEFAULT_DISCOVERY_INPUT,
        help=f"Discovery JSON input (default: {DEFAULT_DISCOVERY_INPUT})",
    )
    parser.add_argument(
        "--top-clusters",
        type=int,
        default=40,
        help="Top N clusters to convert to TODOs",
    )
    parser.add_argument(
        "--top-motifs-per-cluster",
        type=int,
        default=8,
        help="Max unresolved motifs kept per cluster in stubs",
    )
    parser.add_argument(
        "--output-json",
        default=DEFAULT_TODOS_JSON,
        help=f"Output TODO JSON path (default: {DEFAULT_TODOS_JSON})",
    )
    parser.add_argument(
        "--output-markdown",
        default=DEFAULT_TODOS_MD,
        help=f"Output TODO markdown path (default: {DEFAULT_TODOS_MD})",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input discovery JSON not found: {input_path}")

    discovery_payload = _load_json(input_path)
    taxonomy_context = load_taxonomy_context()
    todos_payload = build_todo_payload(
        discovery_payload,
        taxonomy_context=taxonomy_context,
        top_clusters=max(1, int(args.top_clusters)),
        top_motifs_per_cluster=max(1, int(args.top_motifs_per_cluster)),
    )

    out_json = Path(args.output_json).expanduser().resolve()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(todos_payload, indent=2, sort_keys=True), encoding="utf-8")

    out_md = Path(args.output_markdown).expanduser().resolve()
    write_markdown(out_md, todos_payload)

    summary = todos_payload.get("summary") or {}
    print(f"TODO JSON: {out_json}")
    print(f"TODO Markdown: {out_md}")
    print(f"Input clusters: {summary.get('input_clusters', 0)}")
    print(f"TODO entries: {summary.get('todo_entries', 0)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
