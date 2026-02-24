#!/usr/bin/env python3
"""
Audit overlap and drift across taxonomy + retrosynthesis JSON knowledge stores.

This script is intentionally read-only. It reports:
1. Which retro/template/transformation `taxonomy_id` values are canonical.
2. Which are alias-resolvable against `reaction_types.v4.0.json`.
3. Which are orphaned (not canonical and not alias-resolvable).
4. Redundant-field signals (e.g. retron `reaction_name` vs `taxonomy_id` drift).

Usage:
    python scripts/audit_taxonomy_retro_overlap.py
    python scripts/audit_taxonomy_retro_overlap.py --json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

# Allow `python scripts/...` execution from repo root without installing package.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.taxonomy import loader as taxonomy_loader
from chemtools.taxonomy.reaction_catalog import load_reaction_catalog, resolve_reaction_type


def _reaction_ids() -> set[str]:
    payload = taxonomy_loader.load_reaction_types_raw()
    return {
        str(entry.get("id") or "").strip()
        for entry in (payload.get("reaction_types") or [])
        if isinstance(entry, dict) and str(entry.get("id") or "").strip()
    }


def _classify_taxonomy_link(
    label: Any,
    canonical_ids: set[str],
) -> tuple[str, Optional[str]]:
    text = str(label or "").strip()
    if not text:
        return "missing", None
    if text in canonical_ids:
        return "canonical", text
    resolved = resolve_reaction_type(text)
    if resolved:
        return "alias_resolvable", resolved
    return "orphan", None


def _audit_items(
    items: Iterable[Mapping[str, Any]],
    *,
    key_name: str,
    label_name: str,
    canonical_ids: set[str],
) -> Dict[str, Any]:
    counts: Counter[str] = Counter()
    samples: Dict[str, List[Dict[str, Any]]] = {
        "alias_resolvable": [],
        "orphan": [],
        "missing": [],
    }
    unique_labels: set[str] = set()
    unique_canonical: set[str] = set()

    for item in items:
        label = str(item.get(label_name) or "").strip()
        status, resolved = _classify_taxonomy_link(label, canonical_ids)
        counts[status] += 1
        if label:
            unique_labels.add(label)
        if resolved:
            unique_canonical.add(resolved)

        if status in samples and len(samples[status]) < 12:
            samples[status].append(
                {
                    "name": str(item.get(key_name) or ""),
                    label_name: label,
                    "canonical": resolved,
                }
            )

    return {
        "counts": dict(counts),
        "unique_labels": len(unique_labels),
        "unique_canonical_resolved": len(unique_canonical),
        "samples": samples,
    }


def _audit_retron_reaction_name_drift(retrons: Iterable[Mapping[str, Any]]) -> Dict[str, Any]:
    counts: Counter[str] = Counter()
    examples: List[Dict[str, str]] = []

    for retron in retrons:
        rn = str(retron.get("reaction_name") or "").strip()
        tid = str(retron.get("taxonomy_id") or "").strip()
        if not rn and not tid:
            counts["both_missing"] += 1
            continue
        if rn == tid:
            counts["identical"] += 1
        elif rn and tid and rn.lower() == tid.lower():
            counts["case_only"] += 1
        else:
            counts["different_semantics"] += 1
            if len(examples) < 12:
                examples.append(
                    {
                        "name": str(retron.get("name") or ""),
                        "reaction_name": rn,
                        "taxonomy_id": tid,
                    }
                )

    return {"counts": dict(counts), "examples": examples}


def _audit_transformation_patterns(
    canonical_ids: set[str],
) -> Dict[str, Any]:
    payload = taxonomy_loader.load_transformation_patterns()
    patterns = payload.get("reaction_patterns") if isinstance(payload, dict) else {}
    if not isinstance(patterns, dict):
        return {"counts": {}, "samples": {"alias_resolvable": [], "orphan": [], "missing": []}}

    records = [{"name": key, "taxonomy_id": key} for key in patterns.keys()]
    return _audit_items(records, key_name="name", label_name="taxonomy_id", canonical_ids=canonical_ids)


def build_report() -> Dict[str, Any]:
    canonical_ids = _reaction_ids()
    reaction_defs, _ = load_reaction_catalog()
    retron_payload = taxonomy_loader.load_retron_patterns()
    hte_payload = taxonomy_loader.load_hte_templates()

    retrons = retron_payload.get("retrons") if isinstance(retron_payload, dict) else []
    templates = hte_payload.get("templates") if isinstance(hte_payload, dict) else []
    retrons = retrons if isinstance(retrons, list) else []
    templates = templates if isinstance(templates, list) else []

    retron_audit = _audit_items(
        (item for item in retrons if isinstance(item, dict)),
        key_name="name",
        label_name="taxonomy_id",
        canonical_ids=canonical_ids,
    )
    template_audit = _audit_items(
        (item for item in templates if isinstance(item, dict)),
        key_name="name",
        label_name="taxonomy_id",
        canonical_ids=canonical_ids,
    )
    transform_audit = _audit_transformation_patterns(canonical_ids)

    retron_tids = {
        str(item.get("taxonomy_id") or "").strip()
        for item in retrons
        if isinstance(item, dict) and str(item.get("taxonomy_id") or "").strip()
    }
    template_tids = {
        str(item.get("taxonomy_id") or "").strip()
        for item in templates
        if isinstance(item, dict) and str(item.get("taxonomy_id") or "").strip()
    }

    overlap = {
        "retrons_only": sorted(retron_tids - template_tids)[:30],
        "templates_only": sorted(template_tids - retron_tids)[:30],
        "shared_count": len(retron_tids & template_tids),
        "retrons_only_count": len(retron_tids - template_tids),
        "templates_only_count": len(template_tids - retron_tids),
    }

    return {
        "counts": {
            "reaction_types": len(reaction_defs),
            "retrons": len(retrons),
            "hte_templates": len(templates),
            "transformation_pattern_reactions": len(
                (taxonomy_loader.load_transformation_patterns().get("reaction_patterns") or {})
                if isinstance(taxonomy_loader.load_transformation_patterns(), dict)
                else {}
            ),
        },
        "taxonomy_link_status": {
            "retron_patterns": retron_audit,
            "hte_templates": template_audit,
            "transformation_patterns": transform_audit,
        },
        "retron_reaction_name_vs_taxonomy_id": _audit_retron_reaction_name_drift(
            (item for item in retrons if isinstance(item, dict))
        ),
        "taxonomy_id_overlap_between_retro_sources": overlap,
        "cleanup_recommendations": [
            "Keep `reaction_types.v4.0.json` as canonical reaction-family IDs and aliases.",
            "Treat retro JSON `taxonomy_id` as a link field; allow aliases temporarily, but audit and normalize periodically.",
            "Rename retron `reaction_name` to `retro_transform_id` in a future schema revision to avoid confusion with taxonomy IDs.",
            "Move retro-only metadata ownership into a dedicated namespace (e.g. `taxonomy/data/retro/*.json`) or a manifest to clarify boundaries.",
            "Derive display category/description from taxonomy where possible; keep retro-specific descriptions only when they add disconnection detail.",
        ],
    }


def _print_human(report: Mapping[str, Any]) -> None:
    counts = report.get("counts") or {}
    print("Taxonomy / Retro Overlap Audit")
    print("=" * 36)
    print(f"Reaction taxonomy entries : {counts.get('reaction_types', 0)}")
    print(f"Retron entries            : {counts.get('retrons', 0)}")
    print(f"HTE template entries      : {counts.get('hte_templates', 0)}")
    print(f"Transformation mappings   : {counts.get('transformation_pattern_reactions', 0)}")

    print("\nTaxonomy Link Status")
    print("-" * 20)
    link_status = report.get("taxonomy_link_status") or {}
    for source_name in ("retron_patterns", "hte_templates", "transformation_patterns"):
        src = link_status.get(source_name) or {}
        src_counts = src.get("counts") or {}
        print(
            f"{source_name}: canonical={src_counts.get('canonical', 0)}, "
            f"alias_resolvable={src_counts.get('alias_resolvable', 0)}, "
            f"orphan={src_counts.get('orphan', 0)}, missing={src_counts.get('missing', 0)}"
        )
        samples = src.get("samples") or {}
        alias_samples = samples.get("alias_resolvable") or []
        orphan_samples = samples.get("orphan") or []
        if alias_samples:
            print("  alias-resolvable samples:")
            for item in alias_samples[:5]:
                print(
                    f"    - {item.get('name')}: {item.get('taxonomy_id')} -> {item.get('canonical')}"
                )
        if orphan_samples:
            print("  orphan samples:")
            for item in orphan_samples[:5]:
                print(f"    - {item.get('name')}: {item.get('taxonomy_id')}")

    rn_drift = report.get("retron_reaction_name_vs_taxonomy_id") or {}
    rn_counts = rn_drift.get("counts") or {}
    print("\nRetron `reaction_name` vs `taxonomy_id`")
    print("-" * 34)
    print(
        f"identical={rn_counts.get('identical', 0)}, "
        f"case_only={rn_counts.get('case_only', 0)}, "
        f"different_semantics={rn_counts.get('different_semantics', 0)}"
    )
    for item in (rn_drift.get("examples") or [])[:6]:
        print(
            f"  - {item.get('name')}: reaction_name={item.get('reaction_name')} | "
            f"taxonomy_id={item.get('taxonomy_id')}"
        )

    overlap = report.get("taxonomy_id_overlap_between_retro_sources") or {}
    print("\nRetro Source Overlap")
    print("-" * 19)
    print(f"shared taxonomy_id labels : {overlap.get('shared_count', 0)}")
    print(f"retrons-only labels       : {overlap.get('retrons_only_count', 0)}")
    print(f"templates-only labels     : {overlap.get('templates_only_count', 0)}")

    print("\nRecommended Simplification")
    print("-" * 26)
    for item in report.get("cleanup_recommendations") or []:
        print(f"- {item}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON report.")
    args = parser.parse_args()

    report = build_report()
    if args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False))
    else:
        _print_human(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
