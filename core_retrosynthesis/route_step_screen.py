"""Attach source-level structural dispositions to exact released route occurrences."""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import sqlite3
from typing import Any

from condition_recommender.processed_catalog import ProcessedCatalog
from reactive_taxonomy.chemistry.source_mapping import SOURCE_MAPPING_SCREEN_VERSION, source_mapping_conflicts


def apply_route_step_screen(db: sqlite3.Connection, evidence_catalog: Path, report: dict[str, Any]) -> None:
    """Annotate individual CbyD memberships, preserving every original step and edge."""
    # Occurrence IDs, unlike the source publisher's reaction IDs, are unique across routes.
    members = db.execute(
        "SELECT expected_obs, reaction_id FROM route_steps WHERE route_id GLOB 'cbyd:*'"
    ).fetchall()
    expected = dict(members)
    if not expected:
        return
    if None in expected or len(members) != len(expected):
        raise ValueError("CbyD route occurrences require unique expected observation IDs")
    db.execute("""CREATE TABLE IF NOT EXISTS route_step_screen(
        observation_id TEXT PRIMARY KEY, reaction_id TEXT NOT NULL UNIQUE,
        disposition TEXT NOT NULL, reasons TEXT NOT NULL) WITHOUT ROWID""")
    db.execute("DELETE FROM route_step_screen")
    counts: Counter[str] = Counter()
    for oid, rid, reasons in ProcessedCatalog(evidence_catalog).source_step_dispositions("cbyd"):
        if expected.pop(oid, None) != rid:
            raise ValueError(f"Canonical CbyD observation does not match its route occurrence: {oid}")
        disposition = "structural_precedent"
        if reasons:
            disposition = "quarantined" if source_mapping_conflicts(reasons) else "source_lookup_only"
        db.execute("INSERT INTO route_step_screen VALUES (?,?,?,?)",
                   (oid, rid, disposition, json.dumps(reasons)))
        counts[disposition] += 1
    if expected:
        raise ValueError(f"Missing canonical observations for {len(expected)} CbyD route occurrences")
    # A complete tree previously accepted by the graph validator can be withdrawn
    # when a step lacks a usable signature. A gap never creates a new edge.
    db.execute("""UPDATE routes SET status='unresolved',reason='structural_step_excluded',tree=NULL
                  WHERE status='validated_tree' AND id IN (
                      SELECT s.route_id FROM route_steps s JOIN route_step_screen x
                      ON x.reaction_id=s.reaction_id AND x.observation_id=s.expected_obs
                      WHERE x.disposition != 'structural_precedent')""")
    report["step_screen"] = {"definition_version": SOURCE_MAPPING_SCREEN_VERSION,
                             "dataset": "cbyd", "counts": dict(counts)}
    report["counts"] = dict(db.execute("SELECT status,count(*) FROM routes GROUP BY status"))
    report["unresolved_reasons"] = dict(db.execute(
        "SELECT reason,count(*) FROM routes WHERE reason IS NOT NULL GROUP BY reason"))
    for bundle, item in enumerate(report["bundles"]):
        item["counts"] = dict(db.execute(
            "SELECT status,count(*) FROM routes WHERE bundle=? GROUP BY status", (bundle,)))
