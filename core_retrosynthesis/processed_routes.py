"""Preserve released routes and qualify structure-derived observed route trees."""

from __future__ import annotations

from collections import deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import closing
from dataclasses import asdict
import gzip
import hashlib
from itertools import islice
from functools import lru_cache
import json
from pathlib import Path
import sqlite3
from typing import Any, Callable, Iterator, Mapping
import zlib

from condition_recommender.corpus_io import file_sha256
from .route_conversion import build_observed_route_tree, ObservedRouteConversionError
from .route_curation import _parse_reaction, _patent_id, _source_reaction_id, RouteQualityError

ROUTE_CATALOG_VERSION = "processed_routes.v1.1"


@lru_cache(maxsize=4096)
def _physical_step(source_id: str, reaction_smiles: str) -> Any:
    return _parse_reaction({"_id": f"physical_{source_id}", "reaction_smiles": reaction_smiles,
                            "abstracted_reaction_smiles": ""}, "physical")


def normalize_released_route(wrapper: Mapping[str, Any]) -> dict[str, Any]:
    """Reconstruct a unique connected tree, ignoring source order and weak abstractions.

    Repeated physical steps across abstraction subtrees are deduplicated only when
    their original IDs and mapped reactions agree. Conflicts remain unresolved.
    The original wrapper, including algorithmic labels, is retained by the catalog.
    """
    route = wrapper.get("raw_route") or {}
    route_id = route.get("route_id")
    patent_id = _patent_id(route_id)
    original = route.get("original_tree") or {}
    expected = original.get("reaction_ids")
    if not isinstance(expected, list) or not expected or len(set(expected)) != len(expected):
        raise RouteQualityError("invalid_original_reaction_ids")
    parsed: dict[str, Any] = {}
    for subtree in route.get("subtrees") or ():
        for raw in subtree.get("reactions") or ():
            # Algorithmic abstraction is not an observed physical reaction.
            source_id = _source_reaction_id(raw.get("_id"), subtree.get("subtree_id"))
            reaction = _physical_step(source_id, raw.get("reaction_smiles"))
            previous = parsed.get(reaction.source_reaction_id)
            if previous and previous.reaction_smiles != reaction.reaction_smiles:
                raise RouteQualityError("conflicting_original_reaction")
            parsed[reaction.source_reaction_id] = reaction
    if set(expected) != set(parsed) or original.get("num_reactions") != len(parsed):
        raise RouteQualityError("source_reaction_ids_do_not_match")
    products = [r.product_smiles for r in parsed.values()]
    if len(set(products)) != len(products):
        raise RouteQualityError("duplicate_step_product")
    consumed = {p for r in parsed.values() for p in r.precursor_smiles}
    roots = sorted(set(products) - consumed)
    if len(roots) != 1:
        raise RouteQualityError("ambiguous_route_target")
    steps = [asdict(parsed[key]) for key in sorted(parsed)]
    for step in steps:
        step["precursor_smiles"] = list(step["precursor_smiles"])
        step["abstraction_archived"] = True
    return {"schema_version": "released_route_structure_adapter.v1", "route_id": route_id,
            "patent_id": patent_id, "target_smiles": roots[0], "split": "unassigned",
            "steps": steps, "original_reaction_count": len(steps)}


def _encode(value: Any) -> bytes:
    return zlib.compress(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode())


def _convert_batch(rows: list[dict[str, Any]]) -> list[tuple[Any, ...]]:
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    result = []
    for wrapper in rows:
        raw = wrapper.get("raw_route") or {}
        route_id = str(raw.get("route_id") or "")
        if not route_id:
            raise ValueError("A source route has no identity")
        tree = None
        reason = None
        try:
            tree = build_observed_route_tree(normalize_released_route(wrapper)).to_dict()
        except (RouteQualityError, ObservedRouteConversionError) as error:
            reason = error.reason
        result.append((route_id, "validated_tree" if tree else "unresolved", reason,
                       _encode(wrapper), _encode(tree) if tree else None,
                       tuple((raw.get("original_tree") or {}).get("reaction_ids") or ())))
    return result


def _rows(path: Path) -> Iterator[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def build_processed_route_catalog(
    route_source: str | Path, step_source: str | Path, destination: str | Path, *,
    workers: int = 1, progress: Callable[[dict[str, Any]], None] | None = None,
) -> dict[str, Any]:
    """Build a resumable indexed catalog containing every route and step membership.

    Qualification failures preserve both raw memberships and a stable reason.
    Missing or ambiguous original-step joins remain review evidence and prevent
    those routes entering the qualified-tree stream. No train/test scope is invented.
    """
    routes, steps, output = Path(route_source).resolve(), Path(step_source).resolve(), Path(destination).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    binding = hashlib.sha256((file_sha256(routes) + file_sha256(steps) + ROUTE_CATALOG_VERSION).encode()).hexdigest()
    if output.is_file():
        with closing(sqlite3.connect(output.as_uri() + "?mode=ro", uri=True)) as existing:
            report = json.loads(existing.execute("SELECT payload FROM metadata").fetchone()[0])
        if report.get("source_binding") == binding and report.get("build_complete"):
            return report
        raise ValueError("Existing route catalog belongs to different sources or definitions")
    checkpoint = output.with_suffix(".building.sqlite")
    with closing(sqlite3.connect(checkpoint)) as db:
        db.executescript("""
            PRAGMA journal_mode=WAL;
            PRAGMA synchronous=NORMAL;
            CREATE TABLE IF NOT EXISTS state(id INTEGER PRIMARY KEY,binding TEXT,position INTEGER,steps_complete INTEGER);
            CREATE TABLE IF NOT EXISTS step_observations(reaction_id TEXT,observation_id TEXT,
                PRIMARY KEY(reaction_id,observation_id)) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS routes(id TEXT PRIMARY KEY,status TEXT,reason TEXT,source BLOB,tree BLOB) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS route_steps(route_id TEXT,reaction_id TEXT,
                PRIMARY KEY(route_id,reaction_id)) WITHOUT ROWID;
            CREATE INDEX IF NOT EXISTS routes_by_step ON route_steps(reaction_id,route_id);
        """)
        state = db.execute("SELECT binding,position,steps_complete FROM state WHERE id=1").fetchone()
        if state and state[0] != binding:
            raise ValueError("Route checkpoint has different sources or definitions")
        if not state:
            db.execute("INSERT INTO state VALUES (1,?,0,0)", (binding,))
            db.commit()
            state = (binding, 0, 0)
        if not state[2]:
            for row in _rows(steps):
                source = row["source"]
                db.execute("INSERT OR IGNORE INTO step_observations VALUES (?,?)",
                           (source["source_record_id"], row["observation_id"]))
            db.execute("UPDATE state SET steps_complete=1 WHERE id=1")
            db.commit()
        position = state[1]
        rows = iter(islice(_rows(routes), position, None))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            pending: deque[Any] = deque()
            exhausted = False
            while pending or not exhausted:
                while not exhausted and len(pending) < max(1, 2 * workers):
                    batch = list(islice(rows, 500))
                    if not batch:
                        exhausted = True
                        break
                    pending.append(pool.submit(_convert_batch, batch))
                if not pending:
                    break
                converted = pending.popleft().result()
                for route_id, status, reason, source, tree, ids in converted:
                    db.execute("INSERT INTO routes VALUES (?,?,?,?,?)", (route_id, status, reason, source, tree))
                    db.executemany("INSERT OR IGNORE INTO route_steps VALUES (?,?)", ((route_id, str(i)) for i in ids))
                position += len(converted)
                db.execute("UPDATE state SET position=? WHERE id=1", (position,))
                db.commit()
                if progress:
                    progress({"phase": "routes", "routes": position})
        missing = db.execute("SELECT count(*) FROM route_steps s WHERE NOT EXISTS (SELECT 1 FROM step_observations o WHERE o.reaction_id=s.reaction_id)").fetchone()[0]
        ambiguous = db.execute("SELECT count(*) FROM (SELECT reaction_id FROM step_observations GROUP BY reaction_id HAVING count(*)>1)").fetchone()[0]
        if missing:
            db.execute("""UPDATE routes SET status='unresolved',reason='missing_original_step_observation'
                          WHERE status='validated_tree' AND id IN
                          (SELECT s.route_id FROM route_steps s WHERE NOT EXISTS
                           (SELECT 1 FROM step_observations o WHERE o.reaction_id=s.reaction_id))""")
        if ambiguous:
            db.execute("""UPDATE routes SET status='unresolved',reason='ambiguous_original_step_observation'
                          WHERE status='validated_tree' AND id IN
                          (SELECT s.route_id FROM route_steps s JOIN
                           (SELECT reaction_id FROM step_observations GROUP BY reaction_id HAVING count(*)>1) o
                           ON o.reaction_id=s.reaction_id)""")
        counts = dict(db.execute("SELECT status,count(*) FROM routes GROUP BY status"))
        reasons = dict(db.execute("SELECT reason,count(*) FROM routes WHERE reason IS NOT NULL GROUP BY reason"))
        report = {"schema_version": ROUTE_CATALOG_VERSION, "build_complete": True,
                  "source_binding": binding, "route_count": position, "counts": counts,
                  "unresolved_reasons": reasons, "missing_memberships": missing,
                  "ambiguous_step_aliases": ambiguous,
                  "composite_promotion": "unavailable_no_declared_training_scope_or_independent_support_review"}
        db.execute("CREATE TABLE IF NOT EXISTS metadata(payload TEXT NOT NULL)")
        db.execute("DELETE FROM metadata")
        db.execute("INSERT INTO metadata VALUES (?)", (json.dumps(report, sort_keys=True),))
        db.commit()
        db.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        db.execute("PRAGMA journal_mode=DELETE")
    checkpoint.replace(output)
    return report


class ProcessedRouteCatalog:
    """Indexed route source, membership and qualified-tree access."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path).resolve()

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path.as_uri() + "?mode=ro&immutable=1", uri=True)

    def route(self, identity: str, *, include_tree: bool = False,
              include_source: bool = False) -> dict[str, Any]:
        """Fetch one route with exact original-step aliases and observation joins."""
        with closing(self._connect()) as db:
            metadata = json.loads(db.execute("SELECT payload FROM metadata").fetchone()[0])
            if metadata.get("schema_version") != ROUTE_CATALOG_VERSION or not metadata.get("build_complete"):
                raise ValueError("Route catalog is incompatible or incomplete")
            row = db.execute("SELECT status,reason,source,tree FROM routes WHERE id=?", (identity,)).fetchone()
            if row is None:
                raise KeyError(identity)
            memberships = []
            for (reaction_id,) in db.execute("SELECT reaction_id FROM route_steps WHERE route_id=? ORDER BY reaction_id", (identity,)):
                ids = [r[0] for r in db.execute("SELECT observation_id FROM step_observations WHERE reaction_id=? ORDER BY observation_id", (reaction_id,))]
                memberships.append({"reaction_id": reaction_id, "observation_ids": ids,
                                    "join_status": "resolved" if len(ids) == 1 else "missing" if not ids else "ambiguous"})
        result = {"route_id": identity, "status": row[0], "unresolved_reason": row[1],
                  "memberships": memberships, "tree_available": row[3] is not None,
                  "connectivity_scope": "inferred_from_source_molecular_identity"}
        if include_tree and row[3] is not None:
            result["tree"] = json.loads(zlib.decompress(row[3]))
        if include_source:
            result["source"] = json.loads(zlib.decompress(row[2]))
        result["omitted_sections"] = [name for name, included in (("tree", include_tree), ("source", include_source)) if not included]
        return result

    def iter_trees(self) -> Iterator[Any]:
        """Stream qualified trees through the existing validated typed contract."""
        from .route_contract import ReactionRouteTree
        with closing(self._connect()) as db:
            metadata = json.loads(db.execute("SELECT payload FROM metadata").fetchone()[0])
            if metadata.get("schema_version") != ROUTE_CATALOG_VERSION or not metadata.get("build_complete"):
                raise ValueError("Route catalog is incompatible or incomplete")
            for (payload,) in db.execute("SELECT tree FROM routes WHERE status='validated_tree' ORDER BY id"):
                yield ReactionRouteTree.from_dict(json.loads(zlib.decompress(payload)))


def main() -> None:
    """Build route evidence independently while canonical chemistry converts."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default="datasets/intermediate_datasets/routes/higher_level_retrosynthesis")
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    root = Path(args.source)
    print(json.dumps(build_processed_route_catalog(root / "routes.source.jsonl.gz",
          root / "route_steps.observations.jsonl.gz", args.output, workers=args.workers,
          progress=lambda e: print(json.dumps(e), flush=True))), flush=True)


if __name__ == "__main__":
    main()
