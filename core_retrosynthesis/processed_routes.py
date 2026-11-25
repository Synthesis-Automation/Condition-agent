"""Preserve released routes and qualify structure-derived observed route trees."""

from __future__ import annotations

from collections import defaultdict, deque
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
from typing import Any, Callable, Iterator, Mapping, Sequence
import zlib

from condition_recommender.corpus_io import file_sha256
from condition_recommender.ingestion.cbyd import CBYD_ROUTE_SCHEMA_VERSION
from reactive_taxonomy.chemistry.source_mapping import assess_source_mapping
from .route_conversion import build_observed_route_tree, ObservedRouteConversionError
from .route_curation import (
    _canonical_unmapped_smiles, _parse_reaction, _patent_id, _source_reaction_id,
    RouteQualityError,
)

LEGACY_ROUTE_CATALOG_VERSION = "processed_routes.v1.1"
ROUTE_CATALOG_VERSION = "processed_routes.v2"


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
    if wrapper.get("schema_version") == CBYD_ROUTE_SCHEMA_VERSION:
        return _normalize_literature_route(wrapper)
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


def _normalize_literature_route(wrapper: Mapping[str, Any]) -> dict[str, Any]:
    """Verify mapped physical steps and exact molecule identities against publisher edges."""
    source = wrapper["source"]
    rows = wrapper["raw_route"]["synthesis"]
    memberships = wrapper["step_memberships"]
    graph = wrapper.get("publisher_graph")
    if len(rows) != len(memberships) or not rows:
        raise RouteQualityError("source_membership_count_mismatch")
    if not isinstance(graph, dict) or not all(isinstance(graph.get(k), list) for k in ("rx", "nodes", "edges")):
        raise RouteQualityError("invalid_publisher_graph")
    if not all(isinstance(item, dict) for item in graph["rx"]) or [
        (item.get("rid"), item.get("m")) for item in graph["rx"]
    ] != [
        (row.get("reaction_id"), row.get("reaction_smiles_mapped")) for row in rows
    ]:
        raise RouteQualityError("source_graph_reaction_mismatch")
    steps = []
    for row, member in zip(rows, memberships, strict=True):
        step = asdict(_physical_step(member["source_record_id"], row.get("reaction_smiles_mapped")))
        step["precursor_smiles"] = list(step["precursor_smiles"])
        steps.append(step)
    for row in rows:
        agreement = assess_source_mapping(row.get("reaction_smiles") or "", row["reaction_smiles_mapped"])
        if agreement.issues:
            raise RouteQualityError(agreement.issues[0])
    for node in graph["nodes"]:
        if (not isinstance(node, dict) or not isinstance(node.get("id"), str)
            or node.get("t") not in ("m", "r")
            or (node["t"] == "r" and type(node.get("i")) is not int)
            or (node["t"] == "m" and not isinstance(node.get("s"), str))):
            raise RouteQualityError("invalid_publisher_graph_node")
    if graph.get("warn") or any(node.get("inv") for node in graph["nodes"]):
        raise RouteQualityError("publisher_graph_flagged")
    nodes = {node["id"]: node for node in graph["nodes"]}
    targets = [node for node in nodes.values() if node.get("r") == "target" and node["t"] == "m"]
    if len(targets) != 1:
        raise RouteQualityError("ambiguous_route_target")
    reaction_nodes = [node for node in nodes.values() if node["t"] == "r"]
    if sorted(node["i"] for node in reaction_nodes) != list(range(len(steps))):
        raise RouteQualityError("source_graph_step_coverage_mismatch")
    incoming: dict[str, set[str]] = defaultdict(set)
    outgoing: dict[str, set[str]] = defaultdict(set)
    for edge in graph["edges"]:
        if not isinstance(edge, list) or len(edge) != 2 or not all(isinstance(key, str) for key in edge):
            raise RouteQualityError("invalid_publisher_graph_edge")
        left, right = edge
        if left not in nodes or right not in nodes or {nodes[left]["t"], nodes[right]["t"]} != {"r", "m"}:
            raise RouteQualityError("invalid_publisher_graph_edge")
        if right in outgoing[left]:
            raise RouteQualityError("publisher_graph_repeated_edge")
        incoming[right].add(left)
        outgoing[left].add(right)
    # The same molecule can appear at several places in a publisher drawing.
    # Connectivity requires edges to the actual producing occurrences.
    products = {step["product_smiles"]: i for i, step in enumerate(steps)}
    if len(products) != len(steps):
        raise RouteQualityError("duplicate_step_product")
    pending = [targets[0]["id"]]
    reachable = set()
    while pending:
        node_id = pending.pop()
        if node_id not in reachable:
            reachable.add(node_id)
            pending.extend(incoming[node_id])
    if len(reachable) != len(graph["nodes"]) or sum(map(len, outgoing.values())) != len(nodes) - 1:
        raise RouteQualityError("publisher_graph_not_single_tree")
    for node in reaction_nodes:
        step = steps[node["i"]]
        node_id = node["id"]
        producer_ids = {nodes[r]["i"] for m in incoming[node_id] for r in incoming[m]}
        expected = {products[p] for p in step["precursor_smiles"] if p in products}
        if producer_ids != expected:
            raise RouteQualityError("source_dependency_mismatch")
        product_nodes = outgoing[node_id]
        if len(product_nodes) != 1 or _canonical_unmapped_smiles(nodes[next(iter(product_nodes))]["s"]) != step["product_smiles"]:
            raise RouteQualityError("source_graph_product_mismatch")
        precursors = sorted(_canonical_unmapped_smiles(nodes[m]["s"]) for m in incoming[node_id])
        if precursors != step["precursor_smiles"]:
            raise RouteQualityError("source_graph_precursor_mismatch")
    return {"schema_version": "literature_route_structure_adapter.v1",
            "source_kind": "literature", "route_id": wrapper["route_id"],
            "source_dataset_id": f"{source['corpus_id']}:{source['release_id']}",
            "source_reference": source.get("reference") or source["source_file"],
            "target_smiles": _canonical_unmapped_smiles(targets[0]["s"]),
            "split": "unassigned", "steps": steps, "original_reaction_count": len(steps)}


def _encode(value: Any) -> bytes:
    return zlib.compress(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode())


def _convert_batch(rows: list[dict[str, Any]]) -> list[tuple[Any, ...]]:
    from rdkit import rdBase
    with rdBase.BlockLogs():
        return [_convert_wrapper(wrapper) for wrapper in rows]


def _convert_wrapper(wrapper: dict[str, Any]) -> tuple[Any, ...]:
    raw = wrapper.get("raw_route") or {}
    literature = wrapper.get("schema_version") == CBYD_ROUTE_SCHEMA_VERSION
    route_id = str((wrapper.get("route_id") if literature else raw.get("route_id")) or "")
    if not route_id:
        raise ValueError("A source route has no identity")
    tree = None
    reason = None
    try:
        tree = build_observed_route_tree(normalize_released_route(wrapper)).to_dict()
    except (RouteQualityError, ObservedRouteConversionError) as error:
        reason = error.reason
    if literature:
        members = wrapper["step_memberships"]
        if (len(members) != len(raw["synthesis"])
            or any(member["source_row_number"] != pos or
                   member["source_reaction_id"] != row.get("reaction_id")
                   for pos, (member, row) in enumerate(zip(members, raw["synthesis"]), 1))):
            raise ValueError(f"Invalid source memberships in {route_id}")
        memberships = [(m["source_record_id"], m["source_row_number"], m["observation_id"])
                       for m in members]
    else:
        memberships = [(i, None, None) for i in dict.fromkeys((raw.get("original_tree") or {}).get("reaction_ids") or ())]
    return (route_id, "validated_tree" if tree else "unresolved", reason,
            _encode(wrapper), _encode(tree) if tree else None,
            memberships, int(literature and isinstance(wrapper.get("publisher_graph"), dict)))


def _rows(path: Path) -> Iterator[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _converted_batches(
    routes: Iterator[dict[str, Any]], workers: int,
) -> Iterator[list[tuple[Any, ...]]]:
    batches = iter(lambda: list(islice(routes, 128)), [])
    if workers == 1:
        for batch in batches:
            yield _convert_batch(batch)
        return
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending: deque[Any] = deque()
        for batch in batches:
            pending.append(pool.submit(_convert_batch, batch))
            if len(pending) >= 2 * workers:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()


def _verified_legacy_catalog(release_root: str | Path) -> tuple[Path, str] | None:
    """Verify the release checksum before importing unchanged patent route payloads."""
    from condition_recommender.processed_release import resolve_processed_release

    release = resolve_processed_release(release_root)
    if "route_catalog" not in release.manifest["artifacts"]:
        return None
    path = release.artifact("route_catalog")
    if file_sha256(path) != release.manifest["artifacts"]["route_catalog"]["sha256"]:
        raise ValueError("Reusable route catalog checksum mismatch")
    with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as db:
        metadata = json.loads(db.execute("SELECT payload FROM metadata").fetchone()[0])
    if metadata["schema_version"] == LEGACY_ROUTE_CATALOG_VERSION and metadata.get("build_complete"):
        return path, metadata["source_binding"]
    return None


def _copy_legacy_bundle(db: sqlite3.Connection, path: Path, bundle: int) -> None:
    db.execute("ATTACH DATABASE ? AS legacy", (path.as_uri() + "?mode=ro&immutable=1",))
    with db:
        db.execute("INSERT INTO routes SELECT id,status,reason,source,tree,?,0 FROM legacy.routes", (bundle,))
        db.execute("INSERT INTO step_observations SELECT ?,reaction_id,observation_id FROM legacy.step_observations", (bundle,))
        db.execute("INSERT INTO route_steps SELECT route_id,reaction_id,NULL,NULL FROM legacy.route_steps")
        db.execute("""UPDATE bundle_state SET position=(SELECT count(*) FROM legacy.routes),
                      steps=(SELECT count(*) FROM legacy.step_observations),done=1 WHERE bundle=?""", (bundle,))
    db.execute("DETACH DATABASE legacy")


def build_processed_route_catalog(
    route_source: str | Path, step_source: str | Path, destination: str | Path, *,
    workers: int = 1, progress: Callable[[dict[str, Any]], None] | None = None,
) -> dict[str, Any]:
    """Build the catalog for a single route and observation pair."""
    return build_route_bundle_catalog([(route_source, step_source)], destination,
                                      workers=workers, progress=progress)


def build_route_bundle_catalog(
    bundles: Sequence[tuple[str | Path, str | Path]], destination: str | Path, *,
    workers: int = 1, progress: Callable[[dict[str, Any]], None] | None = None,
    reuse_release: str | Path | None = None, evidence_catalog: str | Path | None = None,
) -> dict[str, Any]:
    """Build one atomic, resumable catalog over ordered, independently scoped route sources."""
    if not bundles or workers < 1:
        raise ValueError("Route catalog requires at least one bundle and one worker")
    output = Path(destination).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    pairs = [(Path(routes).resolve(), Path(steps).resolve()) for routes, steps in bundles]
    hashes = [(file_sha256(routes), file_sha256(steps)) for routes, steps in pairs]
    inputs = [ROUTE_CATALOG_VERSION, hashes]
    if evidence_catalog is not None:
        from reactive_taxonomy.chemistry.source_mapping import SOURCE_MAPPING_SCREEN_VERSION

        inputs.extend((SOURCE_MAPPING_SCREEN_VERSION, file_sha256(evidence_catalog)))
    binding = hashlib.sha256(json.dumps(inputs).encode()).hexdigest()
    if output.is_file():
        with closing(sqlite3.connect(output.as_uri() + "?mode=ro", uri=True)) as existing:
            report = json.loads(existing.execute("SELECT payload FROM metadata").fetchone()[0])
        if report.get("source_binding") == binding and report.get("build_complete"):
            return report
        raise ValueError("Existing route catalog belongs to different sources or definitions")
    legacy = _verified_legacy_catalog(reuse_release) if reuse_release else None
    checkpoint = output.with_suffix(".building.sqlite")
    with closing(sqlite3.connect(checkpoint, uri=True)) as db:
        db.executescript("""
            PRAGMA journal_mode=WAL;
            PRAGMA synchronous=NORMAL;
            CREATE TABLE IF NOT EXISTS state(id INTEGER PRIMARY KEY,binding TEXT);
            CREATE TABLE IF NOT EXISTS bundle_state(bundle INTEGER PRIMARY KEY,steps INTEGER,position INTEGER,done INTEGER);
            CREATE TABLE IF NOT EXISTS step_observations(bundle INTEGER,reaction_id TEXT,observation_id TEXT,
                PRIMARY KEY(bundle,reaction_id,observation_id)) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS routes(id TEXT PRIMARY KEY,status TEXT,reason TEXT,source BLOB,tree BLOB,
                bundle INTEGER,publisher_graph INTEGER) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS route_steps(route_id TEXT,reaction_id TEXT,position INTEGER,expected_obs TEXT,
                PRIMARY KEY(route_id,reaction_id)) WITHOUT ROWID;
        """)
        state = db.execute("SELECT binding FROM state WHERE id=1").fetchone()
        if state and state[0] != binding:
            raise ValueError("Route checkpoint has different sources or definitions")
        db.execute("INSERT OR IGNORE INTO state VALUES (1,?)", (binding,))
        db.commit()
        for bundle, ((routes, steps), (route_hash, step_hash)) in enumerate(zip(pairs, hashes, strict=True)):
            db.execute("INSERT OR IGNORE INTO bundle_state VALUES (?,0,0,0)", (bundle,))
            db.commit()
            step_count, position, done = db.execute(
                "SELECT steps,position,done FROM bundle_state WHERE bundle=?", (bundle,)).fetchone()
            if done:
                continue
            legacy_binding = hashlib.sha256((route_hash + step_hash + LEGACY_ROUTE_CATALOG_VERSION).encode()).hexdigest()
            if legacy and legacy[1] == legacy_binding and position == step_count == 0:
                if progress:
                    progress({"phase": "importing_verified_routes", "bundle": bundle})
                _copy_legacy_bundle(db, legacy[0], bundle)
                if progress:
                    progress({"phase": "reused_verified_routes", "bundle": bundle})
                continue
            remaining = islice(_rows(steps), step_count, None)
            for batch in iter(lambda: list(islice(remaining, 500)), []):
                with db:
                    db.executemany("INSERT INTO step_observations VALUES (?,?,?)",
                                   ((bundle, row["source"]["source_record_id"], row["observation_id"]) for row in batch))
                    step_count += len(batch)
                    db.execute("UPDATE bundle_state SET steps=? WHERE bundle=?", (step_count, bundle))
            for converted in _converted_batches(islice(_rows(routes), position, None), workers):
                with db:
                    for route_id, status, reason, source, tree, members, graph in converted:
                        db.execute("INSERT INTO routes VALUES (?,?,?,?,?,?,?)",
                                   (route_id, status, reason, source, tree, bundle, graph))
                        db.executemany("INSERT INTO route_steps VALUES (?,?,?,?)",
                                       ((route_id, *member) for member in members))
                    position += len(converted)
                    db.execute("UPDATE bundle_state SET position=? WHERE bundle=?", (position, bundle))
                if progress:
                    progress({"phase": "routes", "bundle": bundle, "routes": position})
            with db:
                db.execute("UPDATE bundle_state SET done=1 WHERE bundle=?", (bundle,))
        db.executescript("""
            CREATE INDEX IF NOT EXISTS routes_by_step ON route_steps(reaction_id,route_id);
            CREATE INDEX IF NOT EXISTS steps_by_observation ON step_observations(observation_id);
            CREATE INDEX IF NOT EXISTS routes_by_bundle ON routes(bundle,status);
            CREATE TEMP TABLE join_issues AS
            SELECT s.route_id, s.reaction_id, CASE WHEN count(o.observation_id)=0 THEN 'missing_original_step_observation'
                 WHEN count(o.observation_id)>1 THEN 'ambiguous_original_step_observation'
                 ELSE 'conflicting_expected_observation' END AS reason
            FROM route_steps s JOIN routes r ON r.id=s.route_id LEFT JOIN step_observations o
              ON o.bundle=r.bundle AND o.reaction_id=s.reaction_id
            GROUP BY s.route_id,s.reaction_id
            HAVING count(o.observation_id)<>1 OR (s.expected_obs IS NOT NULL AND min(o.observation_id)<>s.expected_obs);
            CREATE INDEX join_issues_by_route ON join_issues(route_id);
            UPDATE routes SET status='unresolved',reason=(SELECT min(reason) FROM join_issues j WHERE j.route_id=routes.id)
              WHERE status='validated_tree' AND id IN (SELECT route_id FROM join_issues);
        """)
        missing = db.execute("SELECT count(*) FROM join_issues WHERE reason='missing_original_step_observation'").fetchone()[0]
        conflicts = db.execute("SELECT count(*) FROM join_issues WHERE reason='conflicting_expected_observation'").fetchone()[0]
        ambiguous = db.execute("""SELECT count(*) FROM (SELECT 1 FROM step_observations
                                  GROUP BY bundle,reaction_id HAVING count(*)>1)""").fetchone()[0]
        counts = dict(db.execute("SELECT status,count(*) FROM routes GROUP BY status"))
        reasons = dict(db.execute("SELECT reason,count(*) FROM routes WHERE reason IS NOT NULL GROUP BY reason"))
        bundle_reports = []
        for bundle, (route_hash, step_hash) in enumerate(hashes):
            nsteps, nroute = db.execute("SELECT steps,position FROM bundle_state WHERE bundle=?", (bundle,)).fetchone()
            states = dict(db.execute("SELECT status,count(*) FROM routes WHERE bundle=? GROUP BY status", (bundle,)))
            bundle_reports.append({"route_sha256": route_hash, "step_sha256": step_hash, "route_count": nroute,
                                   "observation_count": nsteps, "counts": states})
        report = {"schema_version": ROUTE_CATALOG_VERSION, "build_complete": True,
                  "source_binding": binding, "route_count": sum(counts.values()), "counts": counts,
                  "bundles": bundle_reports, "unresolved_reasons": reasons, "missing_memberships": missing,
                  "ambiguous_step_aliases": ambiguous, "conflicting_memberships": conflicts,
                  "composite_promotion": "unavailable_no_declared_training_scope_or_independent_support_review"}
        if evidence_catalog is not None:
            from .route_step_screen import apply_route_step_screen

            apply_route_step_screen(db, Path(evidence_catalog), report)
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

    @staticmethod
    def _version(db: sqlite3.Connection) -> str:
        report = json.loads(db.execute("SELECT payload FROM metadata").fetchone()[0])
        version = report.get("schema_version")
        if version not in {ROUTE_CATALOG_VERSION, LEGACY_ROUTE_CATALOG_VERSION} or not report.get("build_complete"):
            raise ValueError("Route catalog is incompatible or incomplete")
        return version

    def route(self, identity: str, *, include_tree: bool = False,
              include_source: bool = False) -> dict[str, Any]:
        """Fetch one route with exact original-step aliases and observation joins."""
        with closing(self._connect()) as db:
            legacy = self._version(db) == LEGACY_ROUTE_CATALOG_VERSION
            columns = "0,0" if legacy else "bundle,publisher_graph"
            row = db.execute(f"SELECT status,reason,source,tree,{columns} FROM routes WHERE id=?", (identity,)).fetchone()
            if row is None:
                raise KeyError(identity)
            memberships = []
            order = "reaction_id" if legacy else "position,reaction_id"
            details = "NULL,NULL" if legacy else "position,expected_obs"
            for reaction_id, position, expected in db.execute(
                f"SELECT reaction_id,{details} FROM route_steps WHERE route_id=? ORDER BY {order}", (identity,)
            ):
                where = "reaction_id=?" if legacy else "reaction_id=? AND bundle=?"
                params = (reaction_id,) if legacy else (reaction_id, row[4])
                ids = [r[0] for r in db.execute(f"SELECT observation_id FROM step_observations WHERE {where} ORDER BY observation_id", params)]
                status = "resolved" if len(ids) == 1 else "missing" if not ids else "ambiguous"
                if status == "resolved" and expected and ids[0] != expected:
                    status = "conflicting_expected_observation"
                member = {"reaction_id": reaction_id, "observation_ids": ids, "join_status": status}
                if position is not None:
                    member.update(source_row_number=position, expected_observation_id=expected)
                memberships.append(member)
            if not legacy and db.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='route_step_screen'"
            ).fetchone():
                screened = {oid: (disposition, json.loads(reasons)) for oid, disposition, reasons in db.execute(
                    """SELECT x.observation_id,x.disposition,x.reasons FROM route_steps s JOIN route_step_screen x
                       ON s.reaction_id=x.reaction_id AND s.expected_obs=x.observation_id WHERE s.route_id=?""",
                    (identity,),
                )}
                for member in memberships:
                    screen = screened.get(member.get("expected_observation_id"))
                    if screen:
                        member.update(structural_disposition=screen[0], structural_exclusion_reasons=screen[1])
        result = {"route_id": identity, "status": row[0], "unresolved_reason": row[1],
                  "memberships": memberships, "tree_available": row[3] is not None,
                  "publisher_graph_available": bool(row[5]),
                  "connectivity_scope": ("publisher_inferred_from_reported_synthesis" if row[5]
                                         else "inferred_from_source_molecular_identity")}
        if include_tree and row[3] is not None:
            result["tree"] = json.loads(zlib.decompress(row[3]))
        if include_source:
            result["source"] = json.loads(zlib.decompress(row[2]))
        result["omitted_sections"] = [name for name, included in (("tree", include_tree), ("source", include_source)) if not included]
        return result

    def routes_for_observation(self, observation_id: str) -> list[str]:
        """Look up route memberships for an exact intermediate observation identity."""
        with closing(self._connect()) as db:
            version = self._version(db)
            scope = "" if version == LEGACY_ROUTE_CATALOG_VERSION else (
                "AND o.bundle=r.bundle AND (s.expected_obs IS NULL OR s.expected_obs=o.observation_id)"
            )
            rows = db.execute(f"""SELECT DISTINCT r.id FROM step_observations o
                JOIN route_steps s ON s.reaction_id=o.reaction_id JOIN routes r ON r.id=s.route_id
                WHERE o.observation_id=? {scope} ORDER BY r.id""", (observation_id,))
            return [row[0] for row in rows]

    def iter_trees(self) -> Iterator[Any]:
        """Stream qualified trees through the existing validated typed contract."""
        from .route_contract import ReactionRouteTree
        with closing(self._connect()) as db:
            self._version(db)
            for (payload,) in db.execute("SELECT tree FROM routes WHERE status='validated_tree' ORDER BY id"):
                yield ReactionRouteTree.from_dict(json.loads(zlib.decompress(payload)))


def main() -> None:
    """Build route evidence independently while canonical chemistry converts."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    from condition_recommender.processed_build import discover_route_bundles

    parser.add_argument("--source", default="datasets/intermediate_datasets")
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--reuse-release", help="Reuse checksum-verified, unchanged patent routes from a processed release")
    args = parser.parse_args()
    root = Path(args.source)
    print(json.dumps(build_route_bundle_catalog(discover_route_bundles(root),
          args.output, workers=args.workers, reuse_release=args.reuse_release,
          progress=lambda e: print(json.dumps(e), flush=True))), flush=True)


if __name__ == "__main__":
    main()
