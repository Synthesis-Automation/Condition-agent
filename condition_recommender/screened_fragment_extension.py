"""Extend a published fragment index using a screened, separately attributable corpus."""

from __future__ import annotations

from collections import Counter
from contextlib import closing
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
from tempfile import TemporaryDirectory
from typing import Any

from rdkit import Chem

from reactive_taxonomy.chemistry.source_mapping import SOURCE_MAPPING_SCREEN_VERSION
from .conversion.atomic import atomic_output_path
from .corpus_io import file_sha256
from .fragment_index import _finish_library, load_fragment_library, open_fragment_index, unpack
from .processed_catalog import CATALOG_SCHEMA_VERSION, ProcessedCatalog
from .record_storage import json_bytes


def _screen_catalog(source: Path, destination: Path) -> tuple[set[str], dict[str, Any]]:
    catalog = ProcessedCatalog(source)
    metadata = catalog._metadata().copy()
    if metadata["schema_version"] != CATALOG_SCHEMA_VERSION:
        raise ValueError("The supplemental evidence catalog requires indexed reaction graphs")
    excluded = []
    accepted: set[str] = set()
    reasons: Counter[str] = Counter()
    for identity, _, issues in catalog.source_step_dispositions():
        if issues:
            excluded.append(identity)
            reasons.update(issues)
        else:
            accepted.add(identity)
    with atomic_output_path(destination) as temporary:
        shutil.copyfile(source, temporary)
        with closing(sqlite3.connect(temporary)) as db:
            db.execute("CREATE TEMP TABLE excluded(id TEXT PRIMARY KEY) WITHOUT ROWID")
            db.executemany("INSERT INTO excluded VALUES (?)", ((i,) for i in excluded))
            db.execute("DELETE FROM experiment_lookup WHERE observation_id IN excluded")
            counts = metadata["counts"].copy()
            for name, where in (("experiments_with_graph", ""), ("experiments_with_local_core",
                                 "WHERE kind <> 'whole_reaction'")):
                counts[name] = db.execute(
                    f"SELECT COUNT(DISTINCT observation_id) FROM experiment_lookup {where}"
                ).fetchone()[0]
            counts["source_screen_excluded"] = len(excluded)
            metadata["input_catalog_id"] = metadata["catalog_id"]
            metadata["catalog_id"] = hashlib.sha256(json_bytes(
                [metadata["input_catalog_id"], SOURCE_MAPPING_SCREEN_VERSION]
            )).hexdigest()
            metadata.update(counts=counts, source_mapping_screen=SOURCE_MAPPING_SCREEN_VERSION,
                            source_screen_issues=dict(sorted(reasons.items())))
            db.execute("UPDATE metadata SET payload=?", (json.dumps(metadata, sort_keys=True),))
            db.commit()
    return accepted, metadata


def extend_fragment_index(
    base: str | Path, supplement: str | Path, output_directory: str | Path,
) -> dict[str, Any]:
    """Publish screened evidence and a joint search index in a new output directory.

    Existing molecular projections and verified RDKit library entries in the
    base remain usable. Supplement evidence must have the same projection
    policy, distinct observation IDs and a bound canonical source catalog.
    """
    base, supplement, directory = (Path(p).resolve() for p in (base, supplement, output_directory))
    output = directory / "fragment_index.sqlite"
    evidence = directory / "catalogs.sqlite"
    if directory.exists():
        raise FileExistsError(f"Choose a new directory for the screened index: {directory}")
    manifests = []
    catalogs = []
    for path in (base, supplement):
        source, manifest = open_fragment_index(path)
        source.close()
        dependency = manifest.get("evidence_catalog")
        if not dependency or not manifest["source_coverage_complete"] or "supplemental_evidence_catalog" in manifest:
            raise ValueError("Expected two complete standalone fragment indexes with evidence catalogs")
        catalogs.append((path.parent / dependency["relative_path"]).resolve())
        manifests.append(manifest)
    if catalogs[0] == catalogs[1] or base == supplement:
        raise ValueError("Base and supplemental corpora must be distinct")
    inputs = {path: path.stat() for path in (base, supplement, *catalogs)}
    directory.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".screened-fragments-", dir=directory.parent) as staging:
        stage = Path(staging) / "artifacts"
        stage.mkdir()
        screened, report = _screen_catalog(catalogs[1], stage / evidence.name)
        merged = _merge_components(base, supplement, stage / output.name, screened, manifests, report,
                                   catalogs[0], stage / evidence.name, directory)
        for path, expected in inputs.items():
            stat = path.stat()
            if (stat.st_ino, stat.st_size, stat.st_mtime_ns) != (
                expected.st_ino, expected.st_size, expected.st_mtime_ns
            ):
                raise ValueError(f"Source changed while building the screened index: {path}")
        stage.rename(directory)
    return merged


def _merge_components(
    base: Path, supplement: Path, destination: Path, accepted: set[str],
    manifests: list[dict[str, Any]], evidence_report: dict[str, Any],
    base_catalog: Path, supplemental_catalog: Path, publish_directory: Path,
) -> dict[str, Any]:
    primary = manifests[0]
    shutil.copyfile(base, destination)
    with closing(sqlite3.connect(destination.as_uri() + "?mode=rw", uri=True)) as db:
        db.execute("PRAGMA journal_mode=OFF")
        db.execute("PRAGMA synchronous=OFF")
        db.execute("ATTACH DATABASE ? AS extra", (supplement.as_uri() + "?mode=ro",))
        if db.execute("SELECT 1 FROM extra.observations s JOIN main.observations o ON o.id=s.id LIMIT 1").fetchone():
            raise ValueError("Fragment sources have overlapping observation IDs")
        db.execute("CREATE TABLE supplemental_observations(id TEXT PRIMARY KEY) WITHOUT ROWID")
        db.executemany("INSERT INTO supplemental_observations VALUES (?)", ((oid,) for oid in sorted(accepted)))
        db.execute("""INSERT INTO main.observations
                    SELECT o.* FROM extra.observations o JOIN supplemental_observations s ON o.id=s.id""")
        raw_ids = db.execute("SELECT COUNT(*) FROM extra.observations").fetchone()[0]
        present = db.execute("SELECT COUNT(*) FROM extra.observations JOIN supplemental_observations USING(id)").fetchone()[0]
        if raw_ids > evidence_report["counts"]["observations"] or present != len(accepted):
            raise ValueError("Supplemental fragment and canonical observation inventories disagree")
        library = load_fragment_library(db, primary)
        counts = Counter(primary["counts"])
        counts["source_observations"] += evidence_report["counts"]["observations"]
        counts["source_screen_excluded"] += evidence_report["counts"]["source_screen_excluded"]
        db.execute("CREATE TEMP TABLE product_mapping(original INTEGER PRIMARY KEY, merged INTEGER NOT NULL)")
        for old_id, smiles in db.execute("""SELECT DISTINCT p.id, p.smiles FROM extra.products p
                JOIN extra.links l ON l.product_id=p.id
                JOIN supplemental_observations s ON s.id=l.observation_id ORDER BY p.id"""):
            row = db.execute("SELECT id FROM main.products WHERE smiles=?", (smiles,)).fetchone()
            if row:
                target = row[0]
            else:
                target = library.AddMol(Chem.MolFromSmiles(smiles))
                db.execute("INSERT INTO main.products VALUES (?,?)", (target, smiles))
            db.execute("INSERT INTO product_mapping VALUES (?,?)", (old_id, target))
        db.execute("""INSERT INTO main.links
            SELECT m.merged,l.observation_id,l.component_index,l.atom_order,l.evidence,l.side
            FROM extra.links l JOIN product_mapping m ON m.original=l.product_id
            JOIN supplemental_observations s ON s.id=l.observation_id ORDER BY l.rowid""")
        for side, payload in db.execute("""SELECT l.side,l.evidence FROM extra.links l
                JOIN supplemental_observations s ON s.id=l.observation_id"""):
            counts[f"indexed_{side}_components"] += 1
            counts["evidence_" + unpack(payload)["evidence_status"]] += 1
        # Unassigned procedures without an observation ID can refer to excluded
        # occurrences. Only globally identified source observations are joined.
        db.execute("""INSERT INTO main.procedures SELECT p.* FROM extra.procedures p
                   JOIN supplemental_observations s ON s.id=p.observation_id ORDER BY p.rowid""")
        counts["distinct_molecules"] = len(library)
        counts["distinct_products"] = db.execute(
            "SELECT COUNT(DISTINCT product_id) FROM main.links WHERE side='product'").fetchone()[0]
        counts["procedures"] += evidence_report["counts"]["procedures"]
        base_binding = {**primary["evidence_catalog"],
                        "relative_path": os.path.relpath(base_catalog, publish_directory)}
        supplemental_binding = {"catalog_id": evidence_report["catalog_id"],
                                "relative_path": supplemental_catalog.name,
                                "size_bytes": supplemental_catalog.stat().st_size,
                                "sha256": file_sha256(supplemental_catalog)}
        manifest = {key: value for key, value in primary.items()
                    if key not in {"index_id", "library_storage", "sources"}}
        manifest.update(
            source_scope="base_and_screened_supplement", counts=dict(counts),
            supplemental_observations=present, source_mapping_screen=SOURCE_MAPPING_SCREEN_VERSION,
            evidence_catalog=base_binding, supplemental_evidence_catalog=supplemental_binding,
            source_indexes=[{"index_id": m["index_id"], "sha256": file_sha256(path)}
                            for m, path in zip(manifests, (base, supplement))],
        )
        return _finish_library(db, library, manifest)


def main(argv: list[str] | None = None) -> int:
    """Create an opt-in combined fragment index and a screened supplemental catalog."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("supplement", type=Path)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args(argv)
    report = extend_fragment_index(args.base, args.supplement, args.output_directory)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
