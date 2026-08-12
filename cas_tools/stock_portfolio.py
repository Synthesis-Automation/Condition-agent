"""Versioned supplier-stock portfolios with exact molecular lookup."""

from __future__ import annotations

import csv
import gzip
import json
import os
import sqlite3
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Literal, Optional, Sequence

from rdkit import Chem, rdBase

from .molecule_index import (
    CanonicalMoleculeIndex,
    MoleculeIdentity,
    MoleculeIndexMatch,
    molecule_identity,
)


STOCK_PORTFOLIO_SCHEMA_VERSION = "1.0"
STOCK_SOURCE_MANIFEST_VERSION = "1.0"

StockEvidence = Literal[
    "physically_available",
    "supplier_in_stock",
    "catalog_listed",
    "make_on_demand",
]
StockFormat = Literal["smi", "csv", "sdf"]


@dataclass(frozen=True)
class StockSourceDefinition:
    """One authorized, dated supplier snapshot to ingest."""

    path: str
    supplier: str
    collection: str
    snapshot_date: str
    availability_status: str
    evidence_level: StockEvidence
    terminal_eligible: bool
    source_url: str
    terms_url: str
    region: str = "global"
    format: StockFormat = "smi"
    smiles_column: str = "smiles"
    identifier_column: str = "id"
    delimiter: str = "\t"

    def validate(self) -> None:
        """Reject incomplete or chemically misleading source definitions."""

        required = {
            "path": self.path,
            "supplier": self.supplier,
            "collection": self.collection,
            "snapshot_date": self.snapshot_date,
            "availability_status": self.availability_status,
            "source_url": self.source_url,
            "terms_url": self.terms_url,
        }
        missing = [name for name, value in required.items() if not value.strip()]
        if missing:
            raise ValueError(
                "stock source requires non-empty fields: " + ", ".join(missing)
            )
        if self.terminal_eligible and self.evidence_level not in {
            "physically_available",
            "supplier_in_stock",
        }:
            raise ValueError(
                "only physical or supplier-in-stock evidence may terminate routes"
            )
        if self.format not in {"smi", "csv", "sdf"}:
            raise ValueError(f"unsupported stock source format: {self.format}")


@dataclass(frozen=True)
class StockPortfolioBuildReport:
    """Auditable outcome of merging supplier snapshots."""

    output_path: str
    source_count: int
    source_rows: int
    accepted_rows: int
    invalid_structure_rows: int
    unique_molecules: int
    offer_count: int
    terminal_eligible_molecules: int
    source_summaries: tuple[dict[str, Any], ...]
    schema_version: str = STOCK_PORTFOLIO_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible build report."""

        return asdict(self)


def load_stock_source_manifest(
    path: str | Path,
) -> tuple[StockSourceDefinition, ...]:
    """Load and validate a declarative supplier-source manifest."""

    manifest_path = Path(path)
    with manifest_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != STOCK_SOURCE_MANIFEST_VERSION:
        raise ValueError("unsupported stock source manifest schema")
    sources = tuple(
        StockSourceDefinition(**raw) for raw in payload.get("sources") or ()
    )
    if not sources:
        raise ValueError("stock source manifest contains no sources")
    for source in sources:
        source.validate()
        if not Path(source.path).is_file():
            raise FileNotFoundError(f"stock source is unavailable: {source.path}")
    return sources


def _open_text(path: Path):
    if path.name.casefold().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def _iter_smi_records(
    source: StockSourceDefinition,
) -> Iterator[tuple[int, str, str]]:
    with _open_text(Path(source.path)) as handle:
        for row_number, line in enumerate(handle, start=1):
            value = line.strip()
            if not value:
                continue
            fields = value.split(source.delimiter)
            smiles = fields[0].strip()
            record_id = (
                fields[1].strip()
                if len(fields) > 1 and fields[1].strip()
                else str(row_number)
            )
            yield row_number, smiles, record_id


def _iter_csv_records(
    source: StockSourceDefinition,
) -> Iterator[tuple[int, str, str]]:
    with _open_text(Path(source.path)) as handle:
        reader = csv.DictReader(handle, delimiter=source.delimiter)
        fields = tuple(reader.fieldnames or ())
        if source.smiles_column not in fields:
            raise ValueError(
                f"stock CSV requires a {source.smiles_column} column"
            )
        for row_number, row in enumerate(reader, start=2):
            smiles = str(row.get(source.smiles_column) or "").strip()
            record_id = str(row.get(source.identifier_column) or row_number).strip()
            yield row_number, smiles, record_id


def _iter_sdf_records(
    source: StockSourceDefinition,
) -> Iterator[tuple[int, str, str]]:
    path = Path(source.path)
    binary_handle = (
        gzip.open(path, "rb")
        if path.name.casefold().endswith(".gz")
        else path.open("rb")
    )
    with binary_handle as handle, rdBase.BlockLogs():
        supplier = Chem.ForwardSDMolSupplier(handle, sanitize=True, removeHs=True)
        for row_number, molecule in enumerate(supplier, start=1):
            if molecule is None:
                yield row_number, "", str(row_number)
                continue
            record_id = (
                molecule.GetProp(source.identifier_column).strip()
                if molecule.HasProp(source.identifier_column)
                else molecule.GetProp("_Name").strip()
                if molecule.HasProp("_Name")
                else str(row_number)
            )
            yield (
                row_number,
                Chem.MolToSmiles(molecule, canonical=False, isomericSmiles=True),
                record_id or str(row_number),
            )


def _iter_source_records(
    source: StockSourceDefinition,
) -> Iterator[tuple[int, str, str]]:
    if source.format == "smi":
        yield from _iter_smi_records(source)
    elif source.format == "csv":
        yield from _iter_csv_records(source)
    else:
        yield from _iter_sdf_records(source)


def _chunks(
    values: Iterable[tuple[int, str, str]],
    size: int,
) -> Iterator[tuple[tuple[int, str, str], ...]]:
    chunk = []
    for value in values:
        chunk.append(value)
        if len(chunk) >= size:
            yield tuple(chunk)
            chunk = []
    if chunk:
        yield tuple(chunk)


def _canonicalize_chunk(
    values: tuple[tuple[int, str, str], ...],
) -> tuple[tuple[int, str, str], ...]:
    normalized = []
    with rdBase.BlockLogs():
        for row_number, smiles, record_id in values:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                normalized.append((row_number, "", record_id))
                continue
            for atom in molecule.GetAtoms():
                atom.SetAtomMapNum(0)
            try:
                Chem.SanitizeMol(molecule)
                canonical = Chem.MolToSmiles(
                    molecule,
                    canonical=True,
                    isomericSmiles=True,
                )
            except Exception:
                canonical = ""
            normalized.append((row_number, canonical, record_id))
    return tuple(normalized)


def _normalized_chunks(
    source: StockSourceDefinition,
    workers: int,
    chunk_size: int,
) -> Iterator[tuple[tuple[int, str, str], ...]]:
    chunks = _chunks(_iter_source_records(source), chunk_size)
    if workers == 1:
        for chunk in chunks:
            yield _canonicalize_chunk(chunk)
        return
    with ProcessPoolExecutor(max_workers=workers) as executor:
        pending = set()
        for chunk in chunks:
            pending.add(executor.submit(_canonicalize_chunk, chunk))
            if len(pending) < workers * 2:
                continue
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()


def _create_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE sources (
            source_id INTEGER PRIMARY KEY,
            supplier TEXT NOT NULL,
            collection_name TEXT NOT NULL,
            snapshot_date TEXT NOT NULL,
            availability_status TEXT NOT NULL,
            evidence_level TEXT NOT NULL,
            terminal_eligible INTEGER NOT NULL,
            region TEXT NOT NULL,
            source_url TEXT NOT NULL,
            terms_url TEXT NOT NULL,
            source_path TEXT NOT NULL,
            UNIQUE(supplier, collection_name, snapshot_date, region)
        );
        CREATE TABLE molecules (
            molecule_id INTEGER PRIMARY KEY,
            canonical_smiles TEXT NOT NULL UNIQUE
        );
        CREATE TABLE offer_staging (
            canonical_smiles TEXT NOT NULL,
            source_id INTEGER NOT NULL,
            source_row INTEGER NOT NULL,
            record_id TEXT NOT NULL
        );
        CREATE TABLE offers (
            molecule_id INTEGER NOT NULL,
            source_id INTEGER NOT NULL,
            source_row INTEGER NOT NULL,
            record_id TEXT NOT NULL,
            PRIMARY KEY(source_id, record_id),
            FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_id),
            FOREIGN KEY(source_id) REFERENCES sources(source_id)
        );
        """
    )


def build_stock_portfolio(
    sources: Sequence[StockSourceDefinition],
    output_path: str | Path,
    *,
    workers: int = 1,
    chunk_size: int = 2_000,
    progress: Optional[Callable[[str, int, int], None]] = None,
) -> StockPortfolioBuildReport:
    """Atomically normalize and merge authorized supplier snapshots."""

    if not sources:
        raise ValueError("at least one stock source is required")
    if workers < 1 or chunk_size < 1:
        raise ValueError("workers and chunk size must be positive")
    for source in sources:
        source.validate()
        if not Path(source.path).is_file():
            raise FileNotFoundError(f"stock source is unavailable: {source.path}")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    if temporary.exists():
        temporary.unlink()
    connection = sqlite3.connect(temporary)
    summaries = []
    total_rows = 0
    total_accepted = 0
    total_invalid = 0
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=MEMORY")
        connection.execute("PRAGMA cache_size=-200000")
        _create_schema(connection)
        for source in sources:
            cursor = connection.execute(
                """
                INSERT INTO sources(
                    supplier, collection_name, snapshot_date,
                    availability_status, evidence_level, terminal_eligible,
                    region, source_url, terms_url, source_path
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    source.supplier,
                    source.collection,
                    source.snapshot_date,
                    source.availability_status,
                    source.evidence_level,
                    int(source.terminal_eligible),
                    source.region,
                    source.source_url,
                    source.terms_url,
                    str(Path(source.path).resolve()),
                ),
            )
            source_id = int(cursor.lastrowid)
            source_rows = 0
            accepted_rows = 0
            invalid_rows = 0
            for normalized in _normalized_chunks(source, workers, chunk_size):
                source_rows += len(normalized)
                valid = tuple(value for value in normalized if value[1])
                accepted_rows += len(valid)
                invalid_rows += len(normalized) - len(valid)
                connection.executemany(
                    "INSERT OR IGNORE INTO molecules(canonical_smiles) VALUES (?)",
                    ((canonical,) for _, canonical, _ in valid),
                )
                connection.executemany(
                    """
                    INSERT INTO offer_staging(
                        canonical_smiles, source_id, source_row, record_id
                    ) VALUES (?, ?, ?, ?)
                    """,
                    (
                        (canonical, source_id, row_number, record_id)
                        for row_number, canonical, record_id in valid
                    ),
                )
                if progress is not None:
                    progress(source.collection, source_rows, accepted_rows)
            connection.commit()
            summaries.append(
                {
                    "supplier": source.supplier,
                    "collection": source.collection,
                    "snapshot_date": source.snapshot_date,
                    "source_rows": source_rows,
                    "accepted_rows": accepted_rows,
                    "invalid_structure_rows": invalid_rows,
                    "terminal_eligible": source.terminal_eligible,
                }
            )
            total_rows += source_rows
            total_accepted += accepted_rows
            total_invalid += invalid_rows
        connection.execute(
            """
            INSERT OR IGNORE INTO offers(
                molecule_id, source_id, source_row, record_id
            )
            SELECT molecules.molecule_id, offer_staging.source_id,
                   offer_staging.source_row, offer_staging.record_id
            FROM offer_staging
            JOIN molecules USING(canonical_smiles)
            """
        )
        connection.execute("DROP TABLE offer_staging")
        connection.executescript(
            """
            CREATE INDEX molecules_smiles_idx
                ON molecules(canonical_smiles);
            CREATE INDEX offers_molecule_idx
                ON offers(molecule_id);
            CREATE INDEX offers_source_idx
                ON offers(source_id);
            """
        )
        unique_molecules = int(
            connection.execute("SELECT COUNT(*) FROM molecules").fetchone()[0]
        )
        offer_count = int(
            connection.execute("SELECT COUNT(*) FROM offers").fetchone()[0]
        )
        terminal_molecules = int(
            connection.execute(
                """
                SELECT COUNT(DISTINCT offers.molecule_id)
                FROM offers
                JOIN sources USING(source_id)
                WHERE sources.terminal_eligible = 1
                """
            ).fetchone()[0]
        )
        metadata = {
            "schema_version": STOCK_PORTFOLIO_SCHEMA_VERSION,
            "identity_policy": "rdkit_map_free_canonical_isomeric_smiles",
            "source_count": len(sources),
            "source_rows": total_rows,
            "accepted_rows": total_accepted,
            "invalid_structure_rows": total_invalid,
            "unique_molecules": unique_molecules,
            "offer_count": offer_count,
            "terminal_eligible_molecules": terminal_molecules,
            "source_summaries": summaries,
        }
        connection.executemany(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            (
                (key, json.dumps(value, sort_keys=True))
                for key, value in metadata.items()
            ),
        )
        connection.commit()
    except Exception:
        connection.close()
        if temporary.exists():
            temporary.unlink()
        raise
    else:
        connection.close()
    temporary.replace(output)
    return StockPortfolioBuildReport(
        output_path=str(output.resolve()),
        source_count=len(sources),
        source_rows=total_rows,
        accepted_rows=total_accepted,
        invalid_structure_rows=total_invalid,
        unique_molecules=unique_molecules,
        offer_count=offer_count,
        terminal_eligible_molecules=terminal_molecules,
        source_summaries=tuple(summaries),
    )


class StockPortfolio:
    """Read-only exact lookup over combined supplier and internal stock offers."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        if not self.path.is_file():
            raise FileNotFoundError(f"stock portfolio is unavailable: {self.path}")
        self._connection: sqlite3.Connection | None = sqlite3.connect(self.path)
        raw = self._connection.execute(
            "SELECT value FROM metadata WHERE key = 'schema_version'"
        ).fetchone()
        version = json.loads(raw[0]) if raw else None
        if version != STOCK_PORTFOLIO_SCHEMA_VERSION:
            self.close()
            raise ValueError("unsupported stock portfolio schema")

    def close(self) -> None:
        """Close the SQLite connection."""

        if self._connection is not None:
            self._connection.close()
            self._connection = None

    def __enter__(self) -> "StockPortfolio":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def metadata(self) -> dict[str, Any]:
        """Return the portfolio build metadata."""

        if self._connection is None:
            raise RuntimeError("stock portfolio is closed")
        return {
            str(key): json.loads(value)
            for key, value in self._connection.execute(
                "SELECT key, value FROM metadata ORDER BY key"
            )
        }

    def lookup(
        self,
        identity: MoleculeIdentity | str,
        *,
        provenance_limit: int = 5,
    ) -> Optional[MoleculeIndexMatch]:
        """Return exact stock offers with explicit terminal evidence."""

        if provenance_limit < 1:
            raise ValueError("provenance limit must be positive")
        resolved = (
            molecule_identity(identity) if isinstance(identity, str) else identity
        )
        if resolved is None:
            return None
        if self._connection is None:
            raise RuntimeError("stock portfolio is closed")
        molecule = self._connection.execute(
            """
            SELECT molecule_id, canonical_smiles
            FROM molecules WHERE canonical_smiles = ?
            """,
            (resolved.canonical_smiles,),
        ).fetchone()
        if molecule is None:
            return None
        offer_count = int(
            self._connection.execute(
                "SELECT COUNT(*) FROM offers WHERE molecule_id = ?",
                (molecule[0],),
            ).fetchone()[0]
        )
        offers = self._connection.execute(
            """
            SELECT sources.supplier, sources.collection_name,
                   sources.snapshot_date, sources.availability_status,
                   sources.evidence_level, sources.terminal_eligible,
                   sources.region, sources.source_url, sources.terms_url,
                   offers.record_id
            FROM offers
            JOIN sources USING(source_id)
            WHERE offers.molecule_id = ?
            ORDER BY sources.terminal_eligible DESC,
                     sources.evidence_level, sources.supplier,
                     sources.collection_name, offers.record_id
            LIMIT ?
            """,
            (molecule[0], provenance_limit),
        ).fetchall()
        records = tuple(
            {
                "supplier": str(row[0]),
                "source_collection": str(row[1]),
                "snapshot_date": str(row[2]),
                "availability_status": str(row[3]),
                "stock_evidence": str(row[4]),
                "terminal_eligible": "true" if row[5] else "false",
                "region": str(row[6]),
                "source_url": str(row[7]),
                "terms_url": str(row[8]),
                "supplier_record_id": str(row[9]),
                "source_role": (
                    "starting_material" if row[5] else "catalog_listing"
                ),
            }
            for row in offers
        )
        return MoleculeIndexMatch(
            canonical_smiles=str(molecule[1]),
            inchi_key=resolved.inchi_key,
            occurrence_count=offer_count,
            source_records=records,
        )


def open_stock_lookup(
    path: str | Path,
) -> StockPortfolio | CanonicalMoleculeIndex:
    """Open a supplier portfolio or compatible legacy molecule index."""

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"stock lookup is unavailable: {source}")
    connection = sqlite3.connect(source)
    try:
        table_names = {
            str(row[0])
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
    finally:
        connection.close()
    if {"sources", "offers"}.issubset(table_names):
        return StockPortfolio(source)
    return CanonicalMoleculeIndex(source)


__all__ = [
    "STOCK_PORTFOLIO_SCHEMA_VERSION",
    "STOCK_SOURCE_MANIFEST_VERSION",
    "StockPortfolio",
    "StockPortfolioBuildReport",
    "StockSourceDefinition",
    "build_stock_portfolio",
    "load_stock_source_manifest",
    "open_stock_lookup",
]
