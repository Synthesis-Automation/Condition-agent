"""Reusable canonical molecular-identity indexes backed by SQLite."""

from __future__ import annotations

import csv
import json
import os
import sqlite3
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional, Sequence

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors


MOLECULE_INDEX_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class MoleculeIdentity:
    """Strict map-free identity shared by index builders and consumers."""

    canonical_smiles: str
    inchi_key: Optional[str]
    molecular_weight: float


@dataclass(frozen=True)
class MoleculeIndexMatch:
    """One exact normalized match and bounded source provenance."""

    canonical_smiles: str
    inchi_key: Optional[str]
    occurrence_count: int
    source_records: tuple[dict[str, str], ...]


@dataclass(frozen=True)
class MoleculeIndexBuildReport:
    """Auditable result of compiling a CSV molecule column into SQLite."""

    source_path: str
    output_path: str
    smiles_column: str
    provenance_columns: tuple[str, ...]
    source_rows: int
    accepted_rows: int
    invalid_smiles_rows: int
    unique_molecules: int
    schema_version: str = MOLECULE_INDEX_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible report."""

        return asdict(self)


@lru_cache(maxsize=100_000)
def molecule_identity(smiles: str) -> Optional[MoleculeIdentity]:
    """Parse and return canonical isomeric SMILES, InChIKey, and MolWt."""

    molecule = Chem.MolFromSmiles(str(smiles or "").strip())
    if molecule is None:
        return None
    try:
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.SanitizeMol(molecule)
        fragments = Chem.GetMolFrags(
            molecule,
            asMols=True,
            sanitizeFrags=True,
        )
        fragment_smiles = sorted(
            Chem.MolToSmiles(
                fragment,
                canonical=True,
                isomericSmiles=True,
            )
            for fragment in fragments
        )
        canonical = ".".join(fragment_smiles)
        normalized = Chem.MolFromSmiles(canonical)
        if normalized is None:
            return None
        inchi_key = str(Chem.MolToInchiKey(normalized) or "").strip() or None
        molecular_weight = float(Descriptors.MolWt(normalized))
    except Exception:
        return None
    return MoleculeIdentity(
        canonical_smiles=canonical,
        inchi_key=inchi_key,
        molecular_weight=round(molecular_weight, 6),
    )


def _create_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE molecules (
            canonical_smiles TEXT PRIMARY KEY,
            inchi_key TEXT,
            molecular_weight REAL NOT NULL,
            occurrence_count INTEGER NOT NULL
        );
        CREATE INDEX molecules_inchi_key_idx ON molecules(inchi_key);
        CREATE TABLE provenance (
            canonical_smiles TEXT NOT NULL,
            source_row INTEGER NOT NULL,
            metadata_json TEXT NOT NULL,
            PRIMARY KEY (canonical_smiles, source_row),
            FOREIGN KEY (canonical_smiles)
                REFERENCES molecules(canonical_smiles)
        );
        """
    )


def build_canonical_molecule_index(
    source_csv: str | Path,
    output_path: str | Path,
    *,
    smiles_column: str = "compound_smiles",
    provenance_columns: Sequence[str] = (),
) -> MoleculeIndexBuildReport:
    """Build an atomic reusable identity index from one CSV SMILES column."""

    source = Path(source_csv)
    output = Path(output_path)
    if not source.is_file():
        raise FileNotFoundError(f"molecule catalog is unavailable: {source}")
    if not smiles_column.strip():
        raise ValueError("SMILES column must not be empty")
    selected_provenance = tuple(
        dict.fromkeys(str(value).strip() for value in provenance_columns)
    )
    if any(not value for value in selected_provenance):
        raise ValueError("provenance columns must not be empty")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    if temporary.exists():
        temporary.unlink()

    source_rows = 0
    accepted_rows = 0
    invalid_rows = 0
    connection = sqlite3.connect(temporary)
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        _create_schema(connection)
        with source.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = tuple(reader.fieldnames or ())
            if smiles_column not in fields:
                raise ValueError(
                    f"molecule catalog requires a {smiles_column} column"
                )
            missing = [
                value for value in selected_provenance if value not in fields
            ]
            if missing:
                raise ValueError(
                    "molecule catalog is missing provenance columns: "
                    + ", ".join(missing)
                )
            with rdBase.BlockLogs():
                for source_row, row in enumerate(reader, start=2):
                    source_rows += 1
                    identity = molecule_identity(row.get(smiles_column) or "")
                    if identity is None:
                        invalid_rows += 1
                        continue
                    accepted_rows += 1
                    connection.execute(
                        """
                        INSERT INTO molecules(
                            canonical_smiles,
                            inchi_key,
                            molecular_weight,
                            occurrence_count
                        ) VALUES (?, ?, ?, 1)
                        ON CONFLICT(canonical_smiles) DO UPDATE SET
                            occurrence_count = occurrence_count + 1
                        """,
                        (
                            identity.canonical_smiles,
                            identity.inchi_key,
                            identity.molecular_weight,
                        ),
                    )
                    provenance = {
                        column: str(row.get(column) or "").strip()
                        for column in selected_provenance
                    }
                    connection.execute(
                        """
                        INSERT INTO provenance(
                            canonical_smiles,
                            source_row,
                            metadata_json
                        ) VALUES (?, ?, ?)
                        """,
                        (
                            identity.canonical_smiles,
                            source_row,
                            json.dumps(
                                provenance,
                                sort_keys=True,
                                separators=(",", ":"),
                            ),
                        ),
                    )
        unique_molecules = int(
            connection.execute("SELECT COUNT(*) FROM molecules").fetchone()[0]
        )
        metadata = {
            "schema_version": MOLECULE_INDEX_SCHEMA_VERSION,
            "source_path": str(source.resolve()),
            "smiles_column": smiles_column,
            "provenance_columns": selected_provenance,
            "source_rows": source_rows,
            "accepted_rows": accepted_rows,
            "invalid_smiles_rows": invalid_rows,
            "unique_molecules": unique_molecules,
            "identity_policy": (
                "map_free_canonical_isomeric_smiles_or_full_inchi_key"
            ),
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
    return MoleculeIndexBuildReport(
        source_path=str(source.resolve()),
        output_path=str(output.resolve()),
        smiles_column=smiles_column,
        provenance_columns=selected_provenance,
        source_rows=source_rows,
        accepted_rows=accepted_rows,
        invalid_smiles_rows=invalid_rows,
        unique_molecules=unique_molecules,
    )


class CanonicalMoleculeIndex:
    """Read-only exact identity lookup over a compiled molecule catalog."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        if not self.path.is_file():
            raise FileNotFoundError(f"molecule index is unavailable: {self.path}")
        self._connection: sqlite3.Connection | None = sqlite3.connect(self.path)
        raw_version = self._connection.execute(
            "SELECT value FROM metadata WHERE key = 'schema_version'"
        ).fetchone()
        version = json.loads(raw_version[0]) if raw_version else None
        if version != MOLECULE_INDEX_SCHEMA_VERSION:
            self.close()
            raise ValueError("unsupported molecule index schema")

    def close(self) -> None:
        """Close the underlying SQLite connection."""

        if self._connection is not None:
            self._connection.close()
            self._connection = None

    def __enter__(self) -> "CanonicalMoleculeIndex":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def lookup(
        self,
        identity: MoleculeIdentity | str,
        *,
        provenance_limit: int = 5,
    ) -> Optional[MoleculeIndexMatch]:
        """Return an exact canonical-SMILES or full-InChIKey match."""

        if provenance_limit < 1:
            raise ValueError("provenance limit must be positive")
        resolved = (
            molecule_identity(identity) if isinstance(identity, str) else identity
        )
        if resolved is None:
            return None
        if self._connection is None:
            raise RuntimeError("molecule index is closed")
        row = self._connection.execute(
            """
            SELECT canonical_smiles, inchi_key, occurrence_count
            FROM molecules
            WHERE canonical_smiles = ?
               OR (inchi_key IS NOT NULL AND inchi_key = ?)
            ORDER BY canonical_smiles = ? DESC, canonical_smiles
            LIMIT 1
            """,
            (
                resolved.canonical_smiles,
                resolved.inchi_key,
                resolved.canonical_smiles,
            ),
        ).fetchone()
        if row is None:
            return None
        records = self._connection.execute(
            """
            SELECT metadata_json
            FROM provenance
            WHERE canonical_smiles = ?
            ORDER BY source_row
            LIMIT ?
            """,
            (row[0], provenance_limit),
        ).fetchall()
        return MoleculeIndexMatch(
            canonical_smiles=str(row[0]),
            inchi_key=str(row[1]) if row[1] else None,
            occurrence_count=int(row[2]),
            source_records=tuple(json.loads(value[0]) for value in records),
        )


__all__ = [
    "MOLECULE_INDEX_SCHEMA_VERSION",
    "CanonicalMoleculeIndex",
    "MoleculeIdentity",
    "MoleculeIndexBuildReport",
    "MoleculeIndexMatch",
    "build_canonical_molecule_index",
    "molecule_identity",
]
