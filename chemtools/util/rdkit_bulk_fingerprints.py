"""
Bulk molecule fingerprint ingestion utilities backed by RDKit 2025 APIs.

Primary fast path:
- rdkit.Chem.rdMolProcessing.GetFingerprintsForMolsInFile()

Fallback path:
- Python supplier loop + fingerprint generator

The API returns deterministic numpy matrices and row indices so ETL pipelines
can save compact NPZ indexes and preserve source-row traceability.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

import numpy as np


FingerprintType = Literal["morgan", "rdkit", "atom_pair", "topological_torsion"]


@dataclass(frozen=True)
class FingerprintBatch:
    """Container for fingerprinting output."""

    fingerprints: np.ndarray
    row_indices: np.ndarray
    fp_type: FingerprintType
    fp_size: int
    engine: str
    total_records: int
    valid_records: int
    invalid_records: int

    def to_summary(self) -> dict[str, object]:
        return {
            "fp_type": self.fp_type,
            "fp_size": self.fp_size,
            "engine": self.engine,
            "total_records": self.total_records,
            "valid_records": self.valid_records,
            "invalid_records": self.invalid_records,
            "shape": tuple(int(x) for x in self.fingerprints.shape),
        }


def _suppress_rdkit_parse_logs():
    """Best-effort context manager to silence noisy parse warnings."""
    try:
        from rdkit import rdBase  # type: ignore
    except Exception:
        class _NoOp:
            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        return _NoOp()
    blocker = getattr(rdBase, "BlockLogs", None)
    if blocker is None:
        class _NoOp:
            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        return _NoOp()
    return blocker()


def rdkit_bulk_fingerprinting_available() -> bool:
    """Return True when rdMolProcessing bulk API is available."""
    try:
        from rdkit.Chem import rdMolProcessing  # type: ignore

        return hasattr(rdMolProcessing, "GetFingerprintsForMolsInFile")
    except Exception:
        return False


def _build_generator(
    *,
    fp_type: FingerprintType,
    fp_size: int,
    radius: int,
    include_chirality: bool,
):
    from rdkit.Chem import rdFingerprintGenerator  # type: ignore

    if fp_type == "morgan":
        return rdFingerprintGenerator.GetMorganGenerator(
            radius=int(radius),
            fpSize=int(fp_size),
            includeChirality=bool(include_chirality),
        )
    if fp_type == "rdkit":
        return rdFingerprintGenerator.GetRDKitFPGenerator(
            fpSize=int(fp_size),
        )
    if fp_type == "atom_pair":
        return rdFingerprintGenerator.GetAtomPairGenerator(
            fpSize=int(fp_size),
            includeChirality=bool(include_chirality),
        )
    if fp_type == "topological_torsion":
        return rdFingerprintGenerator.GetTopologicalTorsionGenerator(
            fpSize=int(fp_size),
            includeChirality=bool(include_chirality),
        )
    raise ValueError(f"Unsupported fingerprint type: {fp_type}")


def _bitvect_to_array(bitvect, *, fp_size: int) -> np.ndarray:
    from rdkit import DataStructs  # type: ignore

    arr = np.zeros((int(fp_size),), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(bitvect, arr)
    return arr


def _empty_matrix(fp_size: int) -> np.ndarray:
    return np.zeros((0, int(fp_size)), dtype=np.uint8)


def _is_sdf_like(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith(".sdf") or name.endswith(".sdf.gz") or name.endswith(".sd")


def _fingerprints_python_fallback(
    *,
    input_path: Path,
    generator,
    fp_type: FingerprintType,
    fp_size: int,
    delimiter: str,
    smiles_column: int,
    name_column: int,
    title_line: bool,
    sanitize: bool,
    remove_hs: bool,
    skip_invalid: bool,
) -> FingerprintBatch:
    from rdkit import Chem  # type: ignore

    if _is_sdf_like(input_path):
        supplier = Chem.SDMolSupplier(
            str(input_path),
            sanitize=bool(sanitize),
            removeHs=bool(remove_hs),
            strictParsing=True,
        )
    else:
        supplier = Chem.SmilesMolSupplier(
            str(input_path),
            delimiter=str(delimiter),
            smilesColumn=int(smiles_column),
            nameColumn=int(name_column),
            titleLine=bool(title_line),
            sanitize=bool(sanitize),
        )

    rows: list[np.ndarray] = []
    row_indices: list[int] = []
    total = 0

    with _suppress_rdkit_parse_logs():
        for idx, mol in enumerate(supplier):
            total += 1
            if mol is None:
                if skip_invalid:
                    continue
                rows.append(np.zeros((int(fp_size),), dtype=np.uint8))
                row_indices.append(int(idx))
                continue
            bitvect = generator.GetFingerprint(mol)
            rows.append(_bitvect_to_array(bitvect, fp_size=fp_size))
            row_indices.append(int(idx))

    fps = np.vstack(rows) if rows else _empty_matrix(fp_size)
    row_arr = np.array(row_indices, dtype=np.int64)
    valid = int(fps.shape[0])
    invalid = int(total - valid)
    return FingerprintBatch(
        fingerprints=fps,
        row_indices=row_arr,
        fp_type=fp_type,
        fp_size=int(fp_size),
        engine="python_supplier_loop",
        total_records=int(total),
        valid_records=valid,
        invalid_records=invalid,
    )


def fingerprints_from_molecule_file(
    input_path: str | Path,
    *,
    fp_type: FingerprintType = "morgan",
    fp_size: int = 2048,
    radius: int = 2,
    include_chirality: bool = False,
    num_threads: int = 0,
    delimiter: str = "\t",
    smiles_column: int = 0,
    name_column: int = 1,
    title_line: bool = True,
    sanitize: bool = True,
    remove_hs: bool = True,
    strict_parsing: bool = True,
    skip_invalid: bool = True,
    prefer_bulk: bool = True,
) -> FingerprintBatch:
    """
    Fingerprint molecules from an input file into a uint8 matrix.

    Supported inputs include `.smi`, `.smiles`, `.txt`, `.sdf`, `.sdf.gz`.
    """
    path = Path(input_path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if int(fp_size) <= 0:
        raise ValueError("fp_size must be > 0")
    if int(radius) < 0:
        raise ValueError("radius must be >= 0")

    generator = _build_generator(
        fp_type=fp_type,
        fp_size=int(fp_size),
        radius=int(radius),
        include_chirality=bool(include_chirality),
    )

    if not prefer_bulk or not rdkit_bulk_fingerprinting_available():
        return _fingerprints_python_fallback(
            input_path=path,
            generator=generator,
            fp_type=fp_type,
            fp_size=int(fp_size),
            delimiter=str(delimiter),
            smiles_column=int(smiles_column),
            name_column=int(name_column),
            title_line=bool(title_line),
            sanitize=bool(sanitize),
            remove_hs=bool(remove_hs),
            skip_invalid=bool(skip_invalid),
        )

    from rdkit.Chem import rdMolProcessing  # type: ignore

    opts = rdMolProcessing.SupplierOptions()
    opts.numThreads = int(num_threads)
    opts.delimiter = str(delimiter)
    opts.smilesColumn = int(smiles_column)
    opts.nameColumn = int(name_column)
    opts.titleLine = bool(title_line)
    opts.sanitize = bool(sanitize)
    opts.removeHs = bool(remove_hs)
    opts.strictParsing = bool(strict_parsing)

    with _suppress_rdkit_parse_logs():
        raw_fps = rdMolProcessing.GetFingerprintsForMolsInFile(str(path), generator, opts)

    rows: list[np.ndarray] = []
    row_indices: list[int] = []
    for idx, fp in enumerate(raw_fps):
        if fp is None:
            if skip_invalid:
                continue
            rows.append(np.zeros((int(fp_size),), dtype=np.uint8))
            row_indices.append(int(idx))
            continue
        rows.append(_bitvect_to_array(fp, fp_size=int(fp_size)))
        row_indices.append(int(idx))

    fps = np.vstack(rows) if rows else _empty_matrix(int(fp_size))
    row_arr = np.array(row_indices, dtype=np.int64)
    total = int(len(raw_fps))
    valid = int(fps.shape[0])
    invalid = int(total - valid)
    return FingerprintBatch(
        fingerprints=fps,
        row_indices=row_arr,
        fp_type=fp_type,
        fp_size=int(fp_size),
        engine="rdMolProcessing",
        total_records=total,
        valid_records=valid,
        invalid_records=invalid,
    )


def save_fingerprint_batch(
    batch: FingerprintBatch,
    output_path: str | Path,
    *,
    row_id_prefix: str = "row_",
) -> None:
    """Persist a `FingerprintBatch` to compressed NPZ."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    row_ids = np.array([f"{row_id_prefix}{int(i)}" for i in batch.row_indices], dtype=object)
    np.savez_compressed(
        str(out),
        fps=batch.fingerprints.astype("uint8"),
        row_indices=batch.row_indices.astype("int64"),
        row_ids=row_ids,
        fp_type=np.array(batch.fp_type, dtype=object),
        fp_size=np.array(int(batch.fp_size), dtype="int32"),
        engine=np.array(batch.engine, dtype=object),
        total_records=np.array(int(batch.total_records), dtype="int64"),
        valid_records=np.array(int(batch.valid_records), dtype="int64"),
        invalid_records=np.array(int(batch.invalid_records), dtype="int64"),
    )
