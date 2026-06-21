from pathlib import Path

import numpy as np
import pytest

from chemtools.core.rdkit import rdkit_available


pytestmark = pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")


def _write_smi(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_bulk_and_fallback_paths_match(tmp_path: Path):
    from chemtools.util.rdkit_bulk_fingerprints import fingerprints_from_molecule_file

    smi_path = tmp_path / "mols.smi"
    _write_smi(
        smi_path,
        [
            "smiles\tname",
            "CC\tethane",
            "CCC\tpropane",
            "c1ccccc1\tbenzene",
        ],
    )

    bulk = fingerprints_from_molecule_file(
        smi_path,
        fp_type="morgan",
        fp_size=256,
        title_line=True,
        delimiter="\t",
        prefer_bulk=True,
    )
    fallback = fingerprints_from_molecule_file(
        smi_path,
        fp_type="morgan",
        fp_size=256,
        title_line=True,
        delimiter="\t",
        prefer_bulk=False,
    )

    assert bulk.fingerprints.shape == (3, 256)
    assert fallback.fingerprints.shape == (3, 256)
    assert np.array_equal(bulk.row_indices, np.array([0, 1, 2], dtype=np.int64))
    assert np.array_equal(fallback.row_indices, np.array([0, 1, 2], dtype=np.int64))
    assert np.array_equal(bulk.fingerprints, fallback.fingerprints)


def test_invalid_rows_are_filtered_when_requested(tmp_path: Path):
    from chemtools.util.rdkit_bulk_fingerprints import fingerprints_from_molecule_file

    smi_path = tmp_path / "mols_invalid.smi"
    _write_smi(
        smi_path,
        [
            "CC good1",
            "INVALID bad",
            "CCC good2",
        ],
    )

    batch = fingerprints_from_molecule_file(
        smi_path,
        fp_type="morgan",
        fp_size=128,
        title_line=False,
        delimiter=" ",
        skip_invalid=True,
    )
    assert batch.total_records == 3
    assert batch.valid_records == 2
    assert batch.invalid_records == 1
    assert batch.fingerprints.shape == (2, 128)
    assert np.array_equal(batch.row_indices, np.array([0, 2], dtype=np.int64))


def test_save_fingerprint_batch_npz(tmp_path: Path):
    from chemtools.util.rdkit_bulk_fingerprints import (
        fingerprints_from_molecule_file,
        save_fingerprint_batch,
    )

    smi_path = tmp_path / "save_test.smi"
    out_path = tmp_path / "save_test_fp.npz"
    _write_smi(
        smi_path,
        [
            "CC a",
            "CCC b",
        ],
    )

    batch = fingerprints_from_molecule_file(
        smi_path,
        fp_type="morgan",
        fp_size=64,
        title_line=False,
        delimiter=" ",
    )
    save_fingerprint_batch(batch, out_path, row_id_prefix="rec_")

    payload = np.load(str(out_path), allow_pickle=True)
    assert payload["fps"].shape == (2, 64)
    assert payload["row_indices"].tolist() == [0, 1]
    assert payload["row_ids"].tolist() == ["rec_0", "rec_1"]
    assert int(payload["valid_records"]) == 2
