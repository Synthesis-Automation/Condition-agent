from __future__ import annotations

from pathlib import Path

import numpy as np

from chemtools import reaction_similarity
from chemtools.precedent import loader, search
from chemtools.util.drfp_storage import DRFPLoader


def _fake_fp(rsmi: str) -> np.ndarray:
    seed = sum(ord(ch) for ch in rsmi) % 7
    return np.array([(seed + idx) % 2 for idx in range(16)], dtype="uint8")


def test_ensure_precedent_drfp_cache_writes_npz(monkeypatch, tmp_path: Path) -> None:
    rows = [
        {"reaction_id": "rxn:1", "reaction_smiles": "A.B>>P"},
        {"reaction_id": "rxn:2", "reaction_smiles": "C.D>>Q"},
        {"reaction_id": "rxn:2", "reaction_smiles": "ignored.duplicate>>Q"},
    ]

    monkeypatch.setattr(loader, "_PRECEDENT_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(loader, "_csv_max_mtime", lambda: 0.0)
    monkeypatch.setattr(reaction_similarity, "drfp_available", lambda: True)
    monkeypatch.setattr(
        reaction_similarity,
        "encode_drfp_cached",
        lambda rsmi, n_bits=4096, radius=3: _fake_fp(rsmi),
    )

    family_key = ("Suzuki_miyaura",)
    npz_path = loader.ensure_precedent_drfp_cache(family_key, rows, force=True)

    assert npz_path == loader.get_precedent_drfp_cache_path(family_key)
    assert npz_path is not None
    assert Path(npz_path).exists()

    drfp_loader = DRFPLoader(npz_path)
    assert len(drfp_loader) == 2
    np.testing.assert_array_equal(
        drfp_loader.get_fingerprint("rxn:1"),
        _fake_fp("A.B>>P"),
    )


def test_knn_generates_and_reuses_precedent_drfp_cache(monkeypatch, tmp_path: Path) -> None:
    rows = [
        {
            "reaction_id": "r:1",
            "dataset_reaction_id": "Suzuki_miyaura",
            "rxn_type": "Suzuki_miyaura",
            "reaction_smiles": "ArBr.B(OH)2>>ArAr",
            "features": {"bin": "BIN_TEST", "LG": "Br", "nuc_class": "boronic_acid"},
            "yield_value": 88.0,
            "condition_core": "Pd/SPhos",
            "base_uid": "K3PO4",
            "solvent_uid": "THF",
            "reagents": [],
            "solvents": [],
            "source_file": "Suzuki_miyaura",
            "source_group": "literature",
            "reference": "",
            "conditions": {},
            "catalyst": {},
            "full_system": None,
            "catalytic_system": [],
            "T_C": None,
            "time_h": None,
            "precomputed": {},
        },
        {
            "reaction_id": "r:2",
            "dataset_reaction_id": "Suzuki_miyaura",
            "rxn_type": "Suzuki_miyaura",
            "reaction_smiles": "ArCl.B(OH)2>>ArAr",
            "features": {"bin": "BIN_TEST", "LG": "Br", "nuc_class": "boronic_acid"},
            "yield_value": 84.0,
            "condition_core": "Pd/SPhos",
            "base_uid": "K3PO4",
            "solvent_uid": "THF",
            "reagents": [],
            "solvents": [],
            "source_file": "Suzuki_miyaura",
            "source_group": "literature",
            "reference": "",
            "conditions": {},
            "catalyst": {},
            "full_system": None,
            "catalytic_system": [],
            "T_C": None,
            "time_h": None,
            "precomputed": {},
        },
    ]

    calls = {"count": 0}

    def fake_encode(rsmi, n_bits=4096, radius=3):
        calls["count"] += 1
        return _fake_fp(rsmi)

    monkeypatch.setattr(loader, "_PRECEDENT_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(loader, "_csv_max_mtime", lambda: 0.0)
    monkeypatch.setattr(reaction_similarity, "drfp_available", lambda: True)
    monkeypatch.setattr(reaction_similarity, "encode_drfp_cached", fake_encode)
    monkeypatch.setattr(search, "_load", lambda: rows)
    monkeypatch.setattr(search, "_load_selective", lambda families=None: rows)
    monkeypatch.setattr(search, "_family_text", lambda family: family)
    monkeypatch.setattr(search, "rs", reaction_similarity)
    search._DRFP_LOADER_CACHE.clear()
    search._knn_cached.cache_clear()

    out = search.knn(
        family="Suzuki_miyaura",
        features={"bin": "BIN_TEST", "LG": "Br", "nuc_class": "boronic_acid"},
        k=2,
        relax={
            "filter_by_reagent_database": False,
            "selective_loading": True,
            "use_drfp": True,
            "reaction_smiles": "ArI.B(OH)2>>ArAr",
        },
    )

    assert out["drfp_cache_generated"] is True
    assert out["drfp_loader_available"] is True
    assert out["drfp_load_strategy"]["binary"] == 2
    assert out["drfp_load_strategy"]["computed"] == 0
    assert calls["count"] == 3
    assert Path(loader.get_precedent_drfp_cache_path(("Suzuki_miyaura",))).exists()
