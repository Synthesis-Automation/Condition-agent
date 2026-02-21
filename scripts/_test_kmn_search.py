"""Quick integration test for KMN index-based searching."""
import sys
import os
import time
import numpy as np

# Ensure repo root is on sys.path when running as a script
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


def main():
    print("=" * 60)
    print("KMN INDEX INTEGRATION TEST")
    print("=" * 60)

    # ── 1. is_index_built() ───────────────────────────────────
    print("\n[1] is_index_built() ...")
    from chemtools.util.faiss_router import is_index_built, get_default_router, INDEX_DIR

    assert is_index_built(), "FAIL: is_index_built() returned False"
    print(f"    OK  index_dir = {INDEX_DIR}")

    # ── 2. load FAISSRouter ───────────────────────────────────
    print("\n[2] Loading FAISSRouter ...")
    t0 = time.time()
    router = get_default_router()
    elapsed = time.time() - t0
    assert router.is_ready, "FAIL: router.is_ready is False after load"
    ca = "present" if router._cell_assignments is not None else "absent"
    print(f"    OK  {router}  loaded in {elapsed:.2f}s")
    print(f"    cell_assignments: {ca}  (needed for soft-Voronoi)")

    # ── 3. load KMN weights ───────────────────────────────────
    print("\n[3] Loading KernelMetricNetwork ...")
    from chemtools.util.kernel_metric_net import get_default_kmn

    t0 = time.time()
    kmn = get_default_kmn()
    elapsed = time.time() - t0
    print(f"    OK  backend={kmn._backend}  embed_dim={kmn.embed_dim}  loaded in {elapsed:.2f}s")

    # ── 4. end-to-end: MixFP → embed → FAISS search ──────────
    print("\n[4] End-to-end: MixFP → KMN → FAISS search ...")
    from chemtools.util.mix_fingerprint import create_mix_fp

    TEST_SMILES = [
        ("Acetal formation (in-distribution)",
         "CO.O=Cc1ccc(COC2CCCCO2)cc1>>COC(OC)c1ccc(CO)cc1"),
        ("Suzuki-type (cross-family)",
         "c1ccc(B(O)O)cc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"),
        ("Amide coupling",
         "CC(=O)O.Nc1ccccc1>>CC(=O)Nc1ccccc1"),
    ]

    for label, smi in TEST_SMILES:
        fp = create_mix_fp(smi, fp_size=1024)
        if fp is None:
            print(f"    SKIP  {label}: fingerprint computation failed")
            continue
        emb = kmn.get_embeddings(fp.reshape(1, -1))[0]
        ids, scores = router.search(emb, k=5)
        top_scores = [round(s, 4) for s in scores]
        print(f"    [{label}]")
        print(f"      top-5 scores : {top_scores}")
        print(f"      top-1 rxn_id : {ids[0] if ids else 'none'}")
        # Sanity: scores should be in [-1, 1] (cosine)
        assert all(-1.0 <= s <= 1.0 for s in scores), "FAIL: score out of cosine range"

    # ── 5. adaptive_search (confidence-gated) ────────────────
    print("\n[5] adaptive_search (confidence-gated routing) ...")
    fp = create_mix_fp(
        "CO.O=Cc1ccc(COC2CCCCO2)cc1>>COC(OC)c1ccc(CO)cc1", fp_size=1024
    )
    assert fp is not None, "FAIL: could not compute fingerprint"
    emb = kmn.get_embeddings(fp.reshape(1, -1))[0]

    expected_modes = {0.92: "hard_filter", 0.65: "soft_bleed", 0.25: "full_voronoi"}
    for conf in [0.92, 0.65, 0.25]:
        ids, scores, meta = router.adaptive_search(emb, k=10, family_confidence=conf)
        mode = meta["mode"]
        assert mode == expected_modes[conf], (
            f"FAIL: conf={conf} expected mode={expected_modes[conf]}, got {mode}"
        )
        top = round(scores[0], 4) if scores else 0
        print(
            f"    conf={conf}  mode={mode:12s}  nprobe={meta['nprobe']:3d}"
            f"  n={meta['n_candidates']:3d}  top_score={top}"
        )

    # ── 6. auto-detection in search.knn() ────────────────────
    print("\n[6] Auto-detection in chemtools.precedent.search.knn() ...")
    from chemtools.precedent.search import knn, _MIXFP_AVAILABLE, _is_index_built

    print(f"    _MIXFP_AVAILABLE : {_MIXFP_AVAILABLE}")
    auto_on = (
        _MIXFP_AVAILABLE
        and _is_index_built is not None
        and _is_index_built()
    )
    print(f"    auto-enable MixFP when smiles provided: {auto_on}")
    assert auto_on, "FAIL: auto-detection should be True with a built index"

    # Call knn — MixFP should activate automatically
    t0 = time.time()
    result = knn(
        family="Acetal_ketal_formation",
        features={},
        k=5,
        relax={
            "reaction_smiles": "CO.O=Cc1ccc(COC2CCCCO2)cc1>>COC(OC)c1ccc(CO)cc1",
            "debug_timing": True,
        },
    )
    elapsed = time.time() - t0
    precs = result.get("precedents", [])
    err = result.get("error")
    print(f"    knn returned {len(precs)} precedents  error={err}  elapsed={elapsed:.2f}s")
    if precs:
        p = precs[0]
        print(
            f"    top precedent: id={str(p.get('reaction_id','?'))[:50]}"
            f"  score={p.get('score','?')}"
        )

    # Call knn with use_mixfp=False — should skip the index
    result2 = knn(
        family="Acetal_ketal_formation",
        features={},
        k=5,
        relax={
            "reaction_smiles": "CO.O=Cc1ccc(COC2CCCCO2)cc1>>COC(OC)c1ccc(CO)cc1",
            "use_mixfp": False,
        },
    )
    print(f"    knn with use_mixfp=False: {len(result2.get('precedents',[]))} precedents  (index bypassed)")

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    main()
