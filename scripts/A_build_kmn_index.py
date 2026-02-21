#!/usr/bin/env python
"""Build the MixFP + Kernel Metric Network + FAISS routing index.

This script implements the MOSAIC-style data-processing and index-building
pipeline (Li et al., Nature 2026):

    1. Load all reactions from the HTE precedent database.
    2. Compute the mixed reaction fingerprint [reactant | diff | product]
       for every row that has a valid reaction SMILES.
    3. Assign cluster labels from condition_core strings (optionally refined
       by reaction family).
    4. Train the Kernel Metric Network (KMN) on (MixFP, cluster_label) pairs.
    5. Embed all reactions with the trained KMN.
    6. Build and save a FAISS flat/IVF index from those embeddings.
    7. Write index files to ``data/kmn_index/`` (or CHEMTOOLS_KMN_INDEX_DIR).

Usage::

    python scripts/build_kmn_index.py [--fp-size 1024] [--embed-dim 256]
                                       [--epochs 50] [--batch-size 128]
                                       [--lr 1e-3] [--no-train]
                                       [--output-dir PATH]
                                       [--n-processes N]
                                       [--min-cluster-size 5]

Pass ``--no-train`` to skip KMN training and use the random-projection
fallback (fast, no torch required; useful to just build the FAISS index
for evaluation with random embeddings).
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from collections import Counter
from typing import Dict, List, Any, Optional

import numpy as np

# ── ensure repo root is on sys.path ──────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPT_DIR)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("build_kmn_index")


# ── helpers ────────────────────────────────────────────────────────────────

def _encode_cluster_labels(
    records: List[Dict[str, Any]],
    *,
    min_cluster_size: int = 5,
) -> np.ndarray:
    """Derive integer cluster IDs from condition_core + rxn_type.

    Cores with fewer than ``min_cluster_size`` members are merged into a
    single "rare" cluster so the classifier head has enough signal per class.

    Returns an int64 ndarray of shape (N,) parallel to ``records``.
    """
    raw_keys = []
    for r in records:
        core = (r.get("condition_core") or "").strip()
        family = (r.get("rxn_type") or "").strip()
        # Combine family + core for finer granularity
        key = f"{family}|{core}" if core else (family or "unknown")
        raw_keys.append(key)

    # Count occurrences and reassign rare clusters
    ctr = Counter(raw_keys)
    rare_id_str = "__rare__"
    normalised = [
        (k if ctr[k] >= min_cluster_size else rare_id_str)
        for k in raw_keys
    ]

    # Build string → int map
    unique = sorted(set(normalised))
    str2int = {s: i for i, s in enumerate(unique)}
    labels = np.array([str2int[k] for k in normalised], dtype=np.int64)

    n_classes = len(unique)
    n_rare = sum(1 for k in normalised if k == rare_id_str)
    logger.info(
        "Cluster labels: %d unique classes, %d rare reactions merged.",
        n_classes,
        n_rare,
    )
    return labels, str2int


def build(
    *,
    fp_size: int = 1024,
    embed_dim: int = 256,
    hidden_dim: int = 1024,
    epochs: int = 50,
    batch_size: int = 128,
    lr: float = 1e-3,
    train_kmn: bool = True,
    index_only: bool = False,
    output_dir: Optional[str] = None,
    n_processes: Optional[int] = 1,
    min_cluster_size: int = 5,
    use_ivf: bool = False,
) -> None:
    from chemtools.precedent.loader import _load_all_with_progress, _iter_literature_files
    from chemtools.util.mix_fingerprint import create_mix_fp_batch, mix_fp_dim
    from chemtools.util.kernel_metric_net import KernelMetricNetwork
    from chemtools.util.faiss_router import FAISSRouter

    # ── resolve output dir ────────────────────────────────────────────────
    if output_dir is None:
        _env = os.environ.get("CHEMTOOLS_KMN_INDEX_DIR", "").strip()
        output_dir = (
            os.path.abspath(_env)
            if _env
            else os.path.join(_REPO_ROOT, "data", "kmn_index")
        )
    os.makedirs(output_dir, exist_ok=True)
    weights_path = os.path.join(output_dir, "kmn_weights.pt")
    faiss_path = os.path.join(output_dir, "kmn_faiss.index")
    meta_path = os.path.join(output_dir, "kmn_metadata.npz")

    logger.info("Output directory: %s", output_dir)
    logger.info("fp_size=%d  embed_dim=%d  hidden_dim=%d", fp_size, embed_dim, hidden_dim)

    # ── 1. load all precedent rows ────────────────────────────────────────
    n_files = sum(1 for _ in _iter_literature_files())
    logger.info("Loading HTE precedent database (%d CSV files) …", n_files)
    rows = _load_all_with_progress()
    logger.info("  Loaded %d total rows.", len(rows))

    # Filter to rows with a valid reaction SMILES
    valid_rows = [r for r in rows if r.get("reaction_smiles")]
    logger.info("  %d rows have reaction SMILES.", len(valid_rows))

    if not valid_rows:
        logger.error("No valid reaction SMILES found. Aborting.")
        sys.exit(1)

    # ── 2. compute MixFP fingerprints ─────────────────────────────────────
    logger.info("Computing MixFP fingerprints (fp_size=%d) …", fp_size)
    smiles_list = [r["reaction_smiles"] for r in valid_rows]

    fps, successful_indices = create_mix_fp_batch(
        smiles_list,
        fp_size=fp_size,
        use_chirality=True,
        n_processes=n_processes,
        show_progress=True,
    )

    if fps.shape[0] == 0:
        logger.error("No fingerprints computed. Check RDKit installation.")
        sys.exit(1)

    # Align labels and reaction IDs to successfully computed rows
    computed_rows = [valid_rows[i] for i in successful_indices]
    reaction_ids = [r.get("reaction_id", f"row_{i}") for i, r in enumerate(computed_rows)]
    logger.info(
        "  %d / %d fingerprints computed successfully.",
        len(successful_indices),
        len(valid_rows),
    )
    logger.info("  MixFP matrix shape: %s", fps.shape)

    # ── 3. build cluster labels ────────────────────────────────────────────
    logger.info("Building cluster labels …")
    labels, label_map = _encode_cluster_labels(
        computed_rows, min_cluster_size=min_cluster_size
    )
    n_classes = int(labels.max()) + 1

    # Save label map for reference
    label_map_path = os.path.join(output_dir, "kmn_label_map.npz")
    np.savez_compressed(
        label_map_path,
        keys=np.array(list(label_map.keys()), dtype=object),
        values=np.array(list(label_map.values()), dtype=np.int64),
    )
    logger.info("  Label map saved → %s", label_map_path)

    # ── 4. train / load KMN ──────────────────────────────────────────────
    input_dim = mix_fp_dim(fp_size)  # 6 * fp_size
    kmn = KernelMetricNetwork(
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        embed_dim=embed_dim,
    )

    if index_only:
        # --index-only: reuse existing trained weights; only rebuild FAISS.
        # This is the correct mode for incremental dataset additions.
        if not os.path.exists(weights_path):
            logger.error(
                "--index-only requires existing weights at %s\n"
                "Run without --index-only first to train KMN.",
                weights_path,
            )
            sys.exit(1)
        kmn.load(weights_path)
        logger.info(
            "Loaded existing KMN weights (%s backend) from %s",
            kmn._backend,
            weights_path,
        )
    elif train_kmn:
        if kmn._backend != "torch":
            logger.warning(
                "PyTorch not available – KMN training skipped. "
                "Using random projection instead."
            )
        else:
            logger.info(
                "Training KMN: %d samples, %d classes, epochs=%d …",
                len(fps),
                n_classes,
                epochs,
            )
            history = kmn.train(
                fps,
                labels,
                epochs=epochs,
                lr=lr,
                batch_size=batch_size,
            )
            if history:
                last = history[-1]
                val_str = (
                    f"  val_acc={last['val_acc']:.3f}" if "val_acc" in last else ""
                )
                logger.info(
                    "Training complete. Final loss=%.4f%s",
                    last.get("train_loss", 0.0),
                    val_str,
                )
            kmn.save(weights_path)
            logger.info("KMN weights saved → %s", weights_path)
    else:
        logger.info("Skipping KMN training (--no-train); saving RP placeholder.")
        kmn.save(weights_path)  # saves RP projection matrix

    # ── 5. embed all reactions ────────────────────────────────────────────
    logger.info("Generating KMN embeddings for %d reactions …", len(fps))
    chunk_size = 1024
    chunks = range(0, len(fps), chunk_size)
    try:
        from tqdm import tqdm
        chunks = tqdm(chunks, desc="  Embed", unit="chunk", dynamic_ncols=True)
    except ImportError:
        pass
    embeds_list = []
    for start in chunks:
        chunk = fps[start : start + chunk_size]
        embeds_list.append(kmn.get_embeddings(chunk))
    embeddings = np.vstack(embeds_list).astype(np.float32)
    logger.info("  Embedding matrix shape: %s", embeddings.shape)

    # ── 6. build FAISS index ──────────────────────────────────────────────
    logger.info("Building FAISS index …")
    router = FAISSRouter(embed_dim=embed_dim, use_ivf=use_ivf)
    router.build(embeddings, reaction_ids, use_ivf=use_ivf)
    router.save(faiss_path, meta_path)

    # ── summary ───────────────────────────────────────────────────────────
    logger.info("")
    logger.info("╔══════════════════════════════════════════════╗")
    logger.info("║   KMN index build complete                   ║")
    logger.info("╠══════════════════════════════════════════════╣")
    logger.info("║  Reactions indexed : %-24d║", len(reaction_ids))
    logger.info("║  Embedding dim     : %-24d║", embed_dim)
    logger.info("║  KMN backend       : %-24s║", kmn._backend)
    logger.info("║  FAISS backend     : %-24s║", "faiss" if router._index is not None else "numpy")
    logger.info("╠══════════════════════════════════════════════╣")
    logger.info("║  kmn_weights.pt    → %s", os.path.basename(weights_path))
    logger.info("║  kmn_faiss.index   → %s", os.path.basename(faiss_path))
    logger.info("║  kmn_metadata.npz  → %s", os.path.basename(meta_path))
    logger.info("╚══════════════════════════════════════════════╝")
    logger.info("")
    logger.info(
        "To use: pass relax={'use_mixfp': True, 'reaction_smiles': '...'} "
        "to chemtools.precedent.knn()"
    )


# ── CLI ────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build MixFP + KMN + FAISS routing index from HTE precedents."
    )
    p.add_argument("--fp-size", type=int, default=1024,
                   help="Per-molecule fingerprint half-size (default: 1024). "
                        "Total MixFP dim = 6 * fp_size.")
    p.add_argument("--embed-dim", type=int, default=256,
                   help="KMN embedding dimension (default: 256).")
    p.add_argument("--hidden-dim", type=int, default=1024,
                   help="KMN hidden layer size (default: 1024).")
    p.add_argument("--epochs", type=int, default=50,
                   help="KMN training epochs (default: 50).")
    p.add_argument("--batch-size", type=int, default=128,
                   help="Mini-batch size (default: 128).")
    p.add_argument("--lr", type=float, default=1e-3,
                   help="Learning rate (default: 1e-3).")
    p.add_argument("--no-train", action="store_true",
                   help="Skip KMN training; use random projection fallback.")
    p.add_argument("--index-only", action="store_true",
                   help="""
                   Load existing KMN weights (do NOT retrain) and rebuild FAISS only.
                   Use this for incremental dataset additions when the reaction families
                   and condition landscape have not fundamentally changed.
                   Requires a previous full build to exist in --output-dir.
                   """.strip())
    p.add_argument("--use-ivf", action="store_true",
                   help="Use IVF FAISS index (faster for >10k reactions).")
    p.add_argument("--output-dir", type=str, default=None,
                   help="Directory for index files (default: data/kmn_index/).")
    p.add_argument("--n-processes", type=int, default=1,
                   help="Parallel workers for fingerprint computation (default: 1).")
    p.add_argument("--min-cluster-size", type=int, default=5,
                   help="Min reactions per condition-core cluster (default: 5). "
                        "Smaller clusters get merged into '__rare__'.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build(
        fp_size=args.fp_size,
        embed_dim=args.embed_dim,
        hidden_dim=args.hidden_dim,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        train_kmn=not args.no_train,
        index_only=args.index_only,
        output_dir=args.output_dir,
        n_processes=args.n_processes,
        min_cluster_size=args.min_cluster_size,
        use_ivf=args.use_ivf,
    )
