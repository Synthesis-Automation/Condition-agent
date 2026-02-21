"""FAISS-based routing for fast nearest-neighbour condition lookup.

Combines MixFP + KMN embeddings with a FAISS index for sub-millisecond
routing across the entire HTE precedent database.

Architecture summary:

    reaction_smiles
        ↓  mix_fingerprint.create_mix_fp()          [6144-dim]
        ↓  KernelMetricNetwork.get_embeddings()      [256-dim]
        ↓  FAISSRouter.search()                      top-k reaction IDs
        ↓  precedent/search.py lookup                ranked conditions

Graceful degradation:
- If ``faiss`` is not installed → brute-force cosine (numpy dot)
- If KMN weights absent        → random-projection embedding
- If index not built           → empty results (caller falls back to DRFP/cat)

Index files (default ``data/kmn_index/``):
    kmn_weights.pt        – KMN model weights (from kernel_metric_net.py)
    kmn_faiss.index       – serialised FAISS index
    kmn_metadata.npz      – normalised embeddings + reaction_id strings

Build with ``scripts/build_kmn_index.py``.
"""
from __future__ import annotations

import os
import logging
from typing import Optional, List, Tuple, Dict, Any

import numpy as np

logger = logging.getLogger(__name__)

# ── default paths ──────────────────────────────────────────────────────────
_ENV_INDEX_DIR = os.environ.get("CHEMTOOLS_KMN_INDEX_DIR", "").strip()
_BASE = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
INDEX_DIR: str = (
    os.path.abspath(_ENV_INDEX_DIR)
    if _ENV_INDEX_DIR
    else os.path.join(_BASE, "data", "kmn_index")
)
_FAISS_INDEX_PATH = os.path.join(INDEX_DIR, "kmn_faiss.index")
_META_PATH = os.path.join(INDEX_DIR, "kmn_metadata.npz")
_WEIGHTS_PATH = os.path.join(INDEX_DIR, "kmn_weights.pt")


# ── availability helpers ───────────────────────────────────────────────────

def _faiss_available() -> bool:
    try:
        import faiss  # noqa: F401
        return True
    except ImportError:
        return False


# ── core router class ──────────────────────────────────────────────────────

class FAISSRouter:
    """FAISS index mapping KMN embeddings → reaction IDs.

    Supports both exact (flat inner-product) and approximate (IVF) search.
    Falls back to numpy brute-force when *faiss* is unavailable.
    """

    def __init__(
        self,
        embed_dim: int = 256,
        *,
        use_ivf: bool = False,
        nlist: int = 128,
        nprobe: int = 16,
    ):
        self.embed_dim = embed_dim
        self.use_ivf = use_ivf
        self.nlist = nlist
        self.nprobe = nprobe
        self._index: Any = None              # FAISS index object
        self._reaction_ids: List[str] = []
        self._embeddings: Optional[np.ndarray] = None   # L2-normalised; fallback

    # ── build ──────────────────────────────────────────────────────────────

    def build(
        self,
        embeddings: np.ndarray,
        reaction_ids: List[str],
        *,
        use_ivf: Optional[bool] = None,
    ) -> None:
        """Populate in-memory FAISS index from an embedding matrix.

        Embeddings are L2-normalised internally so inner-product search
        is equivalent to cosine similarity.

        Args:
            embeddings:   (N, embed_dim) float32.
            reaction_ids: Parallel list of N reaction ID strings.
            use_ivf:      Override ``self.use_ivf``; auto-selects IVF for N > 10 000.
        """
        if embeddings.ndim != 2 or embeddings.shape[1] != self.embed_dim:
            raise ValueError(
                f"Expected embeddings shape (N, {self.embed_dim}), "
                f"got {embeddings.shape}"
            )

        emb = np.ascontiguousarray(embeddings, dtype=np.float32)

        # L2-normalise (cosine via inner product)
        norms = np.linalg.norm(emb, axis=1, keepdims=True)
        norms = np.where(norms == 0.0, 1.0, norms)
        emb_norm = emb / norms

        self._reaction_ids = list(reaction_ids)
        self._embeddings = emb_norm  # keep for brute-force fallback

        _use_ivf = use_ivf if use_ivf is not None else (
            self.use_ivf or len(reaction_ids) > 10_000
        )

        if _faiss_available():
            import faiss
            if _use_ivf and len(reaction_ids) >= self.nlist * 2:
                quantizer = faiss.IndexFlatIP(self.embed_dim)
                idx = faiss.IndexIVFFlat(
                    quantizer, self.embed_dim, self.nlist, faiss.METRIC_INNER_PRODUCT
                )
                idx.train(emb_norm)
                idx.add(emb_norm)
                idx.nprobe = self.nprobe
                self._index = idx
                logger.info(
                    "FAISSRouter: built IVF-%d index with %d vectors", self.nlist, len(reaction_ids)
                )
            else:
                idx = faiss.IndexFlatIP(self.embed_dim)
                idx.add(emb_norm)
                self._index = idx
                logger.info(
                    "FAISSRouter: built flat index with %d vectors", len(reaction_ids)
                )
        else:
            logger.warning(
                "faiss not installed – FAISSRouter will use numpy brute-force. "
                "Install faiss-cpu for faster routing."
            )

    # ── search ─────────────────────────────────────────────────────────────

    def search(
        self,
        query_embedding: np.ndarray,
        k: int = 20,
    ) -> Tuple[List[str], List[float]]:
        """Find k nearest reaction IDs (descending similarity).

        Args:
            query_embedding: (embed_dim,) or (1, embed_dim) float32.
            k: Number of results.

        Returns:
            Tuple of (reaction_ids, cosine_scores).
        """
        if not self._reaction_ids:
            return [], []

        q = query_embedding.astype(np.float32).reshape(1, -1)
        norm = float(np.linalg.norm(q))
        if norm > 0.0:
            q = q / norm

        k_clamped = min(k, len(self._reaction_ids))

        if self._index is not None and _faiss_available():
            D, I = self._index.search(q, k_clamped)
            valid = [(i, d) for i, d in zip(I[0], D[0]) if 0 <= i < len(self._reaction_ids)]
            ids = [self._reaction_ids[i] for i, _ in valid]
            scores = [float(d) for _, d in valid]
        elif self._embeddings is not None:
            sims = (self._embeddings @ q.T).ravel()
            top_idx = np.argsort(-sims)[:k_clamped]
            ids = [self._reaction_ids[i] for i in top_idx]
            scores = [float(sims[i]) for i in top_idx]
        else:
            return [], []

        return ids, scores

    # ── persistence ────────────────────────────────────────────────────────

    def save(
        self,
        index_path: str = _FAISS_INDEX_PATH,
        meta_path: str = _META_PATH,
    ) -> None:
        """Save FAISS index and metadata (embeddings + reaction IDs) to disk."""
        os.makedirs(os.path.dirname(os.path.abspath(index_path)), exist_ok=True)
        if self._index is not None and _faiss_available():
            import faiss
            faiss.write_index(self._index, index_path)
            logger.info("FAISSRouter: index saved → %s", index_path)
        if self._embeddings is not None:
            np.savez_compressed(
                meta_path,
                embeddings=self._embeddings,
                reaction_ids=np.array(self._reaction_ids, dtype=object),
                embed_dim=np.array([self.embed_dim]),
            )
            logger.info("FAISSRouter: metadata saved → %s", meta_path)

    @classmethod
    def load(
        cls,
        index_path: str = _FAISS_INDEX_PATH,
        meta_path: str = _META_PATH,
        *,
        embed_dim: int = 256,
        nprobe: int = 16,
    ) -> "FAISSRouter":
        """Load a previously saved router from disk.

        Returns an empty (non-ready) router if files are missing.
        """
        router = cls(embed_dim=embed_dim, nprobe=nprobe)

        if os.path.exists(meta_path):
            try:
                data = np.load(meta_path, allow_pickle=True)
                router._embeddings = data["embeddings"].astype(np.float32)
                router._reaction_ids = list(data["reaction_ids"])
                saved_dim = int(data.get("embed_dim", [embed_dim])[0])
                router.embed_dim = saved_dim
                logger.info(
                    "FAISSRouter: loaded %d embeddings ← %s",
                    len(router._reaction_ids),
                    meta_path,
                )
            except Exception as exc:
                logger.warning("FAISSRouter: could not load metadata from %s: %s", meta_path, exc)

        if os.path.exists(index_path) and _faiss_available():
            try:
                import faiss
                router._index = faiss.read_index(index_path)
                # Restore nprobe for IVF indices
                if hasattr(router._index, "nprobe"):
                    router._index.nprobe = nprobe
                logger.info("FAISSRouter: index loaded ← %s", index_path)
            except Exception as exc:
                logger.warning("FAISSRouter: could not load index from %s: %s", index_path, exc)

        return router

    @property
    def is_ready(self) -> bool:
        """True if the router can answer search queries."""
        return bool(self._reaction_ids) and (
            self._index is not None or self._embeddings is not None
        )

    def __len__(self) -> int:
        return len(self._reaction_ids)

    def __repr__(self) -> str:
        backend = "FAISS" if self._index is not None else "numpy"
        return (
            f"FAISSRouter(n={len(self._reaction_ids)}, "
            f"embed_dim={self.embed_dim}, backend={backend})"
        )


# ── module-level singleton helpers ─────────────────────────────────────────

_DEFAULT_ROUTER: Optional[FAISSRouter] = None


def get_default_router(
    *,
    index_path: Optional[str] = None,
    meta_path: Optional[str] = None,
    reload: bool = False,
) -> FAISSRouter:
    """Return (loading if needed) the module-level default FAISSRouter.

    On first call, loads index and metadata from the default ``data/kmn_index/``
    directory.  Subsequent calls return the cached instance.  Pass
    ``reload=True`` to force a fresh load.

    Returns an empty (non-ready) router if the index has not been built yet;
    callers should check ``router.is_ready`` before using.
    """
    global _DEFAULT_ROUTER

    _idx = index_path or _FAISS_INDEX_PATH
    _meta = meta_path or _META_PATH

    if _DEFAULT_ROUTER is None or reload:
        _DEFAULT_ROUTER = FAISSRouter.load(_idx, _meta)
        if not _DEFAULT_ROUTER.is_ready:
            logger.debug(
                "FAISSRouter not ready – run scripts/build_kmn_index.py to build the index."
            )

    return _DEFAULT_ROUTER


def route_reaction(
    reaction_smiles: str,
    k: int = 50,
    *,
    fp_size: int = 1024,
    kmn: Optional[Any] = None,
    router: Optional[FAISSRouter] = None,
) -> Tuple[List[str], List[float]]:
    """One-shot shortcut: MixFP → KMN embedding → FAISS top-k.

    Args:
        reaction_smiles: Input reaction SMILES.
        k: Number of nearest precedents to return.
        fp_size: MixFP dimension parameter.
        kmn: Optional pre-loaded ``KernelMetricNetwork``.  Default → singleton.
        router: Optional pre-loaded ``FAISSRouter``.  Default → singleton.

    Returns:
        Tuple of (reaction_ids, scores) or ([], []) if unavailable.
    """
    from .mix_fingerprint import create_mix_fp
    from .kernel_metric_net import get_default_kmn

    fp = create_mix_fp(reaction_smiles, fp_size=fp_size)
    if fp is None:
        return [], []

    _kmn = kmn if kmn is not None else get_default_kmn()
    _router = router if router is not None else get_default_router()

    if not _router.is_ready:
        return [], []

    emb = _kmn.get_embeddings(fp.reshape(1, -1))  # (1, embed_dim)
    return _router.search(emb[0], k=k)
