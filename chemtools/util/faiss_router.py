"""FAISS-based routing for fast nearest-neighbour condition lookup.

Combines MixFP + KMN embeddings with a FAISS index for sub-millisecond
routing across the entire HTE precedent database.

Architecture summary:

    reaction_smiles
        ↓  mix_fingerprint.create_mix_fp()          [6144-dim]
        ↓  KernelMetricNetwork.get_embeddings()      [256-dim]
        ↓  FAISSRouter.adaptive_search()             top-k reaction IDs
        ↓  precedent/search.py lookup                ranked conditions

Routing modes (selected automatically via ``adaptive_search``):

    High confidence (≥0.85)  →  hard family filter + flat/IVF, nprobe=4
    Medium (0.5–0.85)        →  IVF with moderate nprobe bleed across cells
    Low / unknown (<0.5)     →  full soft-Voronoi: no family gate, wide nprobe

Soft-Voronoi blending (``search_soft_voronoi``):
    Each result is re-scored as  cell_similarity × reaction_similarity
    so reactions in the closest Voronoi cell are upweighted even if their
    raw cosine score is beaten by a reaction in a weaker cell.

Graceful degradation:
- If ``faiss`` is not installed → brute-force cosine (numpy dot)
- If KMN weights absent        → random-projection embedding
- If index not built           → empty results (caller falls back to DRFP/cat)
- If index is flat (not IVF)   → soft-Voronoi falls back to plain search

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


def is_index_built(
    *,
    index_dir: Optional[str] = None,
) -> bool:
    """Return True if the KMN index files exist on disk and are non-empty.

    Does **not** load anything into memory – cheap filesystem check only.
    Used by :func:`chemtools.precedent.search.knn` to auto-enable MixFP
    routing without requiring callers to pass ``use_mixfp=True`` explicitly.
    """
    d = index_dir or INDEX_DIR
    meta = os.path.join(d, "kmn_metadata.npz")
    weights = os.path.join(d, "kmn_weights.pt")
    return os.path.isfile(meta) and os.path.getsize(meta) > 0 and os.path.isfile(weights)


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
        # IVF-only: cell assignment for each reaction (shape N,), dtype int32
        self._cell_assignments: Optional[np.ndarray] = None

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
                # Store per-reaction cell assignments for soft-Voronoi blending
                _, cell_ids = quantizer.search(emb_norm, 1)  # (N, 1)
                self._cell_assignments = cell_ids.ravel().astype(np.int32)
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

    # ── search (plain) ─────────────────────────────────────────────────────

    def search(
        self,
        query_embedding: np.ndarray,
        k: int = 20,
        *,
        nprobe: Optional[int] = None,
    ) -> Tuple[List[str], List[float]]:
        """Find k nearest reaction IDs (descending similarity).

        Args:
            query_embedding: (embed_dim,) or (1, embed_dim) float32.
            k: Number of results.
            nprobe: Override index nprobe for this query only (IVF only).

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
            # Temporarily override nprobe when requested
            old_nprobe = None
            if nprobe is not None and hasattr(self._index, "nprobe"):
                old_nprobe = self._index.nprobe
                self._index.nprobe = nprobe
            D, I = self._index.search(q, k_clamped)
            if old_nprobe is not None:
                self._index.nprobe = old_nprobe
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

    # ── soft-Voronoi search ────────────────────────────────────────────────

    def search_soft_voronoi(
        self,
        query_embedding: np.ndarray,
        k: int = 20,
        *,
        nprobe: Optional[int] = None,
        cell_weight: float = 0.3,
    ) -> Tuple[List[str], List[float]]:
        """Voronoi-cell-weighted search for IVF indices.

        Each candidate reaction is re-scored as::

            blended = (1 - cell_weight) * reaction_sim
                    +      cell_weight  * cell_centroid_sim

        This upweights reactions from the nearest Voronoi cells, so that
        reactions at family boundaries receive contributions from multiple
        specialist clusters rather than only the closest one.

        Falls back to plain ``search()`` when the index is not IVF or
        cell assignments are unavailable.

        Args:
            query_embedding: (embed_dim,) or (1, embed_dim) float32.
            k: Final number of results to return.
            nprobe: Number of Voronoi cells to probe.  Defaults to
                ``self.nprobe``.  Higher values increase recall at the
                cost of latency.
            cell_weight: Weight given to cell centroid similarity vs raw
                reaction similarity (0 = pure reaction sim, 1 = pure cell
                sim).  Recommended range: 0.2–0.4.

        Returns:
            Tuple of (reaction_ids, blended_scores).
        """
        if not self._reaction_ids:
            return [], []

        # Fall back to plain search if IVF data not available
        is_ivf = (
            self._index is not None
            and _faiss_available()
            and self._cell_assignments is not None
            and hasattr(self._index, "quantizer")
        )
        if not is_ivf:
            return self.search(query_embedding, k=k, nprobe=nprobe)

        import faiss  # guaranteed available (is_ivf=True)

        q = query_embedding.astype(np.float32).reshape(1, -1)
        norm = float(np.linalg.norm(q))
        if norm > 0.0:
            q = q / norm

        _nprobe = nprobe if nprobe is not None else self.nprobe

        # 1. Get centroid similarities for each probed cell
        centroid_D, centroid_I = self._index.quantizer.search(q, _nprobe)  # (1, nprobe)
        cell_sim: Dict[int, float] = {
            int(cid): float(csim)
            for cid, csim in zip(centroid_I[0], centroid_D[0])
            if cid >= 0
        }

        # 2. Retrieve a larger candidate pool (k * 2 or at least 50)
        pool_k = min(max(k * 2, 50), len(self._reaction_ids))
        old_nprobe = self._index.nprobe
        self._index.nprobe = _nprobe
        D, I = self._index.search(q, pool_k)
        self._index.nprobe = old_nprobe

        # 3. Re-score each candidate using its cell's centroid similarity
        candidates: List[Tuple[str, float]] = []
        for raw_idx, raw_sim in zip(I[0], D[0]):
            if raw_idx < 0 or raw_idx >= len(self._reaction_ids):
                continue
            cid = int(self._cell_assignments[raw_idx])
            csim = cell_sim.get(cid, 0.0)
            blended = (1.0 - cell_weight) * float(raw_sim) + cell_weight * csim
            candidates.append((self._reaction_ids[raw_idx], blended))

        # 4. Sort by blended score and return top-k
        candidates.sort(key=lambda x: x[1], reverse=True)
        top = candidates[:k]
        return [r for r, _ in top], [s for _, s in top]

    # ── adaptive search (confidence-gated) ────────────────────────────────

    def adaptive_search(
        self,
        query_embedding: np.ndarray,
        k: int = 20,
        *,
        family_confidence: float = 1.0,
        candidate_ids: Optional[List[int]] = None,
        cell_weight: float = 0.3,
    ) -> Tuple[List[str], List[float], Dict[str, Any]]:
        """Confidence-gated Voronoi search.

        Adjusts ``nprobe`` and blending continuously based on how certain
        the reaction-type detector is:

        +--------------------+----------+---------------------------+
        | family_confidence  | nprobe   | mode                      |
        +====================+==========+===========================+
        | ≥ 0.85             | 4        | hard filter (fast/tight)  |
        | 0.50 – 0.85        | 8 – 32   | soft bleed across cells   |
        | < 0.50             | ≥ 48     | full soft-Voronoi, no gate|
        +--------------------+----------+---------------------------+

        Args:
            query_embedding: (embed_dim,) or (1, embed_dim) float32.
            k: Final number of results to return.
            family_confidence: Confidence score from reaction-type detector
                (0.0 = unknown, 1.0 = certain).
            candidate_ids: Optional list of integer indices into the index to
                restrict the search pool (used for hard family gating when
                confidence is high).  When None, full index is searched.
            cell_weight: Blend weight passed to ``search_soft_voronoi``.

        Returns:
            Tuple of::

                (reaction_ids, scores, meta)

            where ``meta`` is a dict with keys
            ``{"nprobe", "mode", "family_confidence", "n_candidates"}``.
        """
        # Compute adaptive nprobe: high confidence → small nprobe (tight)
        # Low confidence → large nprobe (wide Voronoi bleed)
        max_np = self.nlist if self._cell_assignments is not None else 1
        nprobe = max(4, int(max_np * (1.0 - family_confidence) * 0.5) + 1)
        nprobe = min(nprobe, max_np)

        if family_confidence >= 0.85:
            mode = "hard_filter"
        elif family_confidence >= 0.50:
            mode = "soft_bleed"
        else:
            mode = "full_voronoi"

        # Hard-filter: restrict to provided candidate pool
        if mode == "hard_filter" and candidate_ids is not None:
            # Build a temporary numpy-backed sub-search over candidate_ids only
            ids, scores = self._search_subset(query_embedding, k, candidate_ids)
        else:
            # Soft or full Voronoi
            ids, scores = self.search_soft_voronoi(
                query_embedding, k=k, nprobe=nprobe, cell_weight=cell_weight
            )

        meta: Dict[str, Any] = {
            "nprobe": nprobe,
            "mode": mode,
            "family_confidence": family_confidence,
            "n_candidates": len(ids),
        }
        return ids, scores, meta

    # ── subset search (hard-filter helper) ────────────────────────────────

    def _search_subset(
        self,
        query_embedding: np.ndarray,
        k: int,
        candidate_indices: List[int],
    ) -> Tuple[List[str], List[float]]:
        """Brute-force cosine search restricted to a subset of the index.

        Used when family confidence is high and we want to search only within
        a pre-filtered pool (e.g. Suzuki reactions only).
        """
        if self._embeddings is None or not candidate_indices:
            return [], []

        q = query_embedding.astype(np.float32).reshape(-1)
        norm = float(np.linalg.norm(q))
        if norm > 0.0:
            q = q / norm

        idx_arr = np.array(candidate_indices, dtype=np.int64)
        # Clamp to valid range
        idx_arr = idx_arr[(idx_arr >= 0) & (idx_arr < len(self._reaction_ids))]
        if len(idx_arr) == 0:
            return [], []

        sub_emb = self._embeddings[idx_arr]   # (M, embed_dim)
        sims = (sub_emb @ q).ravel()          # (M,)
        k_clamped = min(k, len(idx_arr))
        top_local = np.argsort(-sims)[:k_clamped]
        top_global = idx_arr[top_local]

        ids = [self._reaction_ids[i] for i in top_global]
        scores = [float(sims[i]) for i in top_local]
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
                # persist cell assignments so load() restores soft-Voronoi capability
                cell_assignments=(
                    self._cell_assignments
                    if self._cell_assignments is not None
                    else np.array([], dtype=np.int32)
                ),
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
                # Restore cell assignments if present
                if "cell_assignments" in data and len(data["cell_assignments"]) > 0:
                    router._cell_assignments = data["cell_assignments"].astype(np.int32)
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
        ivf = (
            f"/IVF{self.nlist}"
            if self._cell_assignments is not None
            else ""
        )
        return (
            f"FAISSRouter(n={len(self._reaction_ids)}, "
            f"embed_dim={self.embed_dim}, backend={backend}{ivf})"
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
    family_confidence: float = 1.0,
    candidate_ids: Optional[List[int]] = None,
    cell_weight: float = 0.3,
    kmn: Optional[Any] = None,
    router: Optional[FAISSRouter] = None,
) -> Tuple[List[str], List[float], Dict[str, Any]]:
    """One-shot shortcut: MixFP → KMN embedding → adaptive Voronoi top-k.

    Args:
        reaction_smiles: Input reaction SMILES.
        k: Number of nearest precedents to return.
        fp_size: MixFP dimension parameter.
        family_confidence: Reaction-type detector confidence (0–1).
            Controls nprobe and hard/soft/full-Voronoi mode.
        candidate_ids: Optional integer indices for hard-filter pool
            (used when family_confidence ≥ 0.85).
        cell_weight: Blend weight for soft-Voronoi scoring (0.2–0.4).
        kmn: Optional pre-loaded ``KernelMetricNetwork``.  Default → singleton.
        router: Optional pre-loaded ``FAISSRouter``.  Default → singleton.

    Returns:
        Tuple of (reaction_ids, scores, meta) or ([], [], {}) if unavailable.
    """
    from .mix_fingerprint import create_mix_fp
    from .kernel_metric_net import get_default_kmn

    fp = create_mix_fp(reaction_smiles, fp_size=fp_size)
    if fp is None:
        return [], [], {}

    _kmn = kmn if kmn is not None else get_default_kmn()
    _router = router if router is not None else get_default_router()

    if not _router.is_ready:
        return [], [], {}

    emb = _kmn.get_embeddings(fp.reshape(1, -1))  # (1, embed_dim)
    return _router.adaptive_search(
        emb[0],
        k=k,
        family_confidence=family_confidence,
        candidate_ids=candidate_ids,
        cell_weight=cell_weight,
    )
