"""Kernel Metric Network (KMN) for reaction embedding.

Implements a lightweight MLP that maps mixed reaction fingerprints (from
``chemtools.util.mix_fingerprint``) to a compact learned metric space
optimised for routing to the correct condition cluster.

Architecture (default):
    Input(6144) → Linear(1024) → BN → ReLU → Dropout
                → Linear(512)  → BN → ReLU → Dropout
                → Linear(256)  [embedding output]

Training uses a classification auxiliary head (discarded at inference) over
condition-core cluster labels from the HTE database.  This acts as a proxy
for the Voronoi cell labels used in MOSAIC (Li et al., Nature 2026).

If PyTorch is unavailable, ``KernelMetricNetwork.get_embeddings`` falls
back to a signed-random-projection (Johnson–Lindenstrauss) for dimensionality
reduction, keeping the rest of the pipeline functional without GPU/torch.
"""
from __future__ import annotations

import os
import logging
from typing import Optional, List, Dict, Any

import numpy as np

logger = logging.getLogger(__name__)

# Default architecture hyper-parameters
INPUT_DIM: int = 6144   # mix_fp_dim(1024)
HIDDEN_DIM: int = 1024
EMBED_DIM: int = 256
DROPOUT: float = 0.10


# ── optional torch imports ─────────────────────────────────────────────────

def _torch_available() -> bool:
    try:
        import torch  # noqa: F401
        return True
    except ImportError:
        return False


# ── Random Projection fallback (no torch required) ─────────────────────────

class _RPFallback:
    """Signed random projection: cheap, deterministic O(1) setup."""

    def __init__(self, input_dim: int = INPUT_DIM, embed_dim: int = EMBED_DIM,
                 seed: int = 42):
        rng = np.random.default_rng(seed)
        self._proj = rng.choice([-1.0, 1.0], size=(input_dim, embed_dim)).astype(np.float32)
        self._proj /= np.sqrt(embed_dim)
        self.embed_dim = embed_dim
        self.is_trained = False

    def get_embeddings(self, x: np.ndarray) -> np.ndarray:
        if x.ndim == 1:
            x = x.reshape(1, -1)
        return (x.astype(np.float32) @ self._proj)

    def save(self, path: str) -> None:
        np.save(path, self._proj)

    def load(self, path: str) -> None:
        self._proj = np.load(path).astype(np.float32)
        self.embed_dim = self._proj.shape[1]
        self.is_trained = True


# ── PyTorch model ──────────────────────────────────────────────────────────

def _build_torch_model(input_dim: int, hidden_dim: int, embed_dim: int,
                       dropout: float) -> "torch.nn.Module":
    import torch.nn as nn
    return nn.Sequential(
        nn.Linear(input_dim, hidden_dim),
        nn.BatchNorm1d(hidden_dim),
        nn.ReLU(),
        nn.Dropout(dropout),
        nn.Linear(hidden_dim, hidden_dim // 2),
        nn.BatchNorm1d(hidden_dim // 2),
        nn.ReLU(),
        nn.Dropout(dropout),
        nn.Linear(hidden_dim // 2, embed_dim),
    )


class KernelMetricNetwork:
    """Thin wrapper for KMN that works with or without PyTorch.

    Usage::

        kmn = KernelMetricNetwork()          # auto-detects torch
        embeds = kmn.get_embeddings(fps)     # (N, 256) float32

    Training::

        history = kmn.train(fps, cluster_labels, epochs=50)
        kmn.save("kmn_weights.pt")

    Inference on untrained (random projection):: 

        kmn = KernelMetricNetwork()
        kmn.get_embeddings(fps)  # uses RP fallback until trained

    """

    def __init__(
        self,
        input_dim: int = INPUT_DIM,
        hidden_dim: int = HIDDEN_DIM,
        embed_dim: int = EMBED_DIM,
        dropout: float = DROPOUT,
    ):
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.embed_dim = embed_dim
        self.dropout = dropout
        self.is_trained: bool = False

        if _torch_available():
            import torch
            self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            self._model = _build_torch_model(input_dim, hidden_dim, embed_dim, dropout)
            self._model = self._model.to(self._device)
            self._backend = "torch"
            logger.debug("KMN: using PyTorch backend on %s", self._device)
        else:
            self._rp = _RPFallback(input_dim, embed_dim)
            self._backend = "rp"
            logger.warning(
                "PyTorch not available; KMN will use random projection. "
                "Install torch to enable learned routing."
            )

    # ── inference ──────────────────────────────────────────────────────────

    def get_embeddings(self, x: np.ndarray) -> np.ndarray:
        """Embed fingerprints into metric space.

        Args:
            x: float32 array of shape (N, input_dim) or (input_dim,).

        Returns:
            float32 array of shape (N, embed_dim).
        """
        if x.ndim == 1:
            x = x.reshape(1, -1)

        if self._backend == "rp":
            return self._rp.get_embeddings(x)

        # Torch path
        import torch
        self._model.eval()
        with torch.no_grad():
            t = torch.from_numpy(x.astype(np.float32)).to(self._device)
            out = self._model(t)
            return out.cpu().numpy()

    # ── training ───────────────────────────────────────────────────────────

    def train(
        self,
        fps: np.ndarray,
        labels: np.ndarray,
        *,
        epochs: int = 50,
        lr: float = 1e-3,
        batch_size: int = 128,
        val_split: float = 0.1,
    ) -> List[Dict[str, float]]:
        """Train KMN with a classification auxiliary head over condition clusters.

        The auxiliary head is discarded after training; only the embedding
        backbone is retained.

        Args:
            fps: (N, input_dim) float32 mixed fingerprints.
            labels: (N,) integer cluster labels (e.g., condition-core cluster IDs).
            epochs: Number of training epochs.
            lr: Learning rate.
            batch_size: Mini-batch size.
            val_split: Fraction held out for validation.

        Returns:
            List of per-epoch dicts with ``train_loss`` (and ``val_acc`` when applicable).
        """
        if self._backend == "rp":
            logger.warning("KMN: no PyTorch – skipping training; using random projection.")
            return []

        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset, random_split

        n_classes = int(labels.max()) + 1
        logger.info("KMN training: %d samples, %d classes, %d epochs", len(fps), n_classes, epochs)

        X = torch.from_numpy(fps.astype(np.float32))
        y = torch.from_numpy(labels.astype(np.int64))
        dataset = TensorDataset(X, y)

        # Validation split
        n_val = max(1, int(len(dataset) * val_split)) if val_split > 0 else 0
        n_train = len(dataset) - n_val
        if n_val > 0:
            train_ds, val_ds = random_split(dataset, [n_train, n_val])
        else:
            train_ds, val_ds = dataset, None

        train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, drop_last=False)

        # Auxiliary classification head (discarded at inference)
        classifier = nn.Linear(self.embed_dim, n_classes).to(self._device)
        criterion = nn.CrossEntropyLoss()

        optimizer = torch.optim.Adam(
            list(self._model.parameters()) + list(classifier.parameters()),
            lr=lr,
            weight_decay=1e-4,
        )
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

        history: List[Dict[str, float]] = []
        best_acc = 0.0

        for epoch in range(epochs):
            # ── train ──
            self._model.train()
            classifier.train()
            epoch_loss = 0.0
            n_batches = 0
            for xb, yb in train_loader:
                xb, yb = xb.to(self._device), yb.to(self._device)
                optimizer.zero_grad()
                embed = self._model(xb)
                logits = classifier(embed)
                loss = criterion(logits, yb)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self._model.parameters(), 1.0)
                optimizer.step()
                epoch_loss += loss.item()
                n_batches += 1
            scheduler.step()

            avg_loss = epoch_loss / max(1, n_batches)
            entry: Dict[str, float] = {"epoch": epoch + 1, "train_loss": avg_loss}

            # ── validate ──
            if val_ds is not None and len(val_ds) > 0:
                self._model.eval()
                classifier.eval()
                with torch.no_grad():
                    xv = torch.stack([val_ds[i][0] for i in range(len(val_ds))]).to(self._device)
                    yv = torch.stack([val_ds[i][1] for i in range(len(val_ds))]).to(self._device)
                    preds = classifier(self._model(xv)).argmax(dim=1)
                    acc = (preds == yv).float().mean().item()
                entry["val_acc"] = acc
                if acc > best_acc:
                    best_acc = acc

            history.append(entry)
            if (epoch + 1) % 10 == 0 or (epoch + 1) == epochs:
                val_str = f"  val_acc={entry.get('val_acc', 0):.3f}" if "val_acc" in entry else ""
                logger.info("Epoch %d/%d  loss=%.4f%s", epoch + 1, epochs, avg_loss, val_str)

        self._model.eval()
        self.is_trained = True
        logger.info("KMN training complete. best_val_acc=%.3f", best_acc)
        return history

    # ── persistence ────────────────────────────────────────────────────────

    def save(self, path: str) -> None:
        """Save model weights to disk.

        Torch models are saved as ``.pt``; RP fallback as ``.npy``.
        """
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        if self._backend == "torch":
            import torch
            torch.save(
                {
                    "state_dict": self._model.state_dict(),
                    "input_dim": self.input_dim,
                    "hidden_dim": self.hidden_dim,
                    "embed_dim": self.embed_dim,
                    "dropout": self.dropout,
                    "is_trained": self.is_trained,
                },
                path,
            )
            logger.info("KMN weights saved → %s", path)
        else:
            self._rp.save(path.replace(".pt", ".npy"))

    def load(self, path: str) -> None:
        """Load model weights from disk."""
        if self._backend == "torch":
            import torch
            ckpt = torch.load(path, map_location=self._device)
            # Re-build model if architecture changed
            state = ckpt.get("state_dict", ckpt)
            try:
                self._model.load_state_dict(state)
            except RuntimeError:
                # Architecture mismatch – rebuild from saved dims
                in_d = ckpt.get("input_dim", self.input_dim)
                h_d = ckpt.get("hidden_dim", self.hidden_dim)
                e_d = ckpt.get("embed_dim", self.embed_dim)
                dr = ckpt.get("dropout", self.dropout)
                self._model = _build_torch_model(in_d, h_d, e_d, dr).to(self._device)
                self._model.load_state_dict(state)
                self.input_dim = in_d
                self.hidden_dim = h_d
                self.embed_dim = e_d
                self.dropout = dr
            self.is_trained = ckpt.get("is_trained", True)
            self._model.eval()
            logger.info("KMN weights loaded ← %s", path)
        else:
            rp_path = path.replace(".pt", ".npy")
            if os.path.exists(rp_path):
                self._rp.load(rp_path)
                self.is_trained = True


# ── module-level singleton helpers ─────────────────────────────────────────

_DEFAULT_KMN: Optional[KernelMetricNetwork] = None
_DEFAULT_KMN_PATH: Optional[str] = None


def get_default_kmn(
    weights_path: Optional[str] = None,
    *,
    reload: bool = False,
) -> KernelMetricNetwork:
    """Return (loading if needed) the module-level default KMN instance.

    Args:
        weights_path: Path to ``.pt`` file.  Defaults to
            ``data/kmn_index/kmn_weights.pt`` relative to the repo root.
        reload: Force reload even if already cached.

    Returns:
        A ``KernelMetricNetwork`` ready for ``get_embeddings()``.
    """
    global _DEFAULT_KMN, _DEFAULT_KMN_PATH

    if weights_path is None:
        _base = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        env_dir = os.environ.get("CHEMTOOLS_KMN_INDEX_DIR", "").strip()
        index_dir = os.path.abspath(env_dir) if env_dir else os.path.join(_base, "data", "kmn_index")
        weights_path = os.path.join(index_dir, "kmn_weights.pt")

    if _DEFAULT_KMN is None or reload or weights_path != _DEFAULT_KMN_PATH:
        _DEFAULT_KMN = KernelMetricNetwork()
        if os.path.exists(weights_path):
            try:
                _DEFAULT_KMN.load(weights_path)
                logger.info("Loaded KMN weights from %s", weights_path)
            except Exception as exc:
                logger.warning("Could not load KMN weights from %s: %s", weights_path, exc)
        else:
            logger.debug("No KMN weights found at %s; using RP fallback.", weights_path)
        _DEFAULT_KMN_PATH = weights_path

    return _DEFAULT_KMN
