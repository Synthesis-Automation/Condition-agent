"""Mixed reaction fingerprint: [reactant_fp | diff_fp | product_fp].

Replicates MOSAIC's ``create_rxn_Mix_FP`` strategy (Li et al., Nature 2026).
Each molecule fingerprint is formed by concatenating:
  - RDKit path fingerprint  (maxPath=4, ``fp_size`` bits)
  - Morgan radius-2 bit-vector fingerprint (``fp_size`` bits)

giving a per-molecule vector of 2 * fp_size floats.

The reaction fingerprint is then:
    [reactant_fp | (product_fp - reactant_fp) | product_fp]

giving a total size of 6 * fp_size (default 6 * 1024 = 6144 dims).

This vector feeds into the Kernel Metric Network (KMN) for a learned
routing to the best condition cluster / FAISS expert cell.

If RDKit is not installed, all functions return None gracefully.
"""
from __future__ import annotations

from typing import Optional, Tuple
import logging

import numpy as np  # always available; used in annotations and function bodies

logger = logging.getLogger(__name__)

# Default per-molecule dimension (each half is this many bits from each FP type)
DEFAULT_FP_SIZE: int = 1024
# Total output dim = 3 segments × 2 FP types × fp_size = 6 * fp_size
# Default total = 6144

# ── SMARTS cache (avoid repeated compilation) ──────────────────────────────
try:
    from .smarts_cache import compile_smarts  # noqa: F401 (optional)
except Exception:
    pass


def _rdkit_available() -> bool:
    try:
        from rdkit import Chem  # noqa: F401
        return True
    except ImportError:
        return False


def _mol_fp(mol, fp_size: int, use_chirality: bool) -> np.ndarray:
    """Compute RDKit-path(maxPath=4) + Morgan(r=2) fingerprint for one molecule.

    Returns float32 numpy array of shape (2 * fp_size,).
    """
    from rdkit.Chem import AllChem, DataStructs

    # RDKit path fingerprint (maxPath=4)
    fpgen_rdk = AllChem.GetRDKitFPGenerator(maxPath=4, fpSize=fp_size)
    rdkit_fp = np.array(fpgen_rdk.GetFingerprint(mol).ToList(), dtype=np.float32)

    # Morgan radius-2 bit-vector (new generator API avoids deprecation warning)
    fpgen_morgan = AllChem.GetMorganGenerator(radius=2, fpSize=fp_size,
                                              includeChirality=use_chirality)
    morgan_fp = np.array(fpgen_morgan.GetFingerprint(mol).ToList(), dtype=np.float32)

    return np.concatenate((rdkit_fp, morgan_fp), axis=-1)  # shape: (2 * fp_size,)


def create_mix_fp(
    reaction_smiles: str,
    fp_size: int = DEFAULT_FP_SIZE,
    use_chirality: bool = True,
) -> Optional[np.ndarray]:
    """Create mixed reaction fingerprint [reactant_fp | diff_fp | product_fp].

    Args:
        reaction_smiles: Reaction SMILES in ``reactants>>products`` or
            ``reactants>agents>products`` format.
        fp_size: Per-FP-type bit length (default 1024). Each molecule's
            fingerprint has 2 * fp_size dimensions; total output is 6 * fp_size.
        use_chirality: Whether to include chirality in Morgan fingerprint.

    Returns:
        float32 numpy array of shape (6 * fp_size,) with order
        [reactant_fp | diff_fp | product_fp], or None on failure.
    """
    if not _rdkit_available():
        return None
    try:
        from rdkit import Chem

        # Split reaction: reactants >> products  or  reactants > agents > products
        parts = (reaction_smiles or "").strip().split(">")
        if len(parts) == 3:
            rsmi, _, psmi = parts
        elif ">>" in reaction_smiles:
            rsmi, psmi = reaction_smiles.split(">>", 1)
        else:
            return None

        # Handle multi-component reactant/product SMILES by joining with '.'
        rsmi = rsmi.strip()
        psmi = psmi.strip()
        if not rsmi or not psmi:
            return None

        rct_mol = Chem.MolFromSmiles(rsmi)
        prd_mol = Chem.MolFromSmiles(psmi)
        if rct_mol is None or prd_mol is None:
            return None

        rfp = _mol_fp(rct_mol, fp_size, use_chirality)   # (2*fp_size,)
        pfp = _mol_fp(prd_mol, fp_size, use_chirality)   # (2*fp_size,)
        diff_fp = pfp - rfp                               # (2*fp_size,) float diff

        # [reactant | diff | product]  — order matches MOSAIC concatenation
        return np.concatenate((rfp, diff_fp, pfp), axis=-1)  # (6*fp_size,)

    except Exception as exc:
        logger.debug("create_mix_fp failed for %r: %s", reaction_smiles, exc)
        return None


def create_mix_fp_batch(
    reaction_smiles_list: list[str],
    fp_size: int = DEFAULT_FP_SIZE,
    use_chirality: bool = True,
    *,
    n_processes: Optional[int] = None,
    show_progress: bool = False,
) -> Tuple[np.ndarray, list[int]]:
    """Create mixed fingerprints for a batch of reaction SMILES.

    Args:
        reaction_smiles_list: List of reaction SMILES strings.
        fp_size: Per-FP-type bit length (default 1024).
        use_chirality: Chirality flag forwarded to ``create_mix_fp``.
        n_processes: Number of worker processes.  None → use all CPU cores - 1.

    Returns:
        Tuple of:
          - ndarray of shape (M, 6 * fp_size) float32 for the M successful rows
          - list of original indices that succeeded (length M)
    """
    import numpy as np

    total = fp_size * 6
    features: list = []
    successful: list[int] = []

    # ── tqdm helper (graceful fallback to periodic logging) ───────────────
    def _iter_with_progress(iterable, total_n: int):
        if show_progress:
            try:
                from tqdm import tqdm
                yield from tqdm(iterable, total=total_n, unit="rxn",
                                desc="  MixFP", dynamic_ncols=True)
                return
            except ImportError:
                pass
            # fallback: log every 10 %
            report_every = max(1, total_n // 10)
            for idx, item in enumerate(iterable):
                if idx % report_every == 0:
                    logger.info("  MixFP  %d / %d (%.0f %%)",
                                idx, total_n, 100.0 * idx / total_n)
                yield item
            return
        yield from iterable

    if n_processes is not None and n_processes > 1:
        # Parallel path via multiprocessing (best-effort; fall back on error)
        try:
            from multiprocessing import Pool, cpu_count
            workers = n_processes if n_processes > 0 else max(1, cpu_count() - 1)
            with Pool(processes=workers) as pool:
                results = pool.starmap(
                    create_mix_fp,
                    [(rsmi, fp_size, use_chirality) for rsmi in reaction_smiles_list],
                )
            for i, fp in enumerate(results):
                if fp is not None and fp.shape == (total,):
                    features.append(fp)
                    successful.append(i)
            return np.array(features, dtype=np.float32), successful
        except Exception as exc:
            logger.warning("Parallel batch failed (%s); falling back to serial.", exc)

    # Serial path
    items = list(enumerate(reaction_smiles_list))
    for i, rsmi in _iter_with_progress(items, len(items)):
        fp = create_mix_fp(rsmi, fp_size=fp_size, use_chirality=use_chirality)
        if fp is not None and fp.shape == (total,):
            features.append(fp)
            successful.append(i)

    if not features:
        return np.empty((0, total), dtype=np.float32), []
    return np.stack(features, axis=0).astype(np.float32), successful


def mix_fp_dim(fp_size: int = DEFAULT_FP_SIZE) -> int:
    """Return the total dimension of a mix fingerprint (6 * fp_size)."""
    return 6 * fp_size
