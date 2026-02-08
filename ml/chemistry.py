"""Chemistry-aware featurization primitives for Chan-Lam ML pipeline."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from chemtools.featurizers.formatters.molecule import build_molecule_bundle
from chemtools.util.functional_groups import get_functional_groups


BORON_TOKEN_RE = re.compile(r"B\(|\[B")
TOKEN_SPLIT_RE = re.compile(r"[|/,;]+")

CONDITION_COLUMNS = [
    "catalyst",
    "ligand",
    "base",
    "acid",
    "oxidant",
    "reductant",
    "additive",
    "condensation_agent",
    "other_reagent",
    "solvent",
]

TEXT_COLUMNS = [
    "formed_motifs_tokens",
    "spectator_groups_tokens",
    "sulf_motif_tokens",
    "bor_motif_tokens",
]

NUMERIC_COLUMNS = [
    "sulf_motif_count",
    "sulf_aryl_steric_max",
    "sulf_alkyl_steric_max",
    "sulf_aryl_electronic_avg",
    "bor_motif_count",
    "bor_aryl_steric_max",
    "bor_alkyl_steric_max",
    "bor_aryl_electronic_avg",
]


def tokenize_text_field(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    tokens = [tok.strip().replace(" ", "_") for tok in TOKEN_SPLIT_RE.split(text) if tok.strip()]
    if not tokens:
        return ""
    return " ".join(sorted(set(tokens)))


def extract_reaction_substrates(reaction_smiles: str) -> Tuple[Optional[str], Optional[str]]:
    if not isinstance(reaction_smiles, str) or ">>" not in reaction_smiles:
        return None, None
    left = reaction_smiles.split(">>", 1)[0]
    reactants = [r.strip() for r in left.split(".") if r.strip()]
    if not reactants:
        return None, None
    boronic = None
    non_boron: List[str] = []
    for token in reactants:
        if BORON_TOKEN_RE.search(token):
            boronic = token if boronic is None else boronic
        else:
            non_boron.append(token)
    sulfonamide = None
    for token in non_boron:
        if "S(=O)(=O)" in token and "N" in token:
            sulfonamide = token
            break
    if sulfonamide is None and non_boron:
        sulfonamide = non_boron[0]
    return sulfonamide, boronic


def collect_motif_ids(bundle: Dict[str, Any]) -> List[str]:
    motif_ids: List[str] = []
    for key in ("motifs", "context_motifs"):
        for entry in bundle.get(key, []) or []:
            if not isinstance(entry, dict):
                continue
            motif_id = str(entry.get("compound_id") or entry.get("id") or "").strip()
            if motif_id:
                motif_ids.append(motif_id)
    return sorted(set(motif_ids))


def extract_score_values(payload: Any) -> List[float]:
    scores: List[float] = []
    if payload is None:
        return scores
    if isinstance(payload, dict):
        val = payload.get("score_0_10")
        if isinstance(val, (int, float)):
            scores.append(float(val))
        return scores
    if isinstance(payload, list):
        for item in payload:
            scores.extend(extract_score_values(item))
    return scores


def safe_mean(values: Sequence[float], default: float = 0.0) -> float:
    return default if not values else float(np.mean(values))


def safe_max(values: Sequence[float], default: float = 0.0) -> float:
    return default if not values else float(np.max(values))


class DescriptorBuilder:
    """Cache-aware descriptor extraction for sulfonamide/boronic partners."""

    def __init__(self) -> None:
        self._cache: Dict[str, Dict[str, Any]] = {}

    def _describe_smiles(self, smiles: Optional[str], prefix: str) -> Dict[str, Any]:
        if not smiles:
            return self._empty_descriptor(prefix)
        if smiles in self._cache:
            cached = self._cache[smiles]
            return {f"{prefix}_{k}": v for k, v in cached.items()}

        bundle = build_molecule_bundle(smiles, options={"include_rdkit": True})
        motif_ids = collect_motif_ids(bundle)

        aryl_steric: List[float] = []
        alkyl_steric: List[float] = []
        for entry in (bundle.get("steric", {}) or {}).get("aryl", []) or []:
            aryl_steric.extend(extract_score_values(entry.get("result")))
        for entry in (bundle.get("steric", {}) or {}).get("alkyl", []) or []:
            alkyl_steric.extend(extract_score_values(entry.get("result")))

        electronic: List[float] = []
        for entry in (bundle.get("electronics", {}) or {}).get("aryl", []) or []:
            electronic.extend(extract_score_values(entry.get("result")))

        descriptor = {
            "motif_tokens": " ".join(motif_ids),
            "fg_tokens": " ".join(sorted(get_functional_groups(smiles))),
            "motif_count": float(len(motif_ids)),
            "aryl_steric_max": safe_max(aryl_steric, default=0.0),
            "alkyl_steric_max": safe_max(alkyl_steric, default=0.0),
            "aryl_electronic_avg": safe_mean(electronic, default=5.0),
        }
        self._cache[smiles] = descriptor
        return {f"{prefix}_{k}": v for k, v in descriptor.items()}

    @staticmethod
    def _empty_descriptor(prefix: str) -> Dict[str, Any]:
        out = {
            f"{prefix}_motif_tokens": "",
            f"{prefix}_fg_tokens": "",
            f"{prefix}_motif_count": 0.0,
            f"{prefix}_aryl_steric_max": 0.0,
            f"{prefix}_alkyl_steric_max": 0.0,
            f"{prefix}_aryl_electronic_avg": 5.0,
        }
        return out

    def build_row_descriptors(
        self,
        sulfonamide_smiles: Optional[str],
        boronic_smiles: Optional[str],
    ) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        out.update(self._describe_smiles(sulfonamide_smiles, "sulf"))
        out.update(self._describe_smiles(boronic_smiles, "bor"))
        return out


def load_chanlam_dataset(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if "yield" not in df.columns:
        raise ValueError("Input CSV must contain a 'yield' column.")
    if "reaction_smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'reaction_smiles' column.")
    parsed = df["reaction_smiles"].astype(str).map(extract_reaction_substrates)
    df["sulfonamide_smiles"] = [s for s, _ in parsed]
    df["boronic_smiles"] = [b for _, b in parsed]
    return df


def build_feature_table(df: pd.DataFrame) -> pd.DataFrame:
    builder = DescriptorBuilder()
    records: List[Dict[str, Any]] = []
    for _, row in df.iterrows():
        rec: Dict[str, Any] = {
            "yield": float(row["yield"]),
            "sulfonamide_smiles": row.get("sulfonamide_smiles"),
            "boronic_smiles": row.get("boronic_smiles"),
        }
        for col in CONDITION_COLUMNS:
            rec[col] = "NA" if pd.isna(row.get(col)) else str(row.get(col))
        rec["formed_motifs_tokens"] = tokenize_text_field(row.get("formed_motifs"))
        rec["spectator_groups_tokens"] = tokenize_text_field(row.get("spectator_groups"))
        rec.update(builder.build_row_descriptors(rec["sulfonamide_smiles"], rec["boronic_smiles"]))
        records.append(rec)
    feat = pd.DataFrame.from_records(records)
    for col in NUMERIC_COLUMNS:
        feat[col] = pd.to_numeric(feat[col], errors="coerce").fillna(0.0)
    for col in TEXT_COLUMNS:
        feat[col] = feat[col].fillna("").astype(str)
    return feat
