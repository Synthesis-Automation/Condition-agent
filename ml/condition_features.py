"""Engineered condition-property features for condition recommendation."""

from __future__ import annotations

from functools import lru_cache
import re
from typing import Any

import pandas as pd

from chemtools.reagent.reagent_v2 import ReagentTaxonomyV2
from ml.chemistry import CONDITION_COLUMNS


ROLE_HINT_BY_COLUMN = {
    "catalyst": "metal_catalyst",
    "base": "base",
    "solvent": "solvent",
}

CONDITION_PROP_CAT_COLS = [
    "catalyst_metal",
    "catalyst_family",
    "base_family",
    "base_strength_band",
    "base_inorganic_band",
    "solvent_family",
    "solvent_class",
    "solvent_proticity",
    "solvent_polarity_band",
]

CONDITION_PROP_NUM_COLS = [
    "is_cu_catalyst",
    "is_strong_base",
    "is_polar_solvent",
]

METAL_PATTERNS = {
    "Cu": re.compile(r"\bcu\b|copper", re.I),
    "Pd": re.compile(r"\bpd\b|palladium", re.I),
    "Ni": re.compile(r"\bni\b|nickel", re.I),
    "Co": re.compile(r"\bco\b|cobalt", re.I),
    "Fe": re.compile(r"\bfe\b|iron", re.I),
    "Rh": re.compile(r"\brh\b|rhodium", re.I),
    "Ir": re.compile(r"\bir\b|iridium", re.I),
    "Ru": re.compile(r"\bru\b|ruthenium", re.I),
}


@lru_cache(maxsize=1)
def taxonomy() -> ReagentTaxonomyV2:
    return ReagentTaxonomyV2.from_path()


def normalize_name(value: Any) -> str:
    if value is None:
        return "NA"
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return "NA"
    return text


def infer_catalyst_metal(catalyst: str) -> str:
    for metal, pattern in METAL_PATTERNS.items():
        if pattern.search(catalyst):
            return metal
    return "other_or_unknown"


def infer_base_strength(base: str) -> str:
    text = base.lower()
    if text == "na":
        return "unknown"
    strong_tokens = ("otbu", "ome", "oet", "oh", "hydride", "amide", "dbu", "dbn", "lithium")
    medium_tokens = ("co3", "po4", "f", "carbonate", "phosphate", "acetate")
    weak_tokens = ("et3n", "dipea", "pyridine", "dmap", "amine", "lutidine")
    if any(tok in text for tok in strong_tokens):
        return "strong"
    if any(tok in text for tok in medium_tokens):
        return "medium"
    if any(tok in text for tok in weak_tokens):
        return "weak"
    return "unknown"


def infer_base_inorganic(base: str) -> str:
    text = base.lower()
    if text == "na":
        return "unknown"
    tokens = ("k", "na", "cs", "li", "carbonate", "phosphate", "hydroxide", "fluoride")
    return "inorganic" if any(tok in text for tok in tokens) else "organic_or_unknown"


def infer_solvent_props(solvent: str) -> tuple[str, str, str]:
    text = solvent.lower()
    if text == "na":
        return "unknown", "unknown", "unknown"
    if text in {"dce", "dcm", "ch2cl2"}:
        return "chlorinated_aprotic", "aprotic", "medium"
    if text in {"mecn", "dmf", "dma", "dmso", "nmp"}:
        return "polar_aprotic", "aprotic", "high"
    if text in {"meoh", "etoh", "ipa", "h2o", "water"}:
        return "protic", "protic", "high"
    if text in {"etOAc".lower(), "thf", "dioxane", "mtbe"}:
        return "medium_aprotic", "aprotic", "medium"
    return "other", "unknown", "unknown"


def classify_family(value: str, role_hint: str) -> str:
    if value == "NA":
        return "NA"
    try:
        result = taxonomy().classify({"name": value, "cas": None, "smiles": None})
    except Exception:
        return "unknown"
    if result is None:
        return "unknown"
    if result.role_id != role_hint:
        return "role_mismatch"
    return result.family_id


def add_condition_property_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in CONDITION_COLUMNS:
        if col not in out.columns:
            out[col] = "NA"
        out[col] = out[col].map(normalize_name)

    catalyst_norm = out["catalyst"].astype(str)
    base_norm = out["base"].astype(str)
    solvent_norm = out["solvent"].astype(str)

    out["catalyst_metal"] = catalyst_norm.map(infer_catalyst_metal)
    out["base_strength_band"] = base_norm.map(infer_base_strength)
    out["base_inorganic_band"] = base_norm.map(infer_base_inorganic)

    solv = solvent_norm.map(infer_solvent_props)
    out["solvent_class"] = [x[0] for x in solv]
    out["solvent_proticity"] = [x[1] for x in solv]
    out["solvent_polarity_band"] = [x[2] for x in solv]

    out["catalyst_family"] = catalyst_norm.map(lambda x: classify_family(x, ROLE_HINT_BY_COLUMN["catalyst"]))
    out["base_family"] = base_norm.map(lambda x: classify_family(x, ROLE_HINT_BY_COLUMN["base"]))
    out["solvent_family"] = solvent_norm.map(lambda x: classify_family(x, ROLE_HINT_BY_COLUMN["solvent"]))

    out["is_cu_catalyst"] = (out["catalyst_metal"] == "Cu").astype(float)
    out["is_strong_base"] = (out["base_strength_band"] == "strong").astype(float)
    out["is_polar_solvent"] = out["solvent_class"].isin(["polar_aprotic", "protic"]).astype(float)

    for col in CONDITION_PROP_NUM_COLS:
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
    for col in CONDITION_PROP_CAT_COLS:
        out[col] = out[col].fillna("unknown").astype(str)
    return out
