from __future__ import annotations

import json
import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple

from chemtools.reaction import get_crk_options
from chemtools.reaction import featurize_reaction

_TRANSITION_OR_HEAVY_METALS: Set[str] = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Ga", "In", "Sn", "Sb", "Bi",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
}
_METAL_TOKEN_RE = re.compile(r"\[([A-Z][a-z]?)(?:[^\]]*)\]")
_SIMPLE_COUNTERION_COMPONENTS = {
    "Cl", "Br", "I", "F",
    "[Cl-]", "[Br-]", "[I-]", "[F-]",
    "[Na+]", "[K+]", "[Li+]", "[Cs+]", "[Rb+]",
    "[NH4+]",
}
_COMPLEX_COUNTERION_RE = (
    re.compile(r"\[B\+3\]\(\[F-\]\)\(\[F-\]\)\[F-\]"),
    re.compile(r"\[P\+5\]\(\[F-\]\)\(\[F-\]\)\(\[F-\]\)\(\[F-\]\)\[F-\]"),
)


def normalize_llm_assist_options(
    llm_assist_options: Optional[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    if not isinstance(llm_assist_options, Mapping):
        return None
    payload = dict(llm_assist_options)
    if not payload.get("enabled", True):
        return None
    payload["enabled"] = True
    return payload


def build_reaction_featurization_options(
    *,
    base_options: Optional[Mapping[str, Any]] = None,
    llm_assist_options: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    options = dict(get_crk_options())
    options.update(
        {
            "include_reaction_type": True,
            "detailed": True,
            "motif_site_filter": "substituent",
        }
    )

    if isinstance(base_options, Mapping):
        for key, value in base_options.items():
            if key == "llm_assist":
                continue
            options[key] = value

    options["include_reaction_type"] = True
    options["detailed"] = True
    options["motif_site_filter"] = "substituent"

    llm_payload = normalize_llm_assist_options(
        llm_assist_options
        if llm_assist_options is not None
        else (
            base_options.get("llm_assist")
            if isinstance(base_options, Mapping)
            else None
        )
    )
    if llm_payload:
        options["llm_assist"] = llm_payload
    else:
        options.pop("llm_assist", None)
    return options


def reaction_options_signature(
    *,
    base_options: Optional[Mapping[str, Any]] = None,
    llm_assist_options: Optional[Mapping[str, Any]] = None,
) -> str:
    options = build_reaction_featurization_options(
        base_options=base_options,
        llm_assist_options=llm_assist_options,
    )
    try:
        return json.dumps(options, sort_keys=True, separators=(",", ":"))
    except Exception:
        return ""


def _load_options_from_signature(options_signature: str) -> Dict[str, Any]:
    if not options_signature:
        return build_reaction_featurization_options()
    try:
        payload = json.loads(options_signature)
    except Exception:
        return build_reaction_featurization_options()
    if not isinstance(payload, dict):
        return build_reaction_featurization_options()
    return build_reaction_featurization_options(
        base_options=payload,
        llm_assist_options=payload.get("llm_assist"),
    )


@lru_cache(maxsize=20000)
def cached_featurize_reaction(
    reaction_smiles: str,
    options_signature: str = "",
) -> Dict[str, Any]:
    if not reaction_smiles:
        return {}
    try:
        return featurize_reaction(
            reaction_smiles,
            options=_load_options_from_signature(options_signature),
        )
    except Exception:
        return {}


def reaction_type_from_bundle(bundle: Mapping[str, Any]) -> str:
    raw_value: Any = ""
    if isinstance(bundle, Mapping):
        raw_value = bundle.get("reaction_type")
    if isinstance(raw_value, Mapping):
        value = str(
            raw_value.get("name") or raw_value.get("reaction_type") or ""
        ).strip()
    else:
        value = str(raw_value or "").strip()
    if not value or value.lower() == "unknown":
        return "Unknown"
    return value


def _split_reaction_sides(reaction_smiles: str) -> Tuple[List[str], List[str]]:
    text = str(reaction_smiles or "").strip()
    if not text or ">>" not in text:
        return [], []
    left, right = text.split(">>", 1)
    reactants = [part for part in left.split(".") if part]
    products = [part for part in right.split(".") if part]
    return reactants, products


def _contains_heavy_metal_token(text: str) -> bool:
    tokens = {m for m in _METAL_TOKEN_RE.findall(str(text or "")) if m}
    return bool(tokens & _TRANSITION_OR_HEAVY_METALS)


def _is_coordination_component(smiles: str) -> bool:
    token = str(smiles or "").strip()
    return bool(token) and ("->" in token or "<-" in token) and _contains_heavy_metal_token(token)


def _is_counterion_component(smiles: str) -> bool:
    token = str(smiles or "").strip()
    if not token:
        return False
    if token in _SIMPLE_COUNTERION_COMPONENTS:
        return True
    return any(pattern.search(token) for pattern in _COMPLEX_COUNTERION_RE)


def cleanup_reaction_smiles_for_featurization(
    reaction_smiles: str,
) -> Tuple[str, Dict[str, int]]:
    text = str(reaction_smiles or "").strip()
    stats = {
        "coordination_removed": 0,
        "counterion_removed": 0,
        "cleanup_applied": 0,
    }
    reactants, products = _split_reaction_sides(text)
    if not reactants or not products:
        return text, stats

    def _filter_side(parts: List[str]) -> Tuple[List[str], int, int]:
        kept: List[str] = []
        removed_coord = 0
        removed_counter = 0
        multi_component = len(parts) > 1
        for token in parts:
            if _is_coordination_component(token):
                removed_coord += 1
                continue
            if multi_component and _is_counterion_component(token):
                removed_counter += 1
                continue
            kept.append(token)
        return kept, removed_coord, removed_counter

    cleaned_reactants, left_coord, left_counter = _filter_side(reactants)
    cleaned_products, right_coord, right_counter = _filter_side(products)
    removed_coord = left_coord + right_coord
    removed_counter = left_counter + right_counter

    if not cleaned_reactants or not cleaned_products:
        return text, stats

    cleaned = f"{'.'.join(cleaned_reactants)}>>{'.'.join(cleaned_products)}"
    if cleaned != text:
        stats["coordination_removed"] = removed_coord
        stats["counterion_removed"] = removed_counter
        stats["cleanup_applied"] = 1
    return cleaned, stats


@dataclass(frozen=True)
class ReactionFeaturizationResult:
    input_smiles: str
    reaction_smiles: str
    cleanup_stats: Dict[str, int]
    options_signature: str
    bundle: Dict[str, Any]
    detected_reaction_type: str


def analyze_reaction_featurization(
    reaction_smiles: str,
    *,
    base_options: Optional[Mapping[str, Any]] = None,
    llm_assist_options: Optional[Mapping[str, Any]] = None,
    cleanup: bool = True,
) -> ReactionFeaturizationResult:
    cleaned_smiles = str(reaction_smiles or "").strip()
    cleanup_stats = {
        "coordination_removed": 0,
        "counterion_removed": 0,
        "cleanup_applied": 0,
    }
    if cleanup:
        cleaned_smiles, cleanup_stats = cleanup_reaction_smiles_for_featurization(
            cleaned_smiles
        )

    options_signature = reaction_options_signature(
        base_options=base_options,
        llm_assist_options=llm_assist_options,
    )
    bundle = cached_featurize_reaction(cleaned_smiles, options_signature)
    return ReactionFeaturizationResult(
        input_smiles=str(reaction_smiles or "").strip(),
        reaction_smiles=cleaned_smiles,
        cleanup_stats=cleanup_stats,
        options_signature=options_signature,
        bundle=bundle,
        detected_reaction_type=reaction_type_from_bundle(bundle),
    )

