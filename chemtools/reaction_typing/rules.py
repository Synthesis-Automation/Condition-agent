"""
Rule-based reaction typing using motifs and handle deltas.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set

from .features import aggregate_side, delta_counts

_HALIDES = ("F", "Cl", "Br", "I")


@lru_cache(maxsize=1)
def load_motif_buckets() -> Dict[str, Set[str]]:
    """
    Load motif buckets from reaction_types.v3.3.json if available.
    """
    path = Path(__file__).resolve().parents[1] / "taxonomy" / "v2_data" / "reaction_types.v3.3.json"
    if path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
        motif_sets = payload.get("motif_sets") or {}
        return _buckets_from_sets(motif_sets)
    return _default_buckets()


def _buckets_from_sets(motif_sets: Mapping[str, Iterable[str]]) -> Dict[str, Set[str]]:
    def as_set(key: str) -> Set[str]:
        values = motif_sets.get(key) or []
        return {str(v) for v in values if isinstance(v, str)}

    buckets = {
        "sp2_electrophile": as_set("sp2_electrophiles_default"),
        "organoboron": as_set("organoboron_default"),
        "terminal_alkyne": as_set("terminal_alkynes_default"),
        "terminal_alkene": as_set("terminal_alkenes_default"),
        "amine": as_set("amines_nh_default"),
        "alcohol": as_set("alcohols_oh_default"),
        "thiol": as_set("thiols_sh_default"),
        "acid": as_set("carboxylic_acids_default"),
        "acyl_halide": as_set("acyl_halides_default"),
        "amide": as_set("amides_default"),
        "alkyl_halide": as_set("alkyl_halides_default"),
        "alkyl_sulfonate": as_set("alkyl_sulfonates_default"),
    }
    buckets["sp3_electrophile"] = buckets["alkyl_halide"] | buckets["alkyl_sulfonate"]
    return buckets


def _default_buckets() -> Dict[str, Set[str]]:
    buckets = {
        "sp2_electrophile": {
            "Ar-Br", "Ar-Cl", "Ar-I", "Ar-F", "Ar-OMs", "Ar-OTf", "Ar-OTs",
            "Arom-Br", "Arom-Cl", "Arom-I", "Arom-F", "Arom-OMs", "Arom-OTf", "Arom-OTs",
            "Vinyl-Br", "Vinyl-Cl", "Vinyl-I", "Vinyl-F", "Vinyl-OMs", "Vinyl-OTf", "Vinyl-OTs",
        },
        "organoboron": {
            "Ar-B(OH)2", "Ar-Bpin", "Arom-B(OH)2", "Arom-Bpin", "Vinyl-B(OH)2", "Vinyl-Bpin",
        },
        "terminal_alkyne": {"Ar-CCH", "Arom-CCH", "R-CCH"},
        "terminal_alkene": {"Ar-CHCH2", "Arom-CHCH2", "R-CHCH2"},
        "amine": {"Any-NH2", "Any-NHR", "Ar-NH2", "Ar-NHR", "Arom-NH2", "Arom-NHR", "Vinyl-NH2"},
        "alcohol": {"Any-OH", "Ar-OH", "Arom-OH", "Bn-OH", "Allyl-OH", "Vinyl-OH"},
        "thiol": {"Any-SH", "Ar-SH", "Arom-SH"},
        "acid": {"Any-CO2H", "Ar-CO2H", "Arom-CO2H", "Vinyl-CO2H"},
        "acyl_halide": {"Acyl-Br", "Acyl-Cl", "Any-COCl", "Ar-COCl", "Arom-COCl"},
        "amide": {"Any-CONHR", "Ar-CONHR", "Arom-CONHR"},
        "alkyl_halide": {"Allyl-Br", "Allyl-Cl", "Allyl-I", "Bn-Br", "Bn-Cl", "Bn-I", "R-Br", "R-Cl", "R-I"},
        "alkyl_sulfonate": {"R-OMs", "R-OTf", "R-OTs"},
    }
    buckets["sp3_electrophile"] = buckets["alkyl_halide"] | buckets["alkyl_sulfonate"]
    return buckets


def has_any(motif_counts: Mapping[str, int], bucket: Iterable[str]) -> bool:
    return any(motif_counts.get(key, 0) > 0 for key in bucket)


def _halide_delta(delta: Mapping[str, int]) -> int:
    return sum(delta.get(sym, 0) for sym in _HALIDES)


@dataclass(frozen=True)
class Rule:
    name: str
    priority: int

    def match(self, reac: Dict[str, Dict[str, int]], prod: Dict[str, Dict[str, int]], delta: Dict[str, Dict[str, int]]) -> Optional[float]:
        return None


class SuzukiRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if not has_any(reac["motifs"], buckets["organoboron"]):
            return None
        score = 0.75
        if delta["elements"].get("B", 0) < 0:
            score += 0.10
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.95)


class StilleRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if reac["elements"].get("Sn", 0) <= 0:
            return None
        score = 0.75
        if delta["elements"].get("Sn", 0) < 0:
            score += 0.10
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.95)


class NegishiRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if reac["elements"].get("Zn", 0) <= 0:
            return None
        score = 0.75
        if delta["elements"].get("Zn", 0) < 0:
            score += 0.10
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.95)


class KumadaRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if reac["elements"].get("Mg", 0) <= 0:
            return None
        score = 0.75
        if delta["elements"].get("Mg", 0) < 0:
            score += 0.10
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.95)


class SonogashiraRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if not has_any(reac["motifs"], buckets["terminal_alkyne"]):
            return None
        score = 0.75
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.9)


class HeckRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if not has_any(reac["motifs"], buckets["terminal_alkene"]):
            return None
        score = 0.75
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.9)


class AmideFormationRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["amine"]):
            return None
        if not (has_any(reac["motifs"], buckets["acid"]) or has_any(reac["motifs"], buckets["acyl_halide"])):
            return None
        if prod["fg"].get("amide", 0) <= 0:
            return None
        score = 0.75
        if delta["fg"].get("acid_or_ester", 0) < 0:
            score += 0.10
        return min(score, 0.95)


class CNRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if not has_any(reac["motifs"], buckets["amine"]):
            return None
        score = 0.7
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.9)


class CORule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if not has_any(reac["motifs"], buckets["alcohol"]):
            return None
        score = 0.7
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.9)


class CSRule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp2_electrophile"]):
            return None
        if not has_any(reac["motifs"], buckets["thiol"]):
            return None
        score = 0.7
        if _halide_delta(delta["elements"]) < 0:
            score += 0.05
        return min(score, 0.9)


class SN2Rule(Rule):
    def match(self, reac, prod, delta):
        buckets = load_motif_buckets()
        if not has_any(reac["motifs"], buckets["sp3_electrophile"]):
            return None
        if not (
            has_any(reac["motifs"], buckets["amine"])
            or has_any(reac["motifs"], buckets["alcohol"])
            or has_any(reac["motifs"], buckets["thiol"])
        ):
            return None
        score = 0.7
        if _halide_delta(delta["elements"]) < 0:
            score += 0.10
        return min(score, 0.9)


RULES: List[Rule] = sorted(
    [
        SuzukiRule(name="Suzuki", priority=10),
        StilleRule(name="Stille", priority=15),
        NegishiRule(name="Negishi", priority=20),
        KumadaRule(name="Kumada", priority=25),
        SonogashiraRule(name="Sonogashira", priority=30),
        HeckRule(name="Heck", priority=35),
        AmideFormationRule(name="Amide formation", priority=40),
        CNRule(name="C-N coupling", priority=50),
        CORule(name="C-O coupling", priority=60),
        CSRule(name="C-S coupling", priority=70),
        SN2Rule(name="SN2 substitution", priority=80),
    ],
    key=lambda r: r.priority,
)


def classify_reaction(
    reactant_mols: Iterable[Any],
    product_mols: Iterable[Any],
    query_index: Mapping[str, Iterable[Any]],
    *,
    rules: Optional[List[Rule]] = None,
) -> Dict[str, Any]:
    """
    Classify reaction based on motifs and handle deltas.
    """
    reac = aggregate_side(reactant_mols, query_index)
    prod = aggregate_side(product_mols, query_index)
    delt = delta_counts(prod, reac)

    rules = rules or RULES
    predicted = "Unknown"
    confidence = 0.0
    for rule in rules:
        conf = rule.match(reac, prod, delt)
        if conf is not None:
            predicted = rule.name
            confidence = conf
            break

    return {
        "predicted": predicted,
        "confidence": confidence,
        "evidence": {
            "reactant_hits": sorted(reac["motifs"].keys()),
            "product_hits": sorted(prod["motifs"].keys()),
            "delta": delt,
        },
    }
