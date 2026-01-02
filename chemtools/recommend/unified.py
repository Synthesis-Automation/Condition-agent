"""
Unified recommendation engine over reaction datasets, protocols, and HTE data.

This recommender combines two layers:
  - DRFP similarity for entries with reaction SMILES
  - Feature-tag similarity for entries with reaction/HTE tags
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
import csv
import json
import time

import numpy as np

from chemtools.featurizers.unified import featurize_reaction

try:
    from drfp import DrfpEncoder
    _DRFP_AVAILABLE = True
except Exception:
    DrfpEncoder = None
    _DRFP_AVAILABLE = False


_DEFAULT_INDEX_DIR = Path(__file__).resolve().parents[2] / "build" / "unified_recommendation_index"
_INDEX_FILENAME = "index.jsonl"
_FP_FILENAME = "fingerprints.npz"
_STATS_FILENAME = "stats.json"

_SOURCE_BONUS = {
    "protocol": 1.05,
    "dataset": 1.0,
    "hte": 0.95,
}


@dataclass
class IndexEntry:
    id: str
    name: str
    source_type: str
    family: str
    tags: List[str]
    source_file: str
    record_index: Optional[int] = None
    drfp_index: Optional[int] = None


@dataclass
class RecommendationResult:
    id: str
    name: str
    source_type: str
    family: str
    similarity: float
    rank: int
    drfp_similarity: Optional[float] = None
    feature_similarity: Optional[float] = None
    source_file: Optional[str] = None
    record_index: Optional[int] = None
    full_data: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "source_type": self.source_type,
            "family": self.family,
            "similarity": self.similarity,
            "rank": self.rank,
            "drfp_similarity": self.drfp_similarity,
            "feature_similarity": self.feature_similarity,
            "source_file": self.source_file,
            "record_index": self.record_index,
        }


def _tag_weight(tag: str) -> float:
    prefix = tag.split(":", 1)[0]
    if prefix == "rxn_type":
        return 3.0
    if prefix == "reactant":
        return 2.0
    if prefix == "role":
        return 1.5
    if prefix == "rxn_cat":
        return 1.5
    if prefix == "reactant_cat":
        return 1.2
    if prefix == "motif":
        return 1.0
    if prefix == "fg":
        return 0.7
    return 1.0


def _weighted_jaccard(query_tags: set[str], entry_tags: set[str], entry_weight: float, query_weight: float) -> float:
    if not query_tags or not entry_tags:
        return 0.0
    intersection = 0.0
    for tag in entry_tags:
        if tag in query_tags:
            intersection += _tag_weight(tag)
    union = query_weight + entry_weight - intersection
    if union <= 0:
        return 0.0
    return intersection / union


def _compute_drfp_similarity(query_fp: np.ndarray, fps: np.ndarray) -> np.ndarray:
    intersection = np.logical_and(query_fp, fps).sum(axis=1)
    count_query = query_fp.sum()
    count_index = fps.sum(axis=1)
    union = count_query + count_index - intersection
    union = np.maximum(union, 1)
    return intersection / union


def _reaction_tags(bundle: Dict[str, Any]) -> List[str]:
    tags: set[str] = set()
    reaction = bundle.get("reaction") or {}
    reaction_type = reaction.get("reaction_type") or {}
    rxn_id = reaction_type.get("reaction_type")
    if rxn_id and rxn_id != "Unknown":
        tags.add(f"rxn_type:{rxn_id}")
    category = reaction_type.get("category")
    if category:
        tags.add(f"rxn_cat:{category}")

    aggregates = reaction.get("aggregates") or {}
    for motif in aggregates.get("motif_ids") or []:
        if motif:
            tags.add(f"motif:{motif}")
    for fg in aggregates.get("functional_group_ids") or []:
        if fg:
            tags.add(f"fg:{fg}")

    roles = reaction.get("roles") or {}
    for reactant in roles.get("reactants") or []:
        category_id = reactant.get("category")
        if category_id:
            tags.add(f"reactant:{category_id}")
        role = reactant.get("role")
        if role and category_id:
            tags.add(f"role:{role}:{category_id}")

    return sorted(tags)


def _hte_tags_from_query(reaction_type: Optional[str], reactant_types: Optional[Iterable[str]]) -> List[str]:
    tags: set[str] = set()
    if reaction_type:
        tags.add(f"rxn_type:{reaction_type}")
    if reactant_types:
        for entry in reactant_types:
            if entry:
                tags.add(f"reactant:{entry}")
    return sorted(tags)


class UnifiedRecommender:
    def __init__(self, index_dir: str | Path | None = None):
        self.index_dir = Path(index_dir) if index_dir else _DEFAULT_INDEX_DIR
        self.entries: List[IndexEntry] = []
        self.entry_by_id: Dict[str, IndexEntry] = {}
        self.tag_sets: List[set[str]] = []
        self.tag_weights: List[float] = []
        self._drfp_ids: List[str] = []
        self._drfp_fps: Optional[np.ndarray] = None
        self._drfp_index: Dict[str, int] = {}
        self._drfp_encoder = DrfpEncoder() if _DRFP_AVAILABLE else None
        self._load_index()
        self._load_drfp()

    def _load_index(self) -> None:
        index_path = self.index_dir / _INDEX_FILENAME
        if not index_path.exists():
            raise FileNotFoundError(
                f"Index file not found: {index_path}. "
                "Build it with: python data-processor/build_unified_recommendation_index.py"
            )

        with index_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                data = json.loads(line)
                entry = IndexEntry(
                    id=str(data.get("id") or ""),
                    name=str(data.get("name") or data.get("id") or ""),
                    source_type=str(data.get("source_type") or "unknown"),
                    family=str(data.get("family") or "Unknown"),
                    tags=list(data.get("tags") or []),
                    source_file=str(data.get("source_file") or ""),
                    record_index=data.get("record_index"),
                )
                self.entries.append(entry)
                self.entry_by_id[entry.id] = entry
                tag_set = set(entry.tags)
                self.tag_sets.append(tag_set)
                self.tag_weights.append(sum(_tag_weight(tag) for tag in tag_set))

    def _load_drfp(self) -> None:
        fp_path = self.index_dir / _FP_FILENAME
        if not fp_path.exists():
            return
        fp_data = np.load(fp_path, allow_pickle=True)
        self._drfp_ids = [str(item) for item in fp_data.get("entry_ids", [])]
        fps = fp_data.get("fps")
        if fps is None or len(self._drfp_ids) == 0:
            return
        self._drfp_fps = np.array(fps)
        self._drfp_index = {entry_id: idx for idx, entry_id in enumerate(self._drfp_ids)}
        for entry_id, idx in self._drfp_index.items():
            entry = self.entry_by_id.get(entry_id)
            if entry:
                entry.drfp_index = idx

    def _query_tags(
        self,
        reaction_smiles: Optional[str],
        reaction_type: Optional[str],
        reactant_types: Optional[Iterable[str]],
    ) -> List[str]:
        tags: List[str] = []
        if reaction_smiles:
            try:
                bundle = featurize_reaction(reaction_smiles)
                tags = _reaction_tags(bundle)
            except Exception:
                tags = []
        if reaction_type or reactant_types:
            tags.extend(_hte_tags_from_query(reaction_type, reactant_types))
        return sorted(set(tags))

    def recommend(
        self,
        reaction_smiles: Optional[str] = None,
        *,
        reaction_type: Optional[str] = None,
        reactant_types: Optional[Iterable[str]] = None,
        top_k: int = 5,
        min_similarity: float = 0.0,
        source_types: Optional[List[str]] = None,
        include_details: bool = False,
        drfp_weight: float = 0.7,
        feature_weight: float = 0.3,
    ) -> List[RecommendationResult]:
        if source_types:
            source_types = [str(s).lower() for s in source_types]
        query_tags = self._query_tags(reaction_smiles, reaction_type, reactant_types)
        query_tag_set = set(query_tags)
        query_weight = sum(_tag_weight(tag) for tag in query_tag_set)

        drfp_sim_map: Dict[str, float] = {}
        if reaction_smiles and self._drfp_encoder and self._drfp_fps is not None:
            try:
                query_fp = self._drfp_encoder.encode([reaction_smiles])[0]
                sims = _compute_drfp_similarity(query_fp, self._drfp_fps)
                drfp_sim_map = {
                    entry_id: float(sim) for entry_id, sim in zip(self._drfp_ids, sims)
                }
            except Exception:
                drfp_sim_map = {}

        results: List[Tuple[IndexEntry, float, Optional[float], Optional[float]]] = []

        for entry, entry_tags, entry_weight in zip(self.entries, self.tag_sets, self.tag_weights):
            if source_types and entry.source_type not in source_types:
                continue

            drfp_sim = drfp_sim_map.get(entry.id)
            feature_sim = None
            if query_tag_set and entry_weight > 0:
                feature_sim = _weighted_jaccard(query_tag_set, entry_tags, entry_weight, query_weight)

            if drfp_sim is None and feature_sim is None:
                continue

            if drfp_sim is not None and feature_sim is not None and query_tag_set:
                combined = (drfp_weight * drfp_sim) + (feature_weight * feature_sim)
            elif drfp_sim is not None:
                combined = drfp_sim
            else:
                combined = feature_sim or 0.0

            combined *= _SOURCE_BONUS.get(entry.source_type, 1.0)

            if combined < min_similarity:
                continue
            results.append((entry, combined, drfp_sim, feature_sim))

        results.sort(key=lambda item: item[1], reverse=True)
        results = results[: max(1, int(top_k or 1))]

        output: List[RecommendationResult] = []
        for rank, (entry, combined, drfp_sim, feature_sim) in enumerate(results, start=1):
            result = RecommendationResult(
                id=entry.id,
                name=entry.name,
                source_type=entry.source_type,
                family=entry.family,
                similarity=float(combined),
                rank=rank,
                drfp_similarity=float(drfp_sim) if drfp_sim is not None else None,
                feature_similarity=float(feature_sim) if feature_sim is not None else None,
                source_file=entry.source_file,
                record_index=entry.record_index,
            )
            if include_details:
                result.full_data = self.get_source_details(entry)
            output.append(result)
        return output

    def get_source_details(self, entry: IndexEntry | str) -> Optional[Dict[str, Any]]:
        if isinstance(entry, str):
            entry = self.entry_by_id.get(entry)
        if entry is None:
            return None
        if entry.source_type == "dataset":
            return _load_jsonl_record(entry.source_file, entry.record_index)
        if entry.source_type == "protocol":
            return _load_protocol_record(entry.source_file, entry.record_index)
        if entry.source_type == "hte":
            return _load_csv_record(entry.source_file, entry.record_index)
        return None

    def get_statistics(self) -> Dict[str, Any]:
        stats_path = self.index_dir / _STATS_FILENAME
        if stats_path.exists():
            with stats_path.open("r", encoding="utf-8") as handle:
                return json.load(handle)
        return {
            "total_entries": len(self.entries),
            "by_source": {
                "protocol": sum(1 for e in self.entries if e.source_type == "protocol"),
                "dataset": sum(1 for e in self.entries if e.source_type == "dataset"),
                "hte": sum(1 for e in self.entries if e.source_type == "hte"),
            },
        }


def _load_jsonl_record(path: str, record_index: Optional[int]) -> Optional[Dict[str, Any]]:
    if record_index is None:
        return None
    try:
        with Path(path).open("r", encoding="utf-8") as handle:
            for idx, line in enumerate(handle):
                if idx == record_index:
                    return json.loads(line)
    except Exception:
        return None
    return None


def _load_protocol_record(path: str, record_index: Optional[int]) -> Optional[Dict[str, Any]]:
    try:
        with Path(path).open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except Exception:
        return None
    if isinstance(data, list):
        if record_index is None:
            return data[0] if data else None
        if 0 <= record_index < len(data):
            return data[record_index]
        return None
    return data


def _load_csv_record(path: str, record_index: Optional[int]) -> Optional[Dict[str, Any]]:
    if record_index is None:
        return None
    try:
        with Path(path).open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for idx, row in enumerate(reader):
                if idx == record_index:
                    return row
    except Exception:
        return None
    return None


def _extract_conditions(source_type: str, data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not data:
        return {}
    if source_type == "dataset":
        return {
            "conditions": data.get("conditions") or {},
            "reference": data.get("reference"),
        }
    if source_type == "protocol":
        return {
            "reaction": data.get("reaction"),
            "reaction_setup": data.get("reaction_setup"),
            "conditions": data.get("conditions"),
        }
    if source_type == "hte":
        return {
            "catalyst": data.get("Catalyst"),
            "ligand": data.get("Ligand"),
            "base": data.get("Base"),
            "solvent": data.get("Solvent"),
            "secondary_solvent": data.get("Secondary Solvent"),
            "additive": data.get("Additive"),
            "coupling_reagent": data.get("Coupling Reagent"),
            "z_score": data.get("z-Score"),
            "area_total_reduced": data.get("AREA_TOTAL_REDUCED"),
        }
    return {}


_DEFAULT_RECOMMENDER: Optional[UnifiedRecommender] = None


def _get_default_recommender() -> UnifiedRecommender:
    global _DEFAULT_RECOMMENDER
    if _DEFAULT_RECOMMENDER is None:
        _DEFAULT_RECOMMENDER = UnifiedRecommender()
    return _DEFAULT_RECOMMENDER


def recommend_from_reaction(
    reaction: str,
    k: int = 5,
    *,
    reaction_type: Optional[str] = None,
    reactant_types: Optional[Iterable[str]] = None,
    min_similarity: float = 0.0,
    source_types: Optional[List[str]] = None,
    include_details: bool = True,
    **_: Any,
) -> Dict[str, Any]:
    start = time.time()
    recommender = _get_default_recommender()
    results = recommender.recommend(
        reaction_smiles=reaction,
        reaction_type=reaction_type,
        reactant_types=reactant_types,
        top_k=k,
        min_similarity=min_similarity,
        source_types=source_types,
        include_details=include_details,
    )

    detection = {}
    try:
        bundle = featurize_reaction(reaction)
        reaction_type_info = (bundle.get("reaction") or {}).get("reaction_type") or {}
        detection = {
            "detected_reaction_type": reaction_type_info.get("reaction_type"),
            "confidence": reaction_type_info.get("confidence"),
            "method": "featurizer",
        }
    except Exception:
        detection = {}

    recommendations = []
    for result in results:
        conditions = _extract_conditions(result.source_type, result.full_data)
        recommendations.append(
            {
                "rank": result.rank,
                "id": result.id,
                "name": result.name,
                "source_type": result.source_type,
                "family": result.family,
                "similarity": result.similarity,
                "drfp_similarity": result.drfp_similarity,
                "feature_similarity": result.feature_similarity,
                "conditions": conditions,
            }
        )

    return {
        "meta": {
            "model": "unified_recommender",
            "status": "success",
            "processing_time_ms": round((time.time() - start) * 1000, 2),
        },
        "input": {
            "reaction_smiles": reaction,
            "requested_reaction_type": reaction_type,
        },
        "detection": detection,
        "recommended_conditions": recommendations,
        "recommendations": recommendations,
    }


def recommend_conditions_structured(
    reaction: str,
    reaction_type: Optional[str] = None,
    *,
    k: int = 5,
    limit: int = 5,
    min_similarity: float = 0.0,
    source_types: Optional[List[str]] = None,
    include_details: bool = True,
    **kwargs: Any,
) -> Dict[str, Any]:
    results = recommend_from_reaction(
        reaction,
        k=k,
        reaction_type=reaction_type,
        min_similarity=min_similarity,
        source_types=source_types,
        include_details=include_details,
        **kwargs,
    )
    results["recommended_conditions"] = results.get("recommended_conditions", [])[:limit]
    results["recommendations"] = results["recommended_conditions"]
    return results


__all__ = [
    "UnifiedRecommender",
    "RecommendationResult",
    "recommend_from_reaction",
    "recommend_conditions_structured",
]
