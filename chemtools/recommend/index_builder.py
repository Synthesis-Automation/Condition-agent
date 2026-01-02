"""
Build a unified recommendation index from reaction datasets, protocols, and HTE data.

This builder consolidates three source types into a single index:
  - reaction_dataset JSONL files (conditions + reaction SMILES)
  - protocol_db_v2 JSON files (detailed procedures + reaction SMILES)
  - HTE CSV files (reactant-type based screens, no reaction SMILES)

Outputs:
  - index.jsonl: one line per entry with tags and source locator
  - fingerprints.npz: DRFP fingerprints for entries with reaction SMILES
  - stats.json: summary counts and build metadata
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
import csv
import json
import time

import numpy as np

from chemtools.featurizers.unified import featurize_reaction
from chemtools.recommend.utils import canonical_family

try:
    from drfp import DrfpEncoder
    _DRFP_AVAILABLE = True
except Exception:
    DrfpEncoder = None
    _DRFP_AVAILABLE = False


_INDEX_FILENAME = "index.jsonl"
_FP_FILENAME = "fingerprints.npz"
_STATS_FILENAME = "stats.json"
_PROGRESS_EVERY = 5000


@dataclass
class BuildConfig:
    reaction_dataset_dir: Path
    protocol_dir: Path
    hte_dir: Path
    output_dir: Path
    include_hte: bool = True
    max_records: Optional[int] = None
    skip_drfp: bool = False


@dataclass
class BuildStats:
    started_at: float = field(default_factory=time.time)
    total_entries: int = 0
    by_source: Dict[str, int] = field(default_factory=dict)
    drfp_computed: int = 0
    drfp_failed: int = 0
    tag_coverage: int = 0

    def record_source(self, source_type: str) -> None:
        self.by_source[source_type] = self.by_source.get(source_type, 0) + 1


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


def _hte_tags(row: Dict[str, str]) -> List[str]:
    tags: set[str] = set()
    rxn_type = (row.get("Reaction_Type_Standardized") or "").strip()
    if rxn_type:
        tags.add(f"rxn_type:{rxn_type}")

    for key in ("Reactant_A_Type", "Reactant_B_Type", "Reactant_A", "Reactant_B"):
        value = (row.get(key) or "").strip()
        if value:
            tags.add(f"reactant:{value}")

    for key in ("Reactant_A_Category", "Reactant_B_Category"):
        value = (row.get(key) or "").strip()
        if value:
            tags.add(f"reactant_cat:{value}")

    return sorted(tags)


def _hte_name(row: Dict[str, str]) -> str:
    catalyst = (row.get("Catalyst") or "").strip()
    ligand = (row.get("Ligand") or "").strip()
    base = (row.get("Base") or "").strip()
    solvent = (row.get("Solvent") or "").strip()
    parts = [part for part in (catalyst, ligand, base, solvent) if part]
    return " / ".join(parts) if parts else "HTE condition"


class UnifiedRecommendationIndexBuilder:
    def __init__(self, config: BuildConfig):
        self.config = config
        self.stats = BuildStats()
        self._drfp_encoder = DrfpEncoder() if (_DRFP_AVAILABLE and not config.skip_drfp) else None

    def build(self) -> BuildStats:
        output_dir = self.config.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        entries: List[Dict[str, Any]] = []
        drfp_ids: List[str] = []
        drfp_fps: List[np.ndarray] = []

        if self._drfp_encoder is None:
            print("DRFP: disabled (skip_drfp or drfp not installed)")
        else:
            print("DRFP: enabled")

        self._load_reaction_dataset(entries, drfp_ids, drfp_fps)
        self._load_protocols(entries, drfp_ids, drfp_fps)
        if self.config.include_hte:
            self._load_hte(entries)

        self.stats.total_entries = len(entries)
        self.stats.tag_coverage = sum(1 for entry in entries if entry.get("tags"))

        index_path = output_dir / _INDEX_FILENAME
        with index_path.open("w", encoding="utf-8") as handle:
            for entry in entries:
                handle.write(json.dumps(entry, ensure_ascii=True) + "\n")

        fp_path = output_dir / _FP_FILENAME
        if drfp_ids and drfp_fps:
            np.savez_compressed(
                fp_path,
                entry_ids=np.array(drfp_ids, dtype=object),
                fps=np.array(drfp_fps, dtype=object),
                n_bits=len(drfp_fps[0]) * 8 if drfp_fps else 2048,
                radius=3,
            )

        stats_path = output_dir / _STATS_FILENAME
        with stats_path.open("w", encoding="utf-8") as handle:
            json.dump(
                {
                    "build_time_s": round(time.time() - self.stats.started_at, 2),
                    "total_entries": self.stats.total_entries,
                    "by_source": self.stats.by_source,
                    "drfp": {
                        "computed": self.stats.drfp_computed,
                        "failed": self.stats.drfp_failed,
                        "available": bool(self._drfp_encoder),
                    },
                    "tag_coverage": self.stats.tag_coverage,
                },
                handle,
                indent=2,
            )

        return self.stats

    def _load_reaction_dataset(
        self,
        entries: List[Dict[str, Any]],
        drfp_ids: List[str],
        drfp_fps: List[np.ndarray],
    ) -> None:
        data_dir = self.config.reaction_dataset_dir
        if not data_dir.exists():
            return

        print(f"Loading reaction datasets from: {data_dir}")
        local_count = 0
        for path in sorted(data_dir.glob("*.jsonl")):
            print(f"  File: {path.name}")
            family_hint = canonical_family(path.stem)
            with path.open("r", encoding="utf-8") as handle:
                for line_idx, line in enumerate(handle):
                    if self.config.max_records and self.stats.total_entries >= self.config.max_records:
                        print("  Reached max_records limit, stopping.")
                        return
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        record = json.loads(line)
                    except json.JSONDecodeError:
                        continue

                    reaction_smiles = record.get("reaction_smiles")
                    entry_id = record.get("reaction_id") or f"{path.stem}:{line_idx}"
                    tags = self._tags_from_reaction_smiles(reaction_smiles, family_hint)

                    entry = {
                        "id": str(entry_id),
                        "name": str(entry_id),
                        "source_type": "dataset",
                        "family": family_hint,
                        "tags": tags,
                        "source_file": str(path),
                        "record_index": line_idx,
                    }

                    if reaction_smiles:
                        self._maybe_add_drfp(entry_id, reaction_smiles, drfp_ids, drfp_fps)

                    entries.append(entry)
                    self.stats.record_source("dataset")
                    self.stats.total_entries += 1
                    local_count += 1
                    if local_count % _PROGRESS_EVERY == 0:
                        print(f"  ...loaded {local_count} dataset entries")

    def _load_protocols(
        self,
        entries: List[Dict[str, Any]],
        drfp_ids: List[str],
        drfp_fps: List[np.ndarray],
    ) -> None:
        protocol_dir = self.config.protocol_dir
        if not protocol_dir.exists():
            return

        print(f"Loading protocols from: {protocol_dir}")
        local_count = 0
        for path in sorted(protocol_dir.glob("*.json")):
            if path.name.startswith("."):
                continue
            if path.name.endswith("_index.json") or path.name.endswith(".protocol_index.json"):
                continue

            print(f"  File: {path.name}")
            with path.open("r", encoding="utf-8") as handle:
                try:
                    data = json.load(handle)
                except json.JSONDecodeError:
                    continue

            records = data if isinstance(data, list) else [data]
            for idx, record in enumerate(records):
                if self.config.max_records and self.stats.total_entries >= self.config.max_records:
                    return

                meta = record.get("metadata") or {}
                reaction = record.get("reaction") or {}
                entry_id = meta.get("id") or f"{path.stem}:{idx}"
                name = meta.get("name") or entry_id
                family = canonical_family(reaction.get("family") or path.stem)
                reaction_smiles = reaction.get("reaction_smiles")

                tags = self._tags_from_reaction_smiles(reaction_smiles, family)
                entry = {
                    "id": str(entry_id),
                    "name": str(name),
                    "source_type": "protocol",
                    "family": family,
                    "tags": tags,
                    "source_file": str(path),
                    "record_index": idx,
                }

                if reaction_smiles:
                    self._maybe_add_drfp(entry_id, reaction_smiles, drfp_ids, drfp_fps)

                entries.append(entry)
                self.stats.record_source("protocol")
                self.stats.total_entries += 1
                local_count += 1
                if local_count % _PROGRESS_EVERY == 0:
                    print(f"  ...loaded {local_count} protocol entries")

    def _load_hte(self, entries: List[Dict[str, Any]]) -> None:
        hte_dir = self.config.hte_dir
        if not hte_dir.exists():
            return

        print(f"Loading HTE data from: {hte_dir}")
        local_count = 0
        for path in sorted(hte_dir.glob("*.csv")):
            print(f"  File: {path.name}")
            with path.open("r", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                for row_idx, row in enumerate(reader):
                    if self.config.max_records and self.stats.total_entries >= self.config.max_records:
                        print("  Reached max_records limit, stopping.")
                        return
                    tags = _hte_tags(row)
                    family = canonical_family((row.get("Reaction_Type_Standardized") or "").strip())
                    entry_id = f"{path.stem}:{row_idx}"
                    entry = {
                        "id": entry_id,
                        "name": _hte_name(row),
                        "source_type": "hte",
                        "family": family or "Unknown",
                        "tags": tags,
                        "source_file": str(path),
                        "record_index": row_idx,
                    }
                    entries.append(entry)
                    self.stats.record_source("hte")
                    self.stats.total_entries += 1
                    local_count += 1
                    if local_count % _PROGRESS_EVERY == 0:
                        print(f"  ...loaded {local_count} HTE entries")

    def _tags_from_reaction_smiles(self, reaction_smiles: Optional[str], family_hint: str) -> List[str]:
        if not reaction_smiles:
            return [f"rxn_type:{family_hint}"] if family_hint else []
        try:
            bundle = featurize_reaction(reaction_smiles)
            tags = _reaction_tags(bundle)
        except Exception:
            tags = []
        if family_hint:
            tags.append(f"rxn_type:{family_hint}")
        return sorted(set(tags))

    def _maybe_add_drfp(
        self,
        entry_id: str,
        reaction_smiles: str,
        drfp_ids: List[str],
        drfp_fps: List[np.ndarray],
    ) -> None:
        if self._drfp_encoder is None:
            return
        try:
            fp = self._drfp_encoder.encode([reaction_smiles])[0]
            drfp_ids.append(str(entry_id))
            drfp_fps.append(np.array(fp, dtype="uint8"))
            self.stats.drfp_computed += 1
        except Exception:
            self.stats.drfp_failed += 1


def build_unified_recommendation_index(
    *,
    reaction_dataset_dir: Path | str = "data/reaction_dataset",
    protocol_dir: Path | str = "data/protocol_db_v2",
    hte_dir: Path | str = "data/HTE_db",
    output_dir: Path | str = "build/unified_recommendation_index",
    include_hte: bool = True,
    max_records: Optional[int] = None,
    skip_drfp: bool = False,
) -> BuildStats:
    config = BuildConfig(
        reaction_dataset_dir=Path(reaction_dataset_dir),
        protocol_dir=Path(protocol_dir),
        hte_dir=Path(hte_dir),
        output_dir=Path(output_dir),
        include_hte=include_hte,
        max_records=max_records,
        skip_drfp=skip_drfp,
    )
    builder = UnifiedRecommendationIndexBuilder(config)
    return builder.build()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build unified recommendation index")
    parser.add_argument("--reaction-dataset", default="data/reaction_dataset")
    parser.add_argument("--protocols", default="data/protocol_db_v2")
    parser.add_argument("--hte", default="data/HTE_db")
    parser.add_argument("--output", default="build/unified_recommendation_index")
    parser.add_argument("--no-hte", action="store_true")
    parser.add_argument("--max-records", type=int, default=None)
    parser.add_argument("--skip-drfp", action="store_true")
    args = parser.parse_args()

    stats = build_unified_recommendation_index(
        reaction_dataset_dir=args.reaction_dataset,
        protocol_dir=args.protocols,
        hte_dir=args.hte,
        output_dir=args.output,
        include_hte=not args.no_hte,
        max_records=args.max_records,
        skip_drfp=args.skip_drfp,
    )

    print("Build complete.")
    print(f"Total entries: {stats.total_entries}")
    print(f"Sources: {stats.by_source}")
