"""Contracts and configuration objects for the rebuilt ML system."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


ModelType = Literal["rf_reg", "lgbm_rank"]
FeatureProfile = Literal["core_full", "base_motif_spectator"]


@dataclass(frozen=True)
class DatasetConfig:
    input_csv: Path = Path("data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv")
    group_col: str = "sulfonamide_smiles"
    reaction_col: str = "reaction_smiles"
    yield_col: str = "yield"


@dataclass(frozen=True)
class FeatureConfig:
    profile: FeatureProfile = "core_full"
    with_condition_props: bool = True


@dataclass(frozen=True)
class Stage1Config:
    shortlist_k: int = 30
    shrinkage_m: float = 20.0


@dataclass(frozen=True)
class Stage2Config:
    model_type: ModelType = "rf_reg"
    n_estimators: int = 220
    random_state: int = 42


@dataclass(frozen=True)
class BlendConfig:
    w_stage1: float = 0.8
    w_stage2: float = 0.2

    def normalized(self) -> tuple[float, float]:
        total = float(self.w_stage1 + self.w_stage2)
        if total <= 0.0:
            return 0.5, 0.5
        return float(self.w_stage1 / total), float(self.w_stage2 / total)


@dataclass(frozen=True)
class RecommenderConfig:
    dataset: DatasetConfig = field(default_factory=DatasetConfig)
    feature: FeatureConfig = field(default_factory=FeatureConfig)
    stage1: Stage1Config = field(default_factory=Stage1Config)
    stage2: Stage2Config = field(default_factory=Stage2Config)
    blend: BlendConfig = field(default_factory=BlendConfig)


@dataclass(frozen=True)
class EvalConfig:
    output_json: Path = Path("results/ml/chanlam_rebuild/evaluation_report.json")
    top_ks: tuple[int, ...] = (1, 3, 5)
    thresholds: tuple[float, ...] = (50.0, 70.0)
    max_folds: int = 0
