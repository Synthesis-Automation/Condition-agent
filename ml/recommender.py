"""Two-stage recommender (prior shortlist + descriptor reranker)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

import numpy as np
import pandas as pd

from chemtools.featurizers.unified import featurize_reaction
from ml.chemistry import CONDITION_COLUMNS, DescriptorBuilder, extract_reaction_substrates, tokenize_text_field
from ml.condition_features import add_condition_property_features
from ml.contracts import RecommenderConfig
from ml.features import (
    FeatureSpec,
    build_condition_library,
    ensure_feature_columns,
    profile_condition_columns,
    resolve_feature_spec,
)
from ml.models import ConditionPriorRanker, Stage2Model


@dataclass
class TwoStageRecommender:
    """Recommender with stage-1 prior and stage-2 ML reranking."""

    config: RecommenderConfig
    feature_spec: FeatureSpec
    condition_cols: list[str]
    condition_library: pd.DataFrame
    stage1: ConditionPriorRanker
    stage2: Stage2Model

    @classmethod
    def create(cls, config: RecommenderConfig) -> "TwoStageRecommender":
        spec = resolve_feature_spec(
            profile=config.feature.profile,
            with_condition_props=config.feature.with_condition_props,
        )
        condition_cols = profile_condition_columns(config.feature.profile)
        stage1 = ConditionPriorRanker(condition_cols=condition_cols, shrinkage_m=config.stage1.shrinkage_m)
        stage2 = Stage2Model(
            model_type=config.stage2.model_type,
            feature_spec=spec,
            n_estimators=config.stage2.n_estimators,
            random_state=config.stage2.random_state,
        )
        return cls(
            config=config,
            feature_spec=spec,
            condition_cols=condition_cols,
            condition_library=pd.DataFrame(columns=condition_cols),
            stage1=stage1,
            stage2=stage2,
        )

    def fit(self, train_df: pd.DataFrame) -> "TwoStageRecommender":
        df = train_df.copy()
        if self.config.feature.with_condition_props:
            for col in CONDITION_COLUMNS:
                if col not in df.columns:
                    df[col] = "NA"
            df = add_condition_property_features(df)
        df = ensure_feature_columns(df, self.feature_spec)
        self.condition_library = build_condition_library(df, condition_cols=self.condition_cols)
        y = pd.to_numeric(df[self.config.dataset.yield_col], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        groups = df[self.config.dataset.group_col].fillna("NA").astype(str).to_numpy()
        self.stage1.fit(df, y_col=self.config.dataset.yield_col)
        self.stage2.fit(df, y, groups)
        return self

    def _blend_scores(self, stage1_scores: np.ndarray, stage2_scores: np.ndarray) -> np.ndarray:
        w1, w2 = self.config.blend.normalized()
        return w1 * np.asarray(stage1_scores, dtype=float) + w2 * np.asarray(stage2_scores, dtype=float)

    def score_rows(self, df: pd.DataFrame) -> pd.DataFrame:
        work = df.copy()
        if self.config.feature.with_condition_props:
            for col in CONDITION_COLUMNS:
                if col not in work.columns:
                    work[col] = "NA"
            work = add_condition_property_features(work)
        work = ensure_feature_columns(work, self.feature_spec)
        s1 = self.stage1.score(work)
        s2 = self.stage2.predict(work)
        out = work.copy()
        out["stage1_score"] = s1
        out["stage2_score"] = s2
        out["final_score"] = self._blend_scores(s1, s2)
        return out

    def shortlist_rows(self, df: pd.DataFrame) -> pd.DataFrame:
        scored = self.score_rows(df)
        k = max(1, int(self.config.stage1.shortlist_k))
        # shortlist inside each reaction to keep reaction-level ranking stable
        rxn_col = self.config.dataset.reaction_col
        if rxn_col in scored.columns:
            return (
                scored.sort_values("stage1_score", ascending=False)
                .groupby(rxn_col, as_index=False, group_keys=False)
                .head(k)
                .reset_index(drop=True)
            )
        return scored.sort_values("stage1_score", ascending=False).head(k).reset_index(drop=True)

    @staticmethod
    def _reaction_text_tokens(reaction_smiles: str) -> Dict[str, str]:
        payload = featurize_reaction(reaction_smiles, options={"detailed": True, "motif_site_filter": "substituent"})
        aggregates = (payload or {}).get("aggregates") or {}
        formed = aggregates.get("formed_motifs_center") or aggregates.get("formed_motifs") or []
        spectators = aggregates.get("spectator_groups_ranked") or aggregates.get("spectator_groups_combined") or []
        return {
            "formed_motifs_tokens": tokenize_text_field(" / ".join(str(v) for v in formed)),
            "spectator_groups_tokens": tokenize_text_field(" / ".join(str(v) for v in spectators)),
        }

    def _candidate_rows_for_reaction(self, reaction_smiles: str) -> pd.DataFrame:
        sulfonamide_smiles, boronic_smiles = extract_reaction_substrates(reaction_smiles)
        if not sulfonamide_smiles or not boronic_smiles:
            raise ValueError("Could not parse sulfonamide/boronic reactants from reaction_smiles.")

        tokens = self._reaction_text_tokens(reaction_smiles)
        desc = DescriptorBuilder().build_row_descriptors(sulfonamide_smiles, boronic_smiles)
        records: list[Dict[str, Any]] = []
        for _, cond in self.condition_library.iterrows():
            row: Dict[str, Any] = {
                "reaction_smiles": reaction_smiles,
                "sulfonamide_smiles": sulfonamide_smiles,
                "boronic_smiles": boronic_smiles,
                "yield": 0.0,
                "catalyst": str(cond["catalyst"]) if "catalyst" in cond else "NA",
                "base": str(cond["base"]) if "base" in cond else "NA",
                "solvent": str(cond["solvent"]) if "solvent" in cond else "NA",
                "formed_motifs_tokens": tokens["formed_motifs_tokens"],
                "spectator_groups_tokens": tokens["spectator_groups_tokens"],
                "sulf_motif_tokens": str(desc.get("sulf_motif_tokens", "")),
                "bor_motif_tokens": str(desc.get("bor_motif_tokens", "")),
                "sulf_motif_count": float(desc.get("sulf_motif_count", 0.0)),
                "sulf_aryl_steric_max": float(desc.get("sulf_aryl_steric_max", 0.0)),
                "sulf_alkyl_steric_max": float(desc.get("sulf_alkyl_steric_max", 0.0)),
                "sulf_aryl_electronic_avg": float(desc.get("sulf_aryl_electronic_avg", 5.0)),
                "bor_motif_count": float(desc.get("bor_motif_count", 0.0)),
                "bor_aryl_steric_max": float(desc.get("bor_aryl_steric_max", 0.0)),
                "bor_alkyl_steric_max": float(desc.get("bor_alkyl_steric_max", 0.0)),
                "bor_aryl_electronic_avg": float(desc.get("bor_aryl_electronic_avg", 5.0)),
            }
            records.append(row)
        return pd.DataFrame.from_records(records)

    def recommend(self, reaction_smiles: str, *, top_k: int) -> pd.DataFrame:
        candidates = self._candidate_rows_for_reaction(reaction_smiles)
        scored = self.score_rows(candidates)
        if self.config.feature.profile == "base_motif_spectator":
            scored = scored.drop(columns=["catalyst", "solvent"], errors="ignore")
        return scored.sort_values("final_score", ascending=False).head(max(1, int(top_k))).reset_index(drop=True)
