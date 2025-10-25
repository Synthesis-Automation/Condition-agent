"""
Lightweight condition recommender backed by the curated HTE z-score dataset.

This module exposes the ``HTESimpleConditionRecommender`` class which mirrors
the behaviour of the prototype found under ``data-processor/HTE_data`` but
ships as part of the main ChemTools package.  It loads the standardized
``HTE_0.csv`` file (stored in ``data/HTE_db``) and applies a three-tiered
matching strategy:

1. Exact substrate match (reaction type + reactant A + reactant B)
2. Category match (reaction type + reactant categories)
3. Reaction-type-only match

The output is a structured dictionary summarising the top-ranked condition
combinations together with supporting statistics.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


DEFAULT_DATASET = Path(__file__).resolve().parents[2] / "data" / "HTE_db" / "HTE_0.csv"


@dataclass
class RecommendationOptions:
    """Configuration for the HTE recommender."""

    top_n: int = 3
    min_precedents: int = 3
    zscore_threshold: float = 1.0


class HTESimpleConditionRecommender:
    """
    Simple yet effective recommender using the standardized HTE dataset.

    Parameters
    ----------
    csv_path:
        Path to the standardized HTE dataset. Defaults to ``data/HTE_db/HTE_0.csv``.
    zscore_threshold:
        Minimum z-score for cases to be considered high-performing. Defaults to ``1.0``.
    """

    def __init__(self, csv_path: Path | str = DEFAULT_DATASET, zscore_threshold: float = 1.0) -> None:
        csv_path = Path(csv_path)
        if not csv_path.exists():
            raise FileNotFoundError(f"HTE dataset not found at: {csv_path}")

        self.dataset_path = csv_path
        self.zscore_threshold = zscore_threshold
        self.df = pd.read_csv(csv_path)
        self.high_performers = self.df[self.df["z-Score"] > zscore_threshold]

    # ------------------------------------------------------------------ public API

    def recommend(
        self,
        reaction_type: str,
        reactant_a: str,
        reactant_b: str,
        *,
        options: RecommendationOptions | None = None,
    ) -> Dict:
        """
        Produce recommendations for the supplied standardized reactants.

        Parameters
        ----------
        reaction_type:
            Standardized reaction type (e.g., ``"Buchwald-Hartwig-C-N"``).
        reactant_a:
            Standardized electrophile identifier (e.g., ``"ArBr"``).
        reactant_b:
            Standardized nucleophile identifier (e.g., ``"RNH2"``).
        options:
            Optional :class:`RecommendationOptions` overrides.

        Returns
        -------
        dict
            Structured recommendation payload with metadata and ranked condition sets.
        """
        opts = options or RecommendationOptions()

        exact_matches = self._exact_match(reaction_type, reactant_a, reactant_b)
        if len(exact_matches) >= opts.min_precedents:
            return self._aggregate_conditions(
                exact_matches,
                match_level="exact",
                reaction_type=reaction_type,
                reactant_a=reactant_a,
                reactant_b=reactant_b,
                top_n=opts.top_n,
            )

        reactant_a_category = self._get_category(reactant_a)
        reactant_b_category = self._get_category(reactant_b)

        if reactant_a_category and reactant_b_category:
            category_matches = self._category_match(reaction_type, reactant_a_category, reactant_b_category)
            if len(category_matches) >= opts.min_precedents:
                return self._aggregate_conditions(
                    category_matches,
                    match_level="category",
                    reaction_type=reaction_type,
                    reactant_a=reactant_a,
                    reactant_b=reactant_b,
                    top_n=opts.top_n,
                )

        rxn_matches = self._reaction_type_match(reaction_type)
        if len(rxn_matches) > 0:
            return self._aggregate_conditions(
                rxn_matches,
                match_level="reaction_type",
                reaction_type=reaction_type,
                reactant_a=reactant_a,
                reactant_b=reactant_b,
                top_n=opts.top_n,
            )

        return {
            "match_level": "none",
            "message": f"No precedents found for reaction type: {reaction_type}",
            "recommendations": [],
        }

    # ------------------------------------------------------------------ helpers

    def _exact_match(self, reaction_type: str, reactant_a: str, reactant_b: str) -> pd.DataFrame:
        return self.high_performers[
            (self.high_performers["Reaction_Type_Standardized"] == reaction_type)
            & (self.high_performers["Reactant_A"] == reactant_a)
            & (self.high_performers["Reactant_B"] == reactant_b)
        ]

    def _category_match(self, reaction_type: str, reactant_a_category: str, reactant_b_category: str) -> pd.DataFrame:
        return self.high_performers[
            (self.high_performers["Reaction_Type_Standardized"] == reaction_type)
            & (self.high_performers["Reactant_A_Category"] == reactant_a_category)
            & (self.high_performers["Reactant_B_Category"] == reactant_b_category)
        ]

    def _reaction_type_match(self, reaction_type: str) -> pd.DataFrame:
        return self.high_performers[self.high_performers["Reaction_Type_Standardized"] == reaction_type]

    def _get_category(self, reactant_type: str) -> Optional[str]:
        match = self.df[self.df["Reactant_A"] == reactant_type]
        if len(match) > 0:
            category = match.iloc[0]["Reactant_A_Category"]
            if pd.notna(category) and category != "":
                return category

        match = self.df[self.df["Reactant_B"] == reactant_type]
        if len(match) > 0:
            category = match.iloc[0]["Reactant_B_Category"]
            if pd.notna(category) and category != "":
                return category
        return None

    def _aggregate_conditions(
        self,
        df: pd.DataFrame,
        *,
        match_level: str,
        reaction_type: str,
        reactant_a: str,
        reactant_b: str,
        top_n: int,
    ) -> Dict:
        if df.empty:
            return {
                "reaction_type": reaction_type,
                "substrate": {"reactant_a": reactant_a, "reactant_b": reactant_b},
                "match_level": match_level,
                "total_precedents": 0,
                "recommendations": [],
                "metadata": {"zscore_threshold": self.zscore_threshold, "match_explanation": self._get_match_explanation(match_level)},
            }

        grouped = (
            df.groupby(["Catalyst", "Ligand", "Base", "Solvent"])
            .agg(
                {
                    "z-Score": ["count", "mean", "std", "min", "max"],
                    "AREA_TOTAL_REDUCED": "mean",
                }
            )
            .reset_index()
        )
        grouped.columns = [
            "Catalyst",
            "Ligand",
            "Base",
            "Solvent",
            "count",
            "avg_zscore",
            "std_zscore",
            "min_zscore",
            "max_zscore",
            "avg_area",
        ]

        total_high = len(df)
        grouped["frequency_weight"] = grouped["count"] / max(total_high, 1)

        max_z = df["z-Score"].max()
        if max_z <= self.zscore_threshold:
            grouped["zscore_weight"] = 0.0
        else:
            grouped["zscore_weight"] = (grouped["avg_zscore"] - self.zscore_threshold) / (max_z - self.zscore_threshold)
            grouped["zscore_weight"] = grouped["zscore_weight"].clip(0, 1)

        grouped["consistency_weight"] = 1 / (1 + grouped["std_zscore"].fillna(0))

        grouped["score"] = (
            0.5 * grouped["frequency_weight"] + 0.3 * grouped["zscore_weight"] + 0.2 * grouped["consistency_weight"]
        )
        grouped = grouped.sort_values("score", ascending=False)
        top_conditions = grouped.head(top_n)

        recommendations = []
        for _, row in top_conditions.iterrows():
            recommendations.append(
                {
                    "rank": len(recommendations) + 1,
                    "catalyst": self._string_or_none(row["Catalyst"]),
                    "ligand": self._string_or_none(row["Ligand"]),
                    "base": self._string_or_none(row["Base"]),
                    "solvent": self._string_or_none(row["Solvent"]),
                    "confidence_score": round(float(row["score"]), 3),
                    "evidence": {
                        "successful_cases": int(row["count"]),
                        "total_precedents": total_high,
                        "success_rate": f"{row['count'] / total_high * 100:.1f}%" if total_high else "0.0%",
                        "avg_zscore": round(float(row["avg_zscore"]), 2),
                        "zscore_range": [
                            round(float(row["min_zscore"]), 2),
                            round(float(row["max_zscore"]), 2),
                        ],
                    },
                }
            )

        return {
            "reaction_type": reaction_type,
            "substrate": {"reactant_a": reactant_a, "reactant_b": reactant_b},
            "match_level": match_level,
            "total_precedents": total_high,
            "recommendations": recommendations,
            "metadata": {
                "zscore_threshold": self.zscore_threshold,
                "match_explanation": self._get_match_explanation(match_level),
                "dataset": str(self.dataset_path),
            },
        }

    # ------------------------------------------------------------------ misc

    @staticmethod
    def _string_or_none(value: object) -> Optional[str]:
        if isinstance(value, str):
            stripped = value.strip()
            return stripped or None
        if pd.isna(value):
            return None
        return str(value)

    @staticmethod
    def _get_match_explanation(match_level: str) -> str:
        explanations = {
            "exact": "Exact match found for your substrate combination. High confidence.",
            "category": "No exact match, but similar substrates found (category match). Medium confidence.",
            "reaction_type": "No substrate match. Using general conditions for this reaction type. Lower confidence - consider screening.",
        }
        return explanations.get(match_level, "Unknown match level")


__all__ = ["HTESimpleConditionRecommender", "RecommendationOptions", "DEFAULT_DATASET"]
