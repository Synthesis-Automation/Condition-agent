"""
HTE-based Condition Recommendation System

This module provides condition recommendations based on High-Throughput Experimentation (HTE) data.
Recommendations are primarily based on reactant types since no reaction SMILES is provided.

Key Features:
- Reactant type-based matching using existing chemtools detection
- Z-score based ranking (primary metric for condition success)
- Success-weighted condition selection (yield > 50% threshold)
- Statistical ranking with confidence scores
- Multi-reaction type support with automatic detection
- Condition component recommendations (catalyst, ligand, base, solvent)
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any, Iterable
from collections import defaultdict, Counter
from functools import lru_cache
import pandas as pd
from pathlib import Path
import json

from chemtools.featurizers.structural import featurize_molecule


@lru_cache(maxsize=4)
def _load_hte_database_cached(
    hte_db_path: str,
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame], Dict[str, Counter]]:
    """Load and index the HTE database once per path (cached)."""
    db_path = Path(hte_db_path)
    if not db_path.exists():
        raise FileNotFoundError(f"HTE database not found: {db_path}")

    if db_path.suffix.lower() == ".jsonl":
        df = _load_hte_jsonl(db_path)
    else:
        df = pd.read_csv(db_path)
        if "Reactant_Types_Key" not in df.columns:
            df["Reactant_Types_Key"] = df.apply(
                lambda row: _reactant_key([row.get("Reactant_A_Type"), row.get("Reactant_B_Type")]),
                axis=1,
            )
    print(f"Loaded HTE database: {len(df)} experiments")

    indexed_data: Dict[str, pd.DataFrame] = {}
    reaction_type_patterns: Dict[str, Counter] = {}

    print("Building reactant type indices using Reactant_Types_Key...")

    df["Reactant_Types_Key"] = df["Reactant_Types_Key"].fillna("")
    grouped = df.groupby(["Reactant_Types_Key"])
    for key, group_df in grouped:
        indexed_data[key] = group_df

        rxn_types = group_df["Reaction_Type_Standardized"].value_counts()
        reaction_type_patterns[key] = Counter(rxn_types.to_dict())

    print(f"Indexed {len(indexed_data)} unique reactant type combinations")
    return df, indexed_data, reaction_type_patterns


def _ensure_list(values: Any) -> List[str]:
    if values is None:
        return []
    if isinstance(values, list):
        return [str(v).strip() for v in values if str(v).strip()]
    if isinstance(values, str):
        text = values.strip()
        return [text] if text else []
    text = str(values).strip()
    return [text] if text else []


def _dedupe_list(values: Iterable[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for value in values:
        if not value:
            continue
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _format_list(values: Any) -> str:
    items = _dedupe_list(_ensure_list(values))
    return " / ".join(items)


def _reactant_key(values: Iterable[Optional[str]]) -> str:
    items = _dedupe_list([str(v).strip() for v in values if v])
    return "|".join(sorted(items))


def _load_hte_jsonl(path: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue

            reactant_types = _ensure_list(record.get("reactant_types"))
            reactant_categories = _ensure_list(record.get("reactant_categories"))
            conditions = record.get("conditions") or {}
            metrics = record.get("metrics") or {}

            row = {
                "Reaction_Type_Standardized": record.get("reaction_type") or "Unknown",
                "Reactant_A_Type": reactant_types[0] if len(reactant_types) > 0 else "",
                "Reactant_B_Type": reactant_types[1] if len(reactant_types) > 1 else "",
                "Reactant_A_Category": reactant_categories[0] if len(reactant_categories) > 0 else "",
                "Reactant_B_Category": reactant_categories[1] if len(reactant_categories) > 1 else "",
                "Reactant_Types_Key": _reactant_key(reactant_types),
                "Catalyst": _format_list(conditions.get("catalyst")),
                "Ligand": _format_list(conditions.get("ligand")),
                "Base": _format_list(conditions.get("base")),
                "Solvent": _format_list(conditions.get("solvent")),
                "Secondary Solvent": _format_list(conditions.get("secondary_solvent")),
                "Additive": _format_list(conditions.get("additive")),
                "Coupling Reagent": _format_list(conditions.get("coupling_reagent")),
                "AREA_TOTAL_REDUCED": metrics.get("area_total_reduced"),
                "z-Score": metrics.get("z_score"),
            }
            rows.append(row)

    return pd.DataFrame(rows)


@dataclass
class ConditionRecommendation:
    """Single condition recommendation with metadata
    
    Recommendations are ranked primarily by avg_z_score, which measures
    the success of a condition relative to all experiments in the database.
    """
    catalyst: str
    ligand: str
    base: str
    solvent: str
    secondary_solvent: Optional[str] = None
    additive: Optional[str] = None
    coupling_reagent: Optional[str] = None
    
    # Statistics
    success_rate: float = 0.0  # % of experiments with yield > 50
    avg_yield: float = 0.0
    median_yield: float = 0.0
    num_experiments: int = 0
    avg_z_score: float = 0.0  # Average z-score (PRIMARY ranking metric for condition success)
    confidence_score: float = 0.0  # Secondary score considering z-score and sample size
    
    # Metadata
    reaction_type: Optional[str] = None
    reactant_types: Tuple[str, str] = ("", "")
    z_score_range: Tuple[float, float] = (0.0, 0.0)


@dataclass
class HTERecommendationResult:
    """Complete recommendation result with ranked conditions"""
    reactant_a_smiles: str
    reactant_b_smiles: Optional[str]
    
    # Detected types
    reactant_a_type: Optional[str] = None
    reactant_b_type: Optional[str] = None
    reactant_a_category: Optional[str] = None
    reactant_b_category: Optional[str] = None
    
    # Predicted reaction type
    predicted_reaction_type: Optional[str] = None
    reaction_type_confidence: float = 0.0
    
    # Recommendations
    recommendations: List[ConditionRecommendation] = field(default_factory=list)
    
    # Metadata
    total_matching_experiments: int = 0
    database_coverage: float = 0.0  # % of database that matches this query
    is_fallback_match: bool = False
    matched_motifs: Optional[Tuple[str, str]] = None


class HTERecommender:
    """
    HTE-based condition recommender using reactant type matching.
    
    Architecture:
    1. Load and index HTE database by reactant motifs (V2 Taxonomy)
    2. Classify input reactant SMILES to get motifs
    3. Match against database using motif combinations (with hierarchical fallback)
    4. Rank conditions by success rate and confidence
    5. Return top-k recommendations
    """
    
    def __init__(self, hte_db_path: str = "data/HTE_db/HTE_0.jsonl"):
        """Initialize recommender with HTE database"""
        self.db_path = Path(hte_db_path)
        self.df: Optional[pd.DataFrame] = None
        self.indexed_data: Dict[Tuple[str, str], pd.DataFrame] = {}
        self.reaction_type_patterns: Dict[Tuple[str, str], Counter] = {}

        df, indexed_data, patterns = _load_hte_database_cached(str(self.db_path))
        self.df = df
        self.indexed_data = dict(indexed_data)
        self.reaction_type_patterns = dict(patterns)
    
    def _detect_reactant_types(self, smiles: str) -> Tuple[List[str], Optional[str]]:
        """
        Detect reactant motifs and categories from SMILES using V2 Taxonomy.
        
        Returns:
            (motifs, category) e.g., (["Ar-Br", "Ar-OH"], "Aryl Halide")
        """
        if not smiles:
            return [], None
            
        analysis = featurize_molecule(smiles)
        
        # Extract motif IDs (e.g., "Ar-Br")
        motifs = _dedupe_list(
            [m.get("compound_id", "") for m in analysis.get("motifs", []) if m.get("compound_id")]
        )
        
        # Use the first motif's category as a general category
        category = analysis.get("motifs", [{}])[0].get("category", "Unknown") if analysis.get("motifs") else "Unknown"
        
        return motifs, category
    
    def _filter_by_catalyst(self, df: pd.DataFrame, catalyst_filter: str) -> pd.DataFrame:
        """
        Filter dataframe by catalyst metal type.
        
        Args:
            df: DataFrame to filter
            catalyst_filter: Metal type or symbol (e.g., 'Pd', 'Cu', 'Ni', 'palladium', 'copper')
        
        Returns:
            Filtered DataFrame
        """
        # Normalize filter term
        filter_lower = catalyst_filter.lower()
        
        # Map common names to symbols
        metal_map = {
            'palladium': 'Pd',
            'copper': 'Cu',
            'nickel': 'Ni',
            'iridium': 'Ir',
            'rhodium': 'Rh',
            'ruthenium': 'Ru',
            'platinum': 'Pt',
            'gold': 'Au',
            'silver': 'Ag',
            'iron': 'Fe',
            'cobalt': 'Co',
            'zinc': 'Zn'
        }
        
        # Get the search term (symbol or name)
        search_term = metal_map.get(filter_lower, catalyst_filter)
        
        # Filter catalyst column (case-insensitive)
        mask = df['Catalyst'].str.contains(search_term, case=False, na=False)
        
        return df[mask]
    
    def _predict_reaction_type(
        self,
        type_a: str,
        type_b: str,
    ) -> Tuple[Optional[str], float]:
        """
        Predict reaction type based on reactant type combination.

        Returns:
            (reaction_type, confidence_score)
        """
        key = _reactant_key([type_a, type_b])

        if key not in self.reaction_type_patterns:
            return None, 0.0

        rxn_counts = self.reaction_type_patterns[key]
        if not rxn_counts:
            return None, 0.0
        
        # Most common reaction type
        top_rxn = rxn_counts.most_common(1)[0]
        reaction_type = top_rxn[0]
        count = top_rxn[1]
        total = sum(rxn_counts.values())
        confidence = count / total
        
        return reaction_type, confidence
    
    def _calculate_confidence_score(
        self,
        avg_z_score: float,
        num_experiments: int,
        avg_yield: float
    ) -> float:
        """
        Calculate confidence score combining multiple factors.
        
        Formula: weighted combination of z-score (primary), sample size, and avg yield
        Z-score is the primary metric as it measures the success of a condition.
        """
        # Normalize sample size (log scale, cap at 100)
        sample_score = min(num_experiments, 100) / 100.0
        
        # Z-score weight (primary factor)
        # Normalize z-score: typical range is -3 to +3, scale to 0-1
        # Values above 3 are capped at 1.0, below -3 at 0.0
        z_score_normalized = max(0.0, min(1.0, (avg_z_score + 3.0) / 6.0))
        
        # Average yield weight (secondary factor)
        yield_weight = avg_yield / 100.0
        
        # Combined score with weights
        confidence = (
            0.6 * z_score_normalized +  # 60% from z-score (primary)
            0.25 * sample_score +        # 25% from sample size
            0.15 * yield_weight          # 15% from avg yield
        ) * 100.0
        
        return confidence
    
    def _aggregate_conditions(
        self,
        matched_df: pd.DataFrame,
        top_k: int = 10,
        min_experiments: int = 1
    ) -> List[ConditionRecommendation]:
        """
        Aggregate and rank condition combinations from matched experiments.
        
        Strategy:
        1. Group by (catalyst, ligand, base, solvent) combination
        2. Calculate z-score statistics (primary ranking metric)
        3. Calculate success rate (yield > 50)
        4. Calculate avg/median yield
        5. Compute confidence score (weighted by z-score)
        6. Rank by average z-score (primary), then confidence score
        """
        recommendations = []
        
        # Define success threshold
        SUCCESS_THRESHOLD = 50.0
        
        # Group by condition combination
        condition_cols = ['Catalyst', 'Ligand', 'Base', 'Solvent']
        optional_cols = ['Secondary Solvent', 'Additive', 'Coupling Reagent']
        
        grouped = matched_df.groupby(condition_cols, dropna=False)
        
        for condition_tuple, group_df in grouped:
            if len(group_df) < min_experiments:
                continue
            
            # Extract condition components
            catalyst, ligand, base, solvent = condition_tuple
            
            # Get optional components (most common values)
            sec_solvent = group_df['Secondary Solvent'].mode().iloc[0] if not group_df['Secondary Solvent'].isna().all() else None
            additive = group_df['Additive'].mode().iloc[0] if not group_df['Additive'].isna().all() else None
            coupling_reagent = group_df['Coupling Reagent'].mode().iloc[0] if not group_df['Coupling Reagent'].isna().all() else None
            
            # Calculate statistics
            yields = group_df['AREA_TOTAL_REDUCED']
            num_exp = len(group_df)
            success_count = (yields > SUCCESS_THRESHOLD).sum()
            success_rate = (success_count / num_exp) * 100.0
            avg_yield = yields.mean()
            median_yield = yields.median()
            
            # Z-score statistics (primary ranking metric)
            z_scores = group_df['z-Score']
            avg_z_score = z_scores.mean()
            z_min = z_scores.min()
            z_max = z_scores.max()
            
            # Confidence score (uses z-score as primary factor)
            confidence = self._calculate_confidence_score(
                avg_z_score, num_exp, avg_yield
            )
            
            # Reaction type (most common)
            reaction_type = group_df['Reaction_Type_Standardized'].mode().iloc[0] if not group_df['Reaction_Type_Standardized'].isna().all() else None
            
            # Reactant types (from first row)
            reactant_types = (
                group_df.iloc[0]['Reactant_A_Type'],
                group_df.iloc[0]['Reactant_B_Type']
            )
            
            rec = ConditionRecommendation(
                catalyst=catalyst if pd.notna(catalyst) else "",
                ligand=ligand if pd.notna(ligand) else "",
                base=base if pd.notna(base) else "",
                solvent=solvent if pd.notna(solvent) else "",
                secondary_solvent=sec_solvent if pd.notna(sec_solvent) else None,
                additive=additive if pd.notna(additive) else None,
                coupling_reagent=coupling_reagent if pd.notna(coupling_reagent) else None,
                success_rate=success_rate,
                avg_yield=avg_yield,
                median_yield=median_yield,
                num_experiments=num_exp,
                avg_z_score=avg_z_score,
                confidence_score=confidence,
                reaction_type=reaction_type,
                reactant_types=reactant_types,
                z_score_range=(z_min, z_max)
            )
            
            recommendations.append(rec)
        
        # Sort by average z-score (primary metric), then confidence score
        recommendations.sort(key=lambda x: (x.avg_z_score, x.confidence_score), reverse=True)
        
        return recommendations[:top_k]
    
    def recommend(
        self,
        reactant_a_smiles: str,
        reactant_b_smiles: Optional[str] = None,
        top_k: int = 10,
        min_experiments: int = 2,
        reaction_type_filter: Optional[str] = None,
        catalyst_filter: Optional[str] = None
    ) -> HTERecommendationResult:
        """
        Recommend conditions based on reactant SMILES.
        
        Args:
            reactant_a_smiles: SMILES of first reactant
            reactant_b_smiles: SMILES of second reactant (optional)
            top_k: Number of recommendations to return
            min_experiments: Minimum experiments for a condition to be recommended
            reaction_type_filter: Optional filter for specific reaction type
            catalyst_filter: Optional filter by metal type (e.g., 'Pd', 'Cu', 'Ni', 'palladium', 'copper')
        
        Returns:
            HTERecommendationResult with ranked condition recommendations
        """
        result = HTERecommendationResult(
            reactant_a_smiles=reactant_a_smiles,
            reactant_b_smiles=reactant_b_smiles
        )
        
        # Step 1: Detect reactant types
        type_a, cat_a = self._detect_reactant_types(reactant_a_smiles)
        result.reactant_a_type = ",".join(type_a) if type_a else ""
        result.reactant_a_category = cat_a
        
        if reactant_b_smiles:
            type_b, cat_b = self._detect_reactant_types(reactant_b_smiles)
            result.reactant_b_type = ",".join(type_b) if type_b else ""
            result.reactant_b_category = cat_b
        else:
            type_b, cat_b = [], ""
            result.reactant_b_type = ""
            result.reactant_b_category = ""
        
        # If no type detected, return empty
        if not type_a:
            return result
        
        # Step 3: Match against database
        key = _reactant_key(list(type_a) + list(type_b))
        
        matched_df = None
        if key in self.indexed_data:
            matched_df = self.indexed_data[key].copy()
            result.matched_motifs = (result.reactant_a_type, result.reactant_b_type)
        else:
            list_a = type_a or [""]
            list_b = type_b or [""]
            for ma in list_a:
                for mb in list_b:
                    if not ma and not mb:
                        continue
                    candidate = _reactant_key([ma, mb])
                    if candidate in self.indexed_data:
                        matched_df = self.indexed_data[candidate].copy()
                        result.matched_motifs = (ma, mb)
                        result.is_fallback_match = True
                        break
                if matched_df is not None:
                    break

        if matched_df is None:
            return result
        
        # Step 2: Predict reaction type (using matched motifs)
        if result.matched_motifs:
            pred_rxn, rxn_conf = self._predict_reaction_type(
                result.reactant_a_type, result.reactant_b_type
            )
            result.predicted_reaction_type = pred_rxn
            result.reaction_type_confidence = rxn_conf
        
        # Apply reaction type filter if specified
        if reaction_type_filter:
            matched_df = matched_df[matched_df['Reaction_Type_Standardized'] == reaction_type_filter]
        
        # Apply catalyst filter if specified
        if catalyst_filter:
            matched_df = self._filter_by_catalyst(matched_df, catalyst_filter)
        
        result.total_matching_experiments = len(matched_df)
        result.database_coverage = (len(matched_df) / len(self.df)) * 100.0
        
        # Step 4: Aggregate and rank conditions
        if len(matched_df) > 0:
            result.recommendations = self._aggregate_conditions(
                matched_df, top_k, min_experiments
            )
        
        return result
    
    def generate_screening_set(
        self,
        reactant_a_smiles: str,
        reactant_b_smiles: Optional[str] = None,
        num_conditions: int = 24,
        min_experiments: int = 1,
        reaction_type_filter: Optional[str] = None,
        catalyst_filter: Optional[str] = None,
        diversity_strategy: str = "balanced"
    ) -> HTERecommendationResult:
        """
        Generate a diverse set of conditions for HTE screening (up to 24 for standard plate).
        
        This is the PRIMARY use case for HTE systems - generating a group of diverse conditions
        to test in parallel on a screening plate.
        
        Args:
            reactant_a_smiles: SMILES of first reactant
            reactant_b_smiles: SMILES of second reactant (optional)
            num_conditions: Number of conditions to generate (default 24 for 4x6 plate)
            min_experiments: Minimum experiments for a condition to be included
            reaction_type_filter: Optional filter for specific reaction type
            catalyst_filter: Optional filter by metal type (e.g., 'Pd', 'Cu', 'Ni')
            diversity_strategy: Strategy for condition selection:
                - "balanced": Mix of top performers + diverse alternatives
                - "top_performers": Focus on highest z-score conditions
                - "diverse": Maximize diversity across reagent space
        
        Returns:
            HTERecommendationResult with up to num_conditions diverse conditions
        """
        # Get initial recommendations with larger pool
        initial_top_k = max(num_conditions * 3, 50)  # Get 3x more for diversity selection
        
        result = self.recommend(
            reactant_a_smiles=reactant_a_smiles,
            reactant_b_smiles=reactant_b_smiles,
            top_k=initial_top_k,
            min_experiments=min_experiments,
            reaction_type_filter=reaction_type_filter,
            catalyst_filter=catalyst_filter
        )
        
        if len(result.recommendations) == 0:
            return result
        
        # Apply diversity strategy
        if diversity_strategy == "top_performers":
            # Simply take top N by z-score
            selected = result.recommendations[:num_conditions]
        
        elif diversity_strategy == "diverse":
            # Maximize diversity across all reagent dimensions
            selected = self._select_diverse_conditions(
                result.recommendations, 
                num_conditions,
                prioritize_performance=False
            )
        
        else:  # "balanced" (default)
            # Take top performers + diverse alternatives
            num_top = min(num_conditions // 3, len(result.recommendations))  # ~1/3 top performers
            num_diverse = num_conditions - num_top
            
            top_performers = result.recommendations[:num_top]
            remaining = result.recommendations[num_top:]
            
            if remaining:
                diverse_picks = self._select_diverse_conditions(
                    remaining, 
                    num_diverse,
                    prioritize_performance=True
                )
                selected = top_performers + diverse_picks
            else:
                selected = top_performers
        
        result.recommendations = selected[:num_conditions]
        return result
    
    def _select_diverse_conditions(
        self, 
        conditions: List[ConditionRecommendation],
        num_to_select: int,
        prioritize_performance: bool = True
    ) -> List[ConditionRecommendation]:
        """
        Select diverse conditions maximizing reagent variation.
        
        Args:
            conditions: Pool of conditions to select from
            num_to_select: Number of conditions to select
            prioritize_performance: If True, weight selection by z-score
        
        Returns:
            List of diverse conditions
        """
        if len(conditions) <= num_to_select:
            return conditions
        
        selected = []
        remaining = list(conditions)
        
        # Track used reagents for diversity
        used_catalysts = set()
        used_ligands = set()
        used_bases = set()
        used_solvents = set()
        
        while len(selected) < num_to_select and remaining:
            best_score = -1
            best_idx = 0
            
            for i, cond in enumerate(remaining):
                # Calculate diversity score (how many new reagents does this add?)
                diversity_score = 0
                if cond.catalyst not in used_catalysts:
                    diversity_score += 1
                if cond.ligand not in used_ligands:
                    diversity_score += 1
                if cond.base not in used_bases:
                    diversity_score += 1
                if cond.solvent not in used_solvents:
                    diversity_score += 1
                
                # Combine diversity with performance if requested
                if prioritize_performance:
                    # Normalize z-score to 0-1 range (assuming typical range -3 to +3)
                    normalized_zscore = (cond.avg_z_score + 3) / 6.0
                    normalized_zscore = max(0, min(1, normalized_zscore))
                    
                    # 60% diversity, 40% performance
                    combined_score = (diversity_score / 4.0) * 0.6 + normalized_zscore * 0.4
                else:
                    combined_score = diversity_score
                
                if combined_score > best_score:
                    best_score = combined_score
                    best_idx = i
            
            # Select best condition
            selected_cond = remaining.pop(best_idx)
            selected.append(selected_cond)
            
            # Update used reagents
            used_catalysts.add(selected_cond.catalyst)
            used_ligands.add(selected_cond.ligand)
            used_bases.add(selected_cond.base)
            used_solvents.add(selected_cond.solvent)
        
        return selected
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics"""
        if self.df is None:
            return {}
        
        return {
            'total_experiments': len(self.df),
            'reaction_types': self.df['Reaction_Type_Standardized'].nunique(),
            'unique_type_combinations': len(self.indexed_data),
            'success_rate_overall': (self.df['AREA_TOTAL_REDUCED'] > 50).mean() * 100,
            'avg_yield': self.df['AREA_TOTAL_REDUCED'].mean(),
            'catalysts': self.df['Catalyst'].nunique(),
            'ligands': self.df['Ligand'].nunique(),
            'bases': self.df['Base'].nunique(),
            'solvents': self.df['Solvent'].nunique()
        }


def format_recommendation(rec: ConditionRecommendation, rank: int = 1) -> str:
    """Format a single recommendation for display"""
    lines = [
        f"\n{'='*80}",
        f"Recommendation #{rank}",
        f"{'='*80}",
        f"⭐ Avg Z-Score: {rec.avg_z_score:.2f} (Primary Ranking Metric)",
        f"Confidence Score: {rec.confidence_score:.1f}/100",
        f"Success Rate: {rec.success_rate:.1f}% ({rec.num_experiments} experiments)",
        f"Avg Yield: {rec.avg_yield:.1f}% | Median: {rec.median_yield:.1f}%",
        f"",
        f"🧪 CONDITIONS:",
        f"  Catalyst: {rec.catalyst}",
        f"  Ligand: {rec.ligand}",
        f"  Base: {rec.base}",
        f"  Solvent: {rec.solvent}"
    ]
    
    if rec.secondary_solvent:
        lines.append(f"  Secondary Solvent: {rec.secondary_solvent}")
    if rec.additive:
        lines.append(f"  Additive: {rec.additive}")
    if rec.coupling_reagent:
        lines.append(f"  Coupling Reagent: {rec.coupling_reagent}")
    
    lines.extend([
        f"",
        f"📊 STATISTICS:",
        f"  Reaction Type: {rec.reaction_type}",
        f"  Z-Score: Avg={rec.avg_z_score:.2f}, Range=[{rec.z_score_range[0]:.2f}, {rec.z_score_range[1]:.2f}]",
        f"  Reactant Types: {rec.reactant_types[0]} + {rec.reactant_types[1]}"
    ])
    
    return "\n".join(lines)


def format_result(result: HTERecommendationResult) -> str:
    """Format complete recommendation result for display"""
    lines = [
        "\n" + "="*80,
        "HTE-BASED CONDITION RECOMMENDATION",
        "="*80,
        f"Reactant A: {result.reactant_a_smiles}",
        f"  Type: {result.reactant_a_type} ({result.reactant_a_category})"
    ]
    
    if result.reactant_b_smiles:
        lines.extend([
            f"Reactant B: {result.reactant_b_smiles}",
            f"  Type: {result.reactant_b_type} ({result.reactant_b_category})"
        ])
    
    lines.extend([
        f"",
        f"🎯 PREDICTED REACTION TYPE: {result.predicted_reaction_type}",
        f"   Confidence: {result.reaction_type_confidence*100:.1f}%",
        f"",
        f"📊 DATABASE MATCH:",
        f"   {result.total_matching_experiments} matching experiments",
        f"   ({result.database_coverage:.2f}% of database)",
        f"",
        f"🏆 TOP RECOMMENDATIONS: {len(result.recommendations)} conditions found",
        f"   (Ranked by Average Z-Score)",
        "="*80
    ])
    
    # Add individual recommendations
    for i, rec in enumerate(result.recommendations, 1):
        lines.append(format_recommendation(rec, i))
    
    return "\n".join(lines)
