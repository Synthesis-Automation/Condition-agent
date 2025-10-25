"""
Simple Condition Recommender - Proof of Concept
Based on z-Score Peaks dataset with hierarchical exact matching.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json

class SimpleConditionRecommender:
    """
    Simple yet reliable condition recommender using hierarchical exact matching.
    
    Recommendation strategy:
    1. Exact match: reaction_type + reactant_a + reactant_b
    2. Category match: reaction_type + reactant_a_category + reactant_b_category  
    3. Reaction type match: reaction_type only
    """
    
    def __init__(self, csv_path: str, zscore_threshold: float = 1.0):
        """
        Initialize recommender with dataset.
        
        Args:
            csv_path: Path to standardized z-Score CSV
            zscore_threshold: Minimum z-score for "high performer" (default 1.0)
        """
        print(f"Loading dataset from {csv_path}...")
        self.df = pd.read_csv(csv_path)
        self.zscore_threshold = zscore_threshold
        self.high_performers = self.df[self.df['z-Score'] > zscore_threshold]
        
        print(f"Total reactions: {len(self.df):,}")
        print(f"High performers (z > {zscore_threshold}): {len(self.high_performers):,} ({len(self.high_performers)/len(self.df)*100:.1f}%)")
    
    def recommend(self, 
                  reaction_type: str,
                  reactant_a: str,
                  reactant_b: str,
                  functional_groups: Optional[List[str]] = None,
                  top_n: int = 3,
                  min_precedents: int = 3) -> Dict:
        """
        Get condition recommendations for a given reaction.
        
        Args:
            reaction_type: Standardized reaction type (e.g., "Buchwald-Hartwig-C-N")
            reactant_a: Standardized reactant A (e.g., "ArBr")
            reactant_b: Standardized reactant B (e.g., "RNH2")
            functional_groups: Optional list of functional groups present
            top_n: Number of recommendations to return
            min_precedents: Minimum number of precedents required
            
        Returns:
            Dictionary with recommendations and metadata
        """
        
        # Level 1: Exact match
        exact_matches = self._exact_match(reaction_type, reactant_a, reactant_b)
        if len(exact_matches) >= min_precedents:
            print(f"✓ Found {len(exact_matches)} exact matches (high performers)")
            return self._aggregate_conditions(
                exact_matches, 
                match_level="exact",
                reaction_type=reaction_type,
                reactant_a=reactant_a,
                reactant_b=reactant_b,
                top_n=top_n
            )
        
        # Level 2: Category match
        reactant_a_category = self._get_category(reactant_a)
        reactant_b_category = self._get_category(reactant_b)
        
        if reactant_a_category and reactant_b_category:
            category_matches = self._category_match(reaction_type, reactant_a_category, reactant_b_category)
            if len(category_matches) >= min_precedents:
                print(f"⚠ No exact match. Found {len(category_matches)} category matches")
                print(f"  Category: {reactant_a_category} + {reactant_b_category}")
                return self._aggregate_conditions(
                    category_matches,
                    match_level="category",
                    reaction_type=reaction_type,
                    reactant_a=reactant_a,
                    reactant_b=reactant_b,
                    top_n=top_n
                )
        
        # Level 3: Reaction type only
        rxn_matches = self._reaction_type_match(reaction_type)
        if len(rxn_matches) > 0:
            print(f"⚠ No substrate match. Using general {reaction_type} conditions ({len(rxn_matches)} precedents)")
            return self._aggregate_conditions(
                rxn_matches,
                match_level="reaction_type",
                reaction_type=reaction_type,
                reactant_a=reactant_a,
                reactant_b=reactant_b,
                top_n=top_n
            )
        
        # No recommendations available
        return {
            "match_level": "none",
            "message": f"No precedents found for reaction type: {reaction_type}",
            "recommendations": []
        }
    
    def _exact_match(self, reaction_type: str, reactant_a: str, reactant_b: str) -> pd.DataFrame:
        """Filter for exact substrate match."""
        return self.high_performers[
            (self.high_performers['Reaction_Type_Standardized'] == reaction_type) &
            (self.high_performers['Reactant_A'] == reactant_a) &
            (self.high_performers['Reactant_B'] == reactant_b)
        ]
    
    def _category_match(self, reaction_type: str, reactant_a_category: str, reactant_b_category: str) -> pd.DataFrame:
        """Filter for category-level match."""
        return self.high_performers[
            (self.high_performers['Reaction_Type_Standardized'] == reaction_type) &
            (self.high_performers['Reactant_A_Category'] == reactant_a_category) &
            (self.high_performers['Reactant_B_Category'] == reactant_b_category)
        ]
    
    def _reaction_type_match(self, reaction_type: str) -> pd.DataFrame:
        """Filter for reaction type only."""
        return self.high_performers[
            self.high_performers['Reaction_Type_Standardized'] == reaction_type
        ]
    
    def _get_category(self, reactant_type: str) -> Optional[str]:
        """Get category for a reactant type by checking the dataset."""
        # Look up category from the dataset
        match = self.df[self.df['Reactant_A'] == reactant_type]
        if len(match) > 0:
            category = match.iloc[0]['Reactant_A_Category']
            if pd.notna(category) and category != '':
                return category
        
        match = self.df[self.df['Reactant_B'] == reactant_type]
        if len(match) > 0:
            category = match.iloc[0]['Reactant_B_Category']
            if pd.notna(category) and category != '':
                return category
        
        return None
    
    def _aggregate_conditions(self, 
                             df: pd.DataFrame, 
                             match_level: str,
                             reaction_type: str,
                             reactant_a: str,
                             reactant_b: str,
                             top_n: int = 3) -> Dict:
        """
        Aggregate and rank condition combinations.
        
        Returns top N conditions ranked by:
        - Frequency in high performers
        - Average z-score
        - Consistency (low std dev)
        """
        
        # Group by condition combination
        grouped = df.groupby(['Catalyst', 'Ligand', 'Base', 'Solvent']).agg({
            'z-Score': ['count', 'mean', 'std', 'min', 'max'],
            'AREA_TOTAL_REDUCED': 'mean'
        }).reset_index()
        
        # Flatten column names
        grouped.columns = ['Catalyst', 'Ligand', 'Base', 'Solvent', 
                          'count', 'avg_zscore', 'std_zscore', 'min_zscore', 'max_zscore',
                          'avg_area']
        
        # Calculate composite score
        total_high = len(df)
        grouped['frequency_weight'] = grouped['count'] / total_high
        
        max_z = df['z-Score'].max()
        grouped['zscore_weight'] = (grouped['avg_zscore'] - self.zscore_threshold) / (max_z - self.zscore_threshold)
        grouped['zscore_weight'] = grouped['zscore_weight'].clip(0, 1)
        
        grouped['consistency_weight'] = 1 / (1 + grouped['std_zscore'].fillna(0))
        
        # Composite score: 50% frequency, 30% z-score, 20% consistency
        grouped['score'] = (
            0.5 * grouped['frequency_weight'] +
            0.3 * grouped['zscore_weight'] +
            0.2 * grouped['consistency_weight']
        )
        
        # Sort by score
        grouped = grouped.sort_values('score', ascending=False)
        
        # Get top N
        top_conditions = grouped.head(top_n)
        
        # Format output
        recommendations = []
        for idx, row in top_conditions.iterrows():
            recommendations.append({
                'rank': len(recommendations) + 1,
                'catalyst': str(row['Catalyst']) if pd.notna(row['Catalyst']) else None,
                'ligand': str(row['Ligand']) if pd.notna(row['Ligand']) else None,
                'base': str(row['Base']) if pd.notna(row['Base']) else None,
                'solvent': str(row['Solvent']) if pd.notna(row['Solvent']) else None,
                'confidence_score': round(float(row['score']), 3),
                'evidence': {
                    'successful_cases': int(row['count']),
                    'total_precedents': total_high,
                    'success_rate': f"{row['count']/total_high*100:.1f}%",
                    'avg_zscore': round(float(row['avg_zscore']), 2),
                    'zscore_range': [round(float(row['min_zscore']), 2), 
                                    round(float(row['max_zscore']), 2)]
                }
            })
        
        return {
            'reaction_type': reaction_type,
            'substrate': {
                'reactant_a': reactant_a,
                'reactant_b': reactant_b
            },
            'match_level': match_level,
            'total_precedents': total_high,
            'recommendations': recommendations,
            'metadata': {
                'zscore_threshold': self.zscore_threshold,
                'match_explanation': self._get_match_explanation(match_level)
            }
        }
    
    def _get_match_explanation(self, match_level: str) -> str:
        """Get explanation for match level."""
        explanations = {
            'exact': 'Exact match found for your substrate combination. High confidence.',
            'category': 'No exact match, but similar substrates found (category match). Medium confidence.',
            'reaction_type': 'No substrate match. Using general conditions for this reaction type. Lower confidence - consider screening.'
        }
        return explanations.get(match_level, 'Unknown match level')


def main():
    """Example usage."""
    
    # Initialize recommender
    csv_path = Path(__file__).parent / "z-Score Peaks CLEANED.csv"
    recommender = SimpleConditionRecommender(csv_path)
    
    print("\n" + "=" * 80)
    print("EXAMPLE 1: Buchwald-Hartwig-C-N with ArBr + RNH2")
    print("=" * 80)
    
    result = recommender.recommend(
        reaction_type="Buchwald-Hartwig-C-N",
        reactant_a="ArBr",
        reactant_b="RNH2",
        top_n=3
    )
    
    print(json.dumps(result, indent=2))
    
    print("\n" + "=" * 80)
    print("EXAMPLE 2: Suzuki-Miyaura with ArCl + ArB(OH)2")
    print("=" * 80)
    
    result = recommender.recommend(
        reaction_type="Suzuki-Miyaura",
        reactant_a="ArCl",
        reactant_b="ArB(OH)2",
        top_n=3
    )
    
    print(json.dumps(result, indent=2))
    
    print("\n" + "=" * 80)
    print("EXAMPLE 3: Less common - C-O Coupling with ArBr + ROH-primary")
    print("=" * 80)
    
    result = recommender.recommend(
        reaction_type="CO-Coupling",
        reactant_a="ArBr",
        reactant_b="ROH-primary",
        top_n=3
    )
    
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
