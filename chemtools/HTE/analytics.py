"""
HTE Database Analytics Tools

Provides analytical functions to explore and understand the HTE database:
- List reactant pairs by reaction type and/or catalyst
- Analyze catalyst distributions for specific reactions
- Explore reactant type combinations
- Statistical summaries and coverage analysis
"""

from typing import List, Dict, Optional, Tuple, Any
from collections import defaultdict, Counter
import json
import pandas as pd
from pathlib import Path


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


def _dedupe_list(values: List[str]) -> List[str]:
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
            catalyst_types = _ensure_list(record.get("catalyst_type"))
            conditions = record.get("conditions") or {}
            metrics = record.get("metrics") or {}

            row = {
                "Reaction_Type_Standardized": record.get("reaction_type") or "Unknown",
                "Reactant_A_Type": reactant_types[0] if len(reactant_types) > 0 else "",
                "Reactant_B_Type": reactant_types[1] if len(reactant_types) > 1 else "",
                "Reactant_A_Category": reactant_categories[0] if len(reactant_categories) > 0 else "",
                "Reactant_B_Category": reactant_categories[1] if len(reactant_categories) > 1 else "",
                "Catalyst_Type": _format_list(catalyst_types),
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


class HTEAnalytics:
    """
    Analytics interface for the HTE database.
    
    Enables exploration and analysis of:
    - Reactant type combinations
    - Catalyst usage patterns
    - Reaction type coverage
    - Statistical distributions
    """
    
    def __init__(self, hte_db_path: str = "data/HTE_db/HTE_canonical.csv"):
        """Initialize analytics with HTE database"""
        self.db_path = Path(hte_db_path)
        self.df: Optional[pd.DataFrame] = None
        self._load_database()
    
    def _load_database(self):
        """Load HTE database"""
        if not self.db_path.exists():
            raise FileNotFoundError(f"HTE database not found: {self.db_path}")
        
        if self.db_path.suffix.lower() == ".jsonl":
            self.df = _load_hte_jsonl(self.db_path)
        else:
            self.df = pd.read_csv(self.db_path)
            
            # Map new canonical CSV columns to internal standard names
            column_mapping = {
                "reaction_type": "Reaction_Type_Standardized",
                "reactant_1": "Reactant_A_Type",
                "reactant_2": "Reactant_B_Type",
                "yield": "AREA_TOTAL_REDUCED",
                "z_score": "z-Score",
                "catalyst": "Catalyst",
                "ligand": "Ligand",
                "base": "Base",
                "solvent": "Solvent",
                "additive": "Additive"
            }
            
            # Rename columns if they exist
            self.df = self.df.rename(columns={k: v for k, v in column_mapping.items() if k in self.df.columns})
            
            # Ensure missing columns are present
            required_cols = [
                "Reaction_Type_Standardized", "Reactant_A_Type", "Reactant_B_Type",
                "Catalyst", "Ligand", "Base", "Solvent", "Additive",
                "Secondary Solvent", "Coupling Reagent", "AREA_TOTAL_REDUCED", "z-Score",
                "Reactant_A_Category", "Reactant_B_Category"
            ]
            for col in required_cols:
                if col not in self.df.columns:
                    self.df[col] = "" if col not in ["AREA_TOTAL_REDUCED", "z-Score"] else 0.0

        print(f"📊 Loaded HTE database: {len(self.df):,} experiments")
    
    def _filter_by_catalyst_type(self, df: pd.DataFrame, catalyst_filter: str) -> pd.DataFrame:
        if not catalyst_filter:
            return df
        filter_lower = catalyst_filter.lower()
        metal_map = {
            "palladium": "Pd",
            "copper": "Cu",
            "nickel": "Ni",
            "iridium": "Ir",
            "rhodium": "Rh",
            "ruthenium": "Ru",
            "platinum": "Pt",
            "gold": "Au",
            "silver": "Ag",
            "iron": "Fe",
            "cobalt": "Co",
            "zinc": "Zn",
            "organocatalyst": "organocatalyst",
            "organic catalyst": "organocatalyst",
        }
        search_term = metal_map.get(filter_lower, catalyst_filter)

        if "Catalyst_Type" in df.columns:
            mask = df["Catalyst_Type"].str.contains(search_term, case=False, na=False)
            if mask.any():
                return df[mask]

        return df[df["Catalyst"].str.contains(search_term, case=False, na=False)]

    def list_reactant_pairs(
        self,
        reaction_type: Optional[str] = None,
        catalyst_filter: Optional[str] = None,
        min_experiments: int = 1,
        sort_by: str = "count"
    ) -> pd.DataFrame:
        """
        List all reactant type pairs in the database with optional filters.
        
        Args:
            reaction_type: Filter by reaction type (e.g., "Suzuki", "C-N Coupling")
            catalyst_filter: Filter by catalyst metal (e.g., "Pd", "Cu", "palladium")
            min_experiments: Minimum number of experiments for a pair to be included
            sort_by: Sort by "count" (num experiments) or "success_rate" (avg yield)
        
        Returns:
            DataFrame with columns:
                - Reactant_A_Type
                - Reactant_B_Type
                - Reaction_Type
                - Num_Experiments
                - Avg_Yield
                - Success_Rate (% yield > 50)
                - Top_Catalyst
        """
        df = self.df.copy()
        
        # Apply filters
        if reaction_type:
            df = df[df['Reaction_Type_Standardized'].str.contains(reaction_type, case=False, na=False)]
        
        if catalyst_filter:
            df = self._filter_by_catalyst_type(df, catalyst_filter)
        
        # Group by reactant pairs
        grouped = df.groupby(['Reactant_A_Type', 'Reactant_B_Type', 'Reaction_Type_Standardized']).agg({
            'AREA_TOTAL_REDUCED': ['count', 'mean', 'median', lambda x: (x > 50).sum() / len(x) * 100],
            'Catalyst': lambda x: x.mode().iloc[0] if not x.mode().empty else None
        }).reset_index()
        
        # Flatten column names
        grouped.columns = [
            'Reactant_A_Type', 'Reactant_B_Type', 'Reaction_Type',
            'Num_Experiments', 'Avg_Yield', 'Median_Yield', 'Success_Rate', 'Top_Catalyst'
        ]
        
        # Filter by minimum experiments
        grouped = grouped[grouped['Num_Experiments'] >= min_experiments]
        
        # Sort
        if sort_by == "count":
            grouped = grouped.sort_values('Num_Experiments', ascending=False)
        elif sort_by == "success_rate":
            grouped = grouped.sort_values('Success_Rate', ascending=False)
        
        return grouped.reset_index(drop=True)
    
    def get_catalyst_stats(
        self,
        reaction_type: Optional[str] = None,
        reactant_a_type: Optional[str] = None,
        reactant_b_type: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Get catalyst usage statistics with optional filters.
        
        Args:
            reaction_type: Filter by reaction type
            reactant_a_type: Filter by reactant A type
            reactant_b_type: Filter by reactant B type
        
        Returns:
            DataFrame with columns:
                - Catalyst
                - Metal (extracted metal symbol)
                - Num_Experiments
                - Avg_Yield
                - Success_Rate
                - Reaction_Types (list of reaction types using this catalyst)
        """
        df = self.df.copy()
        
        # Apply filters
        if reaction_type:
            df = df[df['Reaction_Type_Standardized'].str.contains(reaction_type, case=False, na=False)]
        if reactant_a_type:
            df = df[df['Reactant_A_Type'] == reactant_a_type]
        if reactant_b_type:
            df = df[df['Reactant_B_Type'] == reactant_b_type]
        
        # Extract metal from catalyst name
        def extract_metal(catalyst_name):
            """Extract metal symbol from catalyst name"""
            if pd.isna(catalyst_name):
                return None
            metals = ['Pd', 'Cu', 'Ni', 'Ir', 'Rh', 'Ru', 'Pt', 'Au', 'Ag', 'Fe', 'Co', 'Zn']
            for metal in metals:
                if metal in str(catalyst_name):
                    return metal
            return 'Other'
        
        df['Metal'] = df['Catalyst'].apply(extract_metal)
        
        # Group by catalyst
        grouped = df.groupby('Catalyst').agg({
            'AREA_TOTAL_REDUCED': ['count', 'mean', lambda x: (x > 50).sum() / len(x) * 100],
            'Metal': 'first',
            'Reaction_Type_Standardized': lambda x: ', '.join(sorted(set(x.dropna())))
        }).reset_index()
        
        grouped.columns = [
            'Catalyst', 'Num_Experiments', 'Avg_Yield', 'Success_Rate',
            'Metal', 'Reaction_Types'
        ]
        
        # Reorder columns
        grouped = grouped[['Catalyst', 'Metal', 'Num_Experiments', 'Avg_Yield', 'Success_Rate', 'Reaction_Types']]
        
        return grouped.sort_values('Num_Experiments', ascending=False).reset_index(drop=True)
    
    def get_reaction_type_summary(self) -> pd.DataFrame:
        """
        Get summary statistics for all reaction types in the database.
        
        Returns:
            DataFrame with columns:
                - Reaction_Type
                - Num_Experiments
                - Num_Reactant_Pairs
                - Num_Catalysts
                - Avg_Yield
                - Success_Rate
                - Top_Catalyst
                - Top_Reactant_Pair
        """
        grouped = self.df.groupby('Reaction_Type_Standardized').agg({
            'AREA_TOTAL_REDUCED': ['count', 'mean', lambda x: (x > 50).sum() / len(x) * 100],
            'Catalyst': ['nunique', lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else None],
            'Reactant_A_Type': lambda x: len(set(zip(x, self.df.loc[x.index, 'Reactant_B_Type'])))
        }).reset_index()
        
        # Get top reactant pair for each reaction
        def get_top_pair(reaction_type):
            sub_df = self.df[self.df['Reaction_Type_Standardized'] == reaction_type]
            pair_counts = sub_df.groupby(['Reactant_A_Type', 'Reactant_B_Type']).size()
            if len(pair_counts) > 0:
                top_pair = pair_counts.idxmax()
                return f"{top_pair[0]} + {top_pair[1]}"
            return None
        
        grouped.columns = [
            'Reaction_Type', 'Num_Experiments', 'Avg_Yield', 'Success_Rate',
            'Num_Catalysts', 'Top_Catalyst', 'Num_Reactant_Pairs'
        ]
        
        grouped['Top_Reactant_Pair'] = grouped['Reaction_Type'].apply(get_top_pair)
        
        # Reorder columns
        grouped = grouped[[
            'Reaction_Type', 'Num_Experiments', 'Num_Reactant_Pairs', 'Num_Catalysts',
            'Avg_Yield', 'Success_Rate', 'Top_Catalyst', 'Top_Reactant_Pair'
        ]]
        
        return grouped.sort_values('Num_Experiments', ascending=False).reset_index(drop=True)
    
    def analyze_metal_usage(self) -> Dict[str, Any]:
        """
        Analyze metal usage across the entire database.
        
        Returns:
            Dictionary with:
                - metal_distribution: DataFrame with metal counts and percentages
                - by_reaction_type: Dict[metal, Dict[reaction_type, count]]
                - total_experiments: Total number of experiments
        """
        # Extract metals
        def extract_metal(catalyst_name):
            if pd.isna(catalyst_name):
                return None
            metals = ['Pd', 'Cu', 'Ni', 'Ir', 'Rh', 'Ru', 'Pt', 'Au', 'Ag', 'Fe', 'Co', 'Zn']
            for metal in metals:
                if metal in str(catalyst_name):
                    return metal
            return 'Other'
        
        self.df['Metal'] = self.df['Catalyst'].apply(extract_metal)
        
        # Overall distribution
        metal_counts = self.df['Metal'].value_counts()
        metal_dist = pd.DataFrame({
            'Metal': metal_counts.index,
            'Num_Experiments': metal_counts.values,
            'Percentage': (metal_counts.values / len(self.df) * 100).round(2)
        })
        
        # By reaction type
        by_reaction = defaultdict(lambda: defaultdict(int))
        for _, row in self.df.iterrows():
            if pd.notna(row['Metal']) and pd.notna(row['Reaction_Type_Standardized']):
                by_reaction[row['Metal']][row['Reaction_Type_Standardized']] += 1
        
        return {
            'metal_distribution': metal_dist,
            'by_reaction_type': dict(by_reaction),
            'total_experiments': len(self.df)
        }
    
    def find_similar_pairs(
        self,
        reactant_a_type: str,
        reactant_b_type: str,
        similarity_criteria: str = "reaction_type"
    ) -> pd.DataFrame:
        """
        Find reactant pairs similar to the given pair.
        
        Args:
            reactant_a_type: Reactant A type to search for
            reactant_b_type: Reactant B type to search for
            similarity_criteria: How to define similarity:
                - "reaction_type": Same reaction type
                - "catalyst": Same catalyst metal
                - "both": Same reaction type AND catalyst metal
        
        Returns:
            DataFrame of similar reactant pairs with their statistics
        """
        # Find the query pair
        query_df = self.df[
            (self.df['Reactant_A_Type'] == reactant_a_type) &
            (self.df['Reactant_B_Type'] == reactant_b_type)
        ]
        
        if len(query_df) == 0:
            return pd.DataFrame()
        
        # Get the dominant reaction type and catalyst
        dominant_reaction = query_df['Reaction_Type_Standardized'].mode().iloc[0] if not query_df['Reaction_Type_Standardized'].mode().empty else None
        
        def extract_metal(catalyst_name):
            if pd.isna(catalyst_name):
                return None
            metals = ['Pd', 'Cu', 'Ni', 'Ir', 'Rh', 'Ru', 'Pt', 'Au', 'Ag', 'Fe', 'Co', 'Zn']
            for metal in metals:
                if metal in str(catalyst_name):
                    return metal
            return 'Other'
        
        query_df['Metal'] = query_df['Catalyst'].apply(extract_metal)
        dominant_metal = query_df['Metal'].mode().iloc[0] if not query_df['Metal'].mode().empty else None
        
        # Find similar pairs
        similar_df = self.df.copy()
        similar_df['Metal'] = similar_df['Catalyst'].apply(extract_metal)
        
        if similarity_criteria == "reaction_type":
            similar_df = similar_df[similar_df['Reaction_Type_Standardized'] == dominant_reaction]
        elif similarity_criteria == "catalyst":
            similar_df = similar_df[similar_df['Metal'] == dominant_metal]
        elif similarity_criteria == "both":
            similar_df = similar_df[
                (similar_df['Reaction_Type_Standardized'] == dominant_reaction) &
                (similar_df['Metal'] == dominant_metal)
            ]
        
        # Exclude the query pair itself
        similar_df = similar_df[
            ~((similar_df['Reactant_A_Type'] == reactant_a_type) &
              (similar_df['Reactant_B_Type'] == reactant_b_type))
        ]
        
        # Group by reactant pairs
        grouped = similar_df.groupby(['Reactant_A_Type', 'Reactant_B_Type', 'Reaction_Type_Standardized']).agg({
            'AREA_TOTAL_REDUCED': ['count', 'mean', lambda x: (x > 50).sum() / len(x) * 100],
            'Catalyst': lambda x: x.mode().iloc[0] if not x.mode().empty else None
        }).reset_index()
        
        grouped.columns = [
            'Reactant_A_Type', 'Reactant_B_Type', 'Reaction_Type',
            'Num_Experiments', 'Avg_Yield', 'Success_Rate', 'Top_Catalyst'
        ]
        
        return grouped.sort_values('Num_Experiments', ascending=False).reset_index(drop=True)
    
    def export_subset(
        self,
        output_path: str,
        reaction_type: Optional[str] = None,
        catalyst_filter: Optional[str] = None,
        reactant_a_type: Optional[str] = None,
        reactant_b_type: Optional[str] = None,
        min_yield: Optional[float] = None
    ) -> int:
        """
        Export a filtered subset of the database to CSV.
        
        Args:
            output_path: Path to save the CSV file
            reaction_type: Filter by reaction type
            catalyst_filter: Filter by catalyst metal
            reactant_a_type: Filter by reactant A type
            reactant_b_type: Filter by reactant B type
            min_yield: Minimum yield threshold
        
        Returns:
            Number of experiments exported
        """
        df = self.df.copy()
        
        # Apply filters
        if reaction_type:
            df = df[df['Reaction_Type_Standardized'].str.contains(reaction_type, case=False, na=False)]
        
        if catalyst_filter:
            df = self._filter_by_catalyst_type(df, catalyst_filter)
        
        if reactant_a_type:
            df = df[df['Reactant_A_Type'] == reactant_a_type]
        
        if reactant_b_type:
            df = df[df['Reactant_B_Type'] == reactant_b_type]
        
        if min_yield is not None:
            df = df[df['AREA_TOTAL_REDUCED'] >= min_yield]
        
        # Export
        df.to_csv(output_path, index=False)
        print(f"✅ Exported {len(df):,} experiments to {output_path}")
        
        return len(df)
