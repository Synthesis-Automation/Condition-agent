"""
HTE Database Analytics Tools

Provides analytical functions to explore and understand the HTE database:
- List reactant pairs by reaction type and/or catalyst
- Analyze catalyst distributions for specific reactions
- Explore reactant type combinations
- Statistical summaries and coverage analysis
"""

from typing import List, Dict, Optional, Tuple, Any, Iterable, Set
from collections import defaultdict, Counter
from functools import lru_cache
import json
import pandas as pd
from pathlib import Path
import re


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


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _collect_hte_files(db_path: Path) -> List[Path]:
    if db_path.is_file():
        return [db_path]
    if not db_path.exists():
        return []

    # Recursively collect all HTE tabular datasets under the DB root so new
    # dataset folders (for example `protocols/`) are included automatically.
    candidates: List[Path] = []
    candidates.extend(db_path.rglob("*.csv"))
    candidates.extend(db_path.rglob("*.jsonl"))

    seen = set()
    ordered: List[Path] = []
    for path in sorted(candidates, key=lambda p: str(p)):
        key = str(path.resolve())
        if key in seen:
            continue
        seen.add(key)
        ordered.append(path)
    return ordered


def _infer_source_group(source_path: Optional[Path]) -> str:
    if not source_path:
        return "unknown"
    parts = [part.lower() for part in source_path.parts]
    for part in parts:
        if part in ("literature", "datasets", "dataset"):
            return "literature"
        if part == "rules":
            return "rules"
        if part in ("protocols", "protocol"):
            return "literature"
        if part in ("motif", "motifs", "experiments", "experiment", "experiements"):
            return "experiments"
    return "other"


def _format_source_path(source_path: Optional[Path]) -> str:
    if not source_path:
        return ""
    try:
        return source_path.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(source_path)


_MOTIF_SPLIT_RE = re.compile(r"[|,]")
_COMPOUND_LOGIC_FILE = Path(__file__).resolve().parents[1] / "taxonomy" / "data" / "compound_logic.json"


@lru_cache(maxsize=1)
def _load_motif_sets() -> Dict[str, List[str]]:
    if not _COMPOUND_LOGIC_FILE.exists():
        return {}
    try:
        with _COMPOUND_LOGIC_FILE.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return {}
    raw_sets = payload.get("motif_sets") or {}
    motif_sets: Dict[str, List[str]] = {}
    for name, entry in raw_sets.items():
        members: List[str] = []
        if isinstance(entry, dict):
            members = entry.get("members") or []
        elif isinstance(entry, list):
            members = entry
        motif_sets[name] = [str(m).strip() for m in members if str(m).strip()]
    return motif_sets


def _split_motif_tokens(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, float) and pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    return [token.strip() for token in _MOTIF_SPLIT_RE.split(text) if token.strip()]


def _expand_macro_token(token: str, motif_sets: Dict[str, List[str]]) -> List[str]:
    token = token.strip()
    if token.startswith("@"):
        set_name = token[1:]
        members = motif_sets.get(set_name) or []
        if members:
            return members
    return [token]


def _expand_motif_tokens(tokens: Iterable[str], motif_sets: Dict[str, List[str]]) -> List[str]:
    expanded: List[str] = []
    for token in tokens:
        expanded.extend(_expand_macro_token(token, motif_sets))
    return _dedupe_list(expanded)


def _normalize_query_tokens(value: Optional[str], motif_sets: Dict[str, List[str]]) -> Set[str]:
    tokens = _split_motif_tokens(value)
    return set(_expand_motif_tokens(tokens, motif_sets))


def _field_matches_query(value: Any, query_tokens: Set[str], motif_sets: Dict[str, List[str]]) -> bool:
    if not query_tokens:
        return True
    field_tokens = set(_expand_motif_tokens(_split_motif_tokens(value), motif_sets))
    if not field_tokens:
        return False
    return bool(field_tokens & query_tokens)


def _normalize_hte_dataframe(df: pd.DataFrame, source_path: Optional[Path] = None) -> pd.DataFrame:
    df = df.copy()

    column_mapping = {
        "reaction_type": "Reaction_Type_Standardized",
        "detected_reaction_type": "Reaction_Type_Standardized",
        "Detected_Reaction_Type": "Reaction_Type_Standardized",
        "Reaction_Type": "Reaction_Type_Standardized",
        "reactant_1": "Reactant_A_Type",
        "reactant_2": "Reactant_B_Type",
        "yield": "AREA_TOTAL_REDUCED",
        "z_score": "z-Score",
        "catalyst": "Catalyst",
        "ligand": "Ligand",
        "base": "Base",
        "solvent": "Solvent",
        "additive": "Additive",
    }
    df = df.rename(columns={k: v for k, v in column_mapping.items() if k in df.columns})

    required_cols = [
        "Reaction_Type_Standardized", "Reactant_A_Type", "Reactant_B_Type",
        "Catalyst", "Ligand", "Base", "Solvent", "Additive",
        "Secondary Solvent", "Coupling Reagent", "AREA_TOTAL_REDUCED", "z-Score",
        "Reactant_A_Category", "Reactant_B_Category", "Source_File", "Source_Group",
    ]
    for col in required_cols:
        if col not in df.columns:
            if col in ("Source_File", "Source_Group"):
                df[col] = ""
            else:
                df[col] = "" if col not in ["AREA_TOTAL_REDUCED", "z-Score"] else 0.0

    if source_path is not None:
        df["Source_File"] = _format_source_path(source_path)
        df["Source_Group"] = _infer_source_group(source_path)

    return df


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
    
    def __init__(self, hte_db_path: str = "data/HTE_db"):
        """Initialize analytics with HTE database.

        Rule-derived screens are expected under `data/HTE_db/rules`.
        """
        self.db_path = Path(hte_db_path)
        self.df: Optional[pd.DataFrame] = None
        self._load_database()
    
    def _load_database(self):
        """Load HTE database"""
        if not self.db_path.exists():
            raise FileNotFoundError(f"HTE database not found: {self.db_path}")

        file_paths = _collect_hte_files(self.db_path)
        if not file_paths:
            raise FileNotFoundError(f"No HTE CSV/JSONL files found under: {self.db_path}")

        frames: List[pd.DataFrame] = []
        for path in file_paths:
            if path.suffix.lower() == ".jsonl":
                frame = _load_hte_jsonl(path)
            else:
                frame = pd.read_csv(path)
            frame = _normalize_hte_dataframe(frame, source_path=path)
            frames.append(frame)

        non_empty = [f for f in frames if not f.empty]
        self.df = pd.concat(non_empty, ignore_index=True) if non_empty else pd.DataFrame()

        print(f"馃搳 Loaded HTE database: {len(self.df):,} experiments")

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
        motif_sets = _load_motif_sets()
        
        # Apply filters
        if reaction_type:
            df = df[df['Reaction_Type_Standardized'].str.contains(reaction_type, case=False, na=False)]
        if reactant_a_type:
            query_a = _normalize_query_tokens(reactant_a_type, motif_sets)
            df = df[df['Reactant_A_Type'].apply(lambda value: _field_matches_query(value, query_a, motif_sets))]
        if reactant_b_type:
            query_b = _normalize_query_tokens(reactant_b_type, motif_sets)
            df = df[df['Reactant_B_Type'].apply(lambda value: _field_matches_query(value, query_b, motif_sets))]
        
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
    
    @staticmethod
    def _mode_or_none(series: pd.Series) -> Optional[str]:
        modes = series.mode()
        if len(modes) == 0:
            return None
        value = modes.iloc[0]
        return None if pd.isna(value) else str(value)

    @staticmethod
    def _top_nonempty_text(series: pd.Series) -> Optional[str]:
        if series is None:
            return None
        cleaned = series.dropna().astype(str).map(str.strip)
        cleaned = cleaned[~cleaned.isin(["", "nan", "None"])]
        if cleaned.empty:
            return None
        modes = cleaned.mode()
        if len(modes) > 0:
            return str(modes.iloc[0]).strip() or None
        counts = cleaned.value_counts()
        if counts.empty:
            return None
        return str(counts.index[0]).strip() or None

    def _build_reaction_type_detail_map(
        self,
        reaction_type: str,
        sub_df: pd.DataFrame,
        *,
        detail_top_k: int = 5,
    ) -> Dict[str, Any]:
        """Build a compact structured summary for a single reaction type."""
        top_k = max(1, int(detail_top_k))
        if sub_df.empty:
            return {
                "reaction_type": reaction_type,
                "num_experiments": 0,
                "top_reactant_pairs": [],
                "top_catalysts": [],
            }

        pair_group = sub_df.groupby(["Reactant_A_Type", "Reactant_B_Type"]).agg(
            Num_Experiments=("AREA_TOTAL_REDUCED", "count"),
            Avg_Yield=("AREA_TOTAL_REDUCED", "mean"),
            Success_Rate=("AREA_TOTAL_REDUCED", lambda x: (x > 50).sum() / len(x) * 100),
            Top_Catalyst=("Catalyst", self._mode_or_none),
        ).reset_index()
        pair_group = pair_group.sort_values(
            ["Num_Experiments", "Avg_Yield"],
            ascending=[False, False],
        )

        catalyst_group = sub_df.groupby("Catalyst").agg(
            Num_Experiments=("AREA_TOTAL_REDUCED", "count"),
            Avg_Yield=("AREA_TOTAL_REDUCED", "mean"),
            Success_Rate=("AREA_TOTAL_REDUCED", lambda x: (x > 50).sum() / len(x) * 100),
        ).reset_index()
        catalyst_group = catalyst_group.sort_values(
            ["Num_Experiments", "Avg_Yield"],
            ascending=[False, False],
        )

        top_pairs = []
        for _, row in pair_group.head(top_k).iterrows():
            top_pairs.append(
                {
                    "reactant_a_type": row["Reactant_A_Type"],
                    "reactant_b_type": row["Reactant_B_Type"],
                    "count": int(row["Num_Experiments"]),
                    "avg_yield": round(float(row["Avg_Yield"]), 3) if pd.notna(row["Avg_Yield"]) else None,
                    "success_rate": round(float(row["Success_Rate"]), 3) if pd.notna(row["Success_Rate"]) else None,
                    "top_catalyst": row["Top_Catalyst"],
                }
            )

        top_catalysts = []
        for _, row in catalyst_group.head(top_k).iterrows():
            top_catalysts.append(
                {
                    "catalyst": row["Catalyst"],
                    "count": int(row["Num_Experiments"]),
                    "avg_yield": round(float(row["Avg_Yield"]), 3) if pd.notna(row["Avg_Yield"]) else None,
                    "success_rate": round(float(row["Success_Rate"]), 3) if pd.notna(row["Success_Rate"]) else None,
                }
            )

        return {
            "reaction_type": reaction_type,
            "num_experiments": int(len(sub_df)),
            "top_reactant_pairs": top_pairs,
            "top_catalysts": top_catalysts,
        }

    def get_reaction_type_detailed_map(
        self,
        reaction_type: Optional[str] = None,
        *,
        min_rows: int = 1,
    ) -> pd.DataFrame:
        """
        Build a flat reaction-type detailed map similar to
        `reaction_type_detailed_map.csv`.

        One row corresponds to a reaction type + (FG1, FG2) subtype.
        `n_eln` is approximated as the count of unique non-empty `Source_File`
        values contributing to the subtype.
        """
        df = self.df.copy()
        if reaction_type:
            df = df[df["Reaction_Type_Standardized"].str.contains(reaction_type, case=False, na=False)]

        columns = [
            "Reaction Type",
            "FG1",
            "FG2",
            "n_rows",
            "n_eln",
            "top_coupling_reagent",
            "top_catalyst",
            "top_ligand",
            "top_base",
            "top_solvent",
            "top_additive",
            "ReactionType_Detailed",
            "ReactionType_SubtypeTag",
        ]
        if df.empty:
            return pd.DataFrame(columns=columns)

        grouped = (
            df.groupby(["Reaction_Type_Standardized", "Reactant_A_Type", "Reactant_B_Type"], dropna=False)
            .agg(
                n_rows=("AREA_TOTAL_REDUCED", "count"),
                n_eln=(
                    "Source_File",
                    lambda s: int(
                        pd.Series(s)
                        .dropna()
                        .astype(str)
                        .map(str.strip)
                        .loc[lambda x: x != ""]
                        .nunique()
                    )
                ),
                top_coupling_reagent=("Coupling Reagent", self._top_nonempty_text),
                top_catalyst=("Catalyst", self._top_nonempty_text),
                top_ligand=("Ligand", self._top_nonempty_text),
                top_base=("Base", self._top_nonempty_text),
                top_solvent=("Solvent", self._top_nonempty_text),
                top_additive=("Additive", self._top_nonempty_text),
            )
            .reset_index()
        )

        grouped.columns = [
            "Reaction Type",
            "FG1",
            "FG2",
            "n_rows",
            "n_eln",
            "top_coupling_reagent",
            "top_catalyst",
            "top_ligand",
            "top_base",
            "top_solvent",
            "top_additive",
        ]

        grouped["FG1"] = grouped["FG1"].fillna("").astype(str)
        grouped["FG2"] = grouped["FG2"].fillna("").astype(str)
        grouped["Reaction Type"] = grouped["Reaction Type"].fillna("Unknown").astype(str)

        if min_rows > 1:
            grouped = grouped[grouped["n_rows"] >= int(min_rows)]

        def _detailed_label(row: pd.Series) -> str:
            rxn = str(row["Reaction Type"]).strip() or "Unknown"
            fg1 = str(row["FG1"]).strip()
            fg2 = str(row["FG2"]).strip()
            if fg1 and fg2:
                return f"{rxn}: {fg1} + {fg2}"
            if fg1:
                return f"{rxn}: {fg1}"
            return rxn

        def _subtype_tag(row: pd.Series) -> str:
            rxn = str(row["Reaction Type"]).strip() or "Unknown"
            fg1 = str(row["FG1"]).strip()
            fg2 = str(row["FG2"]).strip()
            if fg1 and fg2:
                return f"{rxn}__{fg1} + {fg2}"
            if fg1:
                return f"{rxn}__{fg1}"
            return rxn

        grouped["ReactionType_Detailed"] = grouped.apply(_detailed_label, axis=1)
        grouped["ReactionType_SubtypeTag"] = grouped.apply(_subtype_tag, axis=1)

        for col in (
            "top_coupling_reagent",
            "top_catalyst",
            "top_ligand",
            "top_base",
            "top_solvent",
            "top_additive",
        ):
            grouped[col] = grouped[col].where(grouped[col].notna(), "")

        grouped = grouped.sort_values(
            ["Reaction Type", "n_rows", "FG1", "FG2"],
            ascending=[True, False, True, True],
        ).reset_index(drop=True)

        return grouped[columns]

    def get_reaction_type_summary(
        self,
        reaction_type: Optional[str] = None,
        *,
        include_detailed_map: bool = False,
        detail_top_k: int = 5,
    ) -> pd.DataFrame:
        """
        Get summary statistics for reaction types in the database.

        Args:
            reaction_type: Optional case-insensitive substring filter applied to
                `Reaction_Type_Standardized`.
            include_detailed_map: Include a serialized JSON `Detailed_Map`
                column with top reactant-pair and catalyst breakdowns.
            detail_top_k: Number of entries to keep in each top list in the
                detailed map.
        
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
                - Detailed_Map (optional JSON string)
        """
        df = self.df.copy()
        if reaction_type:
            df = df[df['Reaction_Type_Standardized'].str.contains(reaction_type, case=False, na=False)]

        base_columns = [
            'Reaction_Type', 'Num_Experiments', 'Num_Reactant_Pairs', 'Num_Catalysts',
            'Avg_Yield', 'Success_Rate', 'Top_Catalyst', 'Top_Reactant_Pair'
        ]
        if include_detailed_map:
            base_columns.append('Detailed_Map')
        if df.empty:
            return pd.DataFrame(columns=base_columns)

        grouped = df.groupby('Reaction_Type_Standardized').agg({
            'AREA_TOTAL_REDUCED': ['count', 'mean', lambda x: (x > 50).sum() / len(x) * 100],
            'Catalyst': ['nunique', self._mode_or_none],
            'Reactant_A_Type': lambda x: len(set(zip(x, df.loc[x.index, 'Reactant_B_Type'])))
        }).reset_index()
        
        # Get top reactant pair for each reaction
        def get_top_pair(reaction_type):
            sub_df = df[df['Reaction_Type_Standardized'] == reaction_type]
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

        if include_detailed_map:
            grouped['Detailed_Map'] = grouped['Reaction_Type'].apply(
                lambda rxn_type: json.dumps(
                    self._build_reaction_type_detail_map(
                        str(rxn_type),
                        df[df['Reaction_Type_Standardized'] == rxn_type],
                        detail_top_k=detail_top_k,
                    ),
                    ensure_ascii=True,
                    separators=(",", ":"),
                )
            )
        
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
        motif_sets = _load_motif_sets()
        query_a = _normalize_query_tokens(reactant_a_type, motif_sets)
        query_b = _normalize_query_tokens(reactant_b_type, motif_sets)
        query_df = self.df[
            self.df['Reactant_A_Type'].apply(lambda value: _field_matches_query(value, query_a, motif_sets))
            & self.df['Reactant_B_Type'].apply(lambda value: _field_matches_query(value, query_b, motif_sets))
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
            ~(
                similar_df['Reactant_A_Type'].apply(lambda value: _field_matches_query(value, query_a, motif_sets))
                & similar_df['Reactant_B_Type'].apply(lambda value: _field_matches_query(value, query_b, motif_sets))
            )
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
        motif_sets = _load_motif_sets()
        
        # Apply filters
        if reaction_type:
            df = df[df['Reaction_Type_Standardized'].str.contains(reaction_type, case=False, na=False)]
        
        if catalyst_filter:
            df = self._filter_by_catalyst_type(df, catalyst_filter)
        
        if reactant_a_type:
            query_a = _normalize_query_tokens(reactant_a_type, motif_sets)
            df = df[df['Reactant_A_Type'].apply(lambda value: _field_matches_query(value, query_a, motif_sets))]
        
        if reactant_b_type:
            query_b = _normalize_query_tokens(reactant_b_type, motif_sets)
            df = df[df['Reactant_B_Type'].apply(lambda value: _field_matches_query(value, query_b, motif_sets))]
        
        if min_yield is not None:
            df = df[df['AREA_TOTAL_REDUCED'] >= min_yield]
        
        # Export
        df.to_csv(output_path, index=False)
        print(f"✅ Exported {len(df):,} experiments to {output_path}")
        
        return len(df)
