"""
ChemTools master class - Unified context for all chemistry operations.

This module provides a clean, object-oriented API for accessing all chemtools
functionality with built-in resource management, lazy loading, and caching.

Usage:
    # Use global singleton (recommended for most cases)
    from chemtools import chem
    result = chem.smiles.normalize("CCO")
    
    # Or create custom instance with specific configuration
    from chemtools import ChemTools
    my_chem = ChemTools(datasets=["C_N_Coupling_Pd"], preload=True)
    precedents = my_chem.precedent.knn(...)

Architecture:
    ChemTools (master class)
    ├── Core Operations (stateless, always available)
    �?  ├── smiles - SMILES parsing and normalization
    �?  ├── router - Reaction family detection
    �?  ├── properties - Compound property lookup
    �?  ├── constraints - Constraint validation
    �?  └── functional_groups - Functional group detection
    �?
    ├── Data Operations (stateful, lazy-loaded through context)
    �?  ├── precedent - Precedent reaction search
    �?  ├── recommend - HTE-based condition recommendations
    �?  └── explain - Explanation generation
    �?
    └── Advanced Operations (optional dependencies)
        ├── featurizers - Molecular featurization
        ├── features - Role-aware features (if available)
        └── integrations - External integrations (MCP, etc.)
"""

from __future__ import annotations
from typing import Any, Dict, List, Optional, Set, Tuple
from dataclasses import dataclass, field
from threading import RLock
from pathlib import Path
import os


# ============================================================================
# Configuration
# ============================================================================

@dataclass
class ResourceConfig:
    """Configuration for ChemTools resource management.
    
    Attributes:
        datasets: List of reaction family datasets to load (e.g., ["C_N_Coupling_Pd"])
                 If None, datasets are loaded on-demand
        ml_models: List of ML models to preload (e.g., ["buchwald", "ullmann"])
                  If None, models are loaded on-demand
        reagent_dbs: List of reagent databases to preload (e.g., ["ligand", "base"])
                    If None, databases are loaded on-demand
        preload: If True, eagerly load configured resources at initialization
        cache_size: Maximum number of cached resources (LRU eviction)
        data_dir: Override data directory path (default: auto-detect)
        enable_rdkit: Enable RDKit operations (default: True)
    """
    datasets: Optional[List[str]] = None
    ml_models: Optional[List[str]] = None
    reagent_dbs: Optional[List[str]] = None
    preload: bool = False
    cache_size: int = 32
    data_dir: Optional[Path] = None
    enable_rdkit: bool = True
    
    def __post_init__(self):
        """Normalize configuration values."""
        if self.datasets is not None:
            self.datasets = [str(d) for d in self.datasets]
        if self.ml_models is not None:
            self.ml_models = [str(m) for m in self.ml_models]
        if self.reagent_dbs is not None:
            self.reagent_dbs = [str(r) for r in self.reagent_dbs]
        if self.data_dir is not None:
            self.data_dir = Path(self.data_dir)


# ============================================================================
# Namespace Wrappers - Core Operations (Stateless)
# ============================================================================

class SmilesNamespace:
    """SMILES parsing, normalization, and validation operations.
    
    These are stateless operations that don't require resource loading.
    """
    
    @staticmethod
    def normalize(smi: str) -> Dict[str, Any]:
        """Normalize a SMILES string.
        
        Args:
            smi: Input SMILES string
            
        Returns:
            Dict with keys: input, fragments, largest_smiles, smiles_norm
        """
        from . import smiles as _smiles
        return _smiles.normalize(smi)
    
    @staticmethod
    def normalize_reaction(rsmi: str) -> Dict[str, Any]:
        """Normalize a reaction SMILES string.
        
        Args:
            rsmi: Reaction SMILES (reactants>agents>products)
            
        Returns:
            Dict with normalized reactants, agents, products
        """
        from . import smiles as _smiles
        return _smiles.normalize_reaction(rsmi)


class RouterNamespace:
    """Reaction family routing and detection.
    
    Stateless operations for determining reaction types.
    """
    
    @staticmethod
    def detect_family(reaction: str, **kwargs) -> Dict[str, Any]:
        """Detect reaction family from SMILES.
        
        Args:
            reaction: Reaction SMILES string
            **kwargs: Additional detection parameters
            
        Returns:
            Dict with detected family and confidence
        """
        from . import router as _router
        return _router.detect_family_from_reaction(reaction, **kwargs)


class PropertiesNamespace:
    """Compound property lookup and calculation.
    
    Stateless operations for property queries.
    """
    
    @staticmethod
    def lookup(query: str, **kwargs) -> Dict[str, Any]:
        """Look up compound properties.
        
        Args:
            query: Compound identifier (name, CAS, SMILES)
            **kwargs: Additional lookup parameters
            
        Returns:
            Dict with compound properties
        """
        from . import properties as _properties
        return _properties.lookup(query, **kwargs)


class ConstraintsNamespace:
    """Constraint validation and filtering.
    
    Stateless operations for applying constraints.
    """
    
    @staticmethod
    def filter(candidates: List[Dict], constraints: Dict, **kwargs) -> List[Dict]:
        """Filter candidates by constraints.
        
        Args:
            candidates: List of candidate conditions
            constraints: Constraint rules to apply
            **kwargs: Additional filter parameters
            
        Returns:
            Filtered list of candidates
        """
        from . import constraints as _constraints
        return _constraints.filter(candidates, constraints, **kwargs)


class FunctionalGroupsNamespace:
    """Functional group detection and analysis.
    
    Stateless operations for detecting functional groups in molecules.
    Uses RDKit SMARTS matching when available, falls back to text patterns.
    """
    
    @staticmethod
    def detect(smiles: Optional[str]) -> Dict[str, bool]:
        """Detect all functional groups in a molecule.
        
        Args:
            smiles: SMILES string to analyze
            
        Returns:
            Dict mapping functional group name to presence (True/False)
            
        Example:
            >>> chem.functional_groups.detect("CC(=O)O")
            {'carboxylic_acid': True, 'carbonyl': True, ...}
        """
        from .util import functional_groups as _fg
        return _fg.detect_all(smiles)
    
    @staticmethod
    def get_groups(smiles: Optional[str]) -> List[str]:
        """Get list of functional groups present in a molecule.
        
        Args:
            smiles: SMILES string to analyze
            
        Returns:
            List of functional group names
            
        Example:
            >>> chem.functional_groups.get_groups("c1ccc(Br)cc1N")
            ['aniline', 'amine_primary', 'aryl_bromide', 'aromatic']
        """
        from .util import functional_groups as _fg
        return _fg.get_functional_groups(smiles)
    
    @staticmethod
    def has(smiles: Optional[str], group_name: str) -> bool:
        """Check if molecule has a specific functional group.
        
        Args:
            smiles: SMILES string to analyze
            group_name: Name of functional group to check
            
        Returns:
            True if functional group is present
            
        Example:
            >>> chem.functional_groups.has("CC(=O)O", "carboxylic_acid")
            True
        """
        from .util import functional_groups as _fg
        return _fg.has_functional_group(smiles, group_name)
    
    @staticmethod
    def count(smiles: Optional[str], group_name: str) -> int:
        """Count occurrences of a specific functional group.
        
        Args:
            smiles: SMILES string to analyze
            group_name: Name of functional group to count
            
        Returns:
            Number of times the functional group appears
            
        Example:
            >>> chem.functional_groups.count("O=C(O)CC(=O)O", "carboxylic_acid")
            2
        """
        from .util import functional_groups as _fg
        return _fg.count_functional_groups(smiles, group_name)
    
    @staticmethod
    def categorize(smiles: Optional[str]) -> Dict[str, List[str]]:
        """Organize detected functional groups by chemical category.
        
        Args:
            smiles: SMILES string to analyze
            
        Returns:
            Dict mapping category to list of functional groups
            
        Example:
            >>> chem.functional_groups.categorize("CC(=O)Oc1ccccc1")
            {'oxygen': ['ester', 'carbonyl', 'phenol'], 'aromatic': [...]}
        """
        from .util import functional_groups as _fg
        return _fg.get_group_categories(smiles)
    
    @staticmethod
    def summarize(smiles: Optional[str]) -> str:
        """Get human-readable summary of functional groups.
        
        Args:
            smiles: SMILES string to analyze
            
        Returns:
            Formatted string summary
            
        Example:
            >>> print(chem.functional_groups.summarize("CC(=O)Oc1ccccc1"))
            Oxygen: ester, carbonyl, phenol
            Aromatic: aromatic, phenol
        """
        from .util import functional_groups as _fg
        return _fg.summarize_functional_groups(smiles)
    
    @staticmethod
    def list_available() -> List[str]:
        """Get list of all detectable functional group names.
        
        Returns:
            List of functional group names that can be detected
        """
        from .util import functional_groups as _fg
        return sorted(_fg.FUNCTIONAL_GROUP_SMARTS.keys())


# ============================================================================
# Namespace Wrappers - Data Operations (Stateful, Context-Aware)
# ============================================================================

class PrecedentNamespace:
    """Precedent reaction search and analysis.
    
    Stateful operations that use cached datasets through context.
    """
    
    def __init__(self, context: 'ChemTools'):
        self._ctx = context
    
    def knn(self, family: str, features: Dict, k: int = 5, relax: Optional[Dict] = None, **kwargs) -> List[Dict]:
        """K-nearest neighbor precedent search.
        
        Args:
            family: Reaction family name
            features: Feature vector for similarity search
            k: Number of neighbors to return
            relax: Relaxation rules for constraints
            **kwargs: Additional search parameters
            
        Returns:
            List of k nearest precedent reactions
        """
        # Ensure dataset is loaded
        dataset = self._ctx.get_reaction_dataset(family)
        
        from . import precedent as _precedent
        return _precedent.knn(family, features, k, relax or {}, **kwargs)
    
    def find_reactions_by_core(self, core: str, family: Optional[str] = None, 
                              fuzzy: bool = False, limit: int = 50) -> List[Dict]:
        """Find reactions by core structure.
        
        Args:
            core: Core structure SMILES
            family: Optional reaction family filter
            fuzzy: Enable fuzzy matching
            limit: Maximum results to return
            
        Returns:
            List of matching reactions
        """
        from . import precedent as _precedent
        return _precedent.find_reactions_by_core(core, family, fuzzy, limit)
    
    def list_cores(self, family: Optional[str] = None, top_n: int = 200, 
                   include_counts: bool = True) -> List[Dict]:
        """List available core structures.
        
        Args:
            family: Optional reaction family filter
            top_n: Maximum cores to return
            include_counts: Include reaction counts per core
            
        Returns:
            List of core structures
        """
        from . import precedent as _precedent
        return _precedent.list_cores(family, top_n, include_counts)


class RecommendNamespace:
    """HTE-backed condition recommendations.

    Stateful operations that use the HTE recommender as the primary backend.
    """
    
    def __init__(self, context: 'ChemTools'):
        self._ctx = context
    
    def conditions(self, reaction: str, reaction_type: Optional[str] = None,
                   k: int = 5, limit: int = 10, relax: Optional[Dict] = None,
                   constraints: Optional[Dict] = None, 
                   rerank_strategy: str = 'analytics',
                   filter_unknown_reagents: bool = False,
                   search_all_families: bool = False,
                   **kwargs) -> Dict[str, Any]:
        """Get unified condition recommendations.
        
        Args:
            reaction: Reaction SMILES
            reaction_type: Optional reaction type hint (ignored if search_all_families=True)
            k: Number of precedents for each recommendation
            limit: Maximum recommendations to return
            relax: Relaxation rules
            constraints: Constraint rules
            rerank_strategy: Legacy parameter (ignored by unified recommender)
            filter_unknown_reagents: Legacy parameter (ignored by unified recommender)
            search_all_families: Legacy parameter (ignored by unified recommender)
            **kwargs: Additional parameters
            
        Returns:
            Dict with recommendations and metadata
        """
        from . import recommend as _recommend
        return _recommend.recommend_conditions_structured(
            reaction=reaction,
            reaction_type=reaction_type,
            k=k,
            limit=limit,
            relax=relax or {},
            constraints=constraints or {},
            rerank_strategy=rerank_strategy,
            filter_unknown_reagents=filter_unknown_reagents,
            search_all_families=search_all_families,
        )
    
    def from_reaction(self, reaction: str, k: int = 5, relax: Optional[Dict] = None,
                     constraint_rules: Optional[Dict] = None) -> Dict[str, Any]:
        """Get unified recommendations from reaction SMILES.
        
        Args:
            reaction: Reaction SMILES
            k: Number of precedents per recommendation
            relax: Legacy parameter (ignored by unified recommender)
            constraint_rules: Legacy parameter (ignored by unified recommender)
            
        Returns:
            Dict with recommendations
        """
        from . import recommend as _recommend
        return _recommend.recommend_from_reaction(reaction, k, relax or {}, constraint_rules or {})
    
    def design_plate(self, reaction: str, plate_size: int = 96,
                    relax: Optional[Dict] = None, constraint_rules: Optional[Dict] = None) -> Dict[str, Any]:
        """Design experimental plate layout (not available in unified recommender)."""
        raise NotImplementedError(
            "Plate design is not available in the unified recommendation system."
        )


class ReagentNamespace:
    """Reagent database lookup and enrichment.
    
    Stateless operations for finding and enriching reagent information
    from the data/reagent_db/reagents.csv registry.
    """
    
    @staticmethod
    def find(name: str, reagent_type: str) -> Optional[Dict[str, Any]]:
        """Find a reagent by name in the specified database.
        
        Args:
            name: Reagent name, CAS number, or abbreviation
            reagent_type: Database type ('ligand', 'base', 'solvent', 'metal_catalyst', etc.)
            
        Returns:
            Reagent dictionary if found, None otherwise
        """
        from . import reagent as _reagent
        return _reagent.find_reagent(name, reagent_type)
    
    @staticmethod
    def enrich(name: str, reagent_type: str) -> Dict[str, Any]:
        """Enrich a reagent name with detailed information.
        
        Args:
            name: Reagent name
            reagent_type: Database type
            
        Returns:
            Dictionary with enriched reagent details
        """
        from . import reagent as _reagent
        return _reagent.enrich_reagent_info(name, reagent_type)
    
    @staticmethod
    def enrich_conditions(conditions: Dict[str, Any]) -> Dict[str, Any]:
        """Enrich full condition set with reagent details.
        
        Args:
            conditions: Conditions dictionary from recommendation
            
        Returns:
            Enriched conditions with detailed reagent information
        """
        from . import reagent as _reagent
        return _reagent.enrich_conditions(conditions)
    
    @staticmethod
    def list_types() -> List[str]:
        """Get list of available reagent database types.
        
        Returns:
            List of reagent type names
        """
        from . import reagent as _reagent
        return _reagent.get_all_reagent_types()


class ExplainNamespace:
    """Explanation generation for recommendations.
    
    Stateful operations that generate human-readable explanations.
    """
    
    def __init__(self, context: 'ChemTools'):
        self._ctx = context
    
    def precedents(self, precedents: List[Dict], **kwargs) -> str:
        """Generate explanation from precedents.
        
        Args:
            precedents: List of precedent reactions
            **kwargs: Additional formatting parameters
            
        Returns:
            Human-readable explanation text
        """
        from . import explain as _explain
        return _explain.explain_precedents(precedents, **kwargs)


# ============================================================================
# Namespace Wrappers - Advanced Operations (Optional Dependencies)
# ============================================================================

class FeaturizersNamespace:
    """Molecular and reaction featurization.
    
    Provides access to various featurization methods.
    """
    
    def __init__(self, context: 'ChemTools'):
        self._ctx = context
    
    def molecular(self, **kwargs) -> Any:
        """Access pair featurizer for electrophile/nucleophile substrates."""
        from . import featurizers as _featurizers
        return _featurizers.reaction_pair

    def structural(self, **kwargs) -> Any:
        """Access motif-based steric/electronic analysis for single molecules."""
        from . import featurizers as _featurizers
        return _featurizers.structural

    def reaction(self, reaction_smiles: str, **kwargs) -> Dict[str, Any]:
        """Systematically analyze a reaction: detect type and featurize all reactants."""
        from . import featurizers as _featurizers
        return _featurizers.structural.featurize_reaction(reaction_smiles, **kwargs)


class FeaturesNamespace:
    """Advanced feature extraction (optional dependency).
    
    Role-aware and specialized feature extraction methods.
    """
    
    def __init__(self, context: 'ChemTools'):
        self._ctx = context
        self._available = None
    
    def is_available(self) -> bool:
        """Check if advanced features are available."""
        if self._available is None:
            try:
                from .features import role
                self._available = True
            except ImportError:
                self._available = False
        return self._available
    
    @property
    def role(self):
        """Access role-aware features (if available)."""
        if not self.is_available():
            raise ImportError("Role-aware features not available. Install optional dependencies.")
        from .features import role
        return role


class IntegrationsNamespace:
    """External integrations (MCP, etc.).
    
    Provides access to integration modules.
    """
    
    def __init__(self, context: 'ChemTools'):
        self._ctx = context
    
    @property
    def mcp(self):
        """Access Model Context Protocol integration."""
        from .integrations import mcp
        return mcp


class DatasetAnalyticsNamespace:
    """HTE dataset analytics (CSV/JSONL).

    This replaces the legacy dataset_analytics module with the newer
    chemtools.recommend.analytics.HTEAnalytics backend.
    """

    def __init__(self, hte_db_path: str = "data/HTE_db") -> None:
        self._hte_db_path = hte_db_path
        self._analytics = None

    def _get_analytics(self):
        if self._analytics is None:
            from .recommend.analytics import HTEAnalytics
            self._analytics = HTEAnalytics(self._hte_db_path)
        return self._analytics

    def get_stats(self, reaction_type: str) -> Dict[str, Any]:
        analytics = self._get_analytics()
        summary = analytics.get_reaction_type_summary()
        if summary.empty:
            return {"total_reactions": 0}
        row = summary[summary["Reaction_Type"].str.contains(reaction_type, case=False, na=False)]
        if row.empty:
            return {"total_reactions": 0}
        data = row.iloc[0].to_dict()
        return {
            "total_reactions": int(data.get("Num_Experiments", 0) or 0),
            "avg_yield": float(data.get("Avg_Yield", 0.0) or 0.0),
            "success_rate": float(data.get("Success_Rate", 0.0) or 0.0),
            "num_catalysts": int(data.get("Num_Catalysts", 0) or 0),
            "num_reactant_pairs": int(data.get("Num_Reactant_Pairs", 0) or 0),
            "top_catalyst": data.get("Top_Catalyst") or "",
            "top_reactant_pair": data.get("Top_Reactant_Pair") or "",
        }

    def list_reactant_pairs(
        self,
        reaction_type: Optional[str] = None,
        catalyst_filter: Optional[str] = None,
        min_experiments: int = 1,
        sort_by: str = "count",
    ):
        return self._get_analytics().list_reactant_pairs(
            reaction_type=reaction_type,
            catalyst_filter=catalyst_filter,
            min_experiments=min_experiments,
            sort_by=sort_by,
        )

    def get_catalyst_stats(
        self,
        reaction_type: Optional[str] = None,
        reactant_a_type: Optional[str] = None,
        reactant_b_type: Optional[str] = None,
    ):
        return self._get_analytics().get_catalyst_stats(
            reaction_type=reaction_type,
            reactant_a_type=reactant_a_type,
            reactant_b_type=reactant_b_type,
        )

    def get_reaction_type_summary(self):
        return self._get_analytics().get_reaction_type_summary()

    def analyze_metal_usage(self) -> Dict[str, Any]:
        return self._get_analytics().analyze_metal_usage()

    def export_subset(
        self,
        output_path: str,
        reaction_type: Optional[str] = None,
        catalyst_filter: Optional[str] = None,
        reactant_a_type: Optional[str] = None,
        reactant_b_type: Optional[str] = None,
        min_yield: Optional[float] = None,
    ) -> int:
        return self._get_analytics().export_subset(
            output_path=output_path,
            reaction_type=reaction_type,
            catalyst_filter=catalyst_filter,
            reactant_a_type=reactant_a_type,
            reactant_b_type=reactant_b_type,
            min_yield=min_yield,
        )

    def get_all_families(self) -> List[str]:
        analytics = self._get_analytics()
        df = analytics.df
        if df is None or df.empty:
            return []
        return sorted({str(v) for v in df["Reaction_Type_Standardized"].dropna().unique() if str(v).strip()})

    def print_summary(self, reaction_type: str) -> None:
        stats = self.get_stats(reaction_type)
        print(f"HTE Analytics Summary: {reaction_type}")
        for key, value in stats.items():
            print(f"  {key}: {value}")


# ============================================================================
# Main ChemTools Class
# ============================================================================

class ChemTools:
    """Master class for all ChemTools operations.
    
    Provides unified access to all chemistry tools with resource management,
    lazy loading, and caching.
    
    Example:
        >>> from chemtools import ChemTools
        >>> chem = ChemTools(datasets=["C_N_Coupling_Pd"])
        >>> result = chem.smiles.normalize("CCO")
        >>> precedents = chem.precedent.knn("C_N_Coupling_Pd", features={...})
    """
    
    def __init__(self, 
                 datasets: Optional[List[str]] = None,
                 ml_models: Optional[List[str]] = None,
                 reagent_dbs: Optional[List[str]] = None,
                 preload: bool = False,
                 cache_size: int = 32,
                 data_dir: Optional[Path] = None,
                 enable_rdkit: bool = True,
                 config: Optional[ResourceConfig] = None):
        """Initialize ChemTools context.
        
        Args:
            datasets: List of reaction datasets to load
            ml_models: List of ML models to preload
            reagent_dbs: List of reagent databases to preload
            preload: If True, eagerly load resources at init
            cache_size: Maximum cached resources
            data_dir: Override data directory
            enable_rdkit: Enable RDKit operations
            config: ResourceConfig object (overrides other params)
        """
        # Use provided config or create from parameters
        if config is not None:
            self.config = config
        else:
            self.config = ResourceConfig(
                datasets=datasets,
                ml_models=ml_models,
                reagent_dbs=reagent_dbs,
                preload=preload,
                cache_size=cache_size,
                data_dir=data_dir,
                enable_rdkit=enable_rdkit,
            )
        
        # Thread-safe resource caches
        self._lock = RLock()
        self._datasets: Dict[str, Any] = {}
        self._ml_models: Dict[str, Any] = {}
        self._reagent_dbs: Dict[str, Any] = {}
        self._loaded_families: Set[str] = set()
        
        # Initialize namespace wrappers
        # Core operations (stateless)
        self.smiles = SmilesNamespace()
        self.router = RouterNamespace()
        self.properties = PropertiesNamespace()
        self.constraints = ConstraintsNamespace()
        self.functional_groups = FunctionalGroupsNamespace()  # Functional group detection
        self.reagent = ReagentNamespace()  # Reagent database access
        self.analytics = DatasetAnalyticsNamespace()  # Dataset analytics
        
        # Data operations (stateful, context-aware)
        self.precedent = PrecedentNamespace(self)
        self.recommend = RecommendNamespace(self)
        self.explain = ExplainNamespace(self)
        
        # Advanced operations (optional)
        self.featurizers = FeaturizersNamespace(self)
        self.features = FeaturesNamespace(self)
        self.integrations = IntegrationsNamespace(self)
        
        # Preload resources if requested
        if self.config.preload:
            self._preload_resources()
    
    # ------------------------------------------------------------------------
    # Resource Management - Datasets
    # ------------------------------------------------------------------------
    
    def get_reaction_dataset(self, family: str) -> List[Dict[str, Any]]:
        """Get reaction dataset with lazy loading and caching.
        
        Args:
            family: Reaction family name (e.g., "C_N_Coupling_Pd")
            
        Returns:
            List of reaction dictionaries
        """
        with self._lock:
            # Check if already loaded
            if family in self._datasets:
                return self._datasets[family]
            
            # Check if this family is in allowed list
            if self.config.datasets is not None and family not in self.config.datasets:
                # Load anyway but warn
                import warnings
                warnings.warn(
                    f"Dataset '{family}' not in configured datasets {self.config.datasets}. "
                    f"Loading anyway, but this may be slower."
                )
            
            # Load dataset through precedent module
            from . import precedent as _precedent
            dataset = _precedent._load_selective([family])
            
            # Cache the dataset
            self._datasets[family] = dataset
            self._loaded_families.add(family)
            
            return dataset
    
    def list_loaded_datasets(self) -> List[str]:
        """Get list of currently loaded datasets."""
        with self._lock:
            return sorted(self._loaded_families)
    
    def unload_dataset(self, family: str) -> bool:
        """Unload a dataset to free memory.
        
        Args:
            family: Reaction family name
            
        Returns:
            True if dataset was unloaded, False if not loaded
        """
        with self._lock:
            if family in self._datasets:
                del self._datasets[family]
                self._loaded_families.discard(family)
                return True
            return False
    
    # ------------------------------------------------------------------------
    # Resource Management - ML Models
    # ------------------------------------------------------------------------
    
    def get_ml_model(self, model_name: str) -> Any:
        """Get ML model with lazy loading and caching.
        
        Args:
            model_name: Model identifier (e.g., "buchwald", "ullmann")
            
        Returns:
            Loaded ML model
        """
        with self._lock:
            if model_name in self._ml_models:
                return self._ml_models[model_name]
            
            # Check if this model is in allowed list
            if self.config.ml_models is not None and model_name not in self.config.ml_models:
                import warnings
                warnings.warn(
                    f"ML model '{model_name}' not in configured models {self.config.ml_models}. "
                    f"Loading anyway, but this may be slower."
                )
            
            # Load model (implementation depends on model type)
            # For now, return None and let the caller handle loading
            # This will be implemented in Phase 2
            model = None
            
            # Cache the model
            self._ml_models[model_name] = model
            
            return model
    
    def list_loaded_models(self) -> List[str]:
        """Get list of currently loaded ML models."""
        with self._lock:
            return sorted(self._ml_models.keys())
    
    # ------------------------------------------------------------------------
    # Resource Management - Reagent Databases
    # ------------------------------------------------------------------------
    
    def get_reagent_database(self, reagent_type: str) -> List[Dict[str, Any]]:
        """Get reagent database with lazy loading and caching.
        
        Args:
            reagent_type: Type of reagent (e.g., "ligand", "base", "solvent")
            
        Returns:
            List of reagent dictionaries
        """
        with self._lock:
            if reagent_type in self._reagent_dbs:
                return self._reagent_dbs[reagent_type]
            
            # Load through reagent_lookup module
            from . import reagent as _reagent_lookup
            db = _reagent_lookup.load_reagent_database(reagent_type)
            
            # Cache the database
            self._reagent_dbs[reagent_type] = db
            
            return db
    
    # ------------------------------------------------------------------------
    # Preloading and Cache Management
    # ------------------------------------------------------------------------
    
    def _preload_resources(self) -> None:
        """Preload configured resources."""
        # Preload datasets
        if self.config.datasets:
            for family in self.config.datasets:
                try:
                    self.get_reaction_dataset(family)
                except Exception as e:
                    import warnings
                    warnings.warn(f"Failed to preload dataset '{family}': {e}")
        
        # Preload ML models
        if self.config.ml_models:
            for model in self.config.ml_models:
                try:
                    self.get_ml_model(model)
                except Exception as e:
                    import warnings
                    warnings.warn(f"Failed to preload model '{model}': {e}")
        
        # Preload reagent databases
        if self.config.reagent_dbs:
            for reagent_type in self.config.reagent_dbs:
                try:
                    self.get_reagent_database(reagent_type)
                except Exception as e:
                    import warnings
                    warnings.warn(f"Failed to preload reagent DB '{reagent_type}': {e}")
    
    def clear_cache(self) -> None:
        """Clear all cached resources."""
        with self._lock:
            self._datasets.clear()
            self._ml_models.clear()
            self._reagent_dbs.clear()
            self._loaded_families.clear()
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics.
        
        Returns:
            Dict with cache size and loaded resources
        """
        with self._lock:
            return {
                "datasets_loaded": len(self._datasets),
                "dataset_families": sorted(self._loaded_families),
                "ml_models_loaded": len(self._ml_models),
                "ml_model_names": sorted(self._ml_models.keys()),
                "reagent_dbs_loaded": len(self._reagent_dbs),
                "reagent_db_types": sorted(self._reagent_dbs.keys()),
                "total_resources": len(self._datasets) + len(self._ml_models) + len(self._reagent_dbs),
            }
    
    # ------------------------------------------------------------------------
    # Utility Methods
    # ------------------------------------------------------------------------
    
    def __repr__(self) -> str:
        stats = self.get_cache_stats()
        return (
            f"ChemTools(datasets={stats['datasets_loaded']}, "
            f"models={stats['ml_models_loaded']}, "
            f"reagent_dbs={stats['reagent_dbs_loaded']})"
        )


# ============================================================================
# Global Singleton Instance
# ============================================================================

# Default configuration: minimal preloading, load on-demand
_default_config = ResourceConfig(
    datasets=None,  # Load on-demand
    ml_models=None,  # Load on-demand
    reagent_dbs=None,  # Load on-demand
    preload=False,  # Don't preload anything
    cache_size=32,
    enable_rdkit=True,
)

# Global singleton instance for convenience
chem = ChemTools(config=_default_config)

# Expose for "from chemtools import ChemTools" and "from chemtools import chem"
__all__ = ['ChemTools', 'chem', 'ResourceConfig']
