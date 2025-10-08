# ChemTools Directory Structure & File Summary

**Generated**: October 7, 2025  
**Purpose**: Comprehensive overview of all modules in the `/chemtools` directory

---

## 📁 Directory Structure

```
chemtools/
├── __init__.py                   # Package exports & ChemTools v2.0 API
├── context.py                    # ChemTools master class & namespace wrappers
├── contracts.py                  # Pydantic request/response models for API
│
├── Core Operations (Stateless)
├── smiles.py                     # SMILES normalization & reaction parsing
├── router.py                     # Reaction family detection (SMARTS-based)
├── properties.py                 # Compound property/role lookup
├── constraints.py                # Constraint filters (inventory, environmental)
│
├── Data Operations (Stateful)
├── precedent/                    # ✨ REFACTORED: k-NN precedent search (was 846 lines)
│   ├── __init__.py               # Public API exports
│   ├── search.py                 # Main k-NN search logic
│   ├── loader.py                 # Dataset loading and transformation
│   ├── similarity.py             # Feature similarity calculations
│   ├── catalyst.py               # Catalyst class detection/matching
│   ├── core_utils.py             # Family normalization utilities
│   └── integrations.py           # MolPipeline integration
├── recommend.py                  # ML-based condition recommendations (1454 lines) ⚠️ **NEXT REFACTORING TARGET**
├── recommend_ml.py               # Hybrid ML + k-NN recommender (⚠️ refactoring proposed)
├── reagent_lookup.py             # Reagent database enrichment
├── explain.py                    # Human-readable explanations
│
├── Utilities & Helpers
├── reaction_similarity.py        # DRFP encoding & Tanimoto similarity
├── reaction_type_detector.py     # rxn-insight integration
├── output_formatter.py           # Structured output formatting
├── condition_core.py             # Legacy condition parsing (may be deprecated)
├── selector_payloads.py          # Rule feature builders for amide formation
│
├── Subdirectories
├── featurizers/                  # Molecular featurization
│   ├── molecular.py              # Generic molecular features
│   └── ullmann.py                # Ullmann C-N coupling features
├── util/                         # Utility modules
│   └── rdkit_helpers.py          # RDKit wrapper functions
├── ml/                           # Machine learning models (⚠️ refactoring proposed)
│   ├── drfp_predictor.py         # DRFP-based yield prediction
│   └── evaluation.py             # Model evaluation utilities
├── features/                     # Advanced featurization
│   └── role/                     # Role-aware molecular features
├── integrations/                 # External integrations
│   ├── mcp/                      # Model Context Protocol
│   └── molpipeline.py            # MolPipeline integration
├── rule_scdb_matcher/            # Rule-based scheme matching (727+ lines)
│   ├── matcher.py                # SMARTS-based matcher
│   ├── loader.py                 # Database loader
│   ├── types.py                  # Type definitions
│   ├── ecn.py                    # Essential core normalization
│   └── cli.py                    # Command-line interface
└── agent/                        # Agent-based workflows
    ├── config.py                 # Agent configuration
    └── features/                 # Agent feature extraction
```

---

## 📄 Core Files (Root Level)

### `__init__.py` (57 lines)
**Purpose**: Package initialization and exports for ChemTools v2.0

**Exports**:
- `chem` - Global singleton instance (primary interface)
- `ChemTools` - Master class for custom instances
- `ResourceConfig` - Configuration dataclass
- `scdb_matcher` - Legacy alias for backward compatibility

**Version**: 2.0.0

**Key Features**:
- Unified API with lazy loading
- 50-100x performance improvement
- Thread-safe resource management
- Backward compatibility layer

---

### `context.py` (887 lines) ⭐
**Purpose**: Master ChemTools class and all namespace wrappers

**Classes** (14 total):

#### Configuration
- `ResourceConfig` - Configuration for datasets, ML models, caching

#### Namespace Wrappers (Stateless)
- `SmilesNamespace` - SMILES operations (`normalize`, `normalize_reaction`)
- `RouterNamespace` - Reaction family detection (`detect_family`, `detect_family_from_reaction`)
- `PropertiesNamespace` - Property lookup (`lookup`)
- `ConstraintsNamespace` - Constraint filtering (`filter`)
- `ReagentNamespace` - Reagent database access (`find`, `enrich`, `list_types`)
- `RuleNamespace` - **NEW!** Rule-based matching (`load_database`, `match`, `list_databases`, `clear_cache`)

#### Namespace Wrappers (Stateful)
- `PrecedentNamespace` - Precedent search (`knn`, `list_cores`, `find_reactions_by_core`)
- `RecommendNamespace` - ML recommendations (`hybrid_recommend`, `conditions`, `plate_design`)
- `ExplainNamespace` - Explanation generation (`precedents`)

#### Namespace Wrappers (Advanced/Optional)
- `FeaturizersNamespace` - Featurization access (`ullmann`, `molecular`)
- `FeaturesNamespace` - Role-aware features (`mol`, `reaction`)
- `IntegrationsNamespace` - External integrations (`molpipeline`)

#### Master Class
- `ChemTools` - Main context class with:
  - Resource management (datasets, ML models, reagent DBs)
  - Thread-safe caching
  - Lazy loading
  - Preload support

**Key Methods**:
- `get_reaction_dataset()` - Load reaction dataset with caching
- `_preload_resources()` - Eager loading for API servers

---

### `contracts.py` (70 lines)
**Purpose**: Pydantic models for API request/response validation

**Request Models** (13 total):
- `NormalizeRequest` - SMILES normalization
- `DetectFamilyRequest` - Family detection from reactants
- `DetectTypeRequest` - Reaction type detection
- `FeaturizeUllmannRequest` - Ullmann featurization
- `PropertiesLookupRequest` - Property lookup
- `PrecedentKNNRequest` - k-NN search
- `ConstraintsFilterRequest` - Constraint filtering
- `ExplainPrecedentsRequest` - Explanation generation
- `RecommendFromReactionRequest` - Basic recommendation
- `RecommendConditionsRequest` - Structured recommendation
- `PlateDesignRequest` - Plate design generation
- `SchemeMatchRequest` - Rule-based matching
- `CoreSearchRequest` - Core-based search
- `RoleAwareMolRequest` - Role-aware mol features
- `ConditionCoreValidateRequest` - Validation request

**Data Models**:
- `Reagent` - Reagent representation (uid, role, name, token)

---

## 🔧 Core Operations (Stateless)

### `smiles.py` (133 lines)
**Purpose**: SMILES string normalization and reaction parsing

**Key Functions**:
- `normalize(smi: str) -> Dict[str, Any]` - Normalize SMILES with RDKit fallback
  - Removes salts/counterions
  - Neutralizes charges
  - Canonicalizes structure
  - Picks largest organic fragment
  
- `normalize_reaction(rsmi: str) -> Dict[str, Any]` - Parse reaction SMILES
  - Splits reactants >> products >> agents
  - Normalizes each component
  - Returns structured dict

**Features**:
- Graceful RDKit degradation
- Carboxylate neutralization
- Fragment handling

---

### `router.py` (285 lines)
**Purpose**: Reaction family detection using SMARTS patterns and heuristics

**Key Functions**:
- `detect_family(reactants: List[str]) -> Dict[str, Any]` - Detect from reactant list
- `detect_family_from_reaction(reaction_smiles: str) -> Dict[str, Any]` - **Primary method**
  - Parses reaction SMILES
  - Applies SMARTS patterns
  - Detects metals from agents
  - Returns family with confidence

**Detected Families**:
- Suzuki C-C coupling
- C-N coupling (Pd/Cu/Ni variants)
- Sonogashira coupling
- Amide formation
- Esterification

**SMARTS Patterns**:
- Aryl/vinyl halides
- Boronic acids/esters
- Terminal alkynes
- N/O nucleophiles
- Carboxylic acids

**Metal Detection**:
- From agent formulas (Pd, Cu, Ni, Co)
- Atomic number matching
- Name token matching

---

### `properties.py` (317 lines)
**Purpose**: Compound property and role lookup from taxonomy

**Key Functions**:
- `lookup(query: str) -> Dict[str, Any]` - Lookup by CAS/name/token
  - Searches taxonomy files
  - Alias matching
  - External overrides

**Data Sources**:
- `data/compound_taxonomy/` - JSON taxonomy files
- Environment overrides:
  - `CHEMTOOLS_TAXONOMY_DIR` - Custom taxonomy location
  - `CHEMTOOLS_PROPERTIES_PATH` - External properties JSON

**Returned Properties**:
- `uid` - CAS number
- `role` - BASE, SOLVENT, CATALYST, LIGAND, etc.
- `token` - Shorthand notation
- `name` - Full chemical name
- `pKa_DMSO` / `pKa_water` - Acidity constants
- `KT` - Kamlet-Taft parameters (alpha, beta, pi*)

**Built-in Seed Data**:
- Common bases (K3PO4, KOH)
- Solvents (Water, etc.)
- Catalysts (CuI)
- Ligands (Phenanthroline)

---

### `constraints.py` (142 lines)
**Purpose**: Deterministic constraint filtering for reagents/solvents

**Key Functions**:
- `apply_filter(candidates: List[str], rules: Dict[str, Any]) -> Dict[str, Any]`
  - Returns `{allowed: [...], blocked: [{id, reason}]}`

**Constraint Types**:
- `inventory` - Whitelist (only allow these)
- `blacklist` - Explicit blocklist
- `no_HMPA` - Block hexamethylphosphoramide
- `no_chlorinated` - Block DCM, chloroform, etc.
- `aqueous_only` - Allow only water as solvent
- `min_bp_C` - Minimum boiling point threshold
- `allow_unknown` - Handle unknown properties gracefully

**Environmental Filters**:
- HMPA (CAS 680-31-9)
- Chlorinated solvents (DCM, CHCl3, chlorobenzene, etc.)

---

## 📊 Data Operations (Stateful)

### `precedent/` (Package) ✨ **REFACTORED**
**Purpose**: k-NN precedent search with DRFP similarity and selective loading

**Previous State**: Single 846-line file (`precedent.py`)  
**Current State**: Modular package with 6 focused modules (~150 lines each)

#### Module Breakdown

**`precedent/__init__.py`** - Public API
- Exports: `knn()`, `find_reactions_by_core()`, `list_cores()`
- Internal export: `_load_selective()` (for context.py)

**`precedent/search.py`** (~300 lines) - Main k-NN Logic
- `knn()` - Find k nearest neighbor precedents
- `find_reactions_by_core()` - Core-based reaction lookup
- `list_cores()` - List unique condition cores
- `_candidate_pool()` - Build candidate set with relaxation
- `_knn_impl()` - Core k-NN implementation
- `_knn_cached()` - LRU-cached wrapper

**`precedent/loader.py`** (~280 lines) - Dataset Loading
- `_load()` - Load all datasets (cached)
- `_load_selective()` - Load specific families (50-100x faster)
- `_make_row_from_dataset()` - Transform raw records
- `_dataset_family_map()` - Normalize family names
- `_iter_dataset_files()` - List JSONL files
- `_pick_electrophile_nucleophile()` - Reactant classification

**`precedent/similarity.py`** (~70 lines) - Feature Similarity
- `_similarity()` - Calculate feature similarity (0.0-1.0)
- Weighted categorical matching (LG, nuc_class, ortho, etc.)
- Exponential decay for numeric features (T_C, time_h)

**`precedent/catalyst.py`** (~120 lines) - Catalyst Classification
- `_row_catalyst_class()` - Classify by metal/enzyme/organo
- `_match_catalyst_class()` - Filter by catalyst class
- `_normalize_symbol()` - Metal name → symbol (e.g., "palladium" → "Pd")
- Metal detection from condition_core and full_system

**`precedent/core_utils.py`** (~75 lines) - Utilities
- `_family_text()` - Normalize family names
- `_proto_family_id()` - Generate prototype IDs
- `_parse_bin()` - Parse bin strings (e.g., "LG:Br|NUC:aniline")
- `_parse_core_tokens()` - Parse condition cores (e.g., "Pd/XPhos")
- `_norm_family()` - Normalize family with None handling

**`precedent/integrations.py`** (~90 lines) - MolPipeline
- `_attach_molpipeline_features()` - Add MolPipeline features to results
- Optional advanced featurization for precedents

#### Key Functions (Public API)

**`knn(family: str, features: Dict, k: int = 50) -> Dict`**
- **Main k-NN search** for precedent reactions
- DRFP-based similarity (Tanimoto)
- Categorical feature matching
- Catalyst class filtering
- Relaxation strategies
- Returns: `{prototype_id, support, precedents[]}`

**`find_reactions_by_core(core: str, family: str = None) -> List[Dict]`**
- Find reactions using similar condition cores
- Fuzzy matching on metal/ligand names
- Example: `find_reactions_by_core("Pd/XPhos", family="C_N_Coupling_Pd")`

**`list_cores(family: str = None, top_n: int = None) -> List`**
- List unique condition cores from dataset
- Optional frequency counts
- Filter by reaction family

#### Performance Optimizations
- **Selective loading**: Load only specified families (50-100x faster)
- **DRFP precomputation**: Load from NPZ files (8.5s → instant)
- **LRU caching**: Cache k-NN results (`@lru_cache`)
- **Lazy loading**: Load datasets on-demand

#### Data Sources
- `data/reaction_dataset/*.jsonl` - Per-family datasets
- `artifacts/all_drfp_4096.npz` - Precomputed DRFPs
- Environment: `CHEMTOOLS_DATASET_DIR`

#### Features
- Multi-dimensional similarity (DRFP + categorical)
- Catalyst class matching (Pd/Cu/Ni/enzyme/organo)
- Relaxation strategies (expand search when sparse)
- MolPipeline integration (optional)
- Thread-safe caching

---

### `recommend.py` (1454 lines) ⭐⭐⭐
**Purpose**: Comprehensive ML-based condition recommendation system

**Key Functions**:
- `recommend_from_reaction(reaction: str, k: int = 25) -> Dict[str, Any]` - **Main recommendation**
  - Normalizes reaction
  - Detects family
  - Featurizes substrates
  - Searches precedents
  - Generates recommendations
  
- `recommend_conditions_structured(reaction: str, k: int = 50) -> Dict[str, Any]` - Structured output
  - Multiple recommendation variants
  - Confidence scores
  - Precedent evidence
  - Enriched reagent info
  
- `design_plate_from_reaction(reaction: str, plate_size: int = 24) -> Dict[str, Any]` - Plate design
  - Generates experimental plate layouts
  - Well assignments
  - Condition variations

**Family-Specific Features**:
- **Amide formation**: Acid/amine classification, coupling agent selection
- **C-N coupling**: Electrophile/nucleophile analysis, metal/ligand pairing
- **Suzuki**: Boron/halide matching
- **Sonogashira**: Alkyne/halide pairing

**Rule Feature Builders**:
- `_amide_rule_feature_builder()` - Amide-specific features (538 lines)
  - Acid classification (aromatic, aliphatic, α-amino, etc.)
  - Amine classification (primary, secondary, aniline, etc.)
  - Functional group detection
  - Water tolerance inference
  
- `_default_rule_feature_builder()` - General features (78 lines)
  - Electrophile leaving group
  - Nucleophile class
  - Substituent analysis

**Substrate Analysis**:
- Carboxylic acid types (aromatic, aliphatic, α-amino, β-hydroxy, etc.)
- Amine types (primary, secondary, aniline, benzylic, α-amino, etc.)
- Functional groups (alcohol, phenol, sulfonamide, hydroxylamine)
- Steric/electronic properties

**Role-Aware Featurization** (optional):
- Integration with `chemtools.features.role`
- Mol-level feature extraction

---

### `recommend_ml.py` (226 lines)
**Purpose**: Hybrid ML + k-NN recommender with yield prediction

**Key Functions**:
- `hybrid_recommend(reaction_smiles: str, k: int = 10) -> Dict[str, Any]`
  - Uses k-NN to find precedents
  - Predicts yields with ML if n_precedents >= threshold
  - Re-ranks by predicted yield
  - Falls back to vote-based ranking if sparse
  
- `recommend_with_yield_prediction(reaction: str, precedents: List) -> Dict[str, Any]`
  - DRFP-based yield prediction
  - Confidence scoring
  - Variant re-ranking

**Strategy**:
- **High data**: ML yield prediction (n_precedents >= 50)
- **Low data**: k-NN vote-based system (n_precedents < 50)
- **Always**: Show precedents for interpretability

**ML Integration**:
- Uses `DRFPYieldPredictor` from `chemtools.ml.drfp_predictor`
- Default model: `models/drfp_yield_v1.pkl`

---

### `reagent_lookup.py` (303 lines)
**Purpose**: Reagent database lookup and enrichment

**Key Functions**:
- `find_reagent(name: str, reagent_type: str) -> Optional[Dict]` - Find by name/CAS
- `enrich_reagent_info(name: str, reagent_type: str) -> Dict` - Full enrichment
- `enrich_conditions(conditions: Dict) -> Dict` - Enrich full condition set
- `get_all_reagent_types() -> List[str]` - List available types
- `preload_common_databases()` - Eager loading

**Reagent Types** (13 databases in `data/reagents/`):
- `acid` - Carboxylic acids
- `additive` - Additives
- `base` - Bases
- `ligand` - Ligands
- `metal_precursor` - Metal catalysts
- `oxidant` - Oxidizing agents
- `reductant` - Reducing agents
- `solvent` - Solvents
- `reagent` - General reagents
- `catalyst` - Catalysts
- `nucleophile` - Nucleophiles
- `electrophile` - Electrophiles
- `coupling_reagent` - Coupling agents

**Features**:
- Name/CAS/alias matching
- Fuzzy name normalization
- LRU caching
- Multiple alias support

---

### `explain.py` (110 lines)
**Purpose**: Generate human-readable explanations from recommendations

**Key Functions**:
- `for_pack(pack: Dict, features: Dict) -> Dict[str, Any]` - Generate explanation
  - Bin description
  - Top reagents with vote counts
  - Numeric parameter summaries
  - UID → human-readable names

**Output Format**:
```python
{
    "bin": "LG:Br|NUC:aniline; ortho=1, para_EWG",
    "top_metals": "Pd(OAc)2 (15/20), Pd2(dba)3 (3/20)",
    "top_ligands": "XPhos (12/20), SPhos (5/20)",
    "temperature_C": "80-110°C (median: 100)",
    "precedent_count": 20
}
```

---

## 🔬 Utilities & Helpers

### `reaction_similarity.py` (169 lines)
**Purpose**: DRFP-based reaction similarity calculations

**Key Functions**:
- `drfp_available() -> bool` - Check if DRFP installed
- `encode_drfp(rsmi: str, n_bits: int = 4096) -> Optional[Any]` - Encode reaction
- `tanimoto(a, b) -> float` - Tanimoto similarity (0.0-1.0)
- `encode_drfp_cached(rsmi: str) -> Optional[Any]` - LRU-cached encoding
- `load_precomputed_npz(path: str) -> dict` - Load precomputed fingerprints

**Features**:
- Graceful degradation if DRFP unavailable
- RDKit warning suppression
- NPZ precomputation support
- LRU caching

**Parameters**:
- `n_bits`: 4096 (default)
- `radius`: 3 (default)

---

### `reaction_type_detector.py` (273 lines)
**Purpose**: Optional rxn-insight integration for ML-based reaction classification

**Key Functions**:
- `is_available() -> bool` - Check if rxn-insight installed
- `detect_reaction_type(reaction_smiles: str) -> Dict` - Classify reaction
  - Returns: rxn_class, rxn_name, confidence, mapped_family, catalysts

**Mapping to ChemTools Families**:
- "C-C Coupling" → Suzuki_CC, Sonogashira, etc.
- "C-N Coupling" → C_N_Coupling_Pd/Cu/Ni
- "Amide/Peptide coupling" → Amide_Formation
- "Esterification" → Esterification

**Features**:
- Graceful fallback if unavailable
- Metal detection from agents
- Family refinement (Pd vs Cu vs Ni)
- Confidence scoring

---

### `output_formatter.py` (537 lines)
**Purpose**: Structured output formatting for API responses

**Key Functions**:
- `format_meta()` - Metadata (timestamp, model, status)
- `format_input()` - Input section
- `format_detection()` - Reaction detection info
- `format_conditions()` - Condition formatting
- `format_recommendation()` - Full recommendation
- `format_ml_output()` - ML recommendation format
- `format_rule_output()` - Rule-based format
- `expand_rule_conditions_to_recommendations()` - Expand rule lists

**Output Structure**:
```python
{
    "meta": {...},
    "input": {...},
    "detection": {...},
    "recommendations": [
        {
            "rank": 1,
            "conditions": {...},
            "precedents": [...],
            "confidence": 0.85
        }
    ]
}
```

**Features**:
- Reagent enrichment
- Precedent formatting
- Confidence scoring
- Trace information

---

### `condition_core.py` (302 lines) ⚠️ **MAY BE DEPRECATED**
**Purpose**: Legacy condition parsing from reagent lists

**Key Functions**:
- `parse(reagents: List[Reagent], text: str = None) -> Dict` - Parse conditions
- `parse_core(reagents: List[Reagent], text: str = "") -> Dict` - Core parsing

**Features**:
- Metal/ligand pairing
- Dataset alias resolution
- Normalization

**Status**: ⚠️ **May be obsolete** - Consider deprecating in favor of `reagent_lookup.py`

---

### `selector_payloads.py` (302 lines)
**Purpose**: Build feature payloads for rule-based selector matching (amide formation)

**Key Functions**:
- `build_amide_selector_payload(role_pack: Dict, features: Dict) -> Dict`
  - Analyzes acid (carboxylic acid classification)
  - Analyzes amine (primary, secondary, aniline, etc.)
  - Infers acid/amine classes
  - Determines sterics and nucleophilicity
  - Functional group detection
  - Green chemistry/water tolerance flags

**Helper Functions**:
- `_alpha_amino_context()` - Detect α-amino acids
- `_choose_acid_class()` - Classify acid type
- `_choose_amine_class()` - Classify amine type
- `_amine_sterics()` - Steric hindrance
- `_amine_nucleophilicity()` - Nucleophilicity assessment
- `_functional_group_flags()` - Detect phenol, alcohol, etc.

**Used By**: `recommend.py` for amide formation recommendations

---

## 📂 Subdirectories

### `featurizers/` - Molecular Featurization

#### `molecular.py` (complex)
**Purpose**: Generic molecular feature extraction

**Key Functions**:
- `featurize(elec: str, nuc: str) -> Dict` - Featurize electrophile + nucleophile
- Extracts: aromaticity, ring count, heteroatoms, functional groups, etc.

#### `ullmann.py` (complex)
**Purpose**: Ullmann C-N coupling specific featurization

**Key Functions**:
- `featurize_ullmann(elec: str, nuc: str) -> Dict` - C-N specific features
- Extracts: ortho substituents, leaving group, nucleophile class, etc.

---

### `util/` - Utility Modules

#### `rdkit_helpers.py`
**Purpose**: RDKit wrapper functions with graceful degradation

**Key Functions**:
- `rdkit_available() -> bool` - Check RDKit availability
- `parse_smiles(smi: str) -> Mol | None` - Parse SMILES
- `canonical_smiles(smi: str) -> str` - Canonicalize
- `neutralize_and_standardize(mol: Mol) -> Mol` - Neutralize charges
- `choose_largest_organic_fragment(mol: Mol) -> Mol` - Fragment selection
- `mol_to_canonical_smiles(mol: Mol) -> str` - Mol → SMILES

**Features**:
- Graceful None returns if RDKit unavailable
- Exception handling
- Fragment processing

---

### `ml/` - Machine Learning Models

#### `drfp_predictor.py`
**Purpose**: DRFP-based yield prediction model

**Key Classes**:
- `DRFPYieldPredictor` - Yield prediction model
  - `train()` - Train on reaction data
  - `predict()` - Predict yield
  - `save()` / `load()` - Model persistence

#### `evaluation.py`
**Purpose**: Model evaluation metrics and utilities

---

### `features/role/` - Role-Aware Features
**Purpose**: Advanced role-aware molecular featurization (optional)

**Integration**: Used by `recommend.py` when available

---

### `integrations/` - External Integrations

#### `molpipeline.py`
**Purpose**: MolPipeline integration for advanced features

**Key Functions**:
- `is_available() -> bool` - Check availability
- Feature extraction integration

#### `mcp/` - Model Context Protocol
**Purpose**: MCP integration for tool calling

---

### `rule_scdb_matcher/` - Rule-Based Matching ⭐

#### `matcher.py` (727 lines)
**Purpose**: Deterministic SMARTS-based scheme matching

**Key Functions**:
- `match(db: RuleDB, rxn_smiles: str) -> MatchResult`
  - SMARTS pattern matching
  - Multi-selector support (metal, ligand, base, solvent)
  - Priority and specificity scoring
  - Constraint validation

#### `loader.py`
**Purpose**: Load scheme databases from JSON

**Key Functions**:
- `load_db(path: str) -> SchemeConditionDB | SelectorRuleDB`

#### `types.py`
**Purpose**: Type definitions

**Key Classes**:
- `MatchResult` - Match outcome
- `SchemeEntry` - Scheme definition
- `RuleDB` - Database container
- `SelectorRule` - Selector rule

#### `ecn.py`
**Purpose**: Essential core normalization for reaction matching

#### `cli.py`
**Purpose**: Command-line interface for matching

---

### `scdb_matcher/` - Backward Compatibility Shim (**NEW**)

#### `__init__.py`
**Purpose**: Re-export from `rule_scdb_matcher` for backward compatibility

**Enables**:
```python
# Old code still works
from chemtools.scdb_matcher import load_db, match
```

#### `loader.py`
**Purpose**: Re-export loader module

---

### `agent/` - Agent Workflows

#### `config.py`
**Purpose**: Agent configuration management

#### `features/`
**Purpose**: Agent-specific feature extraction

---

## 📊 File Statistics

### By Size (Lines of Code)
1. **`recommend.py`** - 1,454 lines (largest - next refactoring candidate)
2. **`context.py`** - 887 lines
3. **`matcher.py`** (rule_scdb_matcher) - 727 lines
4. **`output_formatter.py`** - 537 lines
5. **`properties.py`** - 317 lines
6. **`condition_core.py`** - 302 lines (⚠️ may be obsolete)
7. **`selector_payloads.py`** - 302 lines
8. **`precedent/search.py`** - ~300 lines ✅ **REFACTORED**
9. **`precedent/loader.py`** - ~280 lines ✅ **REFACTORED**
10. **`router.py`** - 285 lines

### By Category

| Category | Files | Total Lines | Key Modules |
|----------|-------|-------------|-------------|
| **Core API** | 3 | 1,014 | `__init__.py`, `context.py`, `contracts.py` |
| **Core Operations** | 4 | 1,077 | `smiles.py`, `router.py`, `properties.py`, `constraints.py` |
| **Data Operations** | 6 | ~3,043 | `precedent/` (6 modules), `recommend.py`, `recommend_ml.py`, `reagent_lookup.py`, `explain.py` |
| **Utilities** | 5 | 1,481 | All utility/helper modules |
| **Rule-Based** | 5 | 1,000+ | `rule_scdb_matcher/` modules |
| **Featurizers** | 2 | ? | `molecular.py`, `ullmann.py` |
| **ML** | 2 | ? | `drfp_predictor.py`, `evaluation.py` |

### Refactoring Summary

**✅ Completed Refactorings:**
1. **`precedent.py`** (846 lines) → **`precedent/`** package (6 modules, ~150 lines each)
   - Improved maintainability: Single responsibility per module
   - Better code navigation: Logical separation of concerns
   - Easier testing: Isolated functionality
   - Preserved API: 100% backward compatible

**🔄 Recommended Next Refactorings:**
1. **`recommend.py`** (1,454 lines) - Consider splitting by reaction family
2. **`output_formatter.py`** (537 lines) - Consider splitting by output type

---

## 🔍 Potential Cleanup Candidates

### ⚠️ Files That May Be Obsolete

1. **`condition_core.py`** (302 lines)
   - **Reason**: Replaced by `reagent_lookup.py`
   - **Status**: Still imported by some code, but functionality is duplicated
   - **Recommendation**: Audit usage and consider deprecation
   - **Current imports**: Used in some legacy parsing flows

2. **`selector_payloads.py`** (302 lines)
   - **Reason**: Highly specific to amide formation rules
   - **Status**: Used by `recommend.py` but could be integrated
   - **Recommendation**: Keep if rule-based amide matching is important, otherwise consider merging into `recommend.py`

### 🔄 Files That Need Refactoring

1. **`recommend.py`** (1,454 lines) ⚠️ **HIGH PRIORITY**
   - **Issue**: Very large, multiple responsibilities
   - **Recommendation**: Split into family-specific modules:
     - `recommend/core.py` - Main recommendation logic
     - `recommend/amide.py` - Amide-specific features (538 lines)
     - `recommend/cn_coupling.py` - C-N coupling features
     - `recommend/families.py` - Family mappings
     - `recommend/utils.py` - Helper functions

2. **`output_formatter.py`** (537 lines)
   - **Issue**: Large with multiple formatting responsibilities
   - **Recommendation**: Consider splitting by output type (meta, detection, conditions, etc.)

### ✅ Files That Are Well-Structured

1. **`precedent/`** - ✨ **RECENTLY REFACTORED** - Clean modular structure
2. **`context.py`** - Clean namespace pattern
3. **`smiles.py`** - Focused responsibility
4. **`router.py`** - Clear detection logic
5. **`reagent_lookup.py`** - Well-organized enrichment
6. **`rule_scdb_matcher/`** - Good separation of concerns

---

## 🎯 Recommendations

### Immediate Actions
1. ✅ **DONE**: `precedent.py` refactored into modular `precedent/` package
2. ✅ **Keep**: `rule_scdb_matcher/` and new `scdb_matcher/` shim (integrated)
3. ⚠️ **Audit**: `condition_core.py` - Check if still needed after reagent_lookup integration
4. 📝 **Document**: Add docstrings to larger modules (recommend.py)
5. 📋 **PROPOSED**: Refactor `recommend.py`, `recommend_ml.py`, and `/ml` folder
   - **See**: `RECOMMEND_ML_REFACTORING_PROPOSAL.md` for detailed plan
   - **Goal**: Split 1,454-line monolith into focused `recommend/` package
   - **Benefits**: Better maintainability, family-specific isolation, clearer ML integration

### Medium-Term Improvements
1. **Split `recommend.py`**: Break into family-specific modules (next priority)
2. **Standardize**: Ensure all modules follow same pattern (stateless vs stateful)
3. **Type hints**: Add comprehensive type hints to all public functions
4. **Tests**: Ensure test coverage for all core modules

### Long-Term Architecture
1. **Plugin system**: Make family-specific logic pluggable
2. **Configuration**: Centralize configuration (currently scattered across env vars)
3. **Caching**: Unified caching strategy across all modules
4. **Logging**: Structured logging throughout

---

## 📝 Summary

The `/chemtools` directory is **well-organized** with clear separation between:
- Core stateless operations (SMILES, router, properties, constraints)
- Data-driven stateful operations (precedent, recommend, reagent lookup)
- Utilities and helpers
- Advanced features (ML, integrations, rule-based)

**Strengths**:
- ✅ Clean ChemTools v2.0 API via `context.py`
- ✅ Good separation of concerns
- ✅ Modular architecture (`precedent/` refactored, `rule_scdb_matcher/` well-structured)
- ✅ Performance optimizations (DRFP precomputation, selective loading)
- ✅ Graceful degradation (RDKit, DRFP, rxn-insight)
- ✅ Thread-safe caching and resource management

**Recent Improvements** ✨:
- ✅ `precedent.py` (846 lines) refactored into `precedent/` package (6 focused modules)
- ✅ Improved maintainability and testability
- ✅ 100% backward compatible

**Areas for Improvement**:
- ⚠️ `recommend.py` still large (1,454 lines) - **next refactoring priority**
- ⚠️ Potential duplicate functionality (condition_core.py vs reagent_lookup.py)
- 📝 Documentation could be more comprehensive
- 🧪 Test coverage varies

**Overall Assessment**: 🌟🌟🌟🌟½ (4.5/5 stars)
- Solid, well-refactored architecture
- Well-integrated ChemTools v2.0
- Continuous improvement demonstrated
- Production-ready
