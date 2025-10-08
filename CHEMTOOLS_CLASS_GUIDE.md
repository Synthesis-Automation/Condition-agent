# ChemTools Master Class - User Guide

## Overview

The `ChemTools` master class provides a unified, object-oriented API for all chemistry operations with built-in resource management, lazy loading, and intelligent caching.

## Quick Start

### Using the Global Instance (Recommended)

The simplest way to use ChemTools is through the global singleton:

```python
from chemtools import chem

# SMILES operations
result = chem.smiles.normalize("CCO")
reaction = chem.smiles.normalize_reaction("c1ccccc1Br.c1cccnc1>>c1ccccc1-c1cccnc1")

# Reaction family detection
family = chem.router.detect_family("c1ccccc1Br.c1cccnc1>>c1ccccc1-c1cccnc1")

# Property lookup
props = chem.properties.lookup("water")

# Precedent search
precedents = chem.precedent.knn(
    family="C_N_Coupling_Pd",
    features={"ligand": "PPh3", "base": "K2CO3"},
    k=5
)

# ML recommendations
recommendations = chem.recommend.conditions(
    reaction="c1ccccc1Br.Nc1cccnc1>>c1ccccc1Nc1cccnc1",
    k=5,
    limit=10
)
```

### Creating Custom Instances

For specific use cases, create configured instances:

```python
from chemtools import ChemTools

# Configure for Buchwald reactions only
buchwald_chem = ChemTools(
    datasets=["C_N_Coupling_Pd"],  # Only load this dataset
    reagent_dbs=["ligand", "base"],  # Only these reagent types
    preload=True  # Load everything at startup
)

# Use the configured instance
precedents = buchwald_chem.precedent.knn(
    family="C_N_Coupling_Pd",
    features={...},
    k=10
)

# For API servers - preload everything
api_chem = ChemTools(
    datasets=["C_N_Coupling_Pd", "C_N_Coupling_Cu", "Suzuki"],
    preload=True  # Fast first requests
)

# For CLI tools - minimal config
cli_chem = ChemTools(
    datasets=None,  # Load on-demand only
    preload=False  # Fast startup
)

# For Jupyter notebooks - balanced config
notebook_chem = ChemTools(
    datasets=["C_N_Coupling_Pd"],  # Your current work
    preload=False  # Load when first used
)
```

## Architecture

```
ChemTools (master class)
│
├── Core Operations (stateless, always available)
│   ├── smiles - SMILES parsing and normalization
│   │   ├── normalize(smi: str) -> Dict
│   │   └── normalize_reaction(rsmi: str) -> Dict
│   │
│   ├── router - Reaction family detection
│   │   └── detect_family(reaction: str) -> Dict
│   │
│   ├── properties - Compound property lookup
│   │   └── lookup(query: str) -> Dict
│   │
│   └── constraints - Constraint validation
│       └── filter(candidates: List, constraints: Dict) -> List
│
├── Data Operations (stateful, lazy-loaded through context)
│   ├── precedent - Precedent reaction search
│   │   ├── knn(family, features, k) -> List
│   │   ├── find_reactions_by_core(core, family) -> List
│   │   └── list_cores(family) -> List
│   │
│   ├── recommend - ML-based recommendations
│   │   ├── conditions(reaction, k, limit) -> Dict
│   │   ├── from_reaction(reaction, k) -> Dict
│   │   └── design_plate(reaction, plate_size) -> Dict
│   │
│   └── explain - Explanation generation
│       └── precedents(precedents) -> str
│
└── Advanced Operations (optional dependencies)
    ├── featurizers - Molecular featurization
    │   └── ullmann - Ullmann C-N coupling featurizer
    │
    ├── features - Role-aware features
    │   ├── is_available() -> bool
    │   └── role - Role-aware feature extraction
    │
    └── integrations - External integrations
        └── mcp - Model Context Protocol
```

## API Reference

### Core Operations

#### SMILES Namespace

```python
# Normalize a SMILES string
result = chem.smiles.normalize("CCO")
# Returns: {
#     "input": "CCO",
#     "fragments": ["CCO"],
#     "largest_smiles": "CCO",
#     "smiles_norm": "CCO"
# }

# Normalize a reaction SMILES
result = chem.smiles.normalize_reaction("Br.N>>BrN")
# Returns: {
#     "input": "Br.N>>BrN",
#     "reactants": [...],
#     "agents": [...],
#     "products": [...],
#     "normalized": "Br.N>>BrN"
# }
```

#### Router Namespace

```python
# Detect reaction family
result = chem.router.detect_family("c1ccccc1Br.Nc1cccnc1>>c1ccccc1Nc1cccnc1")
# Returns: {
#     "family": "C_N_Coupling_Pd",
#     "confidence": 0.95,
#     "method": "rule-based"
# }
```

#### Properties Namespace

```python
# Look up compound properties
result = chem.properties.lookup("water")
# Returns: {
#     "name": "Water",
#     "formula": "H2O",
#     "cas": "7732-18-5",
#     "smiles": "O"
# }
```

#### Constraints Namespace

```python
# Filter candidates by constraints
filtered = chem.constraints.filter(
    candidates=[...],
    constraints={"temperature": {"max": 100}}
)
```

### Data Operations

#### Precedent Namespace

```python
# K-nearest neighbor search
precedents = chem.precedent.knn(
    family="C_N_Coupling_Pd",
    features={
        "ligand": "PPh3",
        "base": "K2CO3",
        "solvent": "toluene"
    },
    k=5,
    relax={"ligand": True}
)

# Find reactions by core structure
reactions = chem.precedent.find_reactions_by_core(
    core="c1ccccc1N",
    family="C_N_Coupling_Pd",
    fuzzy=False,
    limit=50
)

# List available cores
cores = chem.precedent.list_cores(
    family="C_N_Coupling_Pd",
    top_n=200,
    include_counts=True
)
```

#### Recommend Namespace

```python
# Get ML-based recommendations
recommendations = chem.recommend.conditions(
    reaction="c1ccccc1Br.Nc1cccnc1>>c1ccccc1Nc1cccnc1",
    reaction_type="C_N_Coupling_Pd",  # Optional hint
    k=5,  # Precedents per recommendation
    limit=10,  # Max recommendations
    relax={"ligand": True},
    constraints={"temperature": {"max": 100}}
)

# Simpler recommendation API
recommendations = chem.recommend.from_reaction(
    reaction="c1ccccc1Br.Nc1cccnc1>>c1ccccc1Nc1cccnc1",
    k=5
)

# Design experimental plate
plate = chem.recommend.design_plate(
    reaction="c1ccccc1Br.Nc1cccnc1>>c1ccccc1Nc1cccnc1",
    plate_size=96
)
```

#### Explain Namespace

```python
# Generate explanation from precedents
explanation = chem.explain.precedents(
    precedents=[...],
    format="markdown"
)
```

### Advanced Operations

#### Featurizers Namespace

```python
# Access Ullmann featurizer
ullmann = chem.featurizers.ullmann
```

#### Features Namespace

```python
# Check if role-aware features are available
if chem.features.is_available():
    role_features = chem.features.role
```

#### Integrations Namespace

```python
# Access MCP integration
mcp = chem.integrations.mcp
```

## Resource Management

### Configuration

```python
from chemtools import ChemTools, ResourceConfig

# Create configuration
config = ResourceConfig(
    datasets=["C_N_Coupling_Pd"],
    ml_models=["buchwald"],
    reagent_dbs=["ligand", "base"],
    preload=True,
    cache_size=32,
    enable_rdkit=True
)

# Create instance with config
chem = ChemTools(config=config)

# Or use constructor parameters
chem = ChemTools(
    datasets=["C_N_Coupling_Pd"],
    preload=True
)
```

### Cache Management

```python
# Get cache statistics
stats = chem.get_cache_stats()
# Returns: {
#     "datasets_loaded": 1,
#     "dataset_families": ["C_N_Coupling_Pd"],
#     "ml_models_loaded": 0,
#     "ml_model_names": [],
#     "reagent_dbs_loaded": 2,
#     "reagent_db_types": ["ligand", "base"],
#     "total_resources": 3
# }

# List loaded datasets
families = chem.list_loaded_datasets()
# Returns: ["C_N_Coupling_Pd"]

# Unload a dataset to free memory
chem.unload_dataset("C_N_Coupling_Pd")

# Clear all caches
chem.clear_cache()
```

### Manual Resource Loading

```python
# Explicitly load a dataset
dataset = chem.get_reaction_dataset("C_N_Coupling_Pd")

# Load a reagent database
ligands = chem.get_reagent_database("ligand")

# Load an ML model (Phase 2)
model = chem.get_ml_model("buchwald")
```

## Performance Characteristics

### Global Instance (Default)

```python
from chemtools import chem

# First call - loads dataset (~1-2s for C_N_Coupling_Pd with precomputed DRFP)
precedents = chem.precedent.knn("C_N_Coupling_Pd", features={...}, k=5)

# Second call - uses cached dataset (~10ms)
precedents = chem.precedent.knn("C_N_Coupling_Pd", features={...}, k=5)
```

**Performance**: 
- First search: ~1-2s (loads 1,343 reactions with precomputed DRFP)
- Subsequent searches: ~10ms (cached dataset)
- Memory: ~50MB for C_N_Coupling_Pd dataset

### Preloaded Instance

```python
from chemtools import ChemTools

# Startup time: ~2-3s (loads all configured resources)
chem = ChemTools(
    datasets=["C_N_Coupling_Pd", "C_N_Coupling_Cu"],
    preload=True
)

# All searches are fast (~10ms)
precedents = chem.precedent.knn("C_N_Coupling_Pd", features={...}, k=5)
```

**Performance**:
- Startup: ~2-3s (loads 2 datasets)
- All searches: ~10ms (everything cached)
- Memory: ~100MB (2 datasets)

### Selective Loading

```python
# Only load specific datasets - 50-100x faster than loading all
chem = ChemTools(datasets=["C_N_Coupling_Pd"])  # 1,343 reactions

# vs loading everything
# chem = ChemTools(datasets=None)  # 99,668 reactions
```

**Performance**:
- Selective (1 family): ~1-2s load, ~50MB memory
- All families: ~85s load, ~2GB memory
- **Speedup**: 50-100x faster with selective loading

## Use Cases

### 1. API Server (FastAPI/Flask)

```python
# app/main.py
from fastapi import FastAPI
from chemtools import ChemTools

# Create preloaded instance at startup
app_chem = ChemTools(
    datasets=["C_N_Coupling_Pd", "C_N_Coupling_Cu", "Suzuki"],
    reagent_dbs=["ligand", "base", "solvent"],
    preload=True  # Fast first requests
)

app = FastAPI()

@app.get("/api/v1/precedent/knn")
def get_precedents(family: str, k: int = 5):
    # Uses cached resources - very fast
    return app_chem.precedent.knn(family, features={...}, k=k)
```

### 2. Command-Line Tool

```python
# cli.py
import click
from chemtools import ChemTools

@click.command()
@click.option('--dataset', default='C_N_Coupling_Pd')
@click.option('--reaction', required=True)
def recommend(dataset, reaction):
    # Create minimal instance - fast startup
    chem = ChemTools(
        datasets=[dataset],
        preload=False  # Load on-demand
    )
    
    result = chem.recommend.conditions(reaction)
    click.echo(result)

if __name__ == '__main__':
    recommend()
```

### 3. Jupyter Notebook

```python
# notebook.ipynb
from chemtools import ChemTools

# Create instance for your workflow
chem = ChemTools(
    datasets=["C_N_Coupling_Pd"],  # Your current research
    preload=False  # Load when first used
)

# Explore precedents
precedents = chem.precedent.knn(
    family="C_N_Coupling_Pd",
    features={"ligand": "PPh3"},
    k=10
)

# Get recommendations
recs = chem.recommend.conditions(
    reaction="c1ccccc1Br.Nc1cccnc1>>c1ccccc1Nc1cccnc1"
)
```

### 4. Testing

```python
# tests/test_recommendations.py
import pytest
from chemtools import ChemTools

@pytest.fixture
def test_chem():
    """Create isolated instance for testing."""
    return ChemTools(
        datasets=["C_N_Coupling_Pd"],
        preload=True
    )

def test_precedent_search(test_chem):
    result = test_chem.precedent.knn(
        family="C_N_Coupling_Pd",
        features={...},
        k=5
    )
    assert len(result) == 5
```

## Best Practices

### 1. Use the Global Instance for Simple Scripts

```python
from chemtools import chem

# Quick and easy
result = chem.smiles.normalize("CCO")
```

### 2. Create Custom Instances for Applications

```python
from chemtools import ChemTools

# Configure for your use case
app_chem = ChemTools(datasets=["C_N_Coupling_Pd"], preload=True)
```

### 3. Selective Loading for Performance

```python
# Load only what you need
chem = ChemTools(
    datasets=["C_N_Coupling_Pd"],  # Not all 10+ families
    reagent_dbs=["ligand", "base"]  # Not all 7+ types
)
```

### 4. Preload for Production

```python
# Fast first requests
production_chem = ChemTools(
    datasets=["C_N_Coupling_Pd", "Suzuki"],
    preload=True
)
```

### 5. On-Demand for Development

```python
# Fast startup during development
dev_chem = ChemTools(
    datasets=None,  # Load as needed
    preload=False
)
```

## Troubleshooting

### Issue: Slow First Search

**Symptom**: First precedent search takes 85+ seconds

**Solution**: Use selective loading
```python
# Instead of loading all datasets
chem = ChemTools(datasets=["C_N_Coupling_Pd"])  # 50-100x faster
```

### Issue: High Memory Usage

**Symptom**: Process using 2+ GB of memory

**Solution**: Load only needed datasets
```python
# Unload unused datasets
chem.unload_dataset("Suzuki")

# Or create instance with specific datasets
chem = ChemTools(datasets=["C_N_Coupling_Pd"])  # ~50MB vs 2GB
```

### Issue: Import Errors

**Symptom**: `ImportError: cannot import name 'chem'`

**Solution**: Update your imports
```python
# Make sure you're using the latest version
from chemtools import chem  # Global instance
from chemtools import ChemTools  # Class for custom instances
```

## Future Enhancements

See `RESOURCE_MANAGER_IMPLEMENTATION_PLAN.md` for planned features:

- **Phase 2**: ML model integration and caching
- **Phase 3**: Rule-based system integration
- **Phase 4**: Observable pattern for cache monitoring
- **Phase 5**: Full API integration with lifespan management
- **Phase 6**: Comprehensive testing and benchmarking

## Related Documentation

- `RESOURCE_MANAGER_IMPLEMENTATION_PLAN.md` - Full architecture plan
- `RESOURCE_MANAGER_QUICK_REF.md` - Quick reference guide
- `RESOURCE_MANAGER_SUMMARY.md` - Executive summary
- `DRFP_PRECOMPUTATION_OPTIMIZATION.md` - DRFP performance optimization
- `REGISTRY_REMOVAL_SUMMARY.md` - Registry migration guide
