# HTE-Based Condition Recommendation System

## Overview

The HTE (High-Throughput Experimentation) module provides condition recommendations based on reactant types using a database of 66,308 experimental results covering 41 reaction types.

**Key Innovation**: Since no reaction SMILES is provided in the HTE database, recommendations are made by:
1. Detecting reactant types from input SMILES using existing `chemtools` reactant classification
2. Matching against experimental conditions indexed by reactant type combinations
3. Ranking conditions by success rate, sample size, and average yield

## Architecture

```
chemtools/HTE/
├── __init__.py          # Package exports
└── recommender.py       # Core recommendation engine

Key Components:
- HTERecommender: Main recommendation class
- HTERecommendationResult: Result dataclass with ranked conditions
- ConditionRecommendation: Individual condition with statistics
```

## Database Structure

**HTE_0.jsonl** (66,308 experiments):
- **reaction_type**: 41 types (c_n_cross_coupling, suzuki_miyaura, c_h_activation, etc.)
- **reactant_types**: Detected member types (Ar-Br, Any-NH2, Ar-B(OH)2, etc.)
- **catalyst_type**: Catalyst type tags (Pd, Cu, Ni, organocatalyst, etc.)
- **conditions**: catalyst, ligand, base, solvent, secondary_solvent, additive, coupling_reagent
- **metrics**: area_total_reduced, z_score

**Statistics**:
- Success rate (yield > 50%): 18.8% overall
- Top reactions: C_N_Coupling (24,012), Suzuki (11,588), Arylation-acidic-C-H (4,152)
- 229 catalysts, 153 ligands, 132 bases, 67 solvents
- 71 unique reactant type combinations

## How It Works

### 1. Reactant Type Detection

Uses existing `chemtools.featurizers.analysis.reactants.classify_reactant_smiles()`:

```python
from chemtools.featurizers.analysis.reactants import classify_reactant_smiles

# Detect types
result_a = classify_reactant_smiles("c1ccc(Br)cc1")  # ArBr (ArX*)
result_b = classify_reactant_smiles("CCN")  # RNH2 (RNH2/R2NH)
```

### 2. Database Indexing

On initialization, the recommender builds indices:
- **By reactant type combination**: `(ArBr, RNH2)` → DataFrame subset
- **Reaction type patterns**: `(ArBr, RNH2)` → `{C_N_Coupling: 1080, ...}`

This enables O(1) lookup by reactant types.

### 3. Condition Aggregation

For matched experiments:
1. Group by `(catalyst, ligand, base, solvent)` combination
2. Calculate statistics per combination:
   - **Success rate**: % with yield > 50%
   - **Average/median yield**
   - **Sample size**: number of experiments
3. Compute **confidence score**:
   ```
   confidence = 0.5 * (success_rate/100)
              + 0.3 * min(num_experiments, 100)/100
              + 0.2 * (avg_yield/100)
   ```

### 4. Ranking & Selection

Conditions ranked by confidence score (descending), return top-k.

## Usage

### Basic Usage

```python
from chemtools.HTE import HTERecommender, format_result

# Initialize (loads database, builds indices)
recommender = HTERecommender()

# Get recommendations
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",  # Bromobenzene
    reactant_b_smiles="CCN",            # Ethylamine
    top_k=5
)

# Display results
print(format_result(result))
```

**Output**:
```
Reactant A: c1ccc(Br)cc1
  Type: ArBr (ArX*)
Reactant B: CCN
  Type: RNH2 (RNH2/R2NH)

🎯 PREDICTED REACTION TYPE: C_N_Coupling
   Confidence: 100.0%

📊 DATABASE MATCH: 1080 matching experiments (1.63% of database)

🏆 TOP RECOMMENDATIONS: 5 conditions found

Recommendation #1
Confidence Score: 65.3/100
Success Rate: 100.0% (2 experiments)
Avg Yield: 73.3% | Median: 73.3%

🧪 CONDITIONS:
  Catalyst: XantPhos Pd(allyl)Cl
  Ligand: XantPhos
  Base: Cs2CO3
  Solvent: Dioxane
  Secondary Solvent: water
```

### Advanced Usage

#### Filter by Reaction Type

```python
# Only recommend Suzuki conditions
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",     # Chlorobenzene
    reactant_b_smiles="c1ccc(B(O)O)cc1",  # Phenylboronic acid
    reaction_type_filter="Suzuki",
    top_k=5
)
```

#### Adjust Confidence Threshold

```python
# Require at least 5 experiments per condition
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    min_experiments=5,
    top_k=10
)
```

#### Single Reactant (e.g., CH-Activation)

```python
# For reactions with one reactant
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles=None,  # No second reactant
    top_k=5
)
```

### Database Statistics

```python
stats = recommender.get_statistics()
print(stats)
# {
#   'total_experiments': 66308,
#   'reaction_types': 41,
#   'unique_type_combinations': 71,
#   'success_rate_overall': 18.84,
#   'avg_yield': 22.55,
#   'catalysts': 229,
#   'ligands': 153,
#   'bases': 132,
#   'solvents': 67
# }
```

## API Reference

### HTERecommender

```python
class HTERecommender:
    def __init__(self, hte_db_path: str = "data/HTE_db/HTE_0.jsonl")
    
    def recommend(
        self,
        reactant_a_smiles: str,
        reactant_b_smiles: Optional[str] = None,
        top_k: int = 10,
        min_experiments: int = 2,
        reaction_type_filter: Optional[str] = None
    ) -> HTERecommendationResult
    
    def get_statistics(self) -> Dict[str, Any]
```

### HTERecommendationResult

```python
@dataclass
class HTERecommendationResult:
    reactant_a_smiles: str
    reactant_b_smiles: Optional[str]
    
    # Detected types
    reactant_a_type: Optional[str]
    reactant_b_type: Optional[str]
    reactant_a_category: Optional[str]
    reactant_b_category: Optional[str]
    
    # Predicted reaction type
    predicted_reaction_type: Optional[str]
    reaction_type_confidence: float
    
    # Recommendations
    recommendations: List[ConditionRecommendation]
    
    # Metadata
    total_matching_experiments: int
    database_coverage: float
```

### ConditionRecommendation

```python
@dataclass
class ConditionRecommendation:
    # Conditions
    catalyst: str
    ligand: str
    base: str
    solvent: str
    secondary_solvent: Optional[str] = None
    additive: Optional[str] = None
    coupling_reagent: Optional[str] = None
    
    # Statistics
    success_rate: float  # % yield > 50
    avg_yield: float
    median_yield: float
    num_experiments: int
    confidence_score: float  # 0-100
    
    # Metadata
    reaction_type: Optional[str]
    reactant_types: Tuple[str, str]
    z_score_range: Tuple[float, float]
```

## Performance Characteristics

- **Initialization**: ~1-2 seconds (loads 66K rows, builds indices)
- **Query time**: <100ms (O(1) index lookup + aggregation)
- **Memory**: ~50MB (full DataFrame + indices)

## Example Results

### C-N Coupling (ArBr + RNH2)

**Input**: Bromobenzene + Ethylamine
- **Match**: 1,080 experiments (1.63% of DB)
- **Predicted**: C_N_Coupling (100% confidence)
- **Top condition**: XantPhos Pd(allyl)Cl / XantPhos / Cs2CO3 / Dioxane
  - Success rate: 100% (2 experiments)
  - Avg yield: 73.3%
  - Confidence: 65.3/100

### Suzuki Coupling (ArCl + ArB(OH)2)

**Input**: Chlorobenzene + Phenylboronic acid
- **Match**: 960 experiments (1.45% of DB)
- **Predicted**: Suzuki (100% confidence)
- **Top condition**: Pd-PEPPSI-IPent / IPENT Cl / Cs2CO3 / Brij 35
  - Success rate: 100% (2 experiments)
  - Avg yield: 97.5%
  - Confidence: 70.1/100

### C-N Coupling (ArCl + ArNH2)

**Input**: Chlorobenzene + Aniline
- **Match**: 2,009 experiments (3.03% of DB)
- **Predicted**: C_N_Coupling (100% confidence)
- **Top condition**: tBuBrettPhos Pd(allyl)OTf / tBuBrettPhos / K3PO4 / MeCN
  - Success rate: 90.9% (88 experiments)
  - Avg yield: 84.5%
  - Confidence: 88.8/100

## Limitations & Considerations

1. **Reactant Type Dependency**: Recommendations only as good as type detection
   - If reactant type not detected → no recommendations
   - Rare type combinations → limited data

2. **No Reaction SMILES**: Cannot account for:
   - Substrate-specific effects (steric hindrance, electronics)
   - Regioselectivity requirements
   - Functional group compatibility beyond types

3. **Statistical Bias**:
   - HTE data may favor certain catalyst/ligand families
   - Success threshold (50% yield) is arbitrary
   - Small sample sizes for some conditions

4. **Missing Information**:
   - Temperature, time, concentration not in database
   - Exact reagent stoichiometry not specified
   - Work-up conditions not included

## Integration with Existing Systems

### Complementary to Rule-Based Systems

```python
# Hybrid approach: combine HTE + rule-based
from chemtools.recommend.unified import recommend_conditions_unified
from chemtools.HTE import HTERecommender

# Get rule-based recommendations
rule_recs = recommend_conditions_unified(
    reaction_smiles="c1ccc(Br)cc1.CCN>>c1ccc(NCC)cc1"
)

# Get HTE recommendations
hte = HTERecommender()
hte_recs = hte.recommend("c1ccc(Br)cc1", "CCN")

# Compare/combine results
```

### CLI Integration

Create `chemtools/HTE/cli.py` for command-line access (see Implementation section).

## Testing

```bash
# Run demo
python demo_hte_recommender.py

# Test specific cases
python -c "
from chemtools.HTE import HTERecommender
rec = HTERecommender()
result = rec.recommend('c1ccc(Br)cc1', 'CCN', top_k=3)
print(f'Found {len(result.recommendations)} recommendations')
"
```

## Future Enhancements

1. **Similarity-based expansion**: Use molecular fingerprints to match similar (not exact) reactant types
2. **Multi-objective optimization**: Balance yield, cost, sustainability
3. **Active learning**: Suggest experiments to fill data gaps
4. **Uncertainty quantification**: Provide prediction intervals
5. **Reaction template integration**: Combine with SMARTS-based matching
6. **Cost/availability data**: Rank by reagent accessibility

## References

- HTE Database: `data/HTE_db/HTE_0.jsonl`
- Reactant Type System: `chemtools/taxonomy/data/reactant_types.json`
- Classification: `chemtools/analysis/reactants.py`
