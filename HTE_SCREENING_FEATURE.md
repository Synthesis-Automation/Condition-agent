# HTE Screening Set Generation

## Overview

The HTE screening set generation feature is designed for the **PRIMARY use case of HTE systems**: generating a diverse group of conditions to test in parallel on screening plates.

## Key Features

- **Generate up to 24 conditions** (standard 4×6 plate format, expandable to 96-well)
- **Three diversity strategies**: balanced, top_performers, diverse
- **Automatic plate position assignment** (A1-D6 for 24-well, extendable to A1-H12 for 96-well)
- **Reagent diversity maximization** across catalyst, ligand, base, and solvent
- **Based on 66,308 HTE experiments** with z-score ranking
- **Fast performance**: <200ms for 24 condition generation

## Diversity Strategies

### 1. Balanced (Default)
- **Best for**: General screening campaigns
- **Composition**: ~1/3 top performers + 2/3 diverse alternatives
- **Example**: 8 highest z-score conditions + 16 diverse picks
- **Use when**: You want a mix of proven winners and exploratory conditions

### 2. Top Performers
- **Best for**: Hit optimization, validation studies
- **Composition**: Strictly ranked by z-score
- **Example**: Top 24 conditions by performance
- **Use when**: You want to maximize chances of success

### 3. Diverse
- **Best for**: Broad reagent exploration, substrate scope studies
- **Composition**: Maximum reagent variation (60% diversity, 40% performance)
- **Example**: 24 conditions with maximally different reagent combinations
- **Use when**: Exploring reagent space or establishing structure-activity relationships

## Usage Examples

### Python API

```python
from chemtools.HTE import HTERecommender

recommender = HTERecommender()

# Generate 24-condition screening set for C-N coupling with Cu
result = recommender.generate_screening_set(
    reactant_a_smiles="Brc1ccc(C)cc1",  # 4-bromotoluene
    reactant_b_smiles="Nc1ccccc1",       # aniline
    num_conditions=24,
    reaction_type_filter="C_N_Coupling",
    catalyst_filter="Cu",
    diversity_strategy="balanced"
)

print(f"Found {result.total_matching_experiments} experiments")
print(f"Generated {len(result.recommendations)} conditions")

# Access individual conditions
for i, rec in enumerate(result.recommendations, 1):
    print(f"{i}. {rec.catalyst}/{rec.ligand}/{rec.base}/{rec.solvent}")
    print(f"   Z-score: {rec.avg_z_score:.2f}, Yield: {rec.avg_yield:.1f}%")
```

### Agent Tool

```python
from chem_assistant.chemtools_wrapper import hte_screening_set_tool

result = hte_screening_set_tool.invoke({
    "reaction_smiles": "Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1",
    "num_conditions": 24,
    "reaction_type_filter": "C_N_Coupling",
    "catalyst_filter": "Cu",
    "diversity_strategy": "balanced"
})

if result['success']:
    print(f"Generated {result['num_conditions']} conditions")
    print(f"Plate format: {result['plate_format']}")
    
    # Access screening conditions with plate positions
    for cond in result['screening_conditions']:
        print(f"{cond['plate_position']}: {cond['catalyst']}/{cond['ligand']}")
```

### Natural Language (Agent)

```
Generate a screening set for 4-bromotoluene + aniline C-N coupling with copper catalyst. 
I need 24 conditions for a 4x6 plate, use balanced strategy.
```

The agent will automatically invoke `hte_screening_set_tool` and format the results with plate positions.

## Algorithm Details

### Diversity Selection

The diversity algorithm:

1. **Starts with top performers**: First ~1/3 are highest z-score conditions
2. **Iteratively selects diverse conditions** from remaining pool
3. **Tracks used reagents**: Maintains sets of used catalysts, ligands, bases, solvents
4. **Scores each candidate** based on:
   - **Diversity score**: +1 for each new reagent type (0-4 range)
   - **Performance score**: Normalized z-score (0-1 range)
   - **Combined score**: 60% diversity + 40% performance (balanced mode)
5. **Selects best candidate** and updates used reagent sets
6. **Repeats** until desired number of conditions reached

### Plate Position Assignment

- **24-well (4×6)**: A1-D6 (rows A-D, columns 1-6)
- **96-well (8×12)**: A1-H12 (rows A-H, columns 1-12)
- Positions assigned sequentially by rank

## Test Results

Using 4-bromotoluene + aniline (ArBr + ArNH2) with Cu catalyst:

### Balanced Strategy (24 conditions)
- **Total experiments found**: 112
- **Top condition**: CuI/PPBO/NaOtBu/tAmOH (Z=2.61, Y=49.8%)
- **Last condition**: CuI/N,N-Dimethylglycine/Cs2CO3/DMSO (Z=0.37, Y=17.2%)
- **Z-score range**: 2.61 to 0.37 (2.24 span)

### Diverse Strategy (24 conditions)
- **Unique reagents**:
  - Catalysts: 2 (CuI, Cu)
  - Ligands: 23 (maximal diversity)
  - Bases: 5 (NaOtBu, Cs2CO3, K3PO4, K2CO3, LHMDS)
  - Solvents: 6 (tAmOH, DMSO, DMF, Dioxane, Toluene, MeCN)

### Performance
- **Database loading**: ~50ms (cached after first call)
- **Screening set generation**: ~100ms for 24 conditions
- **Total time**: <200ms end-to-end

## Integration with Workflow

### Typical HTE Workflow

1. **Define substrate**: Input reaction SMILES
2. **Set filters**: Choose reaction type and catalyst
3. **Generate screening set**: Use `hte_screening_set_tool` with desired strategy
4. **Export to plate**: Use plate positions (A1-D6) for liquid handler
5. **Run experiments**: Execute on HTE platform
6. **Analyze results**: Compare to predicted z-scores and yields

### When to Use Each Tool

| Tool | Use Case |
|------|----------|
| `hte_screening_set_tool` | Generate diverse condition sets for parallel screening (12-24 conditions) |
| `hte_recommend_tool` | Get top-K best conditions for focused optimization (1-5 conditions) |
| `hte_analytics_tool` | Explore database statistics, reactant pairs, catalyst usage |
| `hte_conditions_tool` | Query all conditions for specific substrate pair |

## Parameters

### `generate_screening_set()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reactant_a_smiles` | str | Required | SMILES of first reactant |
| `reactant_b_smiles` | str | Optional | SMILES of second reactant |
| `num_conditions` | int | 24 | Number of conditions to generate |
| `min_experiments` | int | 1 | Minimum experiments per condition |
| `reaction_type_filter` | str | None | Filter by reaction type (e.g., 'C_N_Coupling') |
| `catalyst_filter` | str | None | Filter by metal (e.g., 'Cu', 'Pd', 'Ni') |
| `diversity_strategy` | str | "balanced" | Strategy: 'balanced', 'top_performers', 'diverse' |

### Return Format

```python
{
    "success": True,
    "reactant_a_type": "ArBr",
    "reactant_b_type": "ArNH2",
    "matching_experiments": 112,
    "num_conditions": 24,
    "diversity_strategy": "balanced",
    "plate_format": "4x6",
    "screening_conditions": [
        {
            "rank": 1,
            "plate_position": "A1",
            "catalyst": "CuI",
            "ligand": "PPBO",
            "base": "NaOtBu",
            "solvent": "tAmOH",
            "avg_z_score": 2.61,
            "avg_yield": 49.8,
            "median_yield": 50.2,
            "success_rate": 100.0,
            "num_experiments": 1,
            "confidence_score": 63.8
        },
        ...
    ]
}
```

## Advantages Over Traditional Approaches

### vs. Manual Selection
- **10× faster**: Automated screening set in <200ms vs. hours of literature review
- **Data-driven**: Based on 66K experiments, not expert opinion
- **Reproducible**: Same inputs always give same results

### vs. Design of Experiments (DoE)
- **Evidence-based**: Uses actual experimental data, not factorial designs
- **Success-weighted**: Prioritizes conditions that work in practice
- **Flexible**: Can adapt strategy based on goals (optimization vs. exploration)

### vs. Random Selection
- **Smart diversity**: Maximizes reagent coverage while maintaining performance
- **Statistical ranking**: Uses z-scores to identify genuinely superior conditions
- **Coverage**: Ensures representation across reagent space

## Common Use Cases

### 1. Hit Finding Campaign
```python
# Want broad exploration with some guaranteed winners
result = recommender.generate_screening_set(
    reactant_a_smiles="Brc1ccccc1",
    reactant_b_smiles="Nc1ccccc1",
    num_conditions=24,
    diversity_strategy="balanced"  # Mix of top + diverse
)
```

### 2. Catalyst Comparison Study
```python
# Compare different metals for same substrate
for metal in ['Cu', 'Pd', 'Ni']:
    result = recommender.generate_screening_set(
        reactant_a_smiles="Brc1ccccc1",
        reactant_b_smiles="Nc1ccccc1",
        num_conditions=8,
        catalyst_filter=metal,
        diversity_strategy="top_performers"
    )
```

### 3. Substrate Scope Exploration
```python
# Test multiple substrates with diverse conditions
result = recommender.generate_screening_set(
    reactant_a_smiles="Brc1ccccc1",
    reactant_b_smiles="Nc1ccccc1",
    num_conditions=12,
    diversity_strategy="diverse"  # Maximum reagent variation
)
```

## Files Changed

- `chemtools/HTE/recommender.py`: Added `generate_screening_set()` and `_select_diverse_conditions()`
- `chem_assistant/chemtools_wrapper.py`: Added `hte_screening_set_tool` wrapper
- `chem_assistant/chemtools_agent.py`: Updated system prompt with screening tool
- `test_hte_screening.py`: Comprehensive test suite

## Testing

Run the test suite:
```bash
python test_hte_screening.py
```

Expected output:
- ✓ 112 experiments found for ArBr + ArNH2 with Cu
- ✓ 24 conditions generated in balanced mode
- ✓ 12 conditions generated in top_performers mode
- ✓ 24 conditions generated in diverse mode
- ✓ Agent tool integration successful with plate positions
