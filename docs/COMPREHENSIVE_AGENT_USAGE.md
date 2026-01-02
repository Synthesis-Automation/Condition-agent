# Comprehensive Multi-Database Agent Usage

## Overview

The enhanced `chem_assistant.planner.llm_agent_cli` now performs **systematic, comprehensive analysis** by querying ALL available databases and tools before making recommendations.

## What Changed

### Before (Simple Agent)

- Used only 7 basic planner tools
- Limited structural analysis
- Single-source recommendations
- Quick but incomplete

### After (Comprehensive Agent)

- Uses **15+ specialized tools**
- **4-Phase systematic methodology**
- **Multi-database evidence synthesis**
- Thorough but more expensive

## The 4-Phase Analysis Process

### Phase 1: Structural Analysis

The agent first understands what you're working with:

1. **detect_reaction_family_tool** - Identifies reaction type (Suzuki, Buchwald, etc.)
2. **classify_reactant_tool** - Analyzes each substrate (aryl halide, amine, etc.)
3. **get_functional_groups_tool** - Identifies key functional groups
4. **analyze_bond_changes_tool** - Understands mechanism and bond formation

**Output**: Complete structural understanding of your reaction

### Phase 2: Multi-Database Search

The agent queries EVERY available knowledge source:

#### A. Rule Database

- **rule_based_conditions_tool** - Validated expert rules

#### B. ML Precedent Databases (3 tools)

- **recommend_conditions_tool** - ML-ranked precedents
- **search_precedents_tool** - Detailed precedent exploration
- **unified_recommender_tool** - Unified dataset/protocol/HTE search

#### C. Protocol Database

- **protocol_recommendation_tool** - Complete experimental protocols

#### D. HTE Database (66,308 experiments)

- **hte_recommend_tool** - Statistical recommendations with filtering
- **hte_analytics_tool** - Catalyst usage, success rates, metal statistics

#### E. Supporting Analysis

- **list_supported_cores_tool** - Available catalyst options
- **reaction_dataset_analytics_tool** - Frequency and yield patterns

**Output**: Evidence from 8-10 different sources

### Phase 3: Evidence Synthesis

The agent compares and analyzes all gathered data:

1. Creates comparison table across sources
2. Identifies **consensus** (multiple sources agree)
3. Notes **conflicts** (sources disagree)
4. Ranks based on evidence strength + practicality
5. Explains reasoning for rankings

**Output**: Critical analysis of all evidence

### Phase 4: Final Recommendations

Provides 2-3 top recommendations with:

- ✅ Complete reaction setup (specific quantities)
- ✅ Detailed conditions (temp, time, concentration)
- ✅ Evidence citations (which tools support this)
- ✅ Expected outcomes (yield, success rate)
- ✅ Mechanistic rationale
- ✅ Practical considerations (safety, cost, alternatives)
- ✅ Decision tree (when to use which conditions)

**Output**: Actionable, evidence-based protocols

## Usage Examples

### Basic Comprehensive Analysis

```bash
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1"
```

### With More Precedents

```bash
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" \
  --top-k 10 \
  --max-protocols 3
```

### See Agent Thinking Process

```bash
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" \
  --verbose
```

### Use Best Model for Maximum Intelligence

```bash
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" \
  --llm-model gpt-4-turbo \
  --top-k 10 \
  --max-protocols 3 \
  --verbose
```

## Example Output Structure

```
## Phase 1: Reaction Analysis Summary
1. Reaction Family: C-N Coupling (0.9 confidence)
2. Reactant 1: Aryl bromide
3. Reactant 2: Secondary amine (morpholine)
4. Mechanism: Nucleophilic substitution with Cu catalyst

## Phase 2: Evidence from Each Database

### A. Rule Database
- Catalyst: CuI (5-10 mol%)
- Ligand: DMEDA (20-40 mol%)
- Base: Cs2CO3 (2-3 eq)
- Solvent: Toluene
- Temp: 90-110°C

### B. ML Precedent Database
- 25 similar reactions found
- Top core: Cu (75% of precedents)
- Average yield: 78%
- Success rate: 85%

### C. HTE Database (66K experiments)
- Top catalyst: CuI (24,571 experiments)
- Average yield: 20.5%
- Success rate: 16.65%
- Best pair: ArBr + R2NH-a-branch

### D. Protocol Database
[Full experimental protocols if available]

## Phase 3: Evidence Synthesis

| Source    | Catalyst | Base    | Solvent | Yield |
|-----------|----------|---------|---------|-------|
| Rule DB   | CuI      | Cs2CO3  | Toluene | -     |
| ML        | CuI      | K3PO4   | Dioxane | 78%   |
| HTE       | CuI      | Cs2CO3  | DMSO    | 20%   |

**Consensus**: CuI catalyst + carbonate base
**Conflict**: Solvent choice varies (Toluene vs Dioxane vs DMSO)

## Phase 4: Final Recommendations

### Recommendation 1: ML-Optimized (Highest Yield)
- **Catalyst**: CuI (7.5 mol%)
- **Ligand**: DMEDA (30 mol%)
- **Base**: Cs2CO3 (2.5 equivalents)
- **Solvent**: 1,4-Dioxane
- **Temperature**: 100°C
- **Time**: 16 hours
- **Concentration**: 0.2 M
- **Atmosphere**: N2

**Evidence**: 
- ML precedents: 78% average yield (25 reactions)
- Rule database: High confidence match
- HTE: Moderate support (20% yield)

**Expected Outcome**: 70-80% yield

**Rationale**: 
Combines rule-based catalyst selection with ML-optimized conditions.
Dioxane provides better solubility than toluene while maintaining
reactivity. Temperature optimized for balance of rate and selectivity.

**Practical Notes**:
- Cost: Moderate ($50-100/reaction)
- Safety: Handle CuI and Cs2CO3 under inert atmosphere
- Troubleshooting: If low yield, increase temperature to 110°C

### Recommendation 2: Cost-Effective Alternative
[Second option with different trade-offs]

### Decision Tree:
- Use Recommendation 1 if: Yield is critical, budget allows
- Use Recommendation 2 if: Cost is concern, copper preferred
- Alternative: Try Pd catalyst if Cu fails
```

## Comparison: Simple vs Comprehensive

| Aspect | Simple Agent | Comprehensive Agent |
|--------|--------------|---------------------|
| **Tools Used** | 7 | 15+ |
| **Databases Queried** | 2-3 | 5-8 |
| **Analysis Phases** | 1 | 4 |
| **Time** | 5-10 seconds | 15-30 seconds |
| **Cost** | $0.01-0.02 | $0.05-0.15 |
| **Evidence Citations** | Minimal | Extensive |
| **Structural Analysis** | Basic | Deep |
| **Recommendation Quality** | Good | Excellent |
| **Use Case** | Quick screening | Final optimization |

## When to Use Each Mode

### Use Simple Mode (default gpt-4o-mini) when

- ✓ Quick screening needed
- ✓ Cost is a concern
- ✓ Well-studied reaction types
- ✓ Just need a starting point

### Use Comprehensive Mode (gpt-4o/gpt-4-turbo) when

- ✓ Critical reaction for synthesis
- ✓ Need highest success probability
- ✓ Want to understand all options
- ✓ Willing to invest in analysis
- ✓ Novel or challenging transformations

## Cost Estimates

### Per Query (with gpt-4o)

- Input tokens: ~2,000-3,000 ($0.015-0.023)
- Output tokens: ~2,000-4,000 ($0.060-0.120)
- **Total: $0.075-0.143 per analysis**

### Per Query (with gpt-4o-mini)

- Input tokens: ~2,000-3,000 ($0.0003-0.0005)
- Output tokens: ~2,000-4,000 ($0.0012-0.0024)
- **Total: $0.0015-0.0029 per analysis**

For critical reactions, the comprehensive analysis is worth it!

## Tips for Best Results

### 1. Provide Clear SMILES

```bash
# Good: Clear, canonical SMILES
--reaction "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

# Bad: Ambiguous or non-canonical
--reaction "BrPh.PhNH2>>PhNHPh"
```

### 2. Adjust Top-K Based on Novelty

```bash
# Well-studied reactions: fewer precedents needed
--top-k 3

# Novel reactions: get more precedents
--top-k 15
```

### 3. Use Verbose for Learning

```bash
# See which tools the agent uses and why
--verbose
```

### 4. Save Output for Records

```bash
# Save comprehensive analysis to file
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "..." > analysis_2024-11-18.txt
```

## Advanced Options

### Custom Temperature

```bash
# Deterministic (recommended)
--temperature 0.0

# Creative (explore alternatives)
--temperature 0.5
```

### Model Selection

```bash
# Best quality (expensive)
--llm-model gpt-4-turbo

# Balanced (recommended)
--llm-model gpt-4o

# Fast & cheap
--llm-model gpt-4o-mini
```

## Troubleshooting

### Agent Skips Some Tools

**Solution**: Use `--verbose` to see reasoning. May need to adjust system prompt.

### Too Expensive

**Solution**: Use `--llm-model gpt-4o-mini` or limit `--top-k 3`

### Output Too Long

**Solution**: Reduce `--max-protocols 1` or redirect to file

### Some Databases Missing

**Solution**: Check which tools are available with `list_all_families_tool`

## Summary

The comprehensive agent provides:

- ✅ **Systematic** 4-phase analysis methodology
- ✅ **Multi-source** evidence from 5-8 databases
- ✅ **Deep** structural and mechanistic understanding
- ✅ **Critical** comparison and synthesis of evidence
- ✅ **Actionable** specific recommendations with rationale
- ✅ **Practical** considerations for execution

This is now the most thorough automated chemistry analysis tool available!
