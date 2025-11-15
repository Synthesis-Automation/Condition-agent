# HTE Tools - Chem Assistant Integration

## Overview

The HTE (High-Throughput Experimentation) recommendation and analytics tools are now fully integrated into the `chem_assistant` LangChain agent system. This enables LLM-powered chemistry assistants to access **66,308 experimental results** across **41 reaction types** for intelligent condition recommendations and database exploration.

## Integration Summary

### What Was Added

Two new LangChain tools were integrated into `chem_assistant/chemtools_wrapper.py`:

1. **`hte_recommend_tool`**: HTE-based condition recommendations with catalyst filtering
2. **`hte_analytics_tool`**: HTE database analytics for exploring reactant pairs, catalysts, and metal usage

Total tools increased from **23 → 25**.

### Files Modified

#### `chem_assistant/chemtools_wrapper.py`
- **Lines ~120**: Added HTE imports with availability flag
  ```python
  try:
      from chemtools.HTE import HTERecommender, HTEAnalytics, format_result
      HTE_AVAILABLE = True
  except ImportError:
      HTE_AVAILABLE = False
  ```

- **Lines ~720-800**: Added Pydantic input schemas
  - `HTERecommendInput`: 6 fields (reactants, filters, top_k, etc.)
  - `HTEAnalyticsInput`: 9 fields (query_type, filters, sorting, etc.)

- **Lines ~3230-3590**: Added tool function implementations (360 lines)
  - `hte_recommend_tool`: Wraps HTERecommender with comprehensive docstring
  - `hte_analytics_tool`: Wraps HTEAnalytics with 5 query types

- **Lines ~3600**: Registered tools in CHEMTOOLS_TOOLS list

- **Lines 1-30**: Updated module docstring to reflect 25 total tools

## Tool Capabilities

### 1. HTE Recommendation Tool (`hte_recommend_tool`)

**Purpose**: Recommend reaction conditions based on 66K experimental results.

**Key Features**:
- ⚡ Fast: <100ms query time
- 📊 Statistical confidence with success rates & sample sizes
- 🔬 Catalyst filtering (Pd, Cu, Ni)
- 🎯 No reaction SMILES needed - works with just reactants
- 📈 Covers 41 reaction types (Suzuki, C-N coupling, Buchwald-Hartwig, etc.)

**Input Parameters**:
```python
{
    "reactant_a_smiles": str,           # Required: First reactant
    "reactant_b_smiles": Optional[str], # Optional: Second reactant
    "top_k": int = 5,                   # Number of recommendations (1-20)
    "min_experiments": int = 2,         # Minimum experiments per condition
    "reaction_type_filter": Optional[str], # Filter by reaction type
    "catalyst_filter": Optional[str]    # Filter by metal (Pd, Cu, Ni)
}
```

**Example Usage**:
```python
from chem_assistant.chemtools_wrapper import hte_recommend_tool

# Suzuki coupling with Pd catalyst filter
result = hte_recommend_tool.invoke({
    "reactant_a_smiles": "Brc1ccccc1",        # Bromobenzene
    "reactant_b_smiles": "c1ccc(B(O)O)cc1",  # Phenylboronic acid
    "top_k": 3,
    "catalyst_filter": "Pd"
})

# Returns:
{
    "success": True,
    "predicted_reaction_type": "Suzuki",
    "matching_experiments": 1908,
    "recommendations": [
        {
            "rank": 1,
            "catalyst": "dtbpfPdCl2",
            "ligand": "dtbpf",
            "base": "Na2CO3",
            "solvent": "iPrOAc",
            "success_rate": 100.0,
            "avg_yield": 65.5,
            "confidence_score": 72.1,
            "num_experiments": 3
        },
        ...
    ]
}
```

### 2. HTE Analytics Tool (`hte_analytics_tool`)

**Purpose**: Analyze HTE database to explore patterns, catalysts, and reaction statistics.

**Query Types**:
1. **`list_pairs`**: List reactant type pairs with statistics
2. **`catalysts`**: Analyze catalyst usage and performance
3. **`reactions`**: Summarize all reaction types in database
4. **`metals`**: Analyze metal distribution (Pd, Cu, Ni, etc.)
5. **`similar_pairs`**: Find similar reactant combinations

**Input Parameters**:
```python
{
    "query_type": str,                      # Required: Query type (see above)
    "reaction_type": Optional[str],         # Filter by reaction
    "catalyst_filter": Optional[str],       # Filter by metal
    "reactant_a_type": Optional[str],       # Filter by reactant A type
    "reactant_b_type": Optional[str],       # Filter by reactant B type
    "min_experiments": int = 10,            # Minimum experiments threshold
    "top_n": int = 20,                      # Max results to return
    "sort_by": str = "count",              # Sort by 'count' or 'success_rate'
    "similarity_criteria": Optional[str]    # For similar_pairs: 'reaction_type', 'catalyst', 'both'
}
```

**Example Usage**:

#### List Reactant Pairs for Suzuki with Pd Catalysts
```python
result = hte_analytics_tool.invoke({
    "query_type": "list_pairs",
    "reaction_type": "Suzuki",
    "catalyst_filter": "Pd",
    "min_experiments": 50,
    "top_n": 10
})

# Returns:
{
    "success": True,
    "query_type": "list_pairs",
    "total_results": 14,
    "results": [
        {
            "reactant_a": "ArCl",
            "reactant_b": "ArB(OR)2",
            "reaction": "Suzuki",
            "experiments": 2528,
            "avg_yield": 30.4,
            "success_rate": 25.9,
            "top_catalyst": "dtbpfPdCl2"
        },
        ...
    ]
}
```

#### Analyze Metal Usage Across Database
```python
result = hte_analytics_tool.invoke({
    "query_type": "metals",
    "top_n": 5
})

# Returns:
{
    "success": True,
    "query_type": "metals",
    "total_experiments": 66308,
    "metals": [
        {
            "metal": "Pd",
            "experiments": 45205,
            "percentage": 68.2,
            "top_reactions": [
                {"reaction": "Suzuki", "count": 6524},
                {"reaction": "C_N_Coupling", "count": 4215},
                ...
            ]
        },
        ...
    ]
}
```

#### Compare Catalysts for Suzuki Reactions
```python
result = hte_analytics_tool.invoke({
    "query_type": "catalysts",
    "reaction_type": "Suzuki",
    "top_n": 5
})

# Returns:
{
    "success": True,
    "query_type": "catalysts",
    "results": [
        {
            "catalyst": "dtbpfPdCl2",
            "metal": "Pd",
            "experiments": 1635,
            "avg_yield": 35.2,
            "success_rate": 33.5,
            "reaction_types": 12
        },
        ...
    ]
}
```

## Using with LangChain Agent

### Basic Setup

```python
from chem_assistant.chemtools_wrapper import CHEMTOOLS_TOOLS
from langgraph.prebuilt import create_agent
from langchain_openai import ChatOpenAI

# Create LLM
llm = ChatOpenAI(model="gpt-4", temperature=0)

# Create agent with all 25 chemtools
agent = create_agent(llm, CHEMTOOLS_TOOLS)
```

### Example Queries

#### Query 1: Simple Condition Recommendation
```python
result = agent.invoke({
    "messages": [
        "Recommend conditions for Suzuki coupling of bromobenzene and phenylboronic acid using palladium catalysts"
    ]
})
```

**Agent Response**:
```
Based on 1,908 experimental results, I recommend these top conditions for your Suzuki coupling:

1. **dtbpfPdCl2 + Na2CO3 in iPrOAc** (Rank 1)
   - Success rate: 100%, Avg yield: 65.5%
   - Confidence score: 72.1 (3 experiments)

2. **XantPhos Pd(allyl)Cl + Cs2CO3 in Dioxane** (Rank 2)
   - Success rate: 100%, Avg yield: 88.0%
   - Confidence score: 71.0 (2 experiments)

All conditions use palladium catalysts as requested.
```

#### Query 2: Database Exploration
```python
result = agent.invoke({
    "messages": [
        "What are the most commonly used palladium catalysts for Suzuki reactions? Show me the top 5 with their success rates."
    ]
})
```

**Agent Response**:
```
Based on the HTE database analysis, here are the top 5 palladium catalysts for Suzuki reactions:

1. **dtbpfPdCl2**: 1,635 experiments, 33.5% success rate
2. **P(tBu)3 Pd(crotyl)Cl**: 999 experiments, 19.3% success rate
3. **SPhos Pd(crotyl)Cl**: 815 experiments, 34.0% success rate
4. **XPhos Pd(allyl)Cl**: 682 experiments, 28.4% success rate
5. **RuPhos Pd(allyl)Cl**: 574 experiments, 31.2% success rate

dtbpfPdCl2 is the most extensively tested with over 1,600 experiments.
```

#### Query 3: Catalyst Comparison
```python
result = agent.invoke({
    "messages": [
        "Compare palladium vs copper catalysts for C-N coupling reactions. Which one has better success rates?"
    ]
})
```

**Agent Response**:
```
Analyzing C-N coupling reactions in the HTE database:

**Palladium Catalysts:**
- Total experiments: 4,215
- Average success rate: 28.6%
- Top catalyst: XantPhos Pd(allyl)Cl (45.2% success)

**Copper Catalysts:**
- Total experiments: 1,892
- Average success rate: 31.4%
- Top catalyst: Cu(OAc)2 (38.7% success)

Copper catalysts show slightly better success rates for C-N coupling (31.4% vs 28.6%), though palladium has been tested more extensively. For challenging substrates like aryl chlorides, palladium with XantPhos ligand tends to perform better.
```

#### Query 4: Substrate Scope
```python
result = agent.invoke({
    "messages": [
        "What aryl halide + amine combinations have been tested for C-N coupling? Show me combinations with at least 100 experiments."
    ]
})
```

**Agent Response**:
```
Here are the most extensively studied aryl halide + amine combinations for C-N coupling:

1. **ArBr + RNH2**: 1,284 experiments
   - Success rate: 35.2%, Top catalyst: XantPhos Pd(allyl)Cl

2. **ArCl + RNH2**: 892 experiments  
   - Success rate: 22.8%, Top catalyst: XantPhos Pd(allyl)Cl

3. **ArBr + Ar2NH**: 467 experiments
   - Success rate: 41.6%, Top catalyst: P(tBu)3 Pd(crotyl)Cl

4. **ArCl + Ar2NH**: 285 experiments
   - Success rate: 18.3%, Top catalyst: dtbpfPdCl2

Secondary amines (Ar2NH) show higher success rates than primary amines (RNH2), and aryl bromides are generally more reactive than aryl chlorides.
```

## Error Handling

Both tools include graceful error handling:

```python
# If HTE not available
{
    "success": False,
    "error": "HTE tools not available. Install with: pip install chemtools[hte]",
    "recommendation": "Install HTE dependencies"
}

# If invalid query type
{
    "success": False,
    "error": "Unknown query_type 'invalid'. Valid types: list_pairs, catalysts, reactions, metals, similar_pairs"
}

# If missing required parameters
{
    "success": False,
    "error": "similar_pairs query requires both reactant_a_type and reactant_b_type"
}
```

## Testing

A comprehensive integration test is available in `test_hte_integration.py`:

```bash
python test_hte_integration.py
```

**Test Coverage**:
- ✅ Import verification
- ✅ Tool registration check
- ✅ HTE recommendation (Suzuki coupling)
- ✅ HTE analytics (list_pairs, metals, catalysts queries)
- ✅ C-N coupling recommendation
- ✅ Catalyst comparison analytics

**Example Output**:
```
================================================================================
Testing HTE Tools Integration with Chem Assistant
================================================================================

1. Testing imports...
   ✓ Successfully imported HTE tools
   ✓ HTE_AVAILABLE = True
   ✓ Total tools: 25

2. Checking tool registration...
   ✓ hte_recommend_tool registered
   ✓ hte_analytics_tool registered

3. Testing HTE recommendation tool...
   ✓ Recommendation succeeded
   - Predicted reaction: Suzuki
   - Matching experiments: 1908
   - Top condition: dtbpfPdCl2 + Na2CO3 in iPrOAc
     Success rate: 100.0%, Confidence: 72.1

...

✓ All HTE tools successfully integrated into chem_assistant!
```

## Database Statistics

The HTE database powering these tools contains:

- **Total experiments**: 66,308
- **Reaction types**: 41 (Suzuki, C-N coupling, Buchwald-Hartwig, etc.)
- **Reactant combinations**: 71 unique pairs
- **Catalysts**: 
  - Palladium: 45,205 experiments (68.2%), 191 catalyst types
  - Copper: 5,478 experiments (8.3%), 66 catalyst types
  - Nickel: 328 experiments (0.5%)
  - Other metals: 716 experiments (1.1%)

## Performance

- **Recommendation queries**: <100ms
- **Analytics queries**: <200ms
- **Database loading**: ~500ms (cached after first load)
- **Memory footprint**: ~150MB for full database

## Related Documentation

- **HTE Analytics API**: `docs/HTE_ANALYTICS.md` - Complete API reference
- **HTE Quick Reference**: `HTE_ANALYTICS_QUICKREF.md` - Quick start guide
- **Implementation Summary**: `HTE_ANALYTICS_IMPLEMENTATION.md` - Development notes
- **CLI Usage**: `chemtools/HTE/README.md` - Command-line interface

## Future Enhancements

Potential improvements for future versions:

1. **Batch recommendations**: Process multiple reactions in one call
2. **Custom scoring**: User-defined weighting for success rate vs yield
3. **Yield prediction**: ML-based yield estimation for new substrates
4. **Condition clustering**: Group similar conditions for easier comparison
5. **Export to automation format**: Direct output for HTE platforms
6. **Interactive visualization**: Charts and graphs for analytics results

## Support

For issues or questions:
- Check `docs/HTE_ANALYTICS.md` for detailed API documentation
- Run `pytest tests/test_hte_*.py` to verify HTE functionality
- See `demo_hte_analytics.py` for usage examples

---

**Integration Completed**: 2024
**Tools Added**: 2 (hte_recommend_tool, hte_analytics_tool)
**Total Chem Assistant Tools**: 25
**Database Size**: 66,308 experiments across 41 reaction types
