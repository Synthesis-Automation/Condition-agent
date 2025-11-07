# Agent ChemTools Access Audit

## Overview
This document audits which chemtools capabilities are exposed to the LangChain agent vs what exists in the chemtools library but is not yet wrapped.

---

## ✅ Currently Exposed Tools (17 total)

The agent has access to the following chemtools functionality via `CHEMTOOLS_TOOLS`:

### 1. SMILES & Normalization
- ✅ `normalize_smiles_tool` - Canonicalize SMILES strings
- ✅ `normalize_reaction_tool` - Canonicalize reaction SMILES

### 2. Reaction Analysis
- ✅ `detect_reaction_family_tool` - Detect reaction family/type
- ✅ `analyze_bond_changes_tool` - Analyze bond breaking/formation (NEW)

### 3. Molecule Analysis
- ✅ `classify_reactant_tool` - Classify reactant type (aryl halide, amine, etc.)
- ✅ `get_functional_groups_tool` - Detect functional groups (80+ groups)
- ✅ `calculable_features_tool` - Evaluate 251 curated features
- ✅ `molpipeline_featurize_tool` - Generate molecular features with optional MolPipeline vectors

### 4. Recommendation Tools
- ✅ `recommend_conditions_tool` - ML-based condition recommendations
- ✅ `rule_based_conditions_tool` - Deterministic rule-engine guidance
- ✅ `enhanced_cross_family_recommend_tool` - Cross-family precedent search with mechanism filters
- ✅ `search_precedents_tool` - Search for similar precedent reactions

### 5. Database & Analytics
- ✅ `reaction_dataset_analytics_tool` - Analyze dataset frequency/yield statistics
- ✅ `find_reagent_tool` - Look up reagent information from database
- ✅ `reagent_database_analytics_tool` - Summarize reagent registry statistics
- ✅ `list_supported_cores_tool` - Enumerate catalyst cores in precedents
- ✅ `add_reagent_tool` - Insert or preview reagent taxonomy entries

---

## ❌ Missing Tools (Not Yet Exposed)

The following chemtools capabilities exist but are **NOT** exposed to the agent:

### 1. Protocol-Based Recommendations ⚠️ HIGH PRIORITY
**Module**: `chemtools.protocol`

**What it does**:
- Matches user reactions to standard experimental protocols using DRFP similarity
- Returns detailed procedure information (not just conditions)
- Includes step-by-step instructions, workup procedures, purification methods
- Can filter by reaction family, tags, and SMARTS patterns

**Why it's useful**:
- Provides complete experimental procedures, not just reagents/conditions
- Complementary to ML-based recommendations (protocols vs conditions)
- Users can get actual literature procedures for similar reactions

**Key functions to wrap**:
```python
from chemtools.protocol import ProtocolRecommender, recommend_protocol

recommender = ProtocolRecommender()
results = recommender.recommend(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=5,
    use_smarts_filter=True,
    family_filter='Suzuki'
)
```

**Suggested tool name**: `protocol_recommendation_tool`

---

### 2. Reaction Similarity & Comparison
**Module**: `chemtools.reaction_similarity`

**What it does**:
- Compute DRFP Tanimoto similarity between two reactions
- Compare reaction fingerprints
- Batch similarity calculations

**Why it's useful**:
- Users can ask "How similar are these two reactions?"
- Can compare their reaction to literature precedents
- Useful for analog analysis

**Suggested tool name**: `reaction_similarity_tool`

---

### 3. Advanced Constraint Filtering
**Module**: `chemtools.constraints`

**What it does**:
- Green chemistry filters
- Inventory-based filtering
- Custom constraint rules beyond what's in recommendation tools

**Why it's useful**:
- Users can ask "Avoid HMPA" or "Only use reagents we have"
- More sophisticated filtering than current constraint parser

**Note**: Partially exposed via `constraint_rules` parameter in recommendation tools, but not as standalone tool

---

### 4. Dataset Statistics & Analytics
**Module**: `chemtools.dataset_analytics`

**Exposed functions**:
- ✅ Via `reaction_dataset_analytics_tool`

**Missing functions**:
- ❌ `get_condition_cores()` - Get all catalyst cores for a family
- ❌ `get_all_families()` - List all reaction families in dataset
- ❌ Direct family-specific queries

**Suggested enhancement**: Expand `reaction_dataset_analytics_tool` or create separate tools

---

### 5. Taxonomy & Registry Deep Queries
**Module**: `chemtools.taxonomy`, `chemtools.registry`

**What it does**:
- Deep queries into compound taxonomy
- CAS number resolution
- Reagent property lookups beyond simple find

**Why it's useful**:
- Users can ask detailed questions about reagent properties
- Cross-reference multiple databases
- Validate chemical identifiers

**Note**: Partially exposed via `find_reagent_tool`, but limited

---

### 6. Explanation & Reasoning
**Module**: `chemtools.explain`

**What it does**:
- Generate human-readable justifications for recommendations
- Explain why certain conditions were chosen
- Provide reasoning based on precedents

**Why it's useful**:
- Transparency in recommendations
- Educational for users
- Better trust in AI suggestions

**Note**: Currently used internally by `recommend_from_reaction()` but not directly accessible

**Suggested tool name**: `explain_recommendation_tool`

---

### 7. MCP (Model Context Protocol) Integration
**Module**: `chemtools.integrations.mcp`

**What it does**:
- Integration utilities for MCP
- Client/server communication
- Schema validation

**Why it's useful**:
- For advanced integrations
- Standardized protocol communication

**Note**: May not need direct tool exposure

---

## 📊 Coverage Assessment

### Core Functionality Coverage
| Category | Exposed | Available | Coverage |
|----------|---------|-----------|----------|
| SMILES/Normalization | 2 | 2 | 100% ✅ |
| Molecule Analysis | 4 | 5 | 80% ✅ |
| Reaction Analysis | 2 | 3 | 67% ⚠️ |
| Recommendations | 4 | 5 | 80% ⚠️ |
| Database/Analytics | 5 | 6 | 83% ✅ |
| Protocols | 0 | 1 | 0% ❌ |
| Utilities | 0 | 3 | 0% ⚠️ |

**Overall Coverage**: ~70% of major functionality

---

## 🎯 Recommendations for Enhancement

### Priority 1: Add Protocol Recommendation (HIGH IMPACT)
This is the biggest missing piece. Protocol recommendations provide complete experimental procedures, which is what users often need.

**Estimated effort**: 2-3 hours
**Impact**: HIGH - Adds major new capability

### Priority 2: Add Reaction Similarity Tool (MEDIUM IMPACT)
Useful for comparison questions and analog analysis.

**Estimated effort**: 1 hour
**Impact**: MEDIUM - Nice-to-have for advanced users

### Priority 3: Enhance Explanation Capabilities (MEDIUM IMPACT)
Make the explain module directly accessible so agent can provide reasoning.

**Estimated effort**: 1-2 hours
**Impact**: MEDIUM - Improves transparency and trust

### Priority 4: Expand Analytics Tools (LOW IMPACT)
Add missing dataset analytics functions.

**Estimated effort**: 1 hour
**Impact**: LOW - Minor quality-of-life improvement

---

## 🔍 Tool Design Suggestions

### Protocol Recommendation Tool
```python
@tool(args_schema=ProtocolRecommendInput)
def protocol_recommendation_tool(
    reaction_smiles: str,
    k: int = 5,
    family_filter: Optional[str] = None,
    use_smarts_filter: bool = True,
    min_similarity: float = 0.5
) -> Dict[str, Any]:
    """
    Find experimental protocols for reactions similar to the query.
    
    Returns complete procedure information including:
    - Step-by-step instructions
    - Reagent amounts and equivalents
    - Reaction conditions (temperature, time, atmosphere)
    - Workup and purification procedures
    - Literature references
    
    Args:
        reaction_smiles: Reaction SMILES to match
        k: Number of protocols to return (default 5)
        family_filter: Optional reaction family filter (e.g. "Suzuki", "Buchwald_CN")
        use_smarts_filter: Enable SMARTS-based structural filtering
        min_similarity: Minimum DRFP similarity threshold (0.0-1.0)
    
    Returns:
        Dict with ranked protocol recommendations
    """
    from chemtools.protocol import ProtocolRecommender
    
    recommender = ProtocolRecommender()
    return recommender.recommend(
        reaction_smiles=reaction_smiles,
        k=k,
        family_filter=family_filter,
        use_smarts_filter=use_smarts_filter,
        min_similarity=min_similarity
    )
```

### Reaction Similarity Tool
```python
@tool(args_schema=ReactionSimilarityInput)
def reaction_similarity_tool(
    reaction1_smiles: str,
    reaction2_smiles: str
) -> Dict[str, Any]:
    """
    Calculate DRFP-based similarity between two reactions.
    
    Returns Tanimoto coefficient (0.0-1.0) where:
    - 1.0 = identical reactions
    - 0.7+ = very similar
    - 0.5-0.7 = moderately similar
    - <0.5 = different
    
    Args:
        reaction1_smiles: First reaction SMILES
        reaction2_smiles: Second reaction SMILES
    
    Returns:
        Dict with similarity score and interpretation
    """
    from chemtools.reaction_similarity import compute_similarity
    
    similarity = compute_similarity(reaction1_smiles, reaction2_smiles)
    
    if similarity >= 0.7:
        interpretation = "Very similar reactions"
    elif similarity >= 0.5:
        interpretation = "Moderately similar reactions"
    else:
        interpretation = "Different reactions"
    
    return {
        "success": True,
        "reaction1": reaction1_smiles,
        "reaction2": reaction2_smiles,
        "similarity": round(similarity, 4),
        "interpretation": interpretation
    }
```

---

## 📝 Summary

The agent currently has access to **~70% of core chemtools functionality**. The main gaps are:

1. ❌ **Protocol-based recommendations** (complete experimental procedures)
2. ❌ **Reaction similarity comparisons** (DRFP-based)
3. ⚠️ **Limited explanation/reasoning** capabilities
4. ⚠️ **Some advanced analytics** functions

The existing 17 tools cover the most critical workflows:
- ✅ Molecule/reaction analysis
- ✅ ML-based condition recommendations
- ✅ Rule-based recommendations
- ✅ Precedent search
- ✅ Database lookups
- ✅ Feature detection

**Recommendation**: Add protocol recommendation tool as Priority 1 enhancement to significantly expand agent capabilities.

---

**Last Updated**: 2024
**Agent Tools Version**: v2 (17 tools)
**ChemTools Version**: Current main-v2 branch
