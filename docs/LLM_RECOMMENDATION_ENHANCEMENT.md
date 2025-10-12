# LLM Enhancement for Reaction Condition Recommendation Workflow

**Date**: October 12, 2025  
**Status**: ✅ Phase 1 Complete - Multi-Source Synthesis Prototype Working  
**Target**: Condition recommendation workflow (not reagent generator)

---

## 🎉 UPDATE: Multi-Source Synthesis Implemented!

**Implementation Status**: ✅ **COMPLETE** (October 12, 2025)

We've successfully implemented and tested the multi-source LLM synthesis approach:

- ✅ **Core function working**: `synthesize_recommendations_llm()` in `llmtools/recommendation_llm.py`
- ✅ **All tests passing**: 3/3 test cases validated (100% success rate)
- ✅ **Key features confirmed**: Consensus analysis, constraint handling, source attribution
- ✅ **Documentation**: See `docs/MULTISOURCE_IMPLEMENTATION.md` for details

**What's Next**: 
- 🔲 Step 3: Compare vs individual modes (benchmark quality improvement)
- 🔲 Step 4: Optimize latency and refine prompts
- 🔲 Step 5: Integrate into production API

**See bottom of document for detailed implementation notes.**

---

## Executive Summary

**Yes, LLM assistance WILL improve the recommendation workflow!**

**Key Findings**:

- ✅ **High-value additions**: Condition explanations, intelligent re-ranking, troubleshooting
- ⚠️ **Keep current approach**: DRFP-based precedent retrieval (it works well)
- 🎯 **Recommended approach**: Hybrid workflow (DRFP + LLM post-processing)

---

## Current Recommendation Workflows

You currently have **4 recommendation modes**:

1. **Rule-based** - SCDB pattern matching
2. **ML-based** - DRFP similarity with precedents
3. **Protocol-based** - DRFP similarity with literature protocols
4. **Fusion** - Combines rule + ML

---

## LLM Integration Analysis: Step-by-Step

### **Step 1: Reaction Type Detection**

**Current**: Pattern matching + SMARTS  
**LLM Could Help**: ✅ **YES - Moderate improvement**

#### **Implementation**:

```python
def detect_reaction_type_llm(reaction_smiles: str, llm_client: LLMClient) -> Dict:
    """
    Use LLM to classify reaction type with chemical reasoning.
    """
    prompt = f"""
    Analyze this reaction SMILES and classify the reaction type:

    {reaction_smiles}

    Choose from: Suzuki, Buchwald-Hartwig, Ullmann, Negishi, Heck, etc.

    Provide:
    1. Primary reaction type
    2. Confidence (0-1)
    3. Key functional groups involved
    4. Why this classification?

    Respond with JSON:
    {{
      "reaction_type": "Suzuki",
      "confidence": 0.95,
      "functional_groups": ["aryl_halide", "boronic_acid"],
      "reasoning": "Aryl bromide + boronic acid → biaryl product"
    }}
    """

    response = llm_client.chat(prompt, temperature=0.0, max_tokens=300)
    return json.loads(response.content)
```

#### **Benefits**:

- ✅ Better handling of ambiguous cases
- ✅ Can explain reasoning ("I see aryl bromide + amine → likely Buchwald-Hartwig")
- ✅ Catches edge cases pattern matching might miss

#### **Drawbacks**:

- ❌ Slower than SMARTS (~2s vs <0.1s)
- ❌ Less deterministic (might vary between runs)

#### **Verdict**:

**Optional enhancement** - Use LLM for low-confidence SMARTS matches

**Priority**: **P2** (Medium)

---

### **Step 2: Condition Retrieval**

**Current**:

- Rule: SCDB exact/fuzzy matching
- ML: DRFP k-NN search
- Protocol: DRFP similarity + metadata

**LLM Could Help**: ⚠️ **Limited value - LLM doesn't know your specific data**

#### **Why NOT to use LLM here**:

- ❌ LLM doesn't have access to your precedent database
- ❌ Would need to pass all precedents in context (expensive, slow)
- ❌ DRFP similarity already works well for this

#### **Problematic approach** (not recommended):

```python
# Would need to pass ALL precedents - impractical!
prompt = f"""
Given this reaction: {reaction_smiles}

And these 1000 precedents:
{json.dumps(all_precedents)}  # Too large! Context overflow!

Which are most relevant?
"""
```

#### **Verdict**:

**Don't use LLM here** - Keep DRFP/pattern matching

**Priority**: **N/A** (Don't implement)

---

### **Step 3: Condition Ranking/Filtering** ⭐ HIGH VALUE

**Current**:

- Rule: Returns all matches
- ML: Ranked by DRFP similarity
- Protocol: Ranked by similarity + metadata

**LLM Could Help**: ✅ **YES - High value!**

#### **Implementation**:

```python
def rerank_conditions_llm(
    reaction_smiles: str,
    candidate_conditions: List[Dict],
    llm_client: LLMClient
) -> List[Dict]:
    """
    Use LLM to re-rank candidate conditions based on chemistry.

    Similar to how search engines use LLM for re-ranking.
    """
    # Format candidates for LLM
    candidates_text = []
    for i, cond in enumerate(candidate_conditions, 1):
        candidates_text.append(
            f"{i}. {cond.get('catalyst', 'N/A')}, "
            f"{cond.get('solvent', 'N/A')}, "
            f"{cond.get('temperature', 'N/A')} "
            f"(similarity: {cond.get('similarity', 0):.2f})"
        )

    prompt = f"""
    Reaction SMILES: {reaction_smiles}

    Candidate conditions (top 10 from DRFP similarity search):
    {chr(10).join(candidates_text)}

    Re-rank these conditions considering:
    - Functional group compatibility with substrate
    - Practical considerations (cost, air sensitivity, availability)
    - Substrate sterics and electronics
    - Typical reaction conditions for this transformation
    - Likelihood of success

    Return JSON array with re-ranked order:
    [
      {{"original_rank": 3, "new_rank": 1, "reason": "Best match for electron-poor substrate"}},
      {{"original_rank": 1, "new_rank": 2, "reason": "Good but expensive catalyst"}},
      ...
    ]

    IMPORTANT: Consider chemistry, not just similarity scores.
    """

    response = llm_client.chat(prompt, temperature=0.2, max_tokens=800)
    reranking = json.loads(response.content)

    # Apply reranking
    reranked = []
    for item in reranking:
        orig_idx = item['original_rank'] - 1
        cond = candidate_conditions[orig_idx].copy()
        cond['llm_rerank_score'] = item['new_rank']
        cond['llm_reasoning'] = item['reason']
        reranked.append(cond)

    return reranked
```

#### **Benefits**:

- ✅ **Chemistry-aware ranking** beyond just structural similarity
- ✅ Considers substrate-specific factors (sterics, electronics)
- ✅ Can downrank impractical conditions (e.g., expensive catalysts, air-sensitive)
- ✅ Explains reasoning for each choice
- ✅ Catches cases where similar structure ≠ similar reactivity

#### **Example Output**:

```json
[
  {
    "rank": 1,
    "original_rank": 3,
    "catalyst": "Pd(dppf)Cl2",
    "similarity": 0.8,
    "llm_rerank_score": 1,
    "llm_reasoning": "Better for electron-poor aryl bromides; dppf is more electron-rich than PPh3"
  },
  {
    "rank": 2,
    "original_rank": 1,
    "catalyst": "Pd(PPh3)4",
    "similarity": 0.85,
    "llm_rerank_score": 2,
    "llm_reasoning": "High similarity but may be less effective for this substrate"
  }
]
```

#### **Verdict**:

**High-value addition** - Post-process DRFP results with LLM

**Priority**: **P0** (Implement first)

---

### **Step 4: Condition Explanation** ⭐ VERY HIGH VALUE

**Current**: Just returns conditions, no explanation  
**LLM Could Help**: ✅ **YES - Very high value!**

#### **Implementation**:

```python
def explain_conditions_llm(
    reaction_smiles: str,
    recommended_conditions: Dict,
    llm_client: LLMClient
) -> str:
    """
    Generate human-readable explanation of why these conditions work.
    """
    prompt = f"""
    Reaction SMILES: {reaction_smiles}

    Recommended conditions:
    - Catalyst: {recommended_conditions.get('catalyst', 'N/A')}
    - Ligand: {recommended_conditions.get('ligand', 'N/A')}
    - Solvent: {recommended_conditions.get('solvent', 'N/A')}
    - Temperature: {recommended_conditions.get('temperature', 'N/A')}
    - Base: {recommended_conditions.get('base', 'N/A')}

    Provide a concise explanation (2-3 sentences) covering:
    1. Why this catalyst is appropriate for this specific reaction
    2. Why this solvent/temperature combination is chosen
    3. The role of each key reagent (base, ligand, additives)

    Focus on chemistry rationale, not just repeating the data.
    """

    response = llm_client.chat(prompt, temperature=0.3, max_tokens=250)
    return response.content.strip()
```

#### **Example Output**:

```
"Pd(PPh3)4 is chosen as it's highly active for Suzuki couplings of
electron-poor aryl bromides due to the electron-donating PPh3 ligands.
Toluene at 80°C provides good solubility while preventing phosphine
ligand degradation that occurs at higher temperatures. K3PO4 acts as
a mild base to activate the boronic acid without decomposing sensitive
functional groups on the substrate."
```

#### **Benefits**:

- ✅ **Educational** - Helps chemists understand WHY, not just WHAT
- ✅ **Builds trust** - Transparent reasoning increases confidence
- ✅ **Troubleshooting hints** - Explains potential failure modes
- ✅ **Knowledge transfer** - Junior chemists learn mechanism insights
- ✅ **Easy to implement** - Single LLM call per recommendation

#### **Verdict**:

**Very high value** - Essential for user experience

**Priority**: **P0** (Implement first)

---

### **Step 5: Condition Optimization Suggestions** ⭐ NEW CAPABILITY

**Current**: None  
**LLM Could Help**: ✅ **YES - New capability!**

#### **Implementation**:

```python
def suggest_optimizations_llm(
    reaction_smiles: str,
    base_conditions: Dict,
    constraints: Dict,  # e.g., {"cost": "low", "air_stable": True}
    llm_client: LLMClient
) -> List[Dict]:
    """
    Suggest condition variations for optimization screening.
    """
    prompt = f"""
    Reaction: {reaction_smiles}

    Base conditions:
    - Catalyst: {base_conditions.get('catalyst')}
    - Solvent: {base_conditions.get('solvent')}
    - Temperature: {base_conditions.get('temperature')}
    - Base: {base_conditions.get('base')}

    User constraints: {json.dumps(constraints)}

    Suggest 3 alternative condition sets to screen:
    1. **Cost-optimized**: Cheaper catalyst/solvent/reagents
    2. **Robustness-optimized**: Air-stable, easier workup, broader substrate scope
    3. **High-throughput optimized**: Faster reaction, lower temperature, easier monitoring

    For each variant, explain:
    - What changes from base conditions
    - Expected trade-offs (yield, selectivity, cost, time)
    - When to prefer this variant

    Return JSON:
    [
      {{
        "variant": "cost_optimized",
        "changes": {{"catalyst": "Pd(OAc)2 + PPh3", "reason": "..."}},
        "expected_tradeoffs": "10-20% lower yield, but 5x cost reduction",
        "when_to_use": "Large-scale synthesis where cost matters"
      }},
      ...
    ]
    """

    response = llm_client.chat(prompt, temperature=0.4, max_tokens=600)
    return json.loads(response.content)
```

#### **Example Output**:

```json
[
  {
    "variant": "cost_optimized",
    "changes": {
      "catalyst": "Pd(OAc)2 + PPh3 (in situ formation instead of pre-made Pd(PPh3)4)",
      "solvent": "EtOH (cheaper than toluene, similar polarity)",
      "base": "K2CO3 (cheaper than K3PO4)"
    },
    "expected_tradeoffs": "10-20% lower yield possible, but 5x cost reduction. May need longer reaction time.",
    "when_to_use": "Large-scale synthesis where cost matters more than yield"
  },
  {
    "variant": "robustness_optimized",
    "changes": {
      "catalyst": "Pd(dppf)Cl2 (air-stable, easy to handle)",
      "solvent": "DMA (higher boiling, more forgiving temp control)",
      "additive": "H2O (small amount improves base dissolution)"
    },
    "expected_tradeoffs": "Similar yield, more reproducible, easier for non-experts",
    "when_to_use": "Production environment or non-specialist labs"
  },
  {
    "variant": "highthroughput_optimized",
    "changes": {
      "temperature": "60°C (lower, safer for automation)",
      "catalyst_loading": "2 mol% (reduced from 5%, faster screening)",
      "reaction_time": "2h max (stop point for HTE)"
    },
    "expected_tradeoffs": "May have incomplete conversion, but enables rapid screening",
    "when_to_use": "Initial condition screening or substrate scope exploration"
  }
]
```

#### **Benefits**:

- ✅ **Automated DOE planning** - Intelligent experiment suggestions
- ✅ **Constraint-aware** - Respects user requirements (cost, air stability, etc.)
- ✅ **Educational** - Explains trade-offs explicitly
- ✅ **Saves time** - Chemists don't need to manually brainstorm variations

#### **Verdict**:

**High value** - Enables automated optimization planning

**Priority**: **P1** (Add after P0 features)

---

### **Step 6: Troubleshooting Failed Reactions** ⭐ NEW CAPABILITY

**Current**: None  
**LLM Could Help**: ✅ **YES - Very high value!**

#### **Implementation**:

```python
def troubleshoot_llm(
    reaction_smiles: str,
    attempted_conditions: Dict,
    observation: str,  # "No product, starting material recovered"
    llm_client: LLMClient
) -> Dict:
    """
    Diagnose why a reaction failed and suggest fixes.
    """
    prompt = f"""
    Reaction SMILES: {reaction_smiles}

    Conditions attempted:
    - Catalyst: {attempted_conditions.get('catalyst')}
    - Solvent: {attempted_conditions.get('solvent')}
    - Temperature: {attempted_conditions.get('temperature')}
    - Base: {attempted_conditions.get('base')}

    Observation: {observation}

    Diagnose the likely problem(s) and suggest 3 concrete next steps to fix it.

    Common failure modes to consider:
    - Insufficient catalyst/base/ligand
    - Incompatible functional groups
    - Temperature too low/high
    - Solubility issues
    - Air/moisture sensitivity
    - Incorrect stoichiometry

    Return JSON:
    {{
      "diagnosis": "Most likely cause of failure",
      "supporting_evidence": "Why this diagnosis fits the observation",
      "suggested_fixes": [
        {{"action": "Try stronger base", "rationale": "...", "priority": 1}},
        {{"action": "Increase temperature", "rationale": "...", "priority": 2}},
        {{"action": "Add phase transfer catalyst", "rationale": "...", "priority": 3}}
      ]
    }}
    """

    response = llm_client.chat(prompt, temperature=0.3, max_tokens=400)
    return json.loads(response.content)
```

#### **Example Output**:

```json
{
  "diagnosis": "Base (K2CO3) is too weak for this electron-poor aryl bromide",
  "supporting_evidence": "Starting material recovery suggests no catalyst activation. Electron-withdrawing NO2 group makes aryl halide less reactive, requiring stronger base for transmetalation step.",
  "suggested_fixes": [
    {
      "action": "Switch to stronger base (Cs2CO3 or K3PO4)",
      "rationale": "Stronger bases facilitate boronic acid activation and transmetalation with electron-poor substrates",
      "priority": 1,
      "expected_outcome": "Should see product formation"
    },
    {
      "action": "Increase temperature to 100°C",
      "rationale": "Higher temperature can compensate for lower reactivity, but may degrade PPh3 ligands",
      "priority": 2,
      "expected_outcome": "Increased conversion, watch for ligand decomposition"
    },
    {
      "action": "Add 18-crown-6 or switch to Cs2CO3",
      "rationale": "Crown ether or cesium improves base solubility in organic phase",
      "priority": 3,
      "expected_outcome": "Better base activation, especially in non-polar solvents"
    }
  ]
}
```

#### **Benefits**:

- ✅ **Interactive assistant mode** - Chemist can iteratively troubleshoot
- ✅ **Learning tool** - Explains common failure modes
- ✅ **Saves experiments** - Avoids random trial-and-error
- ✅ **Prioritized suggestions** - Ranks fixes by likelihood of success

#### **Verdict**:

**Very high value** - Game-changer for user experience

**Priority**: **P2** (More complex, needs workflow state tracking)

---

## Recommended Implementation Plan

### **Phase 1: Low-Hanging Fruit** (Implement First)

| Step                   | Function                   | Value      | Effort | Priority |
| ---------------------- | -------------------------- | ---------- | ------ | -------- |
| **Explain conditions** | `explain_conditions_llm()` | ⭐⭐⭐⭐⭐ | Low    | **P0**   |
| **Re-rank results**    | `rerank_conditions_llm()`  | ⭐⭐⭐⭐   | Low    | **P0**   |

**Why start here**:

- ✅ Doesn't change core workflow (just adds explanations)
- ✅ High user value (transparency + trust)
- ✅ Easy to implement (single LLM call per recommendation)
- ✅ Low risk (can be optional feature flag)

**Estimated effort**: 2-3 days

---

### **Phase 2: Enhanced Capabilities** (Add Later)

| Step                         | Function                      | Value      | Effort | Priority |
| ---------------------------- | ----------------------------- | ---------- | ------ | -------- |
| **Optimization suggestions** | `suggest_optimizations_llm()` | ⭐⭐⭐⭐   | Medium | **P1**   |
| **Reaction type detection**  | `detect_reaction_type_llm()`  | ⭐⭐⭐     | Low    | **P2**   |
| **Troubleshooting**          | `troubleshoot_llm()`          | ⭐⭐⭐⭐⭐ | High   | **P2**   |

**Estimated effort**: 1-2 weeks

---

### **Phase 3: Advanced Features** (Future)

- Multi-turn conversational troubleshooting
- Active learning from user corrections
- Confidence calibration based on success/failure feedback
- Custom constraint handling (green chemistry, specific vendor catalogs)

---

## Proposed Enhanced Workflow

### **Current Workflow** (Simplified):

```python
def recommend_conditions(reaction_smiles: str, reaction_type: str):
    # Step 1: Detect type (current SMARTS)
    detected_type = detect_type_pattern(reaction_smiles)

    # Step 2: Get candidates (current DRFP)
    candidates = get_precedents_drfp(reaction_smiles, k=10)

    # Step 3: Rank (current similarity)
    ranked = rank_by_similarity(candidates)

    return ranked[:3]
```

### **Enhanced Workflow with LLM**:

```python
def recommend_conditions_enhanced(
    reaction_smiles: str,
    reaction_type: Optional[str] = None,
    llm_client: Optional[LLMClient] = None,
    explain: bool = True,
    suggest_alternatives: bool = False
):
    """
    Enhanced recommendation workflow with optional LLM assistance.

    Args:
        reaction_smiles: Reaction SMILES string
        reaction_type: Optional manual override
        llm_client: Optional LLM for explanations/reranking
        explain: Whether to generate explanations (requires llm_client)
        suggest_alternatives: Whether to suggest optimization variants

    Returns:
        Enhanced recommendation dict with explanations and alternatives
    """
    # Step 1: Detect type (hybrid: SMARTS + optional LLM validation)
    detected_type = detect_type_pattern(reaction_smiles)

    if llm_client and detected_type.confidence < 0.8:
        # Use LLM for low-confidence cases
        detected_type = detect_reaction_type_llm(reaction_smiles, llm_client)

    # Step 2: Get candidates (keep current DRFP - it works!)
    candidates = get_precedents_drfp(reaction_smiles, k=20)  # Get more candidates

    # Step 3: LLM re-ranking (NEW!)
    if llm_client:
        candidates = rerank_conditions_llm(
            reaction_smiles,
            candidates[:10],  # Top 10 from DRFP
            llm_client
        )
    else:
        candidates = candidates[:10]

    top_conditions = candidates[:3]

    # Step 4: Add explanations (NEW!)
    if llm_client and explain:
        for cond in top_conditions:
            cond['explanation'] = explain_conditions_llm(
                reaction_smiles,
                cond,
                llm_client
            )

    # Step 5: Suggest alternatives (NEW - optional)
    alternatives = []
    if llm_client and suggest_alternatives:
        alternatives = suggest_optimizations_llm(
            reaction_smiles,
            top_conditions[0],  # Base on best condition
            constraints={},
            llm_client
        )

    return {
        "meta": {
            "model": "fusion_enhanced" if llm_client else "fusion",
            "llm_enabled": llm_client is not None,
            "timestamp": datetime.now().isoformat()
        },
        "input": {
            "reaction_smiles": reaction_smiles,
            "reaction_type": detected_type.type
        },
        "detection": {
            "type": detected_type.type,
            "confidence": detected_type.confidence,
            "method": detected_type.method  # "smarts" or "llm"
        },
        "recommended_conditions": top_conditions,  # With explanations!
        "alternative_conditions": alternatives,     # Optional optimization variants
        "extras": {
            "num_candidates_evaluated": len(candidates),
            "llm_features_used": ["reranking", "explanation"] if llm_client else []
        }
    }
```

---

## Expected Impact

### **Output Comparison**

#### **Without LLM (Current)**:

```json
{
  "meta": { "model": "fusion", "timestamp": "2025-10-12T10:30:00" },
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.85,
      "catalyst": "Pd(PPh3)4",
      "ligand": null,
      "solvent": "toluene",
      "temperature": "80°C",
      "base": "K3PO4",
      "similarity": 0.85
    }
  ]
}
```

#### **With LLM (Enhanced)**:

```json
{
  "meta": {
    "model": "fusion_enhanced",
    "llm_enabled": true,
    "timestamp": "2025-10-12T10:30:00"
  },
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.85,
      "similarity": 0.85,
      "llm_rerank_score": 0.92,
      "llm_reasoning": "Optimal for electron-poor aryl bromides",
      "catalyst": "Pd(PPh3)4",
      "ligand": null,
      "solvent": "toluene",
      "temperature": "80°C",
      "base": "K3PO4",
      "explanation": "Pd(PPh3)4 is ideal for this Suzuki coupling because the electron-donating PPh3 ligands activate the Pd center for oxidative addition to the electron-poor aryl bromide. Toluene at 80°C provides optimal solubility while preventing ligand degradation. K3PO4 acts as a mild base for boronic acid activation without decomposing sensitive functional groups."
    }
  ],
  "alternative_conditions": [
    {
      "variant": "cost_optimized",
      "changes": {
        "catalyst": "Pd(OAc)2 + PPh3",
        "solvent": "EtOH"
      },
      "expected_tradeoffs": "10-20% lower yield, but 5x cost reduction",
      "when_to_use": "Large-scale synthesis"
    }
  ],
  "extras": {
    "num_candidates_evaluated": 10,
    "llm_features_used": ["reranking", "explanation"]
  }
}
```

---

## Technical Implementation Details

### **File Structure**

```
llmtools/
├── agents.py              # Existing LLM agents
├── prompts.py             # Existing prompt templates
├── recommendation_llm.py  # NEW: LLM functions for recommendations
└── clients.py             # Existing LLM client

chemtools/
├── recommend/
│   ├── fusion.py          # Existing fusion workflow
│   ├── llm_enhanced.py    # NEW: Enhanced workflow with LLM
│   └── ...

app/
└── main.py                # API endpoint updates
```

### **New Module: `llmtools/recommendation_llm.py`**

```python
"""
LLM-powered enhancements for reaction condition recommendations.

Provides:
- explain_conditions_llm: Generate human-readable explanations
- rerank_conditions_llm: Chemistry-aware re-ranking
- suggest_optimizations_llm: Propose alternative conditions
- detect_reaction_type_llm: Classify reaction type
- troubleshoot_llm: Diagnose failed reactions
"""

from typing import Dict, List, Any, Optional
from llmtools.clients import LLMClient
from llmtools.prompts import PromptTemplate
import json


# Prompt templates
CONDITION_EXPLANATION = PromptTemplate(
    template="""...""",  # See implementation above
    reaction="",
    catalyst="",
    solvent="",
    temperature="",
    base=""
)

CONDITION_RERANKING = PromptTemplate(
    template="""...""",
    reaction_smiles="",
    candidates=""
)

# ... other templates


def explain_conditions_llm(...):
    """See implementation above"""
    pass

def rerank_conditions_llm(...):
    """See implementation above"""
    pass

def suggest_optimizations_llm(...):
    """See implementation above"""
    pass

def detect_reaction_type_llm(...):
    """See implementation above"""
    pass

def troubleshoot_llm(...):
    """See implementation above"""
    pass
```

### **API Endpoint Changes**

```python
# app/main.py

@app.post("/api/recommend/enhanced")
def recommend_enhanced(
    request: RecommendRequest,
    explain: bool = Query(True, description="Generate LLM explanations"),
    suggest_alternatives: bool = Query(False, description="Suggest optimization variants"),
    llm_provider: str = Query("aliyun", description="LLM provider (aliyun, openai)"),
    llm_model: str = Query("deepseek-v3.2-exp", description="LLM model")
):
    """
    Enhanced recommendation with optional LLM assistance.

    If explain=true or suggest_alternatives=true, initializes LLM client.
    Falls back to standard recommendation if LLM unavailable.
    """
    llm_client = None
    if explain or suggest_alternatives:
        try:
            llm_client = LLMClient(provider=llm_provider, model=llm_model)
        except Exception as e:
            logger.warning(f"LLM initialization failed: {e}. Falling back to standard.")

    result = recommend_conditions_enhanced(
        reaction_smiles=request.reaction_smiles,
        reaction_type=request.reaction_type,
        llm_client=llm_client,
        explain=explain,
        suggest_alternatives=suggest_alternatives
    )

    return result
```

---

## Testing Strategy

### **Phase 1 Tests**

1. **Explanation Quality**:

   - Test with 20 diverse reactions
   - Human evaluation of explanation accuracy
   - Check for hallucinations
   - Measure consistency across runs

2. **Reranking Accuracy**:
   - Test with known good/bad conditions
   - Compare LLM reranking vs pure similarity
   - Measure correlation with actual success rates
   - A/B test with chemist preferences

### **Phase 2 Tests**

3. **Optimization Suggestions**:

   - Validate chemistry of suggested alternatives
   - Test constraint handling (cost, air stability)
   - Verify trade-off predictions

4. **Troubleshooting**:
   - Test with real failed reaction reports
   - Validate suggested fixes
   - Track success rate of fixes

---

## Performance Considerations

### **Latency Impact**

| Feature                    | LLM Calls | Est. Latency | Mitigation          |
| -------------------------- | --------- | ------------ | ------------------- |
| Explanation (1 condition)  | 1         | +2s          | Async generation    |
| Explanation (3 conditions) | 3         | +6s          | Parallel calls      |
| Reranking (10 candidates)  | 1         | +3s          | Single batch call   |
| Optimization suggestions   | 1         | +3s          | Optional, cached    |
| **Total (all features)**   | **5**     | **~12s**     | **Background jobs** |

### **Cost Impact**

Assuming DeepSeek pricing (~$0.001/1K tokens):

| Feature      | Tokens/call | Cost/call | Calls/day (100 users) | Daily Cost    |
| ------------ | ----------- | --------- | --------------------- | ------------- |
| Explanation  | ~400        | $0.0004   | 300                   | $0.12         |
| Reranking    | ~800        | $0.0008   | 100                   | $0.08         |
| Optimization | ~600        | $0.0006   | 50                    | $0.03         |
| **Total**    |             |           |                       | **$0.23/day** |

**Annual cost**: ~$84 (negligible)

### **Optimization Strategies**

1. **Caching**: Cache explanations for common reactions
2. **Async**: Generate explanations in background
3. **Batch**: Combine multiple LLM calls where possible
4. **Feature flags**: Allow users to disable LLM features
5. **Fallback**: Gracefully degrade if LLM unavailable

---

## Risk Assessment

### **Risks & Mitigations**

| Risk                     | Impact | Probability | Mitigation                          |
| ------------------------ | ------ | ----------- | ----------------------------------- |
| **LLM hallucinations**   | High   | Medium      | Human validation, confidence scores |
| **Latency too high**     | Medium | Low         | Async processing, caching           |
| **Cost escalation**      | Low    | Low         | Rate limiting, monitoring           |
| **API downtime**         | Medium | Low         | Graceful fallback to non-LLM mode   |
| **Inconsistent quality** | Medium | Medium      | Prompt engineering, temperature=0   |

### **Quality Assurance**

1. **Pre-launch**:

   - Test with 100+ diverse reactions
   - Human chemist review of explanations
   - Benchmark against known literature

2. **Post-launch**:
   - User feedback collection
   - Track explanation ratings
   - Monitor for hallucinations
   - A/B testing vs baseline

---

## Success Metrics

### **Phase 1 (Explanations + Reranking)**

| Metric                     | Target                   | Measurement              |
| -------------------------- | ------------------------ | ------------------------ |
| **Explanation accuracy**   | >90% chemist approval    | Human evaluation         |
| **Explanation usefulness** | >4.0/5.0 rating          | User surveys             |
| **Reranking improvement**  | +10% top-1 accuracy      | Benchmarks vs known good |
| **User adoption**          | >50% enable LLM features | Feature flag analytics   |
| **Latency p95**            | <5s for 3 explanations   | API monitoring           |

### **Phase 2 (Optimization + Troubleshooting)**

| Metric                       | Target                 | Measurement   |
| ---------------------------- | ---------------------- | ------------- |
| **Alternative validity**     | >85% chemically sound  | Expert review |
| **Troubleshooting accuracy** | >70% fix success rate  | User reports  |
| **Time saved**               | 30 min/failed reaction | User surveys  |

---

## Conclusion

### **Recommended Approach: Hybrid Workflow**

✅ **Keep DRFP for retrieval** (it's excellent at finding similar reactions)  
✅ **Add LLM for post-processing** (explanations, reranking, optimization)  
✅ **Make LLM features optional** (graceful degradation)  
✅ **Start with explanations** (lowest risk, highest user value)

### **Implementation Timeline**

- **Week 1-2**: Phase 1 (Explanations + Reranking)
- **Week 3-4**: Testing & refinement
- **Week 5-6**: Phase 2 (Optimization suggestions)
- **Week 7-8**: Phase 3 (Troubleshooting - if validated)

### **Expected Benefits**

1. 📚 **Educational value**: Chemists learn WHY conditions work
2. 🎯 **Better recommendations**: Chemistry-aware ranking
3. ⚡ **Faster optimization**: Automated suggestion generation
4. 🔧 **Easier troubleshooting**: Interactive problem-solving
5. 💪 **Competitive advantage**: Most advanced condition recommendation system

---

## Next Steps

1. **Review this proposal** with team
2. **Prototype `explain_conditions_llm()`** with 10 test cases
3. **Gather user feedback** on explanation quality
4. **Implement Phase 1** if validated
5. **Measure impact** and iterate

---

**Status**: Ready for implementation  
**Owner**: TBD  
**Timeline**: 4-8 weeks (phased rollout)  
**Dependencies**: LLM API access (existing), chemtools enhancements (minimal)
