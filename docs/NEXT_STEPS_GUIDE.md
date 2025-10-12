# Next Steps - Quick Reference Guide

**Current Status**: ✅ Step 1 & 2 Complete  
**Next Action**: Step 3 - Compare vs Individual Modes  
**Estimated Time**: 1-2 days

---

## What's Done ✅

- ✅ **Prototype**: `synthesize_recommendations_llm()` working
- ✅ **Tests**: 3/3 passing (100% success rate)
- ✅ **Documentation**: Complete implementation docs
- ✅ **Validation**: All key features confirmed working

---

## Step 3: Compare vs Individual Modes (TODO Next)

### **Goal**
Quantify quality improvement of multi-source synthesis vs single-source recommendations.

### **Tasks**

1. **Expand test suite to 20 diverse reactions**
   - 5x Suzuki (simple, challenging substrates)
   - 5x Buchwald-Hartwig (various amines)
   - 3x Heck reactions
   - 3x Ullmann couplings
   - 4x Edge cases (ambiguous, low precedent)

2. **Run each reaction through all 4 modes**:
   - ML-only (DRFP precedents)
   - Rule-only (SCDB patterns)
   - Protocol-only (Literature)
   - **Multi-source synthesis** (Our new approach)

3. **Collect metrics**:
   - Recommendation quality (expert evaluation)
   - Explanation clarity (1-5 scale)
   - Handling of constraints
   - Handling of edge cases
   - User preference (if applicable)

4. **Analyze results**:
   - Which mode performs best overall?
   - When does synthesis add most value?
   - Are there cases where single-source is better?

### **Code Template**

```python
# test_mode_comparison.py

BENCHMARK_REACTIONS = [
    {
        "name": "Simple Suzuki",
        "reaction_smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>...",
        "expected_difficulty": "easy",
        "known_good_conditions": {...}
    },
    # ... 19 more reactions
]

def compare_all_modes(reaction):
    """Run all 4 modes and compare."""
    # ML-only
    ml_result = recommend_ml_based(reaction['reaction_smiles'])
    
    # Rule-only
    rule_result = recommend_rule_based(reaction['reaction_smiles'])
    
    # Protocol-only
    protocol_result = recommend_protocol_based(reaction['reaction_smiles'])
    
    # Multi-source synthesis
    synthesis_result = synthesize_recommendations_llm(
        reaction['reaction_smiles'],
        ml_results=ml_result,
        rule_results=rule_result,
        protocol_results=protocol_result,
        llm_client=llm_client
    )
    
    return {
        'ml': ml_result,
        'rule': rule_result,
        'protocol': protocol_result,
        'synthesis': synthesis_result
    }

def evaluate_quality(result, ground_truth):
    """
    Score recommendation quality (1-5):
    5 = Perfect match with ground truth
    4 = Very close, minor differences
    3 = Reasonable alternative
    2 = Suboptimal but viable
    1 = Poor recommendation
    """
    # Compare recommended conditions
    # Check explanation quality
    # Verify constraint handling
    return score

# Run comparison
for reaction in BENCHMARK_REACTIONS:
    results = compare_all_modes(reaction)
    
    for mode, result in results.items():
        score = evaluate_quality(result, reaction['known_good_conditions'])
        print(f"{reaction['name']} - {mode}: {score}/5")
```

### **Expected Outcome**

If multi-source synthesis is working well, we should see:
- **Overall quality**: Synthesis ≥ best single-source mode
- **Edge cases**: Synthesis significantly better (handles uncertainty)
- **Constraints**: Synthesis best at respecting user requirements
- **Explanations**: Synthesis provides best rationale

---

## Step 4: Refine Prompt (TODO After Step 3)

### **Goals**

1. **Optimize latency**:
   - Reduce prompt length (fewer tokens)
   - Compress examples
   - Remove redundant instructions

2. **Improve confidence calibration**:
   - Tune thresholds (when to say high/medium/low)
   - More test data needed

3. **Add specific chemistry guidelines**:
   - Functional group compatibility rules
   - Common failure modes
   - Best practices for specific reaction types

### **Code Changes**

```python
# llmtools/prompts.py

# BEFORE (current): ~1500 tokens
MULTI_SOURCE_SYNTHESIS = PromptTemplate(
    template="""... very detailed prompt ..."""
)

# AFTER (optimized): target ~1000 tokens
MULTI_SOURCE_SYNTHESIS_V2 = PromptTemplate(
    template="""... compressed version with same quality ..."""
)
```

### **Testing**

- Run same 20 reactions with V2 prompt
- Compare quality vs V1
- Measure latency improvement
- Ensure no quality degradation

---

## Step 5: Integrate into Production API (Final Step)

### **Goals**

1. **Create new endpoint**: `POST /api/v1/recommend/synthesize`
2. **Add feature flag**: Gradual rollout
3. **Optimize performance**: Parallel execution
4. **Add monitoring**: Track usage, errors, latency

### **Implementation**

**1. Create endpoint**:

```python
# app/main.py

@app.post("/api/v1/recommend/synthesize")
async def recommend_with_synthesis(
    request: RecommendRequest,
    constraints: Optional[Dict] = Body(None),
    enable_synthesis: bool = Query(True, description="Enable LLM synthesis"),
    llm_provider: str = Query("aliyun"),
    llm_model: str = Query("deepseek-v3.2-exp")
):
    """
    Multi-source LLM synthesis endpoint.
    
    Combines ML, Rule, and Protocol recommendations with LLM analysis.
    Falls back to fusion mode if synthesis disabled or fails.
    """
    try:
        # Run all sources in parallel
        ml_task = asyncio.create_task(recommend_ml_async(request.reaction_smiles))
        rule_task = asyncio.create_task(recommend_rule_async(request.reaction_smiles))
        protocol_task = asyncio.create_task(recommend_protocol_async(request.reaction_smiles))
        
        ml_result, rule_result, protocol_result = await asyncio.gather(
            ml_task, rule_task, protocol_task
        )
        
        # LLM synthesis (if enabled)
        if enable_synthesis:
            llm_client = LLMClient(provider=llm_provider, model=llm_model)
            result = synthesize_recommendations_llm(
                reaction_smiles=request.reaction_smiles,
                ml_results=ml_result,
                rule_results=rule_result,
                protocol_results=protocol_result,
                constraints=constraints,
                llm_client=llm_client
            )
            
            if result['status'] == 'success':
                return result
        
        # Fallback: fusion mode
        return recommend_fusion(ml_result, rule_result)
        
    except Exception as e:
        logger.error(f"Synthesis failed: {e}")
        # Fallback
        return recommend_fusion(request.reaction_smiles, request.reaction_type)
```

**2. Add feature flag**:

```python
# config.py

FEATURE_FLAGS = {
    "enable_llm_synthesis": os.getenv("ENABLE_LLM_SYNTHESIS", "false").lower() == "true",
    "synthesis_rollout_percentage": int(os.getenv("SYNTHESIS_ROLLOUT", "0"))  # 0-100
}

# In endpoint
if not FEATURE_FLAGS["enable_llm_synthesis"]:
    return recommend_fusion(...)

# Gradual rollout
import random
if random.randint(0, 100) > FEATURE_FLAGS["synthesis_rollout_percentage"]:
    return recommend_fusion(...)  # Old behavior
else:
    return synthesize_recommendations_llm(...)  # New behavior
```

**3. Add monitoring**:

```python
from prometheus_client import Counter, Histogram

synthesis_requests = Counter('synthesis_requests_total', 'Total synthesis requests')
synthesis_errors = Counter('synthesis_errors_total', 'Total synthesis errors')
synthesis_latency = Histogram('synthesis_latency_seconds', 'Synthesis latency')

@synthesis_latency.time()
def recommend_with_synthesis_monitored(...):
    synthesis_requests.inc()
    try:
        result = synthesize_recommendations_llm(...)
        return result
    except Exception as e:
        synthesis_errors.inc()
        raise
```

**4. Update documentation**:

```markdown
# API Documentation

## POST /api/v1/recommend/synthesize

Multi-source LLM synthesis endpoint.

### Request
{
  "reaction_smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>...",
  "reaction_type": "Suzuki",
  "constraints": {
    "scale": "multigram",
    "cost": "low",
    "air_sensitivity": "avoid"
  }
}

### Response
{
  "status": "success",
  "synthesis": {
    "confidence_level": "high",
    "recommended_condition": {...},
    "backup_conditions": [...],
    "warnings": [...]
  }
}
```

---

## Timeline

| Step | Duration | Start | End |
|------|----------|-------|-----|
| **Step 3**: Compare modes | 1-2 days | Now | +2 days |
| **Step 4**: Refine prompt | 1 day | +2 days | +3 days |
| **Step 5**: Production integration | 2-3 days | +3 days | +6 days |
| **Total** | **4-6 days** | | |

---

## Success Metrics

### **Step 3 Metrics**

- ✅ Multi-source synthesis ≥ best single-source mode (overall quality)
- ✅ Synthesis handles edge cases better (low precedent scenarios)
- ✅ Constraints respected 100% of the time
- ✅ Explanations rated ≥4/5 for clarity

### **Step 4 Metrics**

- ✅ Latency reduced by ≥30% (prompt optimization)
- ✅ No quality degradation vs V1
- ✅ Confidence calibration improved (match expected levels)

### **Step 5 Metrics**

- ✅ API uptime >99.9%
- ✅ Error rate <1%
- ✅ P95 latency <10 seconds
- ✅ User adoption >50% (if feature flag enabled)

---

## Quick Command Reference

### **Run Tests**

```powershell
# Multi-source synthesis tests
cd c:\Git-softwares\Condition-agent
$env:PYTHONIOENCODING="utf-8"
python test_multisource_synthesis.py

# Mode comparison (Step 3 - to be created)
python test_mode_comparison.py

# Performance benchmarks
python benchmark_synthesis.py
```

### **Check Implementation**

```powershell
# View core function
code llmtools/recommendation_llm.py

# View tests
code test_multisource_synthesis.py

# View documentation
code docs/MULTISOURCE_IMPLEMENTATION.md
code docs/SESSION_SUMMARY_MULTISOURCE.md
```

### **Update Prompt**

```powershell
# Edit prompt template
code llmtools/prompts.py
# Search for: MULTI_SOURCE_SYNTHESIS
```

---

## Key Files

| File | Purpose | Lines |
|------|---------|-------|
| `llmtools/recommendation_llm.py` | Core synthesis function | 327 |
| `llmtools/prompts.py` | Prompt templates | 735 |
| `test_multisource_synthesis.py` | Test suite | 290 |
| `docs/MULTISOURCE_IMPLEMENTATION.md` | Technical docs | 350+ |
| `docs/SESSION_SUMMARY_MULTISOURCE.md` | Session summary | 400+ |
| `docs/LLM_RECOMMENDATION_ENHANCEMENT.md` | Overall plan | 1000+ |

---

## Contact Points

**Questions?**

- Implementation details: See `MULTISOURCE_IMPLEMENTATION.md`
- Test results: See `SESSION_SUMMARY_MULTISOURCE.md`
- Overall strategy: See `LLM_RECOMMENDATION_ENHANCEMENT.md`
- Code: `llmtools/recommendation_llm.py`

---

**Status**: Ready for Step 3 (Comparison Testing)  
**Blockers**: None - All dependencies met  
**Recommendation**: Start Step 3 immediately
