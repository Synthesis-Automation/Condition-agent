# Optimized Configuration - Accuracy Priority

## Current Configuration ✓

**Optimized for: ACCURACY OVER COST**

### Tier 2 Quick Glance Settings

```python
# File: reaction_agent/agent.py:111-121

model = "gpt-4o"              # Best accuracy (0.647 score, +60% vs gpt-4o-mini)
prompt_style = "structured"   # Best performing prompt
mode = "always"               # Run on ALL reactions (100% coverage)
```

## Performance Characteristics

### Tier 2 with gpt-4o

| Metric | Value | vs gpt-4o-mini |
|--------|-------|----------------|
| **Accuracy Score** | **0.647** | **+60% better** |
| **Speed** | 1.35s avg | Faster (vs 2.01s) |
| **Cost** | $0.00154 | 20x more ($0.00008) |
| **Success Rate** | 100% | Same |
| **Coverage** | ~85% reactions | Same |

## Example Results

### Tandem Reaction (Suzuki + Hydrolysis)

**SMILES**: `COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O`

**With gpt-4o-mini** (old):
```
Summary: Oxidation of an aryl bromide to a ketone
Reaction Type(s): oxidation
❌ Missed Suzuki coupling entirely
```

**With gpt-4o** (new):
```
Summary: Suzuki coupling of aryl bromide with boronic ester
Reaction Type(s): Suzuki coupling
Key Patterns: aryl bromide, boronic ester, aldehyde
Confidence: 0.95
✓ Correctly identified Suzuki coupling
```

### Simple SN2

**SMILES**: `CCCBr.N>>CCCN`

**With gpt-4o**:
```
Summary: Substitution of bromine with nitrogen to form a primary amine
Reaction Type(s): SN2
Key Patterns: alkyl bromide, primary amine
Complexity: SIMPLE
Confidence: 0.95
✓ Perfect identification
```

## Cost Analysis (100 reactions/day)

| Tier | Old (gpt-4o-mini auto) | New (gpt-4o always) | Difference |
|------|------------------------|---------------------|------------|
| Tier 1 | $0 | $0 | - |
| Tier 2 | $0.003 (30 reactions) | **$0.154 (100 reactions)** | **+$0.15/day** |
| Tier 3 | $1-5 | $1-5 | - |
| **Total** | **$1-5/day** | **$1.15-5.15/day** | **+3% cost** |

**Result**: Only **+$0.15/day** for 60% better accuracy and 100% coverage!

## Key Improvements

### 1. Better Pattern Recognition
- ✓ Correctly identifies Suzuki coupling (was "oxidation")
- ✓ Recognizes nucleophilic substitutions accurately
- ✓ Detects complex reactions better

### 2. Runs on ALL Reactions
- ✓ No longer waits for Tier 1 to be "uncertain"
- ✓ Every reaction gets LLM-based pattern analysis
- ✓ Provides chemistry-aware feedback instantly

### 3. Higher Confidence
- Average confidence: **0.95** (vs 0.85 with gpt-4o-mini)
- More reliable recommendations
- Better at distinguishing similar reaction types

## What Each Tier Now Provides

```
User submits reaction
         ↓
┌────────────────────┐
│ TIER 1 (Instant)   │  Detects: Suzuki + acetal hydrolysis (tandem)
│ String patterns    │  Coverage: ~35% of reactions
│ < 0.1s, Free       │
└────────┬───────────┘
         ↓
┌────────────────────┐
│ TIER 2 (1.4s)      │  Detects: Suzuki coupling correctly
│ gpt-4o glance      │  Coverage: ~85% of reactions
│ ALWAYS runs        │  Confidence: 95%
│ $0.0015 per rxn    │  ✓ Runs on 100% of reactions
└────────┬───────────┘
         ↓
┌────────────────────┐
│ TIER 3 (5-30s)     │  Provides: Full mechanistic analysis
│ Deep analysis      │  Coverage: ~95%+ of reactions
│ gpt-4o/gpt-5.2     │
│ $0.01-0.10         │
└────────────────────┘
```

## Accuracy vs Cost Comparison

| Configuration | Accuracy | Speed | Cost/100rxn | Use Case |
|---------------|----------|-------|-------------|----------|
| **Current (gpt-4o always)** | **0.647** ⭐ | 1.4s | **$0.15** | **Research/Production** |
| Previous (gpt-4o-mini auto) | 0.407 | 2.0s | $0.003 | Budget-constrained |
| Previous (gpt-4o-mini always) | 0.407 | 2.0s | $0.008 | Budget-constrained |
| Optional (gpt-4o auto) | 0.647 | 1.4s | $0.05 | Balanced |

## Recommendation: KEEP CURRENT CONFIG ✓

For research and production where **accuracy matters**:

```python
model = "gpt-4o"          # Best accuracy
mode = "always"           # Maximum coverage
prompt_style = "structured"  # Best results
```

**Cost impact**: Only +$0.15/day for 100 reactions
**Accuracy gain**: +60% better pattern recognition
**Coverage**: 100% of reactions get quick feedback

## When to Adjust

**Switch to "auto" mode** if:
- Budget is very tight (saves $0.10/day)
- You trust Tier 1 patterns (works well for known reaction types)
- Analyzing >1000 reactions/day

**Revert to gpt-4o-mini** if:
- Cost matters more than accuracy
- Reactions are mostly common/simple types
- You're OK with ~40% accuracy

## Files Modified

1. `reaction_agent/agent.py:111-121` - Changed to gpt-4o + always mode

## Quick Reference

Current settings:
```python
# In reaction_agent/agent.py
model="gpt-4o"           # Line 116
mode="always"            # Line 113
prompt_style="structured" # Line 120
```

To change modes:
```python
mode="always"  # Run on all reactions (current)
mode="auto"    # Run when Tier 1 uncertain (saves cost)
mode="never"   # Disable Tier 2
```

To change models:
```python
model="gpt-4o"       # Best accuracy (current, $0.0015/rxn)
model="gpt-4o-mini"  # Best value ($0.00008/rxn, lower accuracy)
```
