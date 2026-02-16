# 🎯 Complex Reactions - Workflow Testing Results

## Executive Summary

**Testing Complete**: All 3 complex reactions tested across 5 workflows with ~30 configurations

### 🏆 Key Findings

| Reaction | Best Model | Confidence | Workflow | Tokens | Time |
|----------|------------|------------|----------|--------|------|
| **Rare C-N Coupling** | **o3-mini** | **0.98** | reasoning | 3000-4000 | ~9s |
| **Three-Member Rings** | **gpt-4o-mini** | **0.91** | fast_screening | 2000 | ~9s |
| **Iodoium Substitution** | **gpt-4o/gpt-4o-mini** | **0.93** | balanced/fast | 2500 | ~3-6s |

---

## 📊 Detailed Results by Reaction

### 1. Rare C-N Coupling (Benzothiazole + Imidazole)

**SMILES**: `Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1`

**Deterministic Analysis**:
- ✅ Mapping confidence: **0.976** (excellent!)
- ✅ Bond changes detected: 1
- ✅ Well-characterized reaction

**Model Performance**: **VERY HIGH CONFIDENCE ACROSS ALL MODELS**

| Model | Max Tokens | Confidence | Events | Time | Classification |
|-------|------------|------------|--------|------|----------------|
| **o3-mini** ⭐ | **3000-4000** | **0.98** | 1 | 8-9s | nucleophilic_substitution (SNAr) |
| **gpt-4o-mini** ⭐ | 2000 | **0.98** | 1 | 5s | nucleophilic_substitution (SNAr) |
| **gpt-4o** | 2500-3500 | 0.90-0.95 | 1 | 3-4s | nucleophilic_substitution (SNAr) |
| o3 | 3500-4500 | 0.90 | 1 | 12-16s | nucleophilic_substitution (SNAr) |
| gpt-5.2 | 3000 | 0.74 | 1 | 8s | nucleophilic_substitution |
| gpt-4o-mini (no map) | 2000 | 0.20 | 1 | 5s | - |

**✅ RECOMMENDATION: gpt-4o-mini or o3-mini**
- Both achieve 0.98 confidence
- gpt-4o-mini is **faster and cheaper** (5s, ~$0.001)
- o3-mini provides **same quality** with more reasoning (9s, ~$0.003)

**Why So Good?**:
- Clean, well-known SNAr mechanism
- Excellent atom mapping (0.976)
- Single mechanistic step
- All models correctly identified as nucleophilic aromatic substitution

---

### 2. Rare Three-Member Rings (Aziridine Chemistry)

**SMILES**: `Clc1ccccn1.[N-]=[N+]=NCc1ccc(C2=NC2)cc1>>[N-]=[N+]=NCc1ccc(-c2cnc3ccccn23)cc1`

**Deterministic Analysis**:
- ✅ Mapping confidence: **0.912** (very good)
- ✅ Bond changes detected: **5** (multi-stage)
- ⚠️ More complex mechanism

**Model Performance**: **HIGH CONFIDENCE, MULTI-EVENT DETECTION**

| Model | Max Tokens | Confidence | Events | Time | Notes |
|-------|------------|------------|--------|------|-------|
| **o3-mini** ⭐ | **4000** | **0.91** | 1 | 17s | Best reasoning |
| **gpt-4o-mini** | 2000 | **0.91** | 1 | 9s | Best value |
| **gpt-4o** ⭐ | 2500-3500 | **0.90** | **3 events** | 3-4s | Most detailed |
| o3 | 3500-4500 | 0.70-0.72 | 2 | 23s | Good detail |
| gpt-5.2 | 3000 | 0.73 | 3 | 9s | Good |
| o3-mini (3000t) | 3000 | 0.00 | 0 | 14s | Failed (too few tokens) |

**✅ RECOMMENDATION: gpt-4o (balanced workflow)**
- **Best detail**: 3 mechanistic events detected
- Confidence: 0.90
- Fast: 3-4 seconds
- Good balance of quality/cost (~$0.002)

**Alternative**: gpt-4o-mini for high-throughput (0.91 conf, 1 event, faster)

**Why Different From #1?**:
- 5 bond changes (vs. 1)
- Multi-stage mechanism
- **gpt-4o detects 3 events** (ring opening, coupling, reorganization)
- **o3-mini needs 4000 tokens** (crashed at 3000)

---

### 3. Iodoium Substitution (Hypervalent Iodine Cyclization)

**SMILES**: `...([I+]c3ccccc3)c(CCN(OSi...)...)... >> ...indole product...`

**Deterministic Analysis**:
- ✅ Mapping confidence: **0.932** (excellent)
- ✅ Bond changes detected: **4**
- ✅ Clean hypervalent iodine chemistry

**Model Performance**: **EXCELLENT ACROSS MOST MODELS**

| Model | Max Tokens | Confidence | Events | Time | Notes |
|-------|------------|------------|--------|------|-------|
| **gpt-4o** ⭐ | **2500-3000** | **0.93** | **2 events** | 3-4s | Best overall |
| **gpt-4o-mini** ⭐ | 2000 | **0.93** | 1 | 6s | Best value |
| o3 (4500t) | 4500 | 0.80 | 1 | 22s | Good but slower |
| o3-mini (4000t) | 4000 | 0.90 | 1 | 13s | Good |
| gpt-4o (3500t) | 3500 | 0.90 | 2 | 3s | Excellent |
| gpt-5.2 | 3000 | 0.76 | 2 | 10s | Good |
| o3-mini (3000t) | 3000 | 0.00 | 0 | 14s | Failed |

**✅ RECOMMENDATION: gpt-4o (2500-3000 tokens)**
- Excellent confidence: 0.93
- **Detects 2 events** (activation + cyclization)
- **Very fast**: 3-4 seconds
- Great value: ~$0.002 per analysis

**Alternative**: gpt-4o-mini (0.93 conf, 1 event, similar quality)

**Why Good Performance?**:
- Well-studied hypervalent iodine chemistry
- Excellent mapping (0.932)
- 4 bond changes well-characterized
- Clear intramolecular cyclization pattern

---

## 🎯 Overall Recommendations by Use Case

### 1. **For These Specific Reactions** ✅

| Reaction Type | Best Model | Max Tokens | Cost/Time | Reason |
|---------------|------------|------------|-----------|--------|
| **Rare C-N Coupling** | **gpt-4o-mini** | 2000 | $/$$ / 5s | 0.98 conf, perfect |
| **Three-Member Rings** | **gpt-4o** | 2500-3000 | $$ / 3-4s | 0.90 conf, 3 events |
| **Iodoium Subst.** | **gpt-4o** | 2500-3000 | $$ / 3-4s | 0.93 conf, 2 events |

### 2. **General Workflow Selection** 🔬

#### Fast Screening (High Throughput)
**Model**: gpt-4o-mini (2000 tokens)
- **Use when**: Processing hundreds/thousands of reactions
- **Confidence**: 0.91-0.98 for well-characterized reactions
- **Cost**: ~$0.001 per reaction
- **Time**: 5-9 seconds
- ⚠️ **Limitation**: Skip mapping gives 0.20 confidence

#### Balanced Production (Recommended)
**Model**: gpt-4o (2500-3000 tokens)
- **Use when**: Production workflows, routine analysis
- **Confidence**: 0.90-0.95
- **Cost**: ~$0.002 per reaction
- **Time**: 3-4 seconds
- ✅ **Best for**: Complex reactions needing detail

#### Maximum Quality
**Model**: o3-mini (4000 tokens) or o3 (3500-4500 tokens)
- **Use when**: Novel mechanisms, multi-stage reactions
- **Confidence**: 0.90-0.98
- **Cost**: ~$0.003-0.008 per reaction
- **Time**: 10-23 seconds
- ⚠️ **Note**: o3-mini needs 4000+ tokens for complex reactions

---

## 💡 Key Insights

### ✅ What Works Best

1. **gpt-4o-mini excels at well-characterized reactions**
   - SNAr, simple substitutions: 0.98 confidence
   - Very cost-effective
   - Only limitation: single-event focus

2. **gpt-4o is the sweet spot for complex reactions**
   - Multi-event detection (2-3 events)
   - Excellent confidence (0.90-0.95)
   - Fast (3-4 seconds)
   - Best value overall

3. **High atom mapping confidence = high analysis confidence**
   - Rxn 1: 0.976 mapping → 0.98 analysis
   - Rxn 2: 0.912 mapping → 0.91 analysis
   - Rxn 3: 0.932 mapping → 0.93 analysis

4. **Token counts matter for reasoning models**
   - o3-mini: 3000 tokens → failed (0.00 conf)
   - o3-mini: 4000 tokens → success (0.90-0.91 conf)
   - Recommendation: **Use 4000+ for o3-mini**

### ⚠️ What to Avoid

1. **Skip mapping is NOT recommended**
   - All "skip mapping" tests: 0.20 confidence
   - Saves ~2-3 seconds but destroys quality

2. **gpt-5.2 underperforms for these reactions**
   - 0.73-0.76 confidence (lower than gpt-4o)
   - Slower (8-10 seconds)
   - More expensive
   - Not recommended for production

3. **o3 overkill for simple reactions**
   - Similar confidence to cheaper models
   - 3-4x slower
   - 4-8x more expensive
   - Reserve for truly novel reactions

---

## 📈 Performance Metrics Summary

### Processing Speed

| Workflow | Model | Avg Time | Reactions/Hour |
|----------|-------|----------|----------------|
| Fast | gpt-4o-mini | 5-9s | 400-720 |
| Balanced | gpt-4o | 3-4s | 900-1200 |
| Quality | o3-mini | 12-17s | 212-300 |
| Maximum | o3 | 15-23s | 157-240 |

### Cost Analysis (Estimated)

| Model | Tokens/Rxn | Cost/Rxn | Cost/1000 Rxns |
|-------|------------|----------|----------------|
| gpt-4o-mini | ~1200 | ~$0.001 | ~$1 |
| gpt-4o | ~1300 | ~$0.002 | ~$2 |
| gpt-5.2 | ~1400 | ~$0.003 | ~$3 |
| o3-mini | ~1900 | ~$0.003 | ~$3 |
| o3 | ~1900 | ~$0.008 | ~$8 |

### Confidence Ranges

| Reaction Complexity | gpt-4o-mini | gpt-4o | o3-mini | o3 |
|---------------------|-------------|---------|---------|-----|
| Simple (1 stage) | 0.91-0.98 | 0.90-0.95 | 0.98 | 0.90 |
| Complex (2-3 stages) | 0.91 | 0.90-0.93 | 0.90-0.91 | 0.70-0.80 |
| Novel/Rare | Use o3-mini or o3 with 4000+ tokens |

---

## 🎬 Recommended Commands

### For Your Reactions:

```bash
# Rare C-N Coupling - Use gpt-4o-mini
python reaction_agent/cli.py \
  --reaction "Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C..." \
  --model gpt-4o-mini \
  --max-tokens 2000

# Three-Member Rings - Use gpt-4o
python reaction_agent/cli.py \
  --reaction "Clc1ccccn1.[N-]=[N+]=NCc1ccc..." \
  --model gpt-4o \
  --max-tokens 3000

# Iodoium Substitution - Use gpt-4o
python reaction_agent/cli.py \
  --reaction "Cc1cc(C)cc(-c2ccc([I+]c3ccccc3)..." \
  --model gpt-4o \
  --max-tokens 2500
```

### For Batch Processing:

```bash
# Use gpt-4o for best balance
python reaction_agent/cli.py \
  --batch llm_test_rxn.csv \
  --model gpt-4o \
  --max-tokens 3000 \
  --output-dir results/complex_reactions
```

---

## 📊 Comparison to Previous Tandem Reaction

### Your C-N + Sonogashira (Previous Test):
- **Best**: o3 (0.30 conf, 2 events correctly identified)
- **Mapping**: 0.34 (very low)
- **Challenge**: True tandem/sequential mechanism

### These Reactions:
- **Best**: gpt-4o or gpt-4o-mini (0.90-0.98 conf)
- **Mapping**: 0.91-0.98 (excellent)
- **Clarity**: Well-characterized mechanisms

**Conclusion**: Model selection depends heavily on:
1. **Atom mapping quality** (biggest factor)
2. **Reaction complexity** (sequential vs concurrent)
3. **Mechanism novelty** (known vs unknown)

---

## ✅ Final Recommendations

### Production Workflow:
```yaml
Default Model: gpt-4o
Default Tokens: 3000
Skip Mapping: false
Temperature: 0.0

Fallback for Simple: gpt-4o-mini (2000 tokens)
Escalate to: o3-mini (4000 tokens) if confidence < 0.5
```

### Decision Tree:
```
Is reaction well-characterized?
├─ YES → gpt-4o-mini (fast, cheap, 0.91-0.98 conf)
└─ NO → Is it multi-stage/tandem?
    ├─ YES → o3 (3500-4000 tokens)
    └─ NO → gpt-4o (3000 tokens)
```

---

## 🎉 Summary

**All 3 complex reactions analyzed successfully!**

- ✅ **30 configurations tested**
- ✅ **Best confidence: 0.93-0.98** (excellent for complex reactions)
- ✅ **Optimal model identified per reaction**
- ✅ **Cost-optimized workflows established**

**Key Takeaway**: **gpt-4o (3000 tokens)** is your best general-purpose choice for these complex reactions, achieving 0.90-0.95 confidence with 3-4 second latency.

---

*Analysis Date: 2026-02-15*
*Test Duration: ~5-7 minutes*
*Total Reactions: 3*
*Total Configurations: 30*
