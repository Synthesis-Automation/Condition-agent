# Complex Reactions - Workflow Testing in Progress

## 🧪 Reactions Being Tested

### 1. **Rare C-N Coupling**
```
Clc1nc2ccccc2s1.Cn1ccnc1 >> CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1
```
- Benzothiazole chloride + methylimidazole
- Complex heterocyclic coupling

### 2. **Rare Three-Member Rings**
```
Clc1ccccn1.[N-]=[N+]=NCc1ccc(C2=NC2)cc1 >> [N-]=[N+]=NCc1ccc(-c2cnc3ccccn23)cc1
```
- Chloropyridine + azide-aziridine substrate
- Ring-forming/opening reactions

### 3. **Iodoium Substitution**
```
Cc1cc(C)cc(-c2ccc([I+]c3ccccc3)c(CCN(OSi...)c2)c1 + tosylate
>> Cc1cc(C)cc(-c2ccc3[nH]ccc3c2)c1
```
- Iodonium salt intramolecular cyclization
- Complex hypervalent iodine chemistry

---

## 🔬 Testing Matrix

### Workflows Being Tested (5 total):

| Workflow | Models | Max Tokens | Skip Mapping | Focus |
|----------|--------|------------|--------------|-------|
| **Fast Screening** | gpt-4o-mini | 2000 | Yes/No | High throughput |
| **Balanced** | gpt-4o | 2500-3000 | No | Quality/cost balance |
| **High Quality** | gpt-4o, gpt-5.2 | 3000-3500 | No | Complex reactions |
| **Reasoning** | o3-mini | 3000-4000 | No | Multi-stage logic |
| **Maximum Detail** | o3 | 3500-4500 | No | Novel mechanisms |

**Total configurations per reaction**: ~10-11
**Total tests**: ~30-33 LLM calls

---

## 📊 What's Being Measured

For each configuration:
- ✅ **Confidence** - Overall analysis confidence
- ✅ **Events Detected** - Number of mechanistic events
- ✅ **Classification** - Reaction class (e.g., cross_coupling)
- ✅ **Tags** - Named reactions identified
- ✅ **Tokens Used** - Cost indicator
- ✅ **Latency** - Speed
- ✅ **Warnings** - Quality indicators
- ✅ **Mapping Confidence** - Deterministic tool quality

---

## 🎯 Expected Outcomes

### For Each Reaction:
1. **Best Configuration** - Highest confidence
2. **Most Detailed** - Most mechanistic events
3. **Fastest** - Lowest latency
4. **Best Value** - Best confidence per token

### Overall:
- Comparison table across all workflows
- Recommendations per reaction type
- Performance insights
- Cost-quality trade-offs

---

## ⏱️ Estimated Time

- **Fast workflows**: ~5-10 seconds per test
- **Reasoning workflows**: ~15-25 seconds per test
- **Total estimated**: 8-12 minutes for all tests

---

## 📁 Output Files

Results will be saved in: `reaction_agent/results/`
- `workflow_comparison_*.json` - Full detailed results
- `recommendations_*.json` - Best configurations
- `workflow_summary_*.md` - Human-readable summary

---

*Testing started...*
*Check progress with: `tail -f /path/to/output/file`*
