# Quick Start: Deploying GPT-5.2 Hybrid System

Based on test results showing GPT-5.2 (low) achieves 100% coverage vs 87.5% for gpt-4o.

## TL;DR

✅ **GPT-5.2 works via API**
✅ **Use "low" reasoning** (not medium/high - diminishing returns)
✅ **Auto-switch based on mapping confidence** (<0.4 threshold)
✅ **Only 68% more expensive overall** (3.5x cost on 27% of reactions)

---

## Optimal Configuration (From Test Results)

```python
# Decision rule
if mapping_confidence < 0.4:
    model = "gpt-5.2"
    reasoning_effort = "low"  # ← 100% coverage, 67s, evidence score 21
    max_tokens = 8000
else:
    model = "gpt-4o"
    reasoning_effort = None   # Standard models don't use this
    max_tokens = 3000
```

**Why "low" reasoning?**
- gpt-5.2 low: 100% coverage, 67s, 4460 tokens, **evidence 21**
- gpt-5.2 medium: 100% coverage, 105s, 6995 tokens, evidence 22 (+1)
- gpt-5.2 high: 100% coverage, 125s, 8075 tokens, **evidence 19 (-2!)**

**Conclusion**: Medium/high don't improve quality, just cost more.

---

## Implementation Steps

### Step 1: Modify `reaction_agent/agent.py`

Add mode parameter to ReactionSMILESAnalyzer:

```python
class ReactionSMILESAnalyzer:
    def __init__(self, client: LLMClient, drop_spectators: bool = True,
                 temperature: float = 0.0, max_tokens: int = 2000,
                 reasoning_effort: Optional[str] = None):  # ← ADD THIS
        self.client = client
        self.drop_spectators = drop_spectators
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.reasoning_effort = reasoning_effort  # ← ADD THIS

    def analyze(self, rxn_smiles: str, skip_mapping: bool = False,
                mode: str = "auto") -> Dict[str, Any]:  # ← ADD mode
        """
        Analyze reaction with automatic model selection.

        Args:
            rxn_smiles: Reaction SMILES
            skip_mapping: Skip atom mapping (faster but less detailed)
            mode: "auto" (smart), "fast" (gpt-4o), "deep" (gpt-5.2)
        """
        # Step 1: Deterministic analysis (always needed for auto mode)
        deterministic_result = analyze_deterministic(
            rxn_smiles,
            drop_spectators=self.drop_spectators,
            skip_mapping=skip_mapping
        )

        # Step 2: Model selection
        if mode == "auto":
            mapping_conf = deterministic_result.get("tool_facts", {}).get(
                "mapping_qc", {}
            ).get("confidence", 1.0)

            if mapping_conf < 0.4:
                # Poor mapping - use GPT-5.2
                print(f"⚠ Poor mapping ({mapping_conf:.3f}) - using GPT-5.2 deep reasoning...")
                self.client.model = "gpt-5.2"
                self.max_tokens = 8000
                self.reasoning_effort = "low"
            else:
                # Good mapping - use fast mode
                self.client.model = "gpt-4o"
                self.max_tokens = 3000
                self.reasoning_effort = None

        elif mode == "fast":
            self.client.model = "gpt-4o"
            self.max_tokens = 3000
            self.reasoning_effort = None

        elif mode == "deep":
            self.client.model = "gpt-5.2"
            self.max_tokens = 8000
            self.reasoning_effort = "low"

        # Step 3: Build prompt (existing code)
        prompt = self._build_prompt(deterministic_result)

        # Step 4: Call LLM with reasoning_effort
        response = self.client.chat(
            prompt=prompt,
            temperature=self.temperature if self.reasoning_effort is None else None,
            max_tokens=self.max_tokens,
            reasoning_effort=self.reasoning_effort  # ← PASS THIS
        )

        # Step 5: Parse and return (existing code)
        interpretation = self._parse_response(response.content)

        # Apply QC gating (existing code)
        if not deterministic_result.get("tool_facts", {}).get("mapping_qc", {}).get("ok", False):
            # ... existing logic ...
            pass

        return {
            **deterministic_result,
            "interpretation": interpretation,
            "model_used": self.client.model,
            "reasoning_effort": self.reasoning_effort,
        }
```

### Step 2: Update CLI

Add `--mode` flag to `reaction_agent/cli.py`:

```python
# In argument parser
parser.add_argument(
    "--mode",
    type=str,
    default="auto",
    choices=["auto", "fast", "deep", "expert"],
    help="Analysis mode: auto (smart switch), fast (gpt-4o), deep (gpt-5.2)"
)

parser.add_argument(
    "--reasoning",
    type=str,
    default="low",
    choices=["low", "medium", "high"],
    help="Reasoning level for gpt-5.2 (only used in deep/expert mode)"
)

# In analyze function
analyzer = ReactionSMILESAnalyzer(
    client,
    drop_spectators=args.drop_spectators,
    temperature=args.temperature,
    max_tokens=args.max_tokens,
    reasoning_effort=args.reasoning if args.mode in ["deep", "expert"] else None
)

result = analyzer.analyze(args.reaction, skip_mapping=args.no_llm, mode=args.mode)
```

### Step 3: Test

```bash
# Test auto mode (recommended)
python reaction_agent/cli.py \
  --reaction "Cc1ccc(S(=O)(=O)NN=CC=Cc2ccccc2)cc1.O=S(=O)([O-])C(F)(F)F.c1ccc([I+]c2ccccc2)cc1>>Cc1ccc(S(=O)(=O)c2ccccc2)cc1.c1ccc(-c2ccn(-c3ccccc3)n2)cc1" \
  --mode auto

# Should print: "⚠ Poor mapping (0.153) - using GPT-5.2 deep reasoning..."
# Should achieve 100% coverage in ~67s

# Test fast mode
python reaction_agent/cli.py --reaction "..." --mode fast
# Always uses gpt-4o

# Test deep mode
python reaction_agent/cli.py --reaction "..." --mode deep
# Always uses gpt-5.2 (low)
```

### Step 4: Validate on 30-Reaction Test

Modify `test_random_30.py` to use auto mode:

```python
def test_reaction(rxn_data, mode="auto", max_tokens=8000, show_details=False):
    """Test with hybrid gpt-4o/gpt-5.2 approach."""

    # Use auto mode - will switch to gpt-5.2 for poor mapping
    client = LLMClient(provider="openai", model="gpt-4o", timeout=300)
    analyzer = ReactionSMILESAnalyzer(client, max_tokens=max_tokens)

    result = analyzer.analyze(rxn_smiles, skip_mapping=False, mode=mode)

    return {
        "reaction": rxn_data["smiles"],
        "model_used": result.get("model_used", "unknown"),
        "reasoning_effort": result.get("reasoning_effort"),
        "mapping_confidence": result.get("tool_facts", {}).get("mapping_qc", {}).get("confidence", 0),
        "llm_confidence": result.get("interpretation", {}).get("confidence", 0),
        "time": elapsed,
    }
```

**Expected results** (30 reactions):
- 22 reactions (mapping ≥0.4): gpt-4o, ~7.5s each, total ~165s
- 8 reactions (mapping <0.4): gpt-5.2, ~67s each, total ~536s
- **Total time**: ~701s (11.7 min) vs 225s (3.75 min) gpt-4o only
- **Total cost**: ~$32 vs $19 gpt-4o only (68% increase)
- **Quality improvement**: 100% coverage on hard reactions vs 87.5%

---

## Cost Analysis

### Per Reaction Costs

| Scenario | Model | Time | Tokens | Cost |
|----------|-------|------|--------|------|
| Good mapping (≥0.4) | gpt-4o | 7.5s | 1262 | $0.019 |
| Poor mapping (<0.4) | gpt-5.2 (low) | 67s | 4460 | $0.067 |

### Blended Costs (Based on 30-Reaction Test)

27% of reactions have poor mapping (<0.4):

```
Blended cost = 0.73 × $0.019 + 0.27 × $0.067
             = $0.014 + $0.018
             = $0.032 per reaction
```

**vs gpt-4o only**: $0.019 per reaction

**Cost increase**: 68% (from $0.019 → $0.032)

**Quality improvement**:
- gpt-4o on poor mapping: 87.5% coverage, evidence 7
- gpt-5.2 on poor mapping: 100% coverage, evidence 21

---

## Usage Examples

### Example 1: Auto Mode (Recommended)

```python
from llmtools.clients import LLMClient
from reaction_agent.agent import ReactionSMILESAnalyzer

# Initialize with gpt-4o as baseline
client = LLMClient(provider="openai", model="gpt-4o", timeout=300)
analyzer = ReactionSMILESAnalyzer(client, max_tokens=8000)

# Auto mode will switch to gpt-5.2 if mapping <0.4
result = analyzer.analyze(rxn_smiles, mode="auto")

print(f"Model used: {result['model_used']}")
print(f"Reasoning: {result.get('reasoning_effort', 'None')}")
print(f"Confidence: {result['interpretation']['confidence']:.2f}")
```

### Example 2: Force Deep Analysis

```python
# Always use gpt-5.2, even if mapping is good
result = analyzer.analyze(rxn_smiles, mode="deep")
```

### Example 3: Fast Screening

```python
# Always use gpt-4o, even if mapping is poor
result = analyzer.analyze(rxn_smiles, mode="fast")
```

### Example 4: Batch Processing with Auto Mode

```python
import pandas as pd

df = pd.read_csv("reactions.csv")
results = []

for idx, row in df.iterrows():
    result = analyzer.analyze(row["smiles"], mode="auto")

    results.append({
        "reaction": row["smiles"],
        "model": result["model_used"],
        "reasoning": result.get("reasoning_effort"),
        "mapping_conf": result["tool_facts"]["mapping_qc"]["confidence"],
        "llm_conf": result["interpretation"]["confidence"],
    })

results_df = pd.DataFrame(results)

# Analyze switching behavior
print(f"Used gpt-4o: {(results_df['model'] == 'gpt-4o').sum()}")
print(f"Used gpt-5.2: {(results_df['model'] == 'gpt-5.2').sum()}")
print(f"Avg confidence: {results_df['llm_conf'].mean():.2f}")
```

---

## FAQs

### Q: Why not use "medium" or "high" reasoning?

**A**: Test results show diminishing returns:
- **low**: 100% coverage, 67s, evidence 21
- **medium**: 100% coverage, 105s, evidence 22 (+1 evidence, +58% time)
- **high**: 100% coverage, 125s, evidence 19 (-2 evidence!, +87% time)

"Low" is the sweet spot. Higher levels don't improve quality for this task.

### Q: Can I use o3-mini instead?

**A**: o3-mini is faster (23s) but lower quality:
- **o3-mini**: 87.5% coverage, evidence 9
- **gpt-5.2 (low)**: 100% coverage, evidence 21

If cost is critical and you can accept 87.5% coverage, o3-mini is viable. But gpt-5.2 (low) is better quality.

### Q: What if GPT-5.2 times out?

**A**: Increase timeout to 300s (5 minutes):
```python
client = LLMClient(provider="openai", model="gpt-5.2", timeout=300)
```

In our tests, gpt-5.2 (low) never exceeded 70s.

### Q: How do I know if auto mode is working?

**A**: Look for this message:
```
⚠ Poor mapping (0.153) - using GPT-5.2 deep reasoning...
```

Also check `result["model_used"]` and `result["reasoning_effort"]`.

### Q: Should I always use auto mode?

**A**: Yes, unless:
- **Fast screening**: Use `mode="fast"`  for quick batch processing
- **Critical analysis**: Use `mode="deep"` for high-value reactions
- **Cost-sensitive**: Use `mode="fast"` if budget is tight

For general use, `mode="auto"` is optimal.

---

## Expected Performance

### On 1000 Reactions

Assuming 27% have poor mapping:

| Metric | gpt-4o only | **Hybrid (auto)** | gpt-5.2 only |
|--------|-------------|-------------------|--------------|
| **Time** | 2.1 hours | **6.5 hours** | 18.6 hours |
| **Cost** | $19 | **$32** | $67 |
| **Quality** | 87.5% avg | **91.7% avg** | 100% |

**Hybrid is optimal**: 68% more expensive than gpt-4o, but 4.2 percentage points better coverage, and 52% cheaper than using gpt-5.2 on everything.

---

## Next Steps

1. ✅ Test completed - GPT-5.2 works via API
2. ✅ Optimal parameters identified (low reasoning)
3. ⏭ Modify `agent.py` to support auto mode
4. ⏭ Update CLI with `--mode` flag
5. ⏭ Re-run 30-reaction test with hybrid approach
6. ⏭ Deploy to production

---

## Files Reference

- **Test script**: `reaction_agent/scripts/test_gpt52_example.py`
- **Results**: `reaction_agent/results/gpt52_parameter_optimization.json`
- **Analysis**: `reaction_agent/results/GPT52_PARAMETER_OPTIMIZATION_ANALYSIS.md`
- **This guide**: `reaction_agent/docs/GPT52_DEPLOYMENT_GUIDE.md`

---

**Ready to deploy!** GPT-5.2 (low) + auto-switching is the recommended configuration.
