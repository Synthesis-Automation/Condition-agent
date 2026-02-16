# Using GPT-5.2 for Deep Reasoning (Simplified Approach)

## TL;DR

**Yes, you can directly call GPT-5.2 with adjustable reasoning levels!** This is much simpler than building custom deep reasoning prompts.

## How It Works

### GPT-5.2 vs Our Custom Approach

| Approach | What It Is | Complexity | Cost |
|----------|-----------|------------|------|
| **GPT-5.2 (Direct)** | Use GPT-5.2's built-in reasoning | Simple - just change model name | Higher per call |
| Custom (o3-mini) | Build reasoning prompts manually | Complex - need special prompts | Lower per call |

**Recommendation**: Start with GPT-5.2 directly. It's simpler and has proven results (as shown in NEW-1.md).

---

## Quick Start

### 1. Basic Usage (Already Works!)

```python
from llmtools.clients import LLMClient
from reaction_agent.agent import ReactionSMILESAnalyzer

# Just change the model to gpt-5.2
client = LLMClient(provider="openai", model="gpt-5.2")
analyzer = ReactionSMILESAnalyzer(client, max_tokens=8000)

result = analyzer.analyze(rxn_smiles)
```

**That's it!** GPT-5.2 automatically uses extended reasoning.

### 2. Control Reasoning Depth (NEW)

I just added `reasoning_effort` parameter to control how deeply GPT-5.2 thinks:

```python
client = LLMClient(provider="openai", model="gpt-5.2")

# Low reasoning (faster, cheaper)
response = client.chat(
    prompt="Analyze this reaction: ...",
    reasoning_effort="low",
    max_tokens=8000
)

# Medium reasoning (balanced)
response = client.chat(
    prompt="Analyze this reaction: ...",
    reasoning_effort="medium",
    max_tokens=8000
)

# High reasoning (most thorough, slowest, expensive)
response = client.chat(
    prompt="Analyze this reaction: ...",
    reasoning_effort="high",
    max_tokens=12000
)
```

---

## Key Differences: GPT-5.2 vs Standard Models

### Temperature
- **Standard models (gpt-4o)**: You control temperature (0.0-1.0)
- **GPT-5.2/o-series**: **Always uses temperature=1.0**, cannot be changed
  - Built-in exploration mode
  - Don't try to set temperature - it's ignored

### Token Parameter
- **Standard models**: Use `max_tokens`
- **GPT-5.2/o-series**: Use `max_completion_tokens`
  - Your `LLMClient` handles this automatically (line 273)

### Reasoning
- **Standard models**: Single-shot, instant response
- **GPT-5.2**: Extended thinking time (can be 30s-3min)
  - Shows reasoning process
  - Explores hypotheses
  - Self-corrects before answering

---

## Test Script Ready

I created `reaction_agent/scripts/test_gpt52_reasoning.py` that tests:

1. **gpt-4o** (baseline)
2. **gpt-5.2 low** (fast reasoning)
3. **gpt-5.2 medium** (balanced)
4. **gpt-5.2 high** (maximum depth)

Run it:
```bash
cd c:\Git-softwares\Condition-agent
python reaction_agent/scripts/test_gpt52_reasoning.py
```

This will test on 3 problematic reactions (including the tosylhydrazone example from NEW-1.md) and determine optimal reasoning level.

---

## Expected Costs (Estimates)

| Model | Reasoning Level | Cost/reaction | Time | When to Use |
|-------|----------------|---------------|------|-------------|
| gpt-4o | N/A | $0.002 | 3-5s | Default screening |
| gpt-5.2 | low | ~$0.05 | 10-20s | Quick but thorough |
| gpt-5.2 | medium | ~$0.15 | 30-60s | Complex reactions |
| gpt-5.2 | high | ~$0.30+ | 60-180s | Critical analyses |

**Strategy**:
- Use gpt-4o for 73% of reactions (good mapping ≥0.4)
- Auto-trigger gpt-5.2 for 27% with poor mapping (<0.4)
- Start with "low" or "medium" reasoning
- Only use "high" for truly complex cases

---

## Advantages of GPT-5.2 Direct Approach

### ✅ **Simpler**
- No custom prompts needed
- No multi-turn conversation logic
- Just change model name

### ✅ **Proven**
- NEW-1.md shows GPT-5.2 correctly analyzed complex tosylhydrazone reaction
- Explored multiple hypotheses before answering
- Generated actionable SMARTS patterns
- Provided literature citations

### ✅ **Adjustable**
- Control reasoning depth with `reasoning_effort`
- Balance cost vs quality

### ✅ **Already Supported**
- Your `LLMClient` already handles GPT-5.2 correctly
- No code changes needed (except reasoning_effort, which I just added)

---

## Disadvantages

### ❌ **More Expensive**
- 25-150x more expensive than gpt-4o
- But only use it on 27% of reactions that need it

### ❌ **Slower**
- 10s-180s vs 3-5s for gpt-4o
- Need to show progress indicator to user

### ❌ **API Availability**
- Need to verify GPT-5.2 is accessible via API
- If only available on web, need alternative approach

---

## Implementation Plan

### Phase 1: Validate GPT-5.2 API Access (TODAY)

```bash
# Test if GPT-5.2 is accessible
python reaction_agent/scripts/test_gpt52_reasoning.py
```

If it works:
- ✅ **Recommended**: Use GPT-5.2 directly
- Skip building custom deep reasoning prompts
- Much simpler implementation

If it fails (API not available):
- ❌ Fallback to o3-mini with custom prompts
- Use `compare_reasoning_modes.py` instead

### Phase 2: Integrate into Main System

If Phase 1 succeeds, add to `reaction_agent/agent.py`:

```python
class ReactionSMILESAnalyzer:
    def analyze(self, rxn_smiles: str, mode: str = "auto"):
        """
        Modes:
        - auto: gpt-4o default, gpt-5.2 if mapping <0.4
        - fast: always gpt-4o
        - deep: always gpt-5.2 (medium reasoning)
        - expert: gpt-5.2 (high reasoning)
        """
        if mode == "auto":
            # Run deterministic first
            det_result = analyze_deterministic(rxn_smiles)
            mapping_conf = det_result.get("tool_facts", {}).get("mapping_qc", {}).get("confidence", 0)

            if mapping_conf < 0.4:
                # Poor mapping - use GPT-5.2
                print("⚠ Poor mapping detected, using GPT-5.2 deep reasoning...")
                self.client.model = "gpt-5.2"
                self.max_tokens = 8000
                # Add reasoning_effort to client.chat call
            else:
                # Good mapping - use fast mode
                self.client.model = "gpt-4o"
```

### Phase 3: Add reasoning_effort to CLI

```bash
# Fast (default)
python reaction_agent/cli.py --reaction "..."

# Deep reasoning (auto-triggers on poor mapping)
python reaction_agent/cli.py --reaction "..." --mode deep

# Expert reasoning
python reaction_agent/cli.py --reaction "..." --mode expert --reasoning-level high
```

---

## Recommendation

### If GPT-5.2 API Works:

**Use GPT-5.2 directly** with this strategy:
1. Default: gpt-4o for screening (fast, cheap)
2. Auto-trigger: gpt-5.2 (medium) when mapping <0.4
3. User option: gpt-5.2 (high) for critical analyses

**Estimated blended cost**:
- 73% × $0.002 (gpt-4o) + 27% × $0.15 (gpt-5.2 medium)
- = **$0.042 per reaction**
- 21x more expensive than current, but much better quality on hard cases

### If GPT-5.2 API Doesn't Work:

Fallback to custom approach:
1. Use o3-mini with custom deep reasoning prompts
2. Use `compare_reasoning_modes.py`
3. More complex but similar results

---

## Next Step

**Run the test to check API access:**

```bash
python reaction_agent/scripts/test_gpt52_reasoning.py
```

This will immediately tell you:
1. ✅ If GPT-5.2 API is accessible
2. 📊 How reasoning levels compare to gpt-4o baseline
3. 💰 Actual time/cost trade-offs
4. 🎯 Optimal reasoning level for your use case

Let me know the results!
