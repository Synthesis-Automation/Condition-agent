# How to Always Use GPT-5.2

There are now **3 ways** to force GPT-5.2 usage in the CLI:

## Option 1: Use `--model gpt-5.2` (NEW - Simplest)

```bash
python reaction_agent/cli.py --model gpt-5.2 --reaction "SMILES"
```

**What it does:**
- Auto-detects GPT-5.2 model
- Automatically uses `mode=deep` (low reasoning effort)
- Forces GPT-5.2 for all reactions

**Example output:**
```
Model: gpt-5.2 (openai)
ℹ Detected reasoning model 'gpt-5.2' - using deep reasoning mode
✓ Initialized: gpt-5.2
Model Selected: gpt-5.2
Reasoning Effort: low
```

## Option 2: Use `--mode deep`

```bash
python reaction_agent/cli.py --mode deep --reaction "SMILES"
```

**What it does:**
- Forces `gpt-5.2` with `reasoning_effort=low`
- Works regardless of `--model` parameter

## Option 3: Use `--mode expert`

```bash
python reaction_agent/cli.py --mode expert --reaction "SMILES"
```

**What it does:**
- Forces `gpt-5.2` with `reasoning_effort=medium`
- Higher quality analysis (slower, more expensive)

---

## Comparison

| Method | Model | Reasoning | Speed | Cost | Use Case |
|--------|-------|-----------|-------|------|----------|
| `--model gpt-5.2` | gpt-5.2 | low | ~30s | $0.04 | **Most common - always use GPT-5.2** |
| `--mode deep` | gpt-5.2 | low | ~30s | $0.04 | Same as above |
| `--mode expert` | gpt-5.2 | medium | ~60s | $0.06 | Critical analyses requiring higher quality |
| `--mode auto` (default) | gpt-4o or gpt-5.2 | auto | 4-30s | $0.03 | Smart switching based on mapping quality |
| `--mode fast` | gpt-4o | none | ~4s | $0.0034 | Quick screening |

---

## Batch Processing with GPT-5.2

Process all reactions with GPT-5.2:

```bash
python reaction_agent/cli.py --model gpt-5.2 --batch reactions.txt
```

---

## Python API

```python
from llmtools.clients import LLMClient
from reaction_agent.agent import ReactionSMILESAnalyzer

# Initialize with gpt-5.2
client = LLMClient(provider="openai", model="gpt-5.2", timeout=300)
analyzer = ReactionSMILESAnalyzer(client, max_tokens=8000)

# Option 1: Use deep mode (forces gpt-5.2)
result = analyzer.analyze(rxn_smiles, mode="deep")

# Option 2: Since client is already gpt-5.2, any mode will use it
result = analyzer.analyze(rxn_smiles, mode="auto")  # Still uses gpt-5.2
```

---

## Supported Model Names

These model names are auto-detected as reasoning models and trigger deep mode:

- `gpt-5.2` ✓
- `gpt-5.2-*` (any variant) ✓
- `gpt-o3-mini` ✓
- `o3-mini` ✓
- `gpt-o3` ✓
- `o3` ✓

---

## Summary

**Simplest way to always use GPT-5.2:**
```bash
python reaction_agent/cli.py --model gpt-5.2 --reaction "YOUR_SMILES"
```

This automatically:
- Uses GPT-5.2 for analysis
- Sets reasoning effort to "low" (optimal quality/cost)
- Provides detailed mechanistic analysis
- Catches missing transformations that gpt-4o misses
