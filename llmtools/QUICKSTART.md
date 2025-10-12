# LLMTools Quick Reference

## Installation & Setup

```bash
pip install openai
export OPENAI_API_KEY="sk-..."
export ALIYUN_API_KEY="sk-..."  # Optional
```

## Basic Usage

```python
from llmtools import LLMClient, ChemistryAgent

# 1. Create LLM client
client = LLMClient(provider="openai", model="gpt-4o-mini")

# 2. Simple chat
response = client.chat("Suggest a solvent for Suzuki coupling")
print(response.content)

# 3. Chemistry agent (with chemtools)
agent = ChemistryAgent(client, use_chemtools=True)

# 4. Get recommendations
result = agent.suggest_conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    reaction_type="Buchwald_CN"
)
```

## Agent Methods

| Method | Purpose | Example |
|--------|---------|---------|
| `suggest_conditions()` | Recommend conditions | `agent.suggest_conditions(rxn, type)` |
| `explain_mechanism()` | Explain mechanism | `agent.explain_mechanism(rxn, type)` |
| `troubleshoot_reaction()` | Debug problems | `agent.troubleshoot_reaction(rxn, type, problem)` |
| `generate_protocol()` | Create protocol | `agent.generate_protocol(rxn, type, scale)` |
| `search_literature()` | Find papers | `agent.search_literature(rxn, type)` |
| `assess_safety()` | Safety evaluation | `agent.assess_safety(rxn, reagents, solvents)` |

## Prompt Templates

```python
from llmtools.prompts import get_template, list_templates

# List all templates
print(list_templates())

# Use a template
template = get_template("condition_recommendation")
prompt = template.format(
    reaction_smiles="...",
    reaction_type="Suzuki",
    context="Scale-up to 100g"
)

response = client.chat(prompt)
```

## Model Selection

```python
# Auto-select by task type
client = LLMClient.from_env(
    provider="openai",
    model_type="reasoning"  # or "fast", "balanced", "advanced"
)

# Manual selection
client = LLMClient(
    provider="aliyun",
    model="deepseek-r1",
    temperature=0.7,
    max_tokens=2000
)
```

## Temperature Guide

- **0.0-0.3**: Factual (protocols, safety)
- **0.5-0.7**: Balanced (troubleshooting)
- **0.7-1.0**: Creative (retrosynthesis)

## Files

- `llmtools/__init__.py` - Exports
- `llmtools/clients.py` - LLM client
- `llmtools/prompts.py` - Templates
- `llmtools/agents.py` - Agents
- `llmtools/examples.py` - Examples
- `llmtools/README.md` - Full docs

## Testing

```bash
# Test connectivity
python tests/test_llm.py openai --model gpt-4o-mini "Test"

# Run examples
python llmtools/examples.py
```

## Complete Example

```python
#!/usr/bin/env python
from llmtools import LLMClient, ChemistryAgent

# Setup
client = LLMClient.from_env(provider="openai", model_type="balanced")
agent = ChemistryAgent(client, use_chemtools=True)

# Recommend conditions
result = agent.suggest_conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    reaction_type="Buchwald_CN",
    context="Need scalable conditions"
)

print(result["llm_response"])
print(f"Tokens: {result['tokens']}, Precedents: {result['precedents_used']}")

# Generate protocol
protocol = agent.generate_protocol(
    reaction=result["reaction"],
    reaction_type=result["reaction_type"],
    scale="100 g"
)

print(protocol["protocol"])
```
