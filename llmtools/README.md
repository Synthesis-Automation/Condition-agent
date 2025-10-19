# LLM Tools

Advanced language model integration for chemistry workflows.

## Overview

`llmtools/` provides intelligent agents that combine large language models (LLMs) with deterministic `chemtools` functions for advanced chemistry operations:

- **Multi-provider LLM support**: OpenAI, Aliyun/DeepSeek, extensible to others
- **Chemistry-specific prompts**: Pre-built templates for common tasks
- **Intelligent agents**: Combine LLM reasoning with chemtools validation
- **Structured outputs**: Parse and validate LLM responses
- **Usage tracking**: Monitor tokens and costs

## Installation

```bash
# Install OpenAI SDK (required)
pip install openai

# Set API keys
export OPENAI_API_KEY="sk-..."
export ALIYUN_API_KEY="sk-..."  # for DeepSeek
```

## Quick Start

### 1. Basic LLM Client

```python
from llmtools import LLMClient

# Create client
client = LLMClient(provider="openai", model="gpt-4o-mini")

# Simple chat
response = client.chat("Suggest a solvent for Suzuki coupling")
print(response.content)
print(f"Tokens: {response.total_tokens}, Latency: {response.latency_ms}ms")
```

### 2. Chemistry Agent (LLM + ChemTools)

```python
from llmtools import LLMClient, ChemistryAgent

# Create agent
client = LLMClient(provider="openai", model="gpt-4o")
agent = ChemistryAgent(client, use_chemtools=True)

# Recommend conditions (uses LLM + precedent database)
result = agent.suggest_conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    reaction_type="Buchwald_CN",
    context="Prefer mild conditions"
)

print(result["llm_response"])
print(f"Precedents used: {result['precedents_used']}")
```

### 3. Using Prompt Templates

```python
from llmtools.prompts import get_template

# Get template
template = get_template("condition_recommendation")

# Format with your data
prompt = template.format(
    reaction_smiles="...",
    reaction_type="Suzuki",
    context="Scale-up to 100g"
)

# Send to LLM
response = client.chat(prompt)
```

## Module Structure

```
llmtools/
├── __init__.py           # Package exports
├── clients.py            # LLM client management
├── prompts.py            # Chemistry prompt templates
├── agents.py             # Intelligent chemistry agents
├── parsers.py            # Output parsing (TODO)
├── reasoning.py          # Multi-step reasoning chains (TODO)
└── cache.py              # Response caching (TODO)
```

## Available Features

### LLM Clients (`clients.py`)

```python
from llmtools import LLMClient

# Auto-select model by task type
client = LLMClient.from_env(
    provider="openai",
    model_type="reasoning"  # or "fast", "balanced", "advanced"
)

# Custom configuration
client = LLMClient(
    provider="aliyun",
    model="deepseek-r1",
    temperature=0.7,
    max_tokens=2000,
    timeout=60
)

# Get usage stats
print(client.get_usage_summary())
```

### Prompt Templates (`prompts.py`)

Available templates:

- `condition_recommendation` - Recommend reaction conditions
- `retrosynthesis` - Plan synthetic routes
- `mechanism` - Explain reaction mechanisms
- `literature` - Search and summarize papers
- `reagent_selection` - Select optimal reagents
- `safety` - Assess safety and hazards
- `troubleshooting` - Debug reactions
- `spectroscopy` - Interpret NMR/MS/IR data
- `protocol` - Generate experimental procedures

```python
from llmtools.prompts import list_templates, get_template

# See all available templates
print(list_templates())

# Use a template
template = get_template("mechanism")
prompt = template(
    reaction_smiles="...",
    reaction_type="Suzuki",
    reagents="Pd(PPh3)4, K2CO3"
)
```

### Chemistry Agents (`agents.py`)

**ChemistryAgent** - General-purpose assistant:

```python
agent = ChemistryAgent(client)

# Recommend conditions
agent.suggest_conditions(reaction, reaction_type)

# Explain mechanism
agent.explain_mechanism(reaction, reaction_type)

# Troubleshoot problems
agent.troubleshoot_reaction(reaction, reaction_type, problem="Low yield")

# Generate protocol
agent.generate_protocol(reaction, reaction_type, scale="10 mmol")

# Search literature
agent.search_literature(reaction, reaction_type)

# Assess safety
agent.assess_safety(reaction, reagents="...", solvents="...")
```

**ConditionOptimizer** - Iterative optimization (TODO):

```python
optimizer = ConditionOptimizer(client)
result = optimizer.optimize(reaction, reaction_type, constraints={...})
```

**RetrosynthesisPlanner** - Route planning:

```python
planner = RetrosynthesisPlanner(client)
route = planner.plan_route(target_smiles, max_steps=5)
```

**ReactionAnalyzer** - Comprehensive analysis:

```python
analyzer = ReactionAnalyzer(client)
analysis = analyzer.analyze(
    reaction,
    reaction_type,
    include_mechanism=True,
    include_literature=True,
    include_safety=True
)
```

## Supported LLM Providers

### OpenAI

```bash
export OPENAI_API_KEY="sk-..."
export OPENAI_BASE_URL="https://api.openai.com/v1"  # optional
```

Models: `gpt-4o`, `gpt-4o-mini`, `gpt-5-mini`, `o3-mini`, etc.

### Aliyun (DeepSeek)

```bash
export ALIYUN_API_KEY="sk-..."
export ALIYUN_BASE_URL="https://dashscope.aliyuncs.com/compatible-mode/v1"  # optional
```

Models: `deepseek-r1`, `deepseek-v3`, `deepseek-r1-distill-qwen-7b`, etc.

## Integration with ChemTools

Agents automatically integrate with `chemtools` when available:

```python
agent = ChemistryAgent(client, use_chemtools=True)

# This will:
# 1. Query chemtools ML precedent database
# 2. Format precedents as context
# 3. Send to LLM with precedent information
# 4. Return LLM recommendations enhanced by precedents

result = agent.suggest_conditions(reaction, reaction_type)
```

If `chemtools` is not available, agents run in LLM-only mode.

## Examples

### Example 1: Condition Recommendation with Precedents

```python
from llmtools import LLMClient, ChemistryAgent

client = LLMClient.from_env(provider="openai", model_type="balanced")
agent = ChemistryAgent(client, use_chemtools=True, verbose=True)

result = agent.suggest_conditions(
    reaction="BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1",
    reaction_type="Suzuki",
    context="Need mild conditions for scale-up"
)

print(result["llm_response"])
```

### Example 2: Mechanism Explanation

```python
result = agent.explain_mechanism(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    reaction_type="Buchwald_CN",
    reagents="Pd(OAc)2, BINAP, NaOtBu",
    detail_level="Detailed with all intermediates"
)

print(result["mechanism_explanation"])
```

### Example 3: Reaction Troubleshooting

```python
result = agent.troubleshoot_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    reaction_type="Buchwald_CN",
    problem="Getting only 10% yield after 24h",
    current_conditions="Pd(PPh3)4 (5 mol%), K2CO3 (2 eq), toluene, 80°C",
    observed_result="10% conversion by TLC, starting material recovered"
)

print(result["troubleshooting_advice"])
```

### Example 4: Generate Experimental Protocol

```python
conditions = {
    "catalyst": "Pd(OAc)2 (2 mol%)",
    "ligand": "BINAP (4 mol%)",
    "base": "NaOtBu (1.5 eq)",
    "solvent": "Toluene",
    "temperature": "100°C",
    "time": "16 h"
}

result = agent.generate_protocol(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    reaction_type="Buchwald_CN",
    conditions=conditions,
    scale="5 mmol"
)

print(result["protocol"])
```

## Roadmap

### Implemented ✅

- [x] Multi-provider LLM client
- [x] Chemistry prompt templates
- [x] Basic chemistry agent
- [x] ChemTools integration
- [x] Usage tracking

### In Progress 🚧

- [ ] Structured output parsing
- [ ] Multi-step reasoning chains
- [ ] Response caching
- [ ] Condition optimizer
- [ ] Batch processing

### Planned 📋

- [ ] Fine-tuned chemistry models
- [ ] RAG (retrieval-augmented generation) for literature
- [ ] Automated experiment design
- [ ] Integration with protocol database
- [ ] Web interface for agents

## Testing

```bash
# Test basic LLM connectivity
python tests/test_llm.py openai --model gpt-4o-mini "Suggest a solvent"

# Test with DeepSeek
python tests/test_llm.py aliyun --model deepseek-r1 "Explain Suzuki coupling"
```

## Best Practices

1. **Choose the right model**:

   - Reasoning tasks: `o3-mini`, `deepseek-r1`
   - Fast responses: `gpt-4o-mini`, `deepseek-r1-distill-qwen-7b`
   - Balanced: `gpt-4o`, `deepseek-v3`

2. **Use appropriate temperature**:

   - Factual tasks (protocols, safety): `0.0-0.3`
   - Creative tasks (retrosynthesis): `0.7-1.0`
   - Balanced (troubleshooting): `0.5-0.7`

3. **Leverage chemtools precedents**:

   - Always use `use_chemtools=True` when available
   - Precedents ground LLM recommendations in real data

4. **Monitor usage**:

   ```python
   print(agent.get_usage_summary())
   ```

5. **Handle errors gracefully**:
   ```python
   try:
       result = agent.suggest_conditions(...)
   except Exception as e:
       print(f"LLM request failed: {e}")
   ```

## License

Same as main Condition-agent project.

## Contributing

See main AGENTS.md for contribution guidelines.
