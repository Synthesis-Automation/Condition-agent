# ChemTools LangChain Integration

This directory contains a **LangChain/LangGraph wrapper** for ChemTools, exposing existing chemistry functions as AI-callable tools and providing an adaptive agent for intelligent chemistry analysis. It now also ships a deterministic **auto-conditions** pipeline (rules → DRFP precedents → HTE summary → fusion → automation formatting) that can run without an LLM when you want fast, reproducible recommendations.

## 🎯 Overview

The integration enables:

- **LangChain Tools**: ChemTools functions wrapped as LangChain `@tool` decorators
- **LangGraph agent**: Automatically selects tool-calling or ReAct planning to solve chemistry problems
- **No Code Changes**: Pure wrapper - doesn't modify existing chemtools code
- **Interactive CLI**: Command-line interface for conversational chemistry queries

## 📁 Structure

```

lang_chain/
├── chemtools_wrapper.py   # LangChain tool wrappers for chemtools
├── chemtools_agent.py     # LangGraph agent wrapper
├── chemtools_cli.py       # Interactive CLI
├── lang_test.py           # Original LangChain test (DataGen example)
├── planner/               # Auto-conditions pipelines (deterministic + LLM-assisted)
├── planner/cli.py         # CLI for auto-conditions (deterministic)
├── planner/llm_agent_cli.py# LLM-assisted auto-conditions (falls back to deterministic)
└── README.md              # This file
```

## 🛠️ Available Tools

The wrapper exposes 10 core ChemTools functions as LangChain tools:

### SMILES & Structure Tools

1. **normalize_smiles_tool** - Canonicalize SMILES strings
2. **normalize_reaction_tool** - Canonicalize reaction SMILES

### Analysis Tools

3. **detect_reaction_family_tool** - Detect reaction type (Suzuki, Buchwald, etc.)
4. **classify_reactant_tool** - Classify reactant types (aryl halide, amine, etc.)
5. **get_functional_groups_tool** - Detect 80+ functional groups

### Recommendation Tools

6. **recommend_conditions_tool** - ML-based condition recommendations
7. **search_precedents_tool** - Find similar precedent reactions
8. **list_supported_cores_tool** - Summarize catalyst cores found in similar precedents

### Database Tools

9. **find_reagent_tool** - Look up reagent information (CAS, properties, roles)
10. **add_reagent_tool** - Add or dry-run reagent entries in the taxonomy registry

## 🚀 Quick Start

### Installation

First, install the required LangChain packages:

```powershell
pip install langchain langchain-openai langgraph python-dotenv
```

### Environment Setup

Create a `.env` file or set environment variables:

```bash
# Required: Choose your provider
LLM_PROVIDER=openai          # or "aliyun"
LLM_MODEL=gpt-4o             # or other model (gpt-4o-mini, gpt-4, etc.)

# OpenAI (if provider=openai)
OPENAI_API_KEY=sk-your-key-here

# Aliyun (if provider=aliyun)
ALIYUN_API_KEY=sk-your-key-here
```

### Usage Examples

#### 0. Deterministic Auto-Conditions (no LLM)

```bash
# From project root
python -m chem_assistant.planner.cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" \
  --top-k 3 \
  --max-protocols 2 \
  --json
```

Outputs the detected family, rule/precedent candidates, optional HTE summary, and formatted protocol additions (per `docs/AUTOMATION_FORMAT.md`). Add `--out result.json` to save the payload.

#### 0b. LLM-Assisted Auto-Conditions (with fallback)

```bash
# Uses OpenAI if OPENAI_API_KEY is set; otherwise falls back to deterministic summary
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" \
  --top-k 3 \
  --max-protocols 2
```

If no LLM/key is present, it prints the deterministic family, counts, and top protocol steps. With an LLM/key, it drives the `auto_conditions_llm_tool` via a simple ReAct agent.

#### 1. Interactive CLI (Recommended)

```powershell
# From project root
python -m lang_chain.chemtools_cli

# Or from lang_chain directory
cd lang_chain
python chemtools_cli.py

# With verbose mode
python chemtools_cli.py --verbose
```

**Example Session:**

```
You: Recommend conditions for Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1

⠋ Agent is thinking...  [animated spinner while processing]

Agent: This is a Buchwald-Hartwig C-N coupling reaction. Based on precedent 
analysis, I recommend:

Catalyst: Pd/XPhos
Base: Cs2CO3
Solvent: 1,4-dioxane
Temperature: 100°C
Time: 12 hours

These conditions show high success rates (85%+ yield) in similar reactions...

You: What is Cs2CO3?

Agent: Cesium carbonate (Cs2CO3) is an inorganic base with:
- CAS: 534-17-8
- Role: base
- SMILES: [Cs+].[Cs+].[O-]C([O-])=O
- Commonly used in cross-coupling reactions...
```

#### 2. Python API

**Basic Usage:**

```python
from lang_chain.chemtools_agent import ChemToolsAgent

# Create agent
agent = ChemToolsAgent()

# Ask a question
response = agent.run("Normalize this SMILES: c1ccccc1")
print(response)

# With conversation history
history = []
response, history = agent.chat("What is the CAS for Cs2CO3?", history)
response, history = agent.chat("What role does it play?", history)
```

**Quick One-Shot Query:**

```python
from lang_chain.chemtools_agent import quick_query

result = quick_query("Recommend conditions for Suzuki coupling")
print(result)
```

**Custom Configuration:**

```python
from lang_chain.chemtools_agent import ChemToolsAgent

agent = ChemToolsAgent(
    provider="openai",
    model="gpt-4",
    temperature=0.1,
    verbose=True
)

response = agent.run("Classify this reactant: Brc1ccccc1")
```

#### 3. Direct Tool Usage (Advanced)

```python
from lang_chain.chemtools_wrapper import (
    normalize_smiles_tool,
    detect_reaction_family_tool,
    recommend_conditions_tool,
    find_reagent_tool
)

# Use tools directly (not through agent)
result = normalize_smiles_tool.invoke({"smiles": "c1ccccc1"})
print(result)

family = detect_reaction_family_tool.invoke({
    "reaction_smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
})
print(family)
```

#### 4. Custom Agent with Subset of Tools

```python
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_agent
from lang_chain.chemtools_wrapper import (
    normalize_smiles_tool,
    recommend_conditions_tool
)

llm = ChatOpenAI(model="gpt-4o", temperature=0)

# Agent with only 2 tools
custom_agent = create_agent(
    llm,
    [normalize_smiles_tool, recommend_conditions_tool],
    prompt="You are a chemistry assistant focused on SMILES and recommendations."
)

result = custom_agent.invoke({
    "messages": [{"role": "user", "content": "Normalize CCO"}]
})
```

Older LangGraph releases only expose `create_react_agent`; the ChemTools wrapper falls back to it automatically when the newer helper is unavailable.

## Key Features (Updated)

- ✅ **No code changes** to existing chemtools
- ✅ **LangGraph agent** combines LLM reasoning with deterministic chemistry
- ✅ **26 chemistry tools** wrapped and ready to use
- ✅ **Conversational interface** with history
- ✅ **Animated spinner** provides visual feedback during processing
- ✅ **Multi-provider support** (OpenAI, Aliyun)
- ✅ **Comprehensive docs** with examples
- ✅ **Interactive CLI** for exploration
- ✅ **Deterministic auto-conditions** pipeline with rule/precedent/HTE fusion and automation formatting

### normalize_smiles_tool

```python
normalize_smiles_tool.invoke({"smiles": "c1ccccc1"})
# Returns: "c1ccccc1"
```

### normalize_reaction_tool

```python
normalize_reaction_tool.invoke({
    "reaction_smiles": "CCBr.CCO>>CCOCC"
})
# Returns: normalized reaction SMILES
```

### Auto-conditions CLI (deterministic)

```bash
python -m chem_assistant.planner.cli --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1"
# Prints family, candidate counts, and top protocol steps
```

### detect_reaction_family_tool

```python
detect_reaction_family_tool.invoke({
    "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
})
# Returns: {"success": True, "family": "...", ...}
```

### classify_reactant_tool

```python
classify_reactant_tool.invoke({"smiles": "Brc1ccccc1"})
# Returns: {"category": "aryl_halide", ...}
```

### get_functional_groups_tool

```python
get_functional_groups_tool.invoke({
    "smiles": "CCO"
})
# Returns: {"success": True, "alcohol": true, ...}
```

### recommend_conditions_tool

```python
recommend_conditions_tool.invoke({
    "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "k": 25,
    "max_variants": 3,
    "rerank_strategy": "rule",
    "constraint_text": "Pd-free, prefer copper",
    "allow_metals": ["Cu", "Ni"],
    "constraint_rules": {"no_chlorinated": true}
})
# Returns: {"success": True, "recommendation": {...}, "alternatives": {...}, ...}
Optional parameters:
- `constraint_text`: Parse natural language requests such as 'Pd-free' or 'prefer copper'.
- `allow_metals` / `exclude_metals` / `prefer_metals`: Structured overrides for catalyst cores.
- `constraint_rules`: Forward deterministic solvent/base rules (e.g., `{"no_chlorinated": true}`).
- `search_all_families`: Set to `true` to perform cross-family precedent search.

Response metadata:
- `cache_hit`: Indicates whether the recommendation came from the in-process cache.
- `timing_ms`: Elapsed time for the tool call (includes cache lookup when relevant).
```

### list_supported_cores_tool

```python
list_supported_cores_tool.invoke({
    "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "k": 25
})
# Returns: {"success": True, "core_candidates": [...], ...}

Response metadata:
- `cache_hit`: True when the catalyst summary reused a cached recommendation pack.
- `timing_ms`: Elapsed time for the underlying recommendation fetch.
```

### add_reagent_tool

```python
add_reagent_tool.invoke({
    "cas": "50-00-0",
    "name": "Formaldehyde",
    "role": "other_reagent",
    "allow_default_family": True,
    "dry_run": True
})
# Returns: {"success": True, "status": "dry_run", ...} (use dry_run=False to persist)
```

Key parameters:

- `allow_default_family`: Permit fallback when the heuristics cannot confidently assign a family.
- `dry_run`: Preview the entry without writing to disk (recommended before persistence).
- `taxonomy_dir`: Optional path to a writable taxonomy directory (defaults to project data).
- `auto_resolve`: Use CAS lookup to fill missing fields when user input is incomplete.

### search_precedents_tool

```python
search_precedents_tool.invoke({
    "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "k": 10,
    "family": "Buchwald_CN"  # optional
})
# Returns: {"success": True, "precedents": [...], ...}
```

### find_reagent_tool

```python
find_reagent_tool.invoke({"query": "Cs2CO3"})
# Returns: {"name": "Cesium carbonate", "cas": "534-17-8", ...}
```

## 🔧 Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `LLM_PROVIDER` | Provider: "openai" or "aliyun" | `openai` |
| `LLM_MODEL` | Model name | `gpt-4o` |
| `OPENAI_API_KEY` | OpenAI API key | Required if provider=openai |
| `OPENAI_BASE_URL` | OpenAI base URL | `https://api.openai.com/v1` |
| `ALIYUN_API_KEY` | Aliyun API key | Required if provider=aliyun |
| `ALIYUN_BASE_URL` | Aliyun base URL | Aliyun default |

### Agent Parameters

```python
ChemToolsAgent(
    provider="openai",      # LLM provider
    model="gpt-4o",         # Model name
    temperature=0,          # Sampling temperature (0=deterministic)
    system_prompt=None,     # Custom system prompt
    verbose=False           # Print debug info
)
```

### Tool Parameters

Most tools accept standard parameters. Key ones:

**recommend_conditions_tool:**

- `k` (int): Number of precedents to retrieve (default: 25)
- `max_variants` (int): Max condition variants (default: 3)
- `rerank_strategy` (str): "rule", "analytics", or "none" (default: "rule")

**search_precedents_tool:**

- `k` (int): Number of precedents to find (default: 10)
- `family` (str, optional): Filter by reaction family

## 🎓 How It Works

### Architecture

```
User Query
    ↓
LangGraph ReAct Agent (reasoning + planning)
    ↓  [Shows animated spinner while working]
    ↓
LangChain Tool Wrappers (@tool decorators)
    ↓
ChemTools Functions (unmodified)
    ↓
Results returned to Agent
    ↓
Agent formats and explains to user
```

### ReAct Pattern

The agent uses the ReAct (Reasoning + Acting) pattern:

1. **Thought**: Reasons about what tools to use
2. **Action**: Calls appropriate ChemTools functions
3. **Observation**: Receives tool results
4. **Thought**: Decides next step or provides answer
5. Repeat until complete

### Example Agent Reasoning

```
User: "Recommend conditions for Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

Agent Thought: I need to analyze this reaction and recommend conditions.
                First, let me detect the reaction family.

Agent Action: detect_reaction_family_tool("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")

Agent Observation: {"family": "Buchwald_CN", "confidence": "high"}

Agent Thought: This is a Buchwald C-N coupling. Now I'll get recommendations.

Agent Action: recommend_conditions_tool(reaction_smiles="...", k=25)

Agent Observation: {
  "recommendation": {"core": "Pd/XPhos", "base": "Cs2CO3", ...},
  "alternatives": {...}
}

Agent Answer: Based on precedent analysis of Buchwald-Hartwig C-N couplings,
              I recommend using Pd/XPhos catalyst with Cs2CO3 base in
              1,4-dioxane at 100°C for 12 hours...
```

## 🧪 Testing

### Test Individual Tools

```python
from lang_chain.chemtools_wrapper import print_tool_summary

# See all available tools
print_tool_summary()
```

### Test Agent

```python
from lang_chain.chemtools_agent import ChemToolsAgent

agent = ChemToolsAgent(verbose=True)

test_queries = [
    "Normalize c1ccccc1",
    "What functional groups are in CCO?",
    "Recommend conditions for Suzuki coupling"
]

for query in test_queries:
    print(f"\nQuery: {query}")
    response = agent.run(query)
    print(f"Response: {response}\n")
```

### Run Agent Tests

```powershell
python lang_chain/chemtools_agent.py
```

## 📝 CLI Commands

When using the interactive CLI:

| Command | Description |
|---------|-------------|
| `help` | Show example queries |
| `tools` | List available tools |
| `clear` or `cls` | Clear conversation history |
| `verbose on/off` | Toggle verbose mode |
| `quit`, `exit`, or `q` | Exit CLI |

## 🔍 Example Use Cases

### 1. Reaction Condition Discovery

```
You: I want to couple bromobenzene with aniline

Agent: [detects Buchwald coupling, recommends Pd/XPhos, Cs2CO3, dioxane, 100°C]
```

### 2. Reagent Information Lookup

```
You: What is the CAS number for cesium carbonate?

Agent: [looks up Cs2CO3, returns CAS: 534-17-8, role: base, etc.]
```

### 3. Structure Analysis

```
You: What functional groups are in c1ccc(O)cc1?

Agent: [normalizes SMILES, detects phenol, aromatic OH, etc.]
```

### 4. Precedent Search

```
You: Find similar reactions to my Suzuki coupling

Agent: [searches precedents, returns top matches with conditions and yields]
```

### 5. Multi-Step Workflows

```
You: Analyze this reaction and recommend conditions: Brc1ccccc1.c1cccnc1B(O)O>>...

Agent: [normalizes → detects family → classifies reactants → recommends conditions → 
        searches precedents → provides comprehensive answer]
```

## 🔗 Integration with Existing Code

This wrapper integrates seamlessly with your existing ChemTools workflows:

```python
# Your existing code
from chemtools.recommend import recommend_from_reaction
result = recommend_from_reaction("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")

# Now with AI assistance
from lang_chain.chemtools_agent import quick_query
ai_result = quick_query(
    "Explain the recommended conditions for: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
)
# AI will call the same function but explain results in natural language
```

## ⚠️ Limitations

1. **API Costs**: Each query uses LLM tokens (keep track of usage)
2. **Latency**: Agent reasoning adds ~2-10s compared to direct function calls
3. **Reliability**: LLM may occasionally misinterpret queries (use verbose mode to debug)
4. **Tool Errors**: If chemtools functions fail, errors are returned to agent

## 🐛 Troubleshooting

### "Missing API key" Error

```powershell
# Set your API key
$env:OPENAI_API_KEY = "sk-your-key-here"
```

### Import Errors

```powershell
# Install missing packages
pip install langchain langchain-openai langgraph python-dotenv
```

### Agent Not Calling Tools

- Increase `recursion_limit` parameter
- Check system prompt is not overridden
- Enable `verbose=True` to see reasoning

### Tool Returns Error

- Test the underlying chemtools function directly first
- Check SMILES format is correct
- Ensure required data files are present

## 📊 Performance

- **Direct chemtools call**: <100ms
- **Agent with tools**: 2-10s (LLM reasoning overhead)
- **Recommendation**: Use agent for exploration, direct calls for production
- **Caching**: Consecutive recommendation + core-list queries reuse cached DRFP results to trim duplicate compute
- **Cache management**: In the CLI, run `cache show` to inspect entries or `cache clear` to flush them
- **Benchmark (Buchwald example)**: Cold recommendation ≈1150 ms, cached follow-up ≈0.01 ms (>99% reduction)

## 🤝 Contributing

To add new tools:

1. Add a new `@tool` function in `chemtools_wrapper.py`
2. Import the chemtools function you want to wrap
3. Add comprehensive docstring (agent uses this!)
4. Add to `CHEMTOOLS_TOOLS` list
5. Update this README

Example:

```python
from pydantic import BaseModel, Field
from langchain_core.tools import tool

class MyToolInput(BaseModel):
    param: str = Field(..., description="Description of parameter")

@tool(args_schema=MyToolInput)
def my_new_tool(param: str) -> dict:
    """
    Detailed description of what this tool does.
    
    Args:
        param: Description of parameter
    
    Returns:
        Structured payload with success flag
    """
    from chemtools.my_module import my_function
    result = my_function(param)
    return {"success": True, "result": result}

# Add to CHEMTOOLS_TOOLS
CHEMTOOLS_TOOLS.append(my_new_tool)
```

## 📖 Related Documentation

- [ChemTools Main README](../docs/README.md)
- [ChemTools API Documentation](../docs/API_DOCUMENTATION.md)
- [LangChain Documentation](https://python.langchain.com/)
- [LangGraph Documentation](https://langchain-ai.github.io/langgraph/)

## 📜 License

Same as parent ChemTools project.

## 🎉 Summary

You now have a powerful AI agent that can:

- ✅ Call ChemTools functions intelligently
- ✅ Reason about chemistry problems
- ✅ Explain results in natural language
- ✅ Maintain conversation context
- ✅ Handle complex multi-step workflows

**All without modifying any existing ChemTools code!**

---

**Quick Start Reminder:**

```powershell
# 1. Install dependencies
pip install langchain langchain-openai langgraph python-dotenv

# 2. Set API key
$env:OPENAI_API_KEY = "sk-your-key-here"

# 3. Run CLI
python -m lang_chain.chemtools_cli
```

Happy chemistry! 🧪✨
