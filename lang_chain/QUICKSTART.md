# ChemTools LangChain Integration - Quick Start

## What is This?

A wrapper that exposes ChemTools functions as **LangChain tools** and provides a **ReAct agent** for intelligent chemistry analysis.

## Installation

```powershell
# Install LangChain dependencies
pip install -r lang_chain/requirements.txt
```

## Setup API Key

```powershell
# PowerShell (Windows)
$env:OPENAI_API_KEY = "sk-your-key-here"

# Or create .env file
echo 'OPENAI_API_KEY=sk-your-key-here' > .env
```

## Usage

### Option 1: Interactive CLI (Recommended)

```powershell
python -m lang_chain.chemtools_cli
```

### Option 2: Python API

```python
from lang_chain.chemtools_agent import quick_query

result = quick_query("Recommend conditions for Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")
print(result)
```

### Option 3: Run Examples

```powershell
python lang_chain/example_usage.py
```

## Available Tools

1. **normalize_smiles_tool** - Canonicalize SMILES
2. **normalize_reaction_tool** - Canonicalize reactions
3. **detect_reaction_family_tool** - Detect reaction type
4. **classify_reactant_tool** - Classify reactants
5. **get_functional_groups_tool** - Detect functional groups
6. **recommend_conditions_tool** - Get ML recommendations
7. **search_precedents_tool** - Find similar reactions
8. **find_reagent_tool** - Look up reagent info

## Example Queries

- "Normalize c1ccccc1"
- "What reaction is Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1?"
- "Recommend conditions for Suzuki coupling"
- "What is the CAS for Cs2CO3?"
- "What functional groups are in CCO?"

## Key Features

✅ No code changes to existing chemtools  
✅ AI reasoning + deterministic chemistry  
✅ Conversational interface with history  
✅ **Animated spinner** shows agent progress  
✅ 8 chemistry tools available  
✅ Supports OpenAI and Aliyun providers  

## Learn More

- Full documentation: [README.md](README.md)
- Example code: [example_usage.py](example_usage.py)
- Original test: [lang_test.py](lang_test.py)

## Architecture

```
User Query
    ↓
LangGraph ReAct Agent (reasoning)
    ↓
LangChain Tools (wrappers)
    ↓
ChemTools Functions (unchanged)
    ↓
AI-formatted Response
```

## Files

- `chemtools_wrapper.py` - LangChain tool wrappers
- `chemtools_agent.py` - ReAct agent
- `chemtools_cli.py` - Interactive CLI
- `example_usage.py` - Example code
- `requirements.txt` - Dependencies
- `README.md` - Full docs (you are here!)
- `QUICKSTART.md` - This file

## Troubleshooting

**Missing API key?**

```powershell
$env:OPENAI_API_KEY = "sk-your-key-here"
```

**Import errors?**

```powershell
pip install langchain langchain-openai langgraph python-dotenv
```

**Agent not working?**

- Enable verbose mode: `ChemToolsAgent(verbose=True)`
- Check API key is set
- Ensure internet connection

## Next Steps

1. Try the CLI: `python -m lang_chain.chemtools_cli`
2. Run examples: `python lang_chain/example_usage.py`
3. Read full docs: `lang_chain/README.md`
4. Build custom agents with subset of tools

Happy chemistry! 🧪✨
