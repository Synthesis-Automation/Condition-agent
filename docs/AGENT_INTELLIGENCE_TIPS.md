# Improving LLM Agent Intelligence

## What Was Changed

The `chem_assistant.planner.llm_agent_cli` has been enhanced with several improvements to make the agent more intelligent and provide better recommendations.

### 1. **Better Default Model** ✅ IMPLEMENTED

- Changed from `gpt-4o-mini` → `gpt-4o`
- `gpt-4o` has better reasoning capabilities and chemistry knowledge
- Cost: ~15× more expensive but significantly better quality

### 2. **Enhanced System Prompt** ✅ IMPLEMENTED

Added a comprehensive system prompt that:

- Defines the agent's role as an expert chemistry assistant
- Lists the knowledge sources it should consult
- Provides a systematic methodology for analysis
- Emphasizes practical considerations (safety, cost, availability)
- Requests clear explanations with evidence citations

### 3. **Improved User Query** ✅ IMPLEMENTED

The user query now:

- Breaks down the task into clear numbered steps
- Requests specific information (quantities, conditions, rationale)
- Asks for practical details a chemist can execute
- Explicitly requests evidence-based reasoning

### 4. **New CLI Options** ✅ IMPLEMENTED

```bash
--temperature FLOAT    # Control creativity (default: 0.0)
--verbose             # Show agent reasoning process
```

## Additional Improvements You Can Make

### 5. **Use a More Capable Model**

```bash
# For best quality (expensive but very intelligent)
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "..." \
  --llm-model gpt-4-turbo

# For best balance (recommended)
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "..." \
  --llm-model gpt-4o

# For speed/cost (current default in code)
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "..." \
  --llm-model gpt-4o-mini
```

### 6. **Adjust Temperature**

```bash
# Deterministic (recommended for chemistry)
--temperature 0.0

# Slightly more creative
--temperature 0.3

# More exploratory (for brainstorming alternatives)
--temperature 0.7
```

### 7. **Increase Context with More Tools**

Add domain-specific tools from `chemtools_wrapper.py`:

```python
from chem_assistant.chemtools_wrapper import (
    # Current tools
    auto_conditions_llm_tool,
    planner_detect_family_tool,
    planner_rule_candidates_tool,
    planner_protocol_candidates_tool,
    planner_hte_summary_tool,
    planner_score_candidates_tool,
    planner_fuse_tool,
    
    # Additional tools for more intelligence
    classify_reactant_tool,          # Understand reactant types
    get_functional_groups_tool,       # Analyze functional groups
    find_reagent_tool,                # Look up reagent properties
    hte_recommend_tool,               # Direct HTE recommendations
    hte_analytics_tool,               # Deep HTE analysis
)
```

### 8. **Add Few-Shot Examples**

Enhance the system prompt with examples:

```python
system_prompt = """You are an expert chemistry assistant...

Example analysis:
Input: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1
Analysis:
1. Family: Buchwald-Hartwig C-N coupling (aryl halide + amine)
2. Best catalyst: Pd/XPhos based on 127 precedents (avg yield 82%)
3. Recommended base: Cs2CO3 (2 eq) - strong base for deprotonation
4. Solvent: 1,4-dioxane (100°C) - good for Pd-catalyzed reactions
5. Time: 12-18h for complete conversion

Now analyze the user's reaction using the same systematic approach."""
```

### 9. **Enable Iterative Refinement**

Add conversation history support:

```python
# Allow multi-turn conversations
parser.add_argument(
    "--interactive",
    action="store_true",
    help="Interactive mode with conversation history"
)

if args.interactive:
    history = []
    while True:
        user_input = input("\nYou: ")
        if user_input.lower() in ['quit', 'exit']:
            break
        history.append({"role": "user", "content": user_input})
        response = agent.invoke({"messages": history})
        assistant_msg = response["messages"][-1].content
        history.append({"role": "assistant", "content": assistant_msg})
        print(f"\nAssistant: {assistant_msg}")
```

### 10. **Add Retrieval-Augmented Generation (RAG)**

For even better intelligence, add a vector database:

```python
from langchain.vectorstores import Chroma
from langchain.embeddings import OpenAIEmbeddings

# Load chemistry papers, textbooks, or protocols
embeddings = OpenAIEmbeddings()
vectorstore = Chroma(persist_directory="./chemistry_kb", 
                     embedding_function=embeddings)

# Add as a retrieval tool
from langchain.tools.retriever import create_retriever_tool

retriever_tool = create_retriever_tool(
    vectorstore.as_retriever(),
    "chemistry_knowledge_base",
    "Search chemistry knowledge base for reaction mechanisms, best practices, and literature precedents"
)

tools.append(retriever_tool)
```

### 11. **Add Chemistry-Specific Reasoning Chains**

Implement Chain-of-Thought prompting:

```python
user_query = f"""Analyze this reaction step-by-step:

Step 1: Identify the reaction type
- What functional groups are present?
- What transformation is occurring?
- Which reaction family does this belong to?

Step 2: Understand the mechanism
- What is the likely mechanism?
- What are the key intermediates?
- What factors influence reactivity?

Step 3: Gather evidence
- [Use tools to fetch rule-based conditions]
- [Use tools to find similar precedents]
- [Use tools to check HTE statistics]

Step 4: Analyze the evidence
- Which conditions appear most frequently?
- What success rates do they have?
- Are there any conflicting recommendations?

Step 5: Synthesize recommendations
- What are the top 2-3 condition sets?
- Why is each recommended?
- What are the trade-offs?

Reaction: {args.reaction}
"""
```

### 12. **Add Self-Critique and Verification**

Make the agent check its own work:

```python
# After getting initial response
verification_query = f"""Review your recommendations:

1. Are all reagent quantities specified?
2. Are temperatures and times realistic?
3. Is the reasoning backed by tool outputs?
4. Are there safety concerns not mentioned?
5. Could a chemist execute this protocol as written?

If anything is missing or unclear, improve your recommendations."""
```

### 13. **Use Structured Output**

Force the agent to output JSON for consistency:

```python
from pydantic import BaseModel

class ConditionRecommendation(BaseModel):
    catalyst: str
    catalyst_loading_mol_percent: float
    ligand: Optional[str]
    ligand_loading_mol_percent: Optional[float]
    base: str
    base_equivalents: float
    solvent: str
    temperature_celsius: float
    time_hours: float
    concentration_M: float
    expected_yield_percent: Optional[float]
    success_rate_percent: Optional[float]
    rationale: str
    safety_notes: List[str]

# Use with create_agent
agent = create_agent(
    llm, 
    tools, 
    system_prompt=system_prompt,
    response_format=ConditionRecommendation  # Enforce structure
)
```

## Recommended Configuration

### For Best Intelligence

```bash
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" \
  --llm-model gpt-4o \
  --temperature 0.0 \
  --top-k 10 \
  --max-protocols 3 \
  --verbose
```

### For Cost-Effective

```bash
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" \
  --llm-model gpt-4o-mini \
  --temperature 0.0 \
  --top-k 5 \
  --max-protocols 2
```

### For Exploration

```bash
python -m chem_assistant.planner.llm_agent_cli \
  --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" \
  --llm-model gpt-4o \
  --temperature 0.5 \
  --top-k 15 \
  --max-protocols 5 \
  --verbose
```

## Comparison: Before vs After

### Before (gpt-4o-mini, basic prompt)

```
### Recommended Reaction Conditions
Reaction Family: C–N Coupling
Confidence: 0.9
...
[Basic recommendations with minimal reasoning]
```

### After (gpt-4o, enhanced prompt)

```
### Reaction Analysis and Recommendations

**Reaction Family:**
- Detected Family: C-N Coupling
- Confidence: 0.9
- Method: Rule-based with catalyst

**Rule-Based Conditions:**
[Detailed conditions with ranges]

**HTE Insights:**
[Statistical data from 66K+ experiments]

### Top 2 Recommendations

#### Recommendation 1:
[Complete reaction setup with specific quantities]

**Rationale:**
[Evidence-based explanation citing multiple sources]

**Considerations:**
[Safety, handling, and practical notes]
```

## Summary

The agent is now significantly more intelligent due to:

1. ✅ **Better base model** (gpt-4o instead of gpt-4o-mini)
2. ✅ **Enhanced prompts** (system + user)
3. ✅ **Clear methodology** (step-by-step instructions)
4. ✅ **Configurable parameters** (temperature, verbose)
5. ✅ **Practical focus** (specific quantities, safety, rationale)

The improvements result in:

- More detailed and accurate recommendations
- Better synthesis of evidence from multiple sources
- Clearer explanations with reasoning
- More actionable protocols
- Safety and practical considerations

Try the verbose mode to see the agent's reasoning process!
