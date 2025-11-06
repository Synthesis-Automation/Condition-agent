# Agent Tool Selection and Calling Workflow

## 🔄 Complete Workflow Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                        USER INPUT                               │
│  "What conditions should I use for Brc1ccccc1.Nc1ccccc1>>?"     │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                    ChemToolsCLI                                 │
│  - Parses commands (help, tools, quit, constraints, etc.)       │
│  - Manages conversation history                                 │
│  - Displays spinner animation                                   │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                  ChemToolsAgent.run()                           │
│  - Merges session + call constraints                            │
│  - Adds constraint prompt to messages                           │
│  - Calls LangGraph ReAct agent                                  │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│            LangGraph create_react_agent()                       │
│  - Reasoning and Acting loop                                    │
│  - LLM decides which tools to call                              │
│  - Executes tool calls                                          │
│  - Returns to LLM for next step                                 │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                  LLM (GPT-4 / DeepSeek)                         │
│  - Reads system prompt with tool descriptions                   │
│  - Analyzes user query                                          │
│  - Decides which tool(s) to call                                │
│  - Generates tool arguments                                     │
└────────────────────────┬────────────────────────────────────────┘
                         │
         ┌───────────────┴───────────────┐
         ▼                               ▼
    [REASONING]                     [TOOL CALL]
         │                               │
         │                               ▼
         │              ┌─────────────────────────────────┐
         │              │   Tool Execution Router         │
         │              │   (LangChain dispatches)        │
         │              └────────┬────────────────────────┘
         │                       │
         │         ┌─────────────┼─────────────┐
         │         ▼             ▼             ▼
         │    Tool A         Tool B        Tool C
         │   (e.g.,       (e.g.,        (e.g.,
         │  normalize)   recommend)    find_reagent)
         │         │             │             │
         │         └─────────────┴─────────────┘
         │                       │
         │                       ▼
         │         ┌────────────────────────────┐
         │         │   ChemTools Functions      │
         │         │   - chemtools.smiles       │
         │         │   - chemtools.recommend    │
         │         │   - chemtools.rule         │
         │         │   - chemtools.reagent      │
         │         └────────────┬───────────────┘
         │                      │
         │                      ▼
         │         ┌────────────────────────────┐
         │         │   Tool Returns JSON        │
         │         └────────────┬───────────────┘
         │                      │
         └──────────────────────┘
                    │
                    ▼
         ┌────────────────────────────┐
         │   LLM Interprets Results   │
         │   Decides: Continue or End │
         └────────────┬───────────────┘
                      │
        ┌─────────────┴─────────────┐
        ▼                           ▼
   [MORE TOOLS]                [FINAL ANSWER]
        │                           │
        └───────────────────────────┘
                    │
                    ▼
         ┌────────────────────────────┐
         │   Return to ChemToolsAgent │
         │   Extract final message    │
         └────────────┬───────────────┘
                      │
                      ▼
         ┌────────────────────────────┐
         │   Display to User (CLI)    │
         └────────────────────────────┘
```

---

## 📋 Detailed Component Breakdown

### 1. **Entry Point: ChemToolsCLI**

**File**: `chem_assistant/chemtools_cli.py`

```python
class ChemToolsCLI:
    def __init__(self, verbose=False):
        self.agent = ChemToolsAgent(verbose=verbose)  # ← Creates agent
        self.history = []
        self.constraint_spec = ConstraintSpec()

    def run(self):
        """Main CLI loop"""
        while True:
            user_input = self.get_user_input_with_indicator()

            # Handle commands
            if user_input in ["quit", "exit", "q"]:
                break
            elif user_input == "tools":
                self.print_tools()
                continue
            elif user_input.startswith("constraints"):
                self.handle_constraints_command(user_input)
                continue

            # Send to agent
            spinner = Spinner("Agent thinking")
            spinner.start()

            response = self.agent.run(
                user_input,
                history=self.history,
                constraints=self.constraint_spec
            )

            spinner.stop()
            print(f"\n{Colors.OKGREEN}{response}{Colors.ENDC}\n")

            # Update history
            self.history.append(HumanMessage(content=user_input))
            self.history.append(AIMessage(content=response))
```

**Key Responsibilities**:

- ✅ User interface (prompts, colors, spinners)
- ✅ Command parsing (help, tools, constraints, quit)
- ✅ Conversation history management
- ✅ Constraint specification management
- ✅ Delegates actual queries to `ChemToolsAgent`

---

### 2. **Agent Orchestrator: ChemToolsAgent**

**File**: `chem_assistant/chemtools_agent.py`

```python
class ChemToolsAgent:
    def __init__(self, provider=None, model=None, temperature=0):
        # Get LLM client (OpenAI or Aliyun)
        self.llm = get_llm_client(provider, model, temperature)

        # System prompt defines tools and behavior
        self.system_prompt = CHEMISTRY_SYSTEM_PROMPT

        # Create ReAct agent with LangGraph
        self.agent = create_react_agent(
            self.llm,
            CHEMTOOLS_TOOLS,  # ← All tools registered here
            prompt=self.system_prompt
        )

    def run(self, query, history=None, recursion_limit=15, constraints=None):
        """Run agent on a query"""
        messages = list(history or [])

        # Merge constraints and add to context
        active_spec = merge_specs(self.session_constraints, constraints)
        constraint_prompt = format_constraints_for_prompt(active_spec)
        if constraint_prompt:
            messages.append(HumanMessage(content=constraint_prompt))

        messages.append(HumanMessage(content=query))

        # Invoke LangGraph agent
        result = self.agent.invoke(
            {"messages": messages},
            config={"recursion_limit": recursion_limit}
        )

        # Extract final answer
        final_message = result["messages"][-1]
        return final_message.content
```

**Key Responsibilities**:

- ✅ LLM client setup (OpenAI/Aliyun)
- ✅ Creates LangGraph ReAct agent
- ✅ Manages conversation context
- ✅ Handles constraint merging
- ✅ Invokes agent execution loop

---

### 3. **Tool Registry: chemtools_wrapper.py**

**File**: `chem_assistant/chemtools_wrapper.py`

This file defines all available tools using LangChain's `@tool` decorator:

```python
@tool(args_schema=NormalizeSmilesInput)
def normalize_smiles_tool(smiles: str) -> str:
    """
    Canonicalize a SMILES string to a standard form.

    Args:
        smiles: SMILES string to normalize

    Returns:
        Canonicalized SMILES string
    """
    try:
        result = normalize(smiles)
        return json.dumps({"smiles": result}, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e)})


@tool(args_schema=RecommendConditionsInput)
def recommend_conditions_tool(
    reaction_smiles: str,
    k: int = 25,
    max_variants: int = 3,
    # ... constraint parameters ...
) -> str:
    """
    Get ML-based condition recommendations for a reaction.

    Args:
        reaction_smiles: Reaction SMILES (reactants>>products)
        k: Number of precedents to retrieve
        max_variants: Number of recommendation variants
        constraint_text: Natural language constraints
        allow_metals: Whitelist of allowed metals
        exclude_metals: Metals to exclude
        # ... more parameters ...

    Returns:
        JSON with recommendations, detection, precedents
    """
    # Build constraint specification
    constraint_spec = build_constraint_spec(
        text=constraint_text,
        allow_metals=allow_metals,
        exclude_metals=exclude_metals,
        # ...
    )

    # Call chemtools recommendation function
    result = recommend_from_reaction(
        reaction_smiles,
        k=k,
        max_variants=max_variants,
        # ...
    )

    return json.dumps(result, indent=2)


# More tools...
@tool(args_schema=RuleRecommendInput)
def rule_based_conditions_tool(...):
    """Deterministic rule-engine condition recommendation"""
    # ...

@tool(args_schema=FindReagentInput)
def find_reagent_tool(...):
    """Look up reagent information"""
    # ...

# All tools collected in a list
CHEMTOOLS_TOOLS = [
    normalize_smiles_tool,
    normalize_reaction_tool,
    detect_reaction_family_tool,
    classify_reactant_tool,
    get_functional_groups_tool,
    recommend_conditions_tool,
    rule_based_conditions_tool,
    enhanced_cross_family_recommend_tool,
    search_precedents_tool,
    reaction_dataset_analytics_tool,
    find_reagent_tool,
    reagent_database_analytics_tool,
    list_supported_cores_tool,
    add_reagent_tool,
    calculable_features_tool,
    molpipeline_featurize_tool,
    analyze_bond_changes_tool,
]
```

**Key Responsibilities**:

- ✅ Wraps ChemTools functions as LangChain tools
- ✅ Defines input schemas with Pydantic
- ✅ Adds docstrings (tool descriptions for LLM)
- ✅ Handles errors and returns JSON
- ✅ Manages caching (for recommendation tool)
- ✅ Exports `CHEMTOOLS_TOOLS` list

---

### 4. **LangGraph ReAct Loop**

**Library**: `langgraph.prebuilt.create_react_agent`

This is a built-in LangGraph function that creates a ReAct (Reasoning + Acting) agent:

```python
# Pseudo-code of what LangGraph does internally
def create_react_agent(llm, tools, prompt):
    graph = StateGraph()

    # Add nodes
    graph.add_node("agent", lambda state: call_llm(llm, state, tools))
    graph.add_node("tools", lambda state: execute_tools(state))

    # Add edges
    graph.add_edge(START, "agent")
    graph.add_conditional_edge(
        "agent",
        should_continue,  # If tool calls → "tools", else → END
        {"continue": "tools", "end": END}
    )
    graph.add_edge("tools", "agent")  # Loop back

    return graph.compile()
```

**Execution Flow**:

1. **Agent Node** (LLM reasoning):

   - LLM receives: system prompt + conversation history + tool descriptions
   - LLM decides: "I need to call tool X with arguments Y"
   - Returns: `AIMessage` with `tool_calls` attribute

2. **Tools Node** (Tool execution):

   - Extract tool name and arguments from `tool_calls`
   - Find matching tool in `CHEMTOOLS_TOOLS`
   - Execute the tool function
   - Returns: `ToolMessage` with results

3. **Loop**:

   - Results go back to Agent Node
   - LLM sees tool results and decides next action
   - Continues until LLM returns final answer (no more tool calls)

4. **Termination**:
   - When LLM returns `AIMessage` without `tool_calls`
   - Or when `recursion_limit` is reached

---

### 5. **LLM Decision Making**

The LLM (e.g., GPT-4, DeepSeek) receives this context:

**System Prompt** (from `CHEMISTRY_SYSTEM_PROMPT`):

```
You are ChemBot, an expert chemistry assistant with access to powerful ChemTools...

You have access to the following tools:
1. normalize_smiles_tool: Canonicalize SMILES strings
2. normalize_reaction_tool: Canonicalize reaction SMILES
3. detect_reaction_family_tool: Identify reaction type
...

How to help users:
- For reaction condition recommendations:
  1. First normalize the reaction SMILES if needed
  2. Detect the reaction family
  3. Use recommend_conditions_tool
  ...
```

**Conversation History**:

```json
[
  { "role": "system", "content": "You are ChemBot..." },
  { "role": "user", "content": "What conditions for Brc1ccccc1.Nc1ccccc1>>?" }
]
```

**LLM Response** (Tool Call):

```json
{
  "role": "assistant",
  "content": null,
  "tool_calls": [
    {
      "id": "call_abc123",
      "type": "function",
      "function": {
        "name": "normalize_reaction_tool",
        "arguments": "{\"reaction_smiles\": \"Brc1ccccc1.Nc1ccccc1>>\"}"
      }
    }
  ]
}
```

---

### 6. **Tool Execution**

LangGraph extracts the tool call and executes:

```python
# LangGraph internally does:
tool_name = "normalize_reaction_tool"
tool_args = {"reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>"}

# Find the tool
tool_func = CHEMTOOLS_TOOLS_MAP[tool_name]

# Execute
result = tool_func(**tool_args)
# Returns: '{"reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"}'
```

**Tool Result** added to conversation:

```json
{
  "role": "tool",
  "tool_call_id": "call_abc123",
  "name": "normalize_reaction_tool",
  "content": "{\"reaction\": \"Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1\"}"
}
```

---

### 7. **Multi-Turn Loop Example**

**Turn 1**: User asks for conditions

```
User: "What conditions for Brc1ccccc1.Nc1ccccc1>>?"

LLM: I'll normalize this reaction first.
     → Calls: normalize_reaction_tool("Brc1ccccc1.Nc1ccccc1>>")

Tool: Returns normalized reaction

LLM: Now I'll detect the reaction type.
     → Calls: detect_reaction_family_tool(normalized_reaction)

Tool: Returns {"family": "Buchwald_CN", "confidence": 0.95}

LLM: This is a Buchwald-Hartwig C-N coupling. Let me get recommendations.
     → Calls: recommend_conditions_tool(reaction, k=25)

Tool: Returns full recommendation JSON

LLM: Based on the analysis, I recommend...
     [Formats final answer]
```

**Turn 2**: Follow-up question

```
User: "Can I use Cu instead?"

LLM: (Has previous context about Buchwald reaction)
     I'll search for Cu-based alternatives.
     → Calls: recommend_conditions_tool(reaction, allow_metals=["Cu"])

Tool: Returns Cu-focused recommendations

LLM: Yes, you can use Cu catalysts. Here are options...
```

---

## 🔍 Tool Selection Process

### How does the LLM choose which tools to call?

1. **System Prompt Analysis**:

   - LLM reads tool descriptions and parameters
   - Understands when each tool is appropriate

2. **Query Parsing**:

   - Extracts intent: "recommend conditions", "find reagent", "normalize SMILES"
   - Maps intent to relevant tools

3. **Sequential Planning**:

   - System prompt provides workflow guidance:
     ```
     For reaction recommendations:
     1. Normalize reaction SMILES
     2. Detect family
     3. Get recommendations
     ```

4. **Context-Aware**:
   - Uses conversation history
   - Avoids redundant tool calls
   - Can chain multiple tools

### Example Tool Selection Logic:

```
Query: "What's the CAS for Cs2CO3?"
→ LLM thinks: This is a reagent lookup
→ Calls: find_reagent_tool(query="Cs2CO3", reagent_type="base")

Query: "Recommend conditions for Suzuki coupling"
→ LLM thinks: Need specific reaction SMILES first
→ Responds: "Please provide the reaction SMILES..."

Query: "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 with Pd-free"
→ LLM thinks: Condition recommendation with constraint
→ Calls: recommend_conditions_tool(
    reaction_smiles="...",
    exclude_metals=["Pd"]
  )
```

---

## 📊 Key Data Structures

### Message Flow:

```python
# Input to agent
messages = [
    HumanMessage(content="Query"),
    AIMessage(content="Previous response"),
    # ...
]

# Agent returns
{
    "messages": [
        ...original_messages,
        AIMessage(
            content=None,
            tool_calls=[{...}]  # If calling tools
        ),
        ToolMessage(
            content="Tool result",
            tool_call_id="..."
        ),
        AIMessage(
            content="Final answer"  # No more tool_calls
        )
    ]
}
```

### Tool Input/Output:

```python
# Tool receives (via Pydantic schema)
{
    "reaction_smiles": "...",
    "k": 25,
    "constraint_text": "Pd-free"
}

# Tool returns (as JSON string)
{
    "meta": {...},
    "recommendations": [...],
    "detection": {...},
    "precedents": {...}
}
```

---

## 🎯 Configuration Points

### LLM Selection:

**Environment variables**:

```bash
# Use OpenAI
export LLM_PROVIDER=openai
export OPENAI_API_KEY=sk-...
export LLM_MODEL=gpt-4o

# Or use Aliyun/DeepSeek
export LLM_PROVIDER=aliyun
export ALIYUN_API_KEY=sk-...
export LLM_MODEL=deepseek-chat
```

### Tool Configuration:

**In `chemtools_wrapper.py`**:

- Add new tools with `@tool` decorator
- Define Pydantic schemas for inputs
- Add to `CHEMTOOLS_TOOLS` list

### System Prompt Customization:

**In `chemtools_agent.py`**:

```python
agent = ChemToolsAgent(
    system_prompt="""
    Custom instructions for the agent...
    """
)
```

---

## 🚀 Summary

**The workflow is**:

1. **CLI** receives user input
2. **ChemToolsAgent** wraps it with context/constraints
3. **LangGraph** creates ReAct loop
4. **LLM** decides which tools to call
5. **Tools** execute and return results
6. **LLM** interprets results and decides next action
7. Loop continues until final answer
8. **Response** displayed to user

**Key innovation**: The LLM acts as the "brain" that decides which ChemTools functions to call and in what order, based on:

- System prompt guidance
- Tool descriptions
- User query
- Conversation history
- Tool results

This creates a flexible, intelligent agent that can handle complex chemistry workflows without hardcoded logic!
