# MCP Server Implementation Guide

**Date:** October 6, 2025  
**Purpose:** Guide for creating an MCP server to expose chemtools functionality to robotic systems

## Executive Summary

✅ **YES, it's relatively easy to create an MCP server** based on your current codebase!

Your codebase is **well-suited** for MCP server implementation because:
- ✅ Clean separation: Core logic in `chemtools/`, APIs in `app/`
- ✅ Already has tool wrappers in `chemtools/integrations/mcp/tools/`
- ✅ Deterministic functions with clear inputs/outputs
- ✅ Pydantic-based contracts for validation
- ✅ Existing FastAPI server as reference architecture

## Use Case: Robotic System Integration

**Your stated need:**
> Accept input (reaction SMILES) → Output recommended conditions → Used by robotic system for execution

**Perfect fit for MCP because:**
- 🤖 **Robotic systems** can call MCP tools directly
- 🔄 **Standardized protocol** - Works with any MCP client
- 📊 **Structured outputs** - JSON with schema validation
- 🔧 **Multiple tools** - Normalize, detect, recommend, etc.

## Current Assets You Already Have

### 1. Core Recommendation Functions ✅

**ML-based recommendations:**
```python
# chemtools/recommend.py
def recommend_conditions_structured(
    reaction: str,
    reaction_type: Optional[str] = None,
    *,
    k: int = 50,
    limit: int = 5,
    relax: Optional[Dict[str, Any]] = None,
    constraints: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Structured condition recommendations for API/UI consumers."""
```

**Rule-based recommendations:**
```python
# chemtools/scdb_matcher/matcher.py
def match(
    db: RuleDB,
    rxn_smiles: str | None = None,
    *,
    features: Mapping[str, Any] | None = None,
) -> MatchResult:
    """Match a reaction to database rules."""
```

### 2. Utility Tools ✅

**Already wrapped in MCP-style interfaces:**

```python
# chemtools/integrations/mcp/tools/
- normalize_reaction()  # SMILES canonicalization
- detect_family()       # Reaction type detection
- featurize_substrates() # Molecular descriptors
```

### 3. FastAPI Reference Implementation ✅

Your `app/main.py` shows how these functions are exposed via HTTP:
- Request/response contracts (Pydantic models)
- Error handling patterns
- Input validation

### 4. Output Formatter ✅

```python
# chemtools/output_formatter.py
def format_ml_output(...)  # Structured ML recommendations
def format_rule_output(...) # Structured rule recommendations
```

## MCP Server Architecture

### Recommended Structure

```
chemtools-mcp-server/
├── server.py                    # MCP server entry point
├── tools/
│   ├── __init__.py
│   ├── normalize.py             # Use existing wrapper
│   ├── detect.py                # Use existing wrapper
│   ├── recommend_ml.py          # NEW: Wrap recommend_conditions_structured
│   ├── recommend_rules.py       # NEW: Wrap scdb_matcher
│   └── format_conditions.py     # NEW: Output formatting for robots
├── schemas/
│   ├── condition_set.json       # Already exists!
│   └── robot_command.json       # NEW: Robot-specific format
├── requirements.txt
└── README.md
```

## Implementation Steps

### Step 1: Install MCP SDK

```bash
pip install mcp
```

### Step 2: Create MCP Server

**Create `server.py`:**

```python
"""MCP Server for ChemTools - Robotic System Interface"""

import asyncio
from typing import Any
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

# Import existing chemtools functionality
from chemtools.integrations.mcp.tools import (
    normalize_reaction,
    detect_family,
    featurize_substrates,
)
from chemtools import recommend
from chemtools.scdb_matcher import load_db, match
from chemtools.output_formatter import format_ml_output, format_rule_output

# Initialize server
app = Server("chemtools-robot-server")

# Database caching (load once)
_RULE_DATABASES = {}

def _get_rule_db(reaction_type: str):
    """Load and cache rule databases."""
    if reaction_type not in _RULE_DATABASES:
        db_paths = {
            "C_N_Coupling_Pd": "data/conditionDB/cn_coupling_pd_db.json",
            "C_N_Coupling_Cu": "data/conditionDB/cn_coupling_cu_db.json",
            "C_N_Coupling_Ni": "data/conditionDB/cn_coupling_ni_db.json",
            "Suzuki_CC": "data/conditionDB/suzuki_db.json",
        }
        if reaction_type in db_paths:
            _RULE_DATABASES[reaction_type] = load_db(db_paths[reaction_type])
    return _RULE_DATABASES.get(reaction_type)


# Tool 1: Normalize Reaction SMILES
@app.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="normalize_reaction_smiles",
            description="Canonicalize and validate reaction SMILES for robotic execution",
            inputSchema={
                "type": "object",
                "properties": {
                    "smiles_rxn": {
                        "type": "string",
                        "description": "Reaction SMILES (reactants>agents>products or reactants>>products)"
                    }
                },
                "required": ["smiles_rxn"]
            }
        ),
        Tool(
            name="detect_reaction_type",
            description="Automatically detect reaction type from reactants",
            inputSchema={
                "type": "object",
                "properties": {
                    "reactants": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of reactant SMILES"
                    }
                },
                "required": ["reactants"]
            }
        ),
        Tool(
            name="recommend_conditions_ml",
            description="Get ML-based condition recommendations for a reaction (for robotic execution)",
            inputSchema={
                "type": "object",
                "properties": {
                    "reaction_smiles": {
                        "type": "string",
                        "description": "Reaction SMILES string"
                    },
                    "reaction_type": {
                        "type": "string",
                        "description": "Reaction type (C_N_Coupling_Pd, Suzuki_CC, etc.) or null for auto-detect"
                    },
                    "top_k": {
                        "type": "integer",
                        "description": "Number of recommendations to return",
                        "default": 3
                    }
                },
                "required": ["reaction_smiles"]
            }
        ),
        Tool(
            name="recommend_conditions_rules",
            description="Get rule-based condition recommendations for a reaction (for robotic execution)",
            inputSchema={
                "type": "object",
                "properties": {
                    "reaction_smiles": {
                        "type": "string",
                        "description": "Reaction SMILES string"
                    },
                    "reaction_type": {
                        "type": "string",
                        "description": "Reaction type (must match available rule databases)"
                    }
                },
                "required": ["reaction_smiles", "reaction_type"]
            }
        ),
        Tool(
            name="format_for_robot",
            description="Format condition recommendations for robotic system execution",
            inputSchema={
                "type": "object",
                "properties": {
                    "recommendations": {
                        "type": "object",
                        "description": "Raw recommendations from ML or rule system"
                    },
                    "robot_format": {
                        "type": "string",
                        "enum": ["json", "opentrons", "tecan"],
                        "description": "Target robot platform format",
                        "default": "json"
                    }
                },
                "required": ["recommendations"]
            }
        ),
    ]


@app.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    """Handle tool calls from MCP clients."""
    
    if name == "normalize_reaction_smiles":
        result = normalize_reaction({
            "smiles_rxn": arguments["smiles_rxn"],
            "include_agents": True
        })
        return [TextContent(
            type="text",
            text=f"Normalized reaction:\n{result['normalized']}\n\n" +
                 f"Reactants: {', '.join(result['reactants'])}\n" +
                 f"Products: {', '.join(result['products'])}\n" +
                 f"Agents: {', '.join(result.get('agents', []))}"
        )]
    
    elif name == "detect_reaction_type":
        result = detect_family({
            "reactants": arguments["reactants"]
        })
        return [TextContent(
            type="text",
            text=f"Detected reaction type: {result['family']}\n" +
                 f"Confidence: {result['confidence']:.2%}\n" +
                 f"Analysis: {result['hits']}"
        )]
    
    elif name == "recommend_conditions_ml":
        # Call the core ML recommendation function
        reaction = arguments["reaction_smiles"]
        reaction_type = arguments.get("reaction_type")
        top_k = arguments.get("top_k", 3)
        
        try:
            recommendations = recommend.recommend_conditions_structured(
                reaction=reaction,
                reaction_type=reaction_type,
                limit=top_k
            )
            
            # Format for robotic consumption
            formatted = format_ml_output(
                reaction_smiles=reaction,
                requested_type=reaction_type or "Auto-detect",
                detected_type=recommendations.get("detected_family", "Unknown"),
                recommendations_data=recommendations.get("recommendations", []),
                processing_time_ms=recommendations.get("processing_time_ms", 0)
            )
            
            return [TextContent(
                type="text",
                text=f"ML Recommendations:\n\n{formatted}"
            )]
        except Exception as e:
            return [TextContent(
                type="text",
                text=f"Error generating ML recommendations: {str(e)}"
            )]
    
    elif name == "recommend_conditions_rules":
        # Call the rule-based matcher
        reaction = arguments["reaction_smiles"]
        reaction_type = arguments["reaction_type"]
        
        try:
            db = _get_rule_db(reaction_type)
            if db is None:
                return [TextContent(
                    type="text",
                    text=f"Error: No rule database available for reaction type '{reaction_type}'"
                )]
            
            result = match(db, rxn_smiles=reaction)
            
            # Format for robotic consumption
            formatted = format_rule_output(
                reaction_smiles=reaction,
                requested_type=reaction_type,
                detected_type=result.reaction_type,
                recommendations_data=[{
                    "conditions": result.conditions,
                    "entry_id": result.entry_id,
                    "match_type": result.match_type
                }],
                database_name=reaction_type,
                processing_time_ms=0
            )
            
            return [TextContent(
                type="text",
                text=f"Rule-based Recommendations:\n\n{formatted}"
            )]
        except Exception as e:
            return [TextContent(
                type="text",
                text=f"Error generating rule-based recommendations: {str(e)}"
            )]
    
    elif name == "format_for_robot":
        # Convert recommendations to robot-specific format
        recommendations = arguments["recommendations"]
        robot_format = arguments.get("robot_format", "json")
        
        if robot_format == "opentrons":
            # Format for Opentrons robot
            robot_output = _format_for_opentrons(recommendations)
        elif robot_format == "tecan":
            # Format for Tecan robot
            robot_output = _format_for_tecan(recommendations)
        else:
            # Default JSON format
            robot_output = recommendations
        
        return [TextContent(
            type="text",
            text=f"Robot-formatted output:\n{robot_output}"
        )]
    
    else:
        raise ValueError(f"Unknown tool: {name}")


def _format_for_opentrons(recommendations: dict) -> dict:
    """Format recommendations for Opentrons liquid handling robot."""
    # Example: Convert to Opentrons protocol format
    return {
        "protocol_version": "1.0",
        "metadata": {
            "reaction": recommendations.get("input", {}).get("reaction_smiles"),
            "reaction_type": recommendations.get("detection", {}).get("detected_reaction_type")
        },
        "commands": [
            # Extract reagents and volumes
            {
                "action": "transfer",
                "source": rec.get("reagents", {}).get("solvent", {}).get("name"),
                "volume_ul": rec.get("reagents", {}).get("solvent", {}).get("volume_ml", 1.0) * 1000,
                "destination": "reaction_vial"
            }
            for rec in recommendations.get("recommendations", [])[:1]  # First recommendation
        ],
        "conditions": recommendations.get("recommendations", [{}])[0].get("conditions", {})
    }


def _format_for_tecan(recommendations: dict) -> dict:
    """Format recommendations for Tecan liquid handling robot."""
    # Example: Convert to Tecan worklist format
    return {
        "worklist_version": "1.0",
        "steps": [
            {
                "step": i + 1,
                "reagent": reagent.get("name"),
                "cas": reagent.get("cas"),
                "volume_ul": reagent.get("volume", 100)
            }
            for i, rec in enumerate(recommendations.get("recommendations", [])[:1])
            for reagent in rec.get("reagents", {}).values()
            if isinstance(reagent, dict)
        ]
    }


async def main():
    """Run the MCP server."""
    async with stdio_server() as (read_stream, write_stream):
        await app.run(
            read_stream,
            write_stream,
            app.create_initialization_options()
        )


if __name__ == "__main__":
    asyncio.run(main())
```

### Step 3: Create Package Structure

**Create `pyproject.toml`:**

```toml
[project]
name = "chemtools-mcp-server"
version = "0.1.0"
description = "MCP Server for ChemTools - Robotic System Interface"
dependencies = [
    "mcp>=0.9.0",
    "chemtools>=0.1.0",  # Your package
]

[project.scripts]
chemtools-mcp-server = "chemtools_mcp_server.server:main"
```

### Step 4: Test Locally

```bash
# Install in development mode
pip install -e .

# Run the server
python server.py

# Test with MCP inspector (if available)
mcp-inspector chemtools-mcp-server
```

### Step 5: Robot Integration

**Example: Python robot client using MCP**

```python
"""Robot client using ChemTools MCP server"""

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def run_reaction_workflow():
    """Example workflow for robotic system."""
    
    # Connect to MCP server
    server_params = StdioServerParameters(
        command="python",
        args=["server.py"]
    )
    
    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            
            # Step 1: Normalize reaction
            reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
            
            normalize_result = await session.call_tool(
                "normalize_reaction_smiles",
                {"smiles_rxn": reaction}
            )
            print("Normalized:", normalize_result)
            
            # Step 2: Get recommendations
            recommendations = await session.call_tool(
                "recommend_conditions_ml",
                {
                    "reaction_smiles": reaction,
                    "reaction_type": "C_N_Coupling_Pd",
                    "top_k": 3
                }
            )
            print("Recommendations:", recommendations)
            
            # Step 3: Format for robot
            robot_commands = await session.call_tool(
                "format_for_robot",
                {
                    "recommendations": recommendations,
                    "robot_format": "opentrons"
                }
            )
            print("Robot commands:", robot_commands)
            
            # Step 4: Execute on robot (robot-specific code)
            # execute_on_robot(robot_commands)


if __name__ == "__main__":
    import asyncio
    asyncio.run(run_reaction_workflow())
```

## Advantages of MCP for Robotic Systems

### 1. **Standardized Protocol** ✅
- Industry-standard way for AI/LLM systems to call tools
- Works with any MCP-compatible client (Claude, custom robots, etc.)
- No need to build custom HTTP APIs

### 2. **Type Safety** ✅
- JSON Schema validation for inputs
- Structured outputs
- Reduces robot programming errors

### 3. **Tool Discovery** ✅
- Robots can query available tools: `list_tools()`
- Self-documenting: Each tool has description and schema
- Easy to add new tools without changing robot code

### 4. **Async/Concurrent** ✅
- MCP is async-native (Python asyncio)
- Can handle multiple robot requests simultaneously
- Good for multi-robot labs

### 5. **Language Agnostic** ✅
- MCP protocol is language-independent
- Robot software in Python, TypeScript, Java, etc. can all connect
- JSON-based communication

## Alternative: Extend Existing FastAPI Server

If MCP overhead is not desired, you can also:

**Option A: Add robot-specific endpoints to `app/main.py`:**

```python
@app.post("/robot/recommend")
async def robot_recommend(
    reaction_smiles: str,
    reaction_type: Optional[str] = None,
    robot_format: str = "opentrons"
):
    """Endpoint specifically for robotic system consumption."""
    
    # Get recommendations
    recommendations = recommend.recommend_conditions_structured(
        reaction=reaction_smiles,
        reaction_type=reaction_type,
        limit=3
    )
    
    # Format for robot
    if robot_format == "opentrons":
        return _format_for_opentrons(recommendations)
    elif robot_format == "tecan":
        return _format_for_tecan(recommendations)
    else:
        return recommendations
```

**Option B: Create dedicated robot API service:**

```python
# app/robot_api.py
from fastapi import FastAPI

robot_app = FastAPI(title="ChemTools Robot API")

@robot_app.post("/execute-reaction")
async def execute_reaction(request: RobotExecutionRequest):
    """Single endpoint: SMILES in → Robot commands out."""
    # Normalize → Detect → Recommend → Format
    ...
```

## Comparison: MCP vs FastAPI vs Direct Integration

| Feature | MCP Server | FastAPI Extension | Direct Python Import |
|---------|-----------|-------------------|---------------------|
| **Protocol** | MCP (stdio/SSE) | HTTP REST | In-process |
| **Discovery** | Built-in (`list_tools()`) | Manual docs | Code inspection |
| **Type Safety** | JSON Schema | Pydantic | Type hints |
| **Language Support** | Any MCP client | HTTP clients | Python only |
| **LLM Integration** | ✅ Native (Claude, etc.) | ⚠️ Manual | ❌ None |
| **Robot Integration** | ✅ Easy | ✅ Easy | ⚠️ Requires Python |
| **Deployment** | Stdio or SSE server | HTTP server | Library import |
| **Complexity** | Medium | Low | Very Low |
| **Best For** | AI/LLM + Robots | Traditional APIs | Python robots |

## Recommendation for Your Use Case

### If your robotic system is:

**1. Python-based and local to codebase:**
→ **Use direct imports** (simplest)
```python
from chemtools import recommend
result = recommend.recommend_conditions_structured(reaction)
```

**2. HTTP-capable (most commercial robots):**
→ **Extend FastAPI** (`app/main.py`)
- Already implemented
- Add robot-specific endpoints
- RESTful and well-understood

**3. LLM-assisted or needs AI integration:**
→ **Build MCP server** (as shown above)
- Claude/GPT can call tools
- Standardized protocol
- Future-proof for AI lab automation

**4. Multi-platform (Python + non-Python robots):**
→ **MCP server OR FastAPI**
- Both provide language-agnostic interfaces
- MCP better for tool discovery
- FastAPI better for traditional REST workflows

## Estimated Implementation Time

| Task | Complexity | Time |
|------|-----------|------|
| **Basic MCP server** (3-5 tools) | Low | 4-8 hours |
| **Add robot formatters** (Opentrons, Tecan) | Medium | 8-16 hours |
| **Testing & integration** | Medium | 8-16 hours |
| **Documentation** | Low | 4-8 hours |
| **Total** | - | **1-3 days** |

## Next Steps

### Immediate (Quick Win):
1. ✅ **Test existing MCP wrappers** - Run `chemtools/integrations/mcp/tools/*` functions
2. ✅ **Create simple MCP server** - Use code template above
3. ✅ **Test with MCP inspector** - Verify tools work

### Short-term (Production):
1. 🔨 **Add ML/Rule recommendation tools** - Wrap existing functions
2. 🔨 **Implement robot formatters** - Opentrons, Tecan, etc.
3. 🔨 **Error handling** - Robust error messages for robots
4. 🔨 **Logging** - Track robot requests for debugging

### Long-term (Scale):
1. 🚀 **Multi-robot support** - Handle concurrent requests
2. 🚀 **Caching** - Speed up repeated calls
3. 🚀 **Monitoring** - Prometheus metrics
4. 🚀 **Authentication** - Secure robot access

## Conclusion

**YES, creating an MCP server is straightforward** because:
- ✅ Your code is already well-structured
- ✅ Core functions are deterministic
- ✅ Tool wrappers already exist
- ✅ Output formatters in place
- ✅ FastAPI server as reference

**For robotic system integration, you have 3 good options:**
1. **Direct Python imports** (simplest, Python-only robots)
2. **FastAPI extension** (traditional REST API, well-understood)
3. **MCP server** (modern, AI-compatible, standardized)

**Recommended path:**
Start with **FastAPI extension** (quickest, most compatible), then add **MCP server** if you need LLM integration or standardized tool discovery.

Would you like me to create a working MCP server implementation or help extend the FastAPI server with robot-specific endpoints?
