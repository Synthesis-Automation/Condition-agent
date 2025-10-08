# MCP Integration Summary

**Date:** October 6, 2025  
**Status:** Legacy MCP Server Removed, Tool Wrappers Remain

## Overview

The `chemtools/integrations/mcp/` folder contains **lightweight tool wrappers** for exposing chemtools functionality through the Model Context Protocol (MCP) interface. This is **NOT an MCP server implementation**.

## What's Inside

### Current Structure

```
chemtools/integrations/
├── mcp/
│   ├── __init__.py              # Version and exports
│   ├── tools/                   # Tool wrappers (current)
│   │   ├── __init__.py
│   │   ├── base.py              # Common helpers (SchemaStamped, validation)
│   │   ├── detect.py            # Reaction family detection wrapper
│   │   ├── featurize.py         # Substrate featurization wrapper
│   │   └── normalize.py         # SMILES normalization wrapper
│   └── resources/
│       └── schemas/
│           └── condition_set.json  # JSON schema for condition recommendations
└── molpipeline.py               # Optional MolPipeline integration
```

### What These Are

**Tool Wrappers (NOT a server):**
- Thin, deterministic functions that wrap chemtools functionality
- Return JSON-serializable payloads with `schema_version` field
- Designed to be consumed by MCP server implementations
- Pydantic-based validation for inputs/outputs

**Purpose:**
- Provide standardized interfaces for external MCP servers to call
- Schema-stamped outputs for consistent data exchange
- Can be used standalone or integrated into an MCP server

## The 3 Tool Wrappers

### 1. `normalize_reaction` (normalize.py)

**Purpose:** Canonicalize reaction SMILES strings

**Input:**
```python
{
    "smiles_rxn": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "include_agents": true
}
```

**Output:**
```python
{
    "schema_version": "0.1.0",
    "input": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "normalized": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "reactants": ["Brc1ccccc1", "Nc1ccccc1"],
    "products": ["c1ccc(Nc2ccccc2)cc1"],
    "agents": [],
    "errors": []
}
```

**Backend:** Uses `chemtools.smiles.normalize_reaction()`

### 2. `detect_family` (detect.py)

**Purpose:** Detect reaction family from reactant SMILES

**Input:**
```python
{
    "reactants": ["Brc1ccccc1", "Nc1ccccc1"]
}
```

**Output:**
```python
{
    "schema_version": "0.1.0",
    "family": "C_N_Coupling_Pd",
    "confidence": 0.8,
    "hits": {
        "has_aryl_halide": true,
        "has_amine": true
    },
    "has_conflict": false
}
```

**Backend:** Uses `chemtools.router.detect_family()`

### 3. `featurize_substrates` (featurize.py)

**Purpose:** Extract molecular descriptors from substrates

**Input:**
```python
{
    "reactants": ["Brc1ccccc1", "Nc1ccccc1"]
}
```

**Output:**
```python
{
    "schema_version": "0.1.0",
    "electrophile": "Brc1ccccc1",
    "nucleophile": "Nc1ccccc1",
    "descriptors": {
        "electrophile_features": {...},
        "nucleophile_features": {...}
    }
}
```

**Backend:** Uses `chemtools.featurizers.molecular.featurize()`

## Common Infrastructure

### Base Module (`base.py`)

**Key Components:**

1. **`SchemaStamped`** - Pydantic base model that auto-injects `schema_version`
2. **`ToolError`** - Custom exception for validation failures
3. **`validate_payload()`** - Pydantic validation wrapper
4. **`pick_first()`** - Helper to select first non-empty string
5. **`flatten_smiles_block()`** - Extract SMILES from normalized blocks

**Example:**
```python
from chemtools.integrations.mcp.tools.base import SchemaStamped

class MyOutput(SchemaStamped):
    result: str
    value: float

output = MyOutput(result="success", value=0.95)
# Automatically includes: {"schema_version": "0.1.0", "result": "success", "value": 0.95}
```

### Schema Versioning

- Default version: `0.1.0` (from `chemtools.integrations.mcp.DEFAULT_SCHEMA_VERSION`)
- All outputs include `schema_version` field
- Enables backward compatibility tracking

## JSON Schema Resources

### `condition_set.json`

Defines the structure for condition recommendations returned by an MCP-compatible system:

**Key Sections:**
- **`core`**: Catalytic core (metal/ligand pair)
- **`components`**: Reagents with roles, identifiers, CAS numbers
- **`conditions`**: Temperature, time, pressure, atmosphere
- **`metrics`**: Confidence scores, support counts
- **`constraints`**: Constraint violations
- **`rationale`**: Explanations and citations
- **`inputs`**: Echo of input parameters

This schema is referenced at:
```
https://synthesis-automation.github.io/condition-mcp/schema/condition_set.json
```

## Legacy MCP Server (REMOVED)

### What Was Removed

The repository **previously** had an MCP server implementation that has been **deprecated and removed**:

**Former structure (no longer exists):**
```
chemtools/integrations/mcp/
├── server/
│   └── server.py        # ❌ REMOVED - Old MCP server
└── client.py            # ❌ REMOVED - Client adapter
```

**Evidence of removal:**

1. **`chemtools/agent/config.py`** - Contains deprecation warning:
   ```python
   warnings.warn(
       "chemtools.agent.config.load_config is deprecated because the rule-based "
       "system has been removed; a placeholder configuration is returned.",
       DeprecationWarning,
   )
   ```

2. **`docs/old/rules_mcp.md`** - Moved to `old/` directory, references removed code

3. **No server implementation files** - `file_search` found no `server.py`

### Why It Was Removed

The legacy MCP server was replaced by **`chemtools.scdb_matcher`** (the rule-based SchemeConditionDB matcher):

**Old approach (MCP server):**
- External subprocess server
- JSON-RPC communication
- Complex client/server setup
- Harder to maintain

**New approach (scdb_matcher):**
- Direct Python module in `chemtools/scdb_matcher/`
- In-process matching
- Deterministic and testable
- Simpler integration

## Current Usage

### How Tool Wrappers Are Used

These tool wrappers are **consumed by the application** but **not deployed as an MCP server**:

**In tests:**
```python
# tests/chemtools/integrations/test_mcp_tools.py
from chemtools.integrations.mcp.tools import detect_family, featurize_substrates, normalize_reaction

# Test tool outputs
result = normalize_reaction({"smiles_rxn": "reactants>>products"})
assert result["schema_version"] == "0.1.0"
```

**Potential external use:**
- Could be wrapped in an actual MCP server implementation
- Could be called by external systems expecting MCP-style interfaces
- Provides standardized data exchange format

### Not Used For

❌ **NOT a running MCP server** - No server process, no JSON-RPC endpoints  
❌ **NOT client/server architecture** - Direct Python function calls  
❌ **NOT deployed independently** - Part of chemtools package  

## Relationship to Main Systems

### ML Recommendations
```
chemtools/recommend.py (ML engine)
    ↓
app/ui_simple.py (UI layer)
    ↓
User interface
```

**MCP tools NOT involved** - ML system uses direct imports

### Rule-Based Recommendations
```
chemtools/scdb_matcher/ (Rule engine)
    ↓
app/ui_simple.py (UI layer)
    ↓
User interface
```

**MCP tools NOT involved** - Rule system uses direct imports

### Tool Wrappers
```
chemtools/integrations/mcp/tools/
    ↓
[EXTERNAL MCP SERVER COULD WRAP THESE]
    ↓
External clients (if implemented)
```

**Currently unused** - Present for potential future MCP server integration

## Key Differences

| Feature | Legacy MCP Server (Removed) | Tool Wrappers (Current) | scdb_matcher (Current) |
|---------|---------------------------|------------------------|----------------------|
| **Type** | External server process | Python functions | Python module |
| **Communication** | JSON-RPC over stdio | Direct function calls | Direct imports |
| **Status** | ❌ Removed | ✅ Present (unused) | ✅ Active |
| **Purpose** | Rule-based matching | MCP interface layer | Rule-based matching |
| **Integration** | Subprocess client | Import and call | Import and call |
| **Maintenance** | ❌ Deprecated | ⚠️ Maintained (dormant) | ✅ Active |

## Summary

### What `chemtools/integrations/mcp/` IS:

✅ **Tool wrapper library** - Thin, schema-stamped functions  
✅ **Data standardization layer** - Consistent JSON outputs  
✅ **Pydantic-validated interfaces** - Type-safe inputs/outputs  
✅ **Potential MCP building blocks** - Ready for server implementation  

### What it IS NOT:

❌ **Not an MCP server** - No running server process  
❌ **Not actively used** - Not called by current UI/API  
❌ **Not the rule-based system** - That's `scdb_matcher/`  
❌ **Not the ML system** - That's `recommend.py`  

### Current State:

- **Legacy server code:** Removed (see `docs/old/rules_mcp.md`)
- **Tool wrappers:** Present but not actively used
- **Schema definitions:** Maintained for potential future use
- **Active recommendation systems:** `chemtools/recommend.py` (ML) and `chemtools/scdb_matcher/` (rules)

## Recommendations

### For External MCP Server Implementation

If you want to create an actual MCP server using these tools:

1. **Install MCP SDK:**
   ```bash
   pip install mcp
   ```

2. **Create server wrapper:**
   ```python
   # hypothetical_mcp_server.py
   from mcp.server import Server
   from chemtools.integrations.mcp.tools import normalize_reaction, detect_family, featurize_substrates
   
   app = Server("chemtools-mcp")
   
   @app.tool()
   def normalize(smiles_rxn: str) -> dict:
       return normalize_reaction({"smiles_rxn": smiles_rxn})
   
   @app.tool()
   def detect(reactants: list[str]) -> dict:
       return detect_family({"reactants": reactants})
   
   if __name__ == "__main__":
       app.run()
   ```

3. **Deploy as MCP server** following MCP protocol specifications

### For Current Codebase

- **Keep tool wrappers** - Useful abstraction layer, low maintenance cost
- **Keep schemas** - Documents data structures for recommendations
- **Remove if unused long-term** - If MCP integration never materializes
- **Focus on active systems** - `scdb_matcher` and `recommend.py` are production code

---

**Conclusion:** The `chemtools/integrations/mcp/` folder contains **tool wrappers** (not a server), the legacy MCP server was removed in favor of the `scdb_matcher` module, and these wrappers are currently dormant infrastructure that could enable future MCP server integration if needed.
