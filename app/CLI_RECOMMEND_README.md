# Interactive CLI for Reaction Condition Recommendations

## Overview

A natural language interface for requesting chemistry condition recommendations. Uses LLM-powered parsing to convert user requirements into structured API calls, with progressive refinement until validation passes.

**Key Feature**: Only a valid reaction SMILES is required to complete - all constraints are optional.

## Features

- ✅ **Natural Language Input**: Describe requirements in plain English
- ✅ **LLM-Powered Parsing**: Automatically extracts constraints from descriptions
- ✅ **Progressive Refinement**: Asks targeted questions until payload validates
- ✅ **Minimal Requirements**: Only valid reaction SMILES needed - constraints optional
- ✅ **Structured JSON Schema**: Ensures API compatibility
- ✅ **Validation Loop**: Iterates until SMILES format is correct
- ✅ **Direct API Integration**: Calls actual recommendation API when available
- ✅ **Draft Saving**: Save incomplete requests for later

## Installation

No additional dependencies required - uses existing `llmtools` package.

## Usage

### Basic Usage

```bash
# Run with default settings (Aliyun/DeepSeek)
python app/cli_recommend.py

# Use OpenAI instead
python app/cli_recommend.py --provider openai --model gpt-4o

# Enable debug logging
python app/cli_recommend.py --debug
```

### Environment Variables

Set your API key before running:

```bash
# For Aliyun (DashScope)
export DASHSCOPE_API_KEY=your_api_key_here

# For OpenAI
export OPENAI_API_KEY=your_api_key_here
```

## Example Sessions

### Session 1: Valid SMILES, No Constraints

```
======================================================================
  🧪 CHEMISTRY CONDITION RECOMMENDATION CLI
  Powered by LLM-assisted natural language parsing
======================================================================

Please provide your reaction details:

Reaction SMILES (reactants>>products format):
Example: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
> Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1

📝 Describe your requirements in natural language (optional):
Examples:
  • no strong base, room temperature
  • avoid expensive catalysts, prefer DMF solvent
  • mild conditions, no air-sensitive reagents
  • (press Enter to skip if no specific requirements)

> 

🤖 Parsing input with LLM...

----------------------------------------------------------------------

📋 PARSED REQUEST:

✅ Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
   Type: Suzuki_Miyaura

✅ Request is VALID and ready to submit!

----------------------------------------------------------------------
✅ Request is ready for submission!

Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Type: Suzuki_Miyaura
Constraints: None (will use default recommendations)

🔍 API Request Preview:
{
  "reaction": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
  "reaction_type": "Suzuki_Miyaura",
  "k": 10,
  "limit": 5,
  "constraints": {},
  "relax": {}
}

======================================================================
Submit this request? (yes/no): yes

======================================================================
  📤 SUBMITTING REQUEST...
======================================================================

⏳ Calling API endpoint: POST /api/v1/recommend/conditions

✅ REQUEST SUCCESSFUL!

======================================================================
RESULTS:
======================================================================
{
  "recommendations": [...],
  "processing_time_ms": 234.5
}

✨ Recommendation complete!
```

### Session 2: Invalid SMILES with Refinement

```
======================================================================
  🧪 CHEMISTRY CONDITION RECOMMENDATION CLI
  Powered by LLM-assisted natural language parsing
======================================================================

Please provide your reaction details:

Reaction SMILES (reactants>>products format):
Example: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
> invalid_smiles

📝 Describe your requirements in natural language (optional):
> no strong base

🤖 Parsing input with LLM...

----------------------------------------------------------------------

📋 PARSED REQUEST:

❌ Reaction: invalid_smiles

📌 Constraints:
  • base_strength: moderate

⚠️  VALIDATION ISSUES:
  • Invalid reaction SMILES format

----------------------------------------------------------------------

❌ VALIDATION ISSUES (must fix):

1. Invalid reaction SMILES format

Please fix the validation issues above:
> Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1

🤖 Refining with your input...

----------------------------------------------------------------------

📋 PARSED REQUEST:

✅ Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
   Type: Suzuki_Miyaura

📌 Constraints:
  • base_strength: moderate

✅ Request is VALID and ready to submit!

----------------------------------------------------------------------
✅ Request is ready for submission!
...
```

### Session 3: With Constraints

```
Please provide your reaction details:

Reaction SMILES (reactants>>products format):
> Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1

📝 Describe your requirements in natural language (optional):
> no strong base, room temperature, avoid expensive palladium catalysts

🤖 Parsing input with LLM...

----------------------------------------------------------------------

📋 PARSED REQUEST:

✅ Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
   Type: Buchwald_Hartwig

📌 Constraints:
  • temperature:
      max: 30
  • base_strength: moderate
  • cost_level: low

⚠️  CLARIFICATIONS NEEDED (optional):

1. You mentioned "room temperature" - should this be strict (max 25°C) or flexible (up to 40°C)?
2. "Avoid expensive palladium catalysts" - are common Pd sources like Pd(OAc)2 acceptable?

Please provide clarifications (or press Enter to skip):
> Flexible up to 40C. Yes common Pd is fine.

🤖 Refining with your input...

----------------------------------------------------------------------

📋 PARSED REQUEST:

✅ Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
   Type: Buchwald_Hartwig

📌 Constraints:
  • temperature:
      max: 40
  • base_strength: moderate
  • cost_level: medium

✅ Request is VALID and ready to submit!
```

## JSON Schema

The CLI uses structured outputs to ensure API compatibility:

```json
{
  "type": "object",
  "properties": {
    "reaction_smiles": {
      "type": "string",
      "description": "Valid SMILES (reactants>>products format)"
    },
    "reaction_smiles_is_valid": {
      "type": "boolean"
    },
    "reaction_type": {
      "type": ["string", "null"],
      "description": "e.g., Suzuki, Buchwald, Ullmann"
    },
    "constraints": {
      "type": "object",
      "properties": {
        "temperature": {
          "type": "object",
          "properties": {
            "max": {"type": "number"},
            "min": {"type": "number"}
          }
        },
        "base_strength": {
          "enum": ["weak", "moderate", "strong", "any", null]
        },
        "required_reagents": {
          "type": "array",
          "items": {"type": "string"},
          "description": "Reagents that must be included"
        },
        "metal_preference": {
          "type": ["string", "null"],
          "enum": ["Pd", "Cu", "Ni", "Fe", "Ru", "Rh", "Ir", "Au", "Ag", "any", null],
          "description": "Preferred metal for catalysis"
        },
        "exclude_reagents": {
          "type": "array",
          "items": {"type": "string"}
        },
        "solvent_preference": {
          "enum": ["polar_aprotic", "polar_protic", "nonpolar", "aqueous", "any", null]
        },
        "air_sensitive": {"type": "boolean"},
        "cost_level": {
          "enum": ["low", "medium", "high", "any", null]
        }
      }
    },
    "validation_issues": {
      "type": "array",
      "items": {"type": "string"}
    },
    "clarification_needed": {
      "type": "array",
      "items": {"type": "string"}
    }
  },
  "required": ["reaction_smiles", "reaction_smiles_is_valid", "validation_issues"]
}
```

## Constraint Vocabulary

The LLM understands common chemical terms and maps them to structured constraints:

### Temperature
- "room temperature" → max: 30°C
- "no high temperature" → max: 80°C
- "reflux" → context-dependent
- "0°C" / "ice bath" → max: 0°C

### Base Strength
- "no strong base" → "moderate"
- "mild base" → "weak"
- "strong base acceptable" → "strong"

### Metal Preference
- "use copper catalyst", "copper" → "Cu"
- "palladium catalyst", "Pd catalyst" → "Pd"
- "nickel", "Ni catalyst" → "Ni"
- "iron", "Fe catalyst" → "Fe"

### Required Reagents (Positive Constraints)
- "use X", "with X", "prefer X" → adds X to required_reagents
- "copper catalyst" → adds "copper" to required list
- "base is required" → adds "base" to required list

### Excluded Reagents (Negative Constraints)
- "no X", "avoid X", "exclude X" → adds X to exclude_reagents
- "no Pd(PPh3)4" → exclude specific compound
- "avoid expensive catalysts" → combines with cost_level

### Solvent
- "DMF", "DMSO", "acetonitrile" → "polar_aprotic"
- "methanol", "ethanol", "water" → "polar_protic"
- "toluene", "hexane" → "nonpolar"

### Cost
- "cheap", "cost-effective" → "low"
- "expensive acceptable" → "high"
- "moderate cost" → "medium"

### Air Sensitivity
- "inert atmosphere", "glovebox", "Schlenk" → air_sensitive: true
- "air stable", "bench-top" → air_sensitive: false

## Validation Logic

### Required for Submission
1. ✅ Valid reaction SMILES format (must have ">>" and non-empty reactants)

### Optional
- All constraints are optional
- Empty constraints object is perfectly acceptable
- LLM will only ask for clarification on ambiguous constraints

### Validation Loop
1. Parse initial input with LLM
2. Check if SMILES is valid
3. If invalid → ask user to fix
4. If valid but constraints ambiguous → optionally clarify
5. Repeat until `validation_issues` array is empty

## API Integration

When run from the main app context, the CLI will automatically call the actual API:

```python
from chemtools.contracts import RecommendConditionsRequest
from app.main import api_recommend_conditions

req = RecommendConditionsRequest(**api_request)
result = api_recommend_conditions(req)
```

## Error Handling

- **Invalid SMILES**: LLM detects format issues, asks for correction
- **LLM parsing failure**: Graceful fallback with error message
- **API unavailable**: Shows request payload that would be sent
- **KeyboardInterrupt**: Saves draft and exits cleanly

## Draft Saving

If the session is interrupted or cancelled, the current state is saved:

```bash
💾 Draft saved to: draft_request.json
```

Load it later or use for debugging.

## Command-Line Options

```
python app/cli_recommend.py [OPTIONS]

Options:
  --provider {openai,aliyun}    LLM provider (default: aliyun)
  --model MODEL                  LLM model name (default: deepseek-v3.2-exp)
  --api-key API_KEY             API key (or set via env var)
  --debug                       Enable debug logging
  -h, --help                    Show help message
```

## Troubleshooting

### "Failed to initialize LLM client"

Make sure your API key is set:
```bash
export DASHSCOPE_API_KEY=your_key  # For Aliyun
export OPENAI_API_KEY=your_key     # For OpenAI
```

### "LLM returned invalid JSON"

Try adding `--debug` flag to see full LLM response. The parser automatically extracts JSON from markdown code blocks.

### "Maximum iterations reached"

The validation loop has a 5-iteration limit. If you reach this, check:
- Is your SMILES format correct? (must have `>>` separator)
- Are you providing clear answers to clarification questions?

## Architecture

```
┌─────────────────────────────────────────────────────┐
│                 InteractiveCLI                      │
│  • User interaction                                 │
│  • Display formatting                               │
│  • Confirmation flows                               │
└─────────────────┬───────────────────────────────────┘
                  │
                  v
┌─────────────────────────────────────────────────────┐
│            NaturalLanguageParser                    │
│  • LLM client integration                           │
│  • Prompt templating                                │
│  • JSON extraction                                  │
└─────────────────┬───────────────────────────────────┘
                  │
                  v
┌─────────────────────────────────────────────────────┐
│              ParsedRequest Model                    │
│  • Data validation                                  │
│  • API request conversion                           │
│  • State management                                 │
└─────────────────┬───────────────────────────────────┘
                  │
                  v
┌─────────────────────────────────────────────────────┐
│          API Integration (main.py)                  │
│  • api_recommend_conditions()                       │
│  • RecommendConditionsRequest                       │
└─────────────────────────────────────────────────────┘
```

## Future Enhancements

- [ ] Batch processing from file
- [ ] Resume from saved draft
- [ ] Multi-reaction support
- [ ] Constraint templates (e.g., "green chemistry", "high throughput")
- [ ] Result visualization
- [ ] Export recommendations to CSV/Excel

## License

Same as main project.
