# CLI Reaction Type Routing Enhancement

## Overview

The CLI now includes **intelligent reaction type determination** before submission, which automatically routes recommendations to the appropriate system (rule-based vs machine learning) based on:

1. Detected reaction family from SMILES structure
2. User-specified metal catalyst preferences
3. Required reagents mentioned in constraints

## Key Features

### 🎯 Automatic Reaction Type Detection

The CLI now performs a **final determination step** before submission that:
- Analyzes the reaction SMILES structure
- Checks for metal catalyst preferences in constraints
- Determines the specific reaction subtype (e.g., Ullmann vs Buchwald for C-N coupling)
- Routes to the appropriate recommendation engine

### 🎨 Color-Coded Output

- **Cyan**: Headers, prompts, informational text
- **Green**: Success messages, valid states
- **Yellow**: Warnings, examples, routing information
- **Red**: Errors, validation issues
- **Magenta**: Main headers

### 🧪 Test Mode

Use `--test` flag to validate input processing without actual API submission:
```bash
python app/cli_recommend.py --test
```

## Reaction Type Routing Logic

### C-N Coupling Reactions

The system intelligently routes C-N coupling reactions based on the metal catalyst:

| Catalyst | Final Type | Reaction Name | System | Description |
|----------|------------|---------------|--------|-------------|
| **Cu** (Copper) | `C_N_Coupling_Cu` | Ullmann | **Rule-based** | Deterministic constraint matching |
| **Pd** (Palladium) | `C_N_Coupling_Pd` | Buchwald-Hartwig | **ML-based** | Similarity search with neural networks |
| **Ni** (Nickel) | `C_N_Coupling_Ni` | Nickel-catalyzed | **ML-based** | Similarity search |
| None specified | Auto-detected | Generic | **Auto-detect** | Uses structural analysis |

### Detection Priority

1. **Explicit metal preference**: `constraints.metal_preference` (Cu, Pd, Ni)
2. **Required reagents**: Searches for metal mentions in `constraints.required_reagents`
3. **Structural analysis**: SMARTS-based pattern matching on reaction SMILES
4. **User-provided type**: Falls back to initially provided reaction type

## Usage Examples

### Example 1: Copper-Catalyzed C-N Coupling (Ullmann)

```bash
python app/cli_recommend.py --test
```

**Input:**
- Reaction: `Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1`
- Requirements: "use copper catalyst, no strong base"

**Output:**
```
🔬 REACTION TYPE DETERMINATION

Initial/Provided Type: None
Detected Family: Ullmann_CN (confidence: 0.90)
Detected Metal Catalyst: Cu

✓ Final Reaction Type: C_N_Coupling_Cu
Routing: metal-specific: Ullmann (rule-based)
🎯 Recommendation System: RULE-BASED (deterministic constraints)

ℹ️  C-N Coupling Routing:
   → Ullmann reaction (copper-catalyzed)
   → Uses rule-based constraint matching
```

### Example 2: Palladium-Catalyzed C-N Coupling (Buchwald)

**Input:**
- Reaction: `Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1`
- Requirements: "prefer palladium catalyst, mild conditions"

**Output:**
```
🔬 REACTION TYPE DETERMINATION

Initial/Provided Type: None
Detected Family: Ullmann_CN (confidence: 0.90)
Detected Metal Catalyst: Pd

✓ Final Reaction Type: C_N_Coupling_Pd
Routing: metal-specific: Buchwald (ML-based)
🎯 Recommendation System: MACHINE LEARNING (similarity-based)

ℹ️  C-N Coupling Routing:
   → Buchwald-Hartwig reaction (palladium-catalyzed)
   → Uses ML-based similarity search
```

### Example 3: Auto-Detection from Structure

**Input:**
- Reaction: `Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1`
- Requirements: "no specific requirements"

**Output:**
```
🔬 REACTION TYPE DETERMINATION

Initial/Provided Type: None
Detected Family: Suzuki_CC (confidence: 0.90)
Detected Metal Catalyst: None

✓ Final Reaction Type: Suzuki
Routing: auto-detected
```

## Command-Line Options

```bash
python app/cli_recommend.py [OPTIONS]
```

### Available Options

| Option | Description | Default |
|--------|-------------|---------|
| `--provider {openai,aliyun}` | LLM provider | `aliyun` |
| `--model MODEL` | LLM model name | `deepseek-v3.2-exp` |
| `--api-key API_KEY` | API key | Environment variable |
| `--test` | Test mode (stop before submission) | `False` |
| `--no-color` | Disable colored output | `False` |
| `--debug` | Enable debug logging | `False` |

### Examples

```bash
# Test mode with default settings
python app/cli_recommend.py --test

# Production mode with OpenAI
python app/cli_recommend.py --provider openai --model gpt-4o

# Test mode without colors (for logging)
python app/cli_recommend.py --test --no-color

# Debug mode to see detailed logs
python app/cli_recommend.py --debug --test
```

## How Metal Detection Works

### 1. Explicit Metal Preference
```python
constraints = {
    "metal_preference": "Cu"  # Directly specifies copper
}
# → Routes to C_N_Coupling_Cu (Ullmann)
```

### 2. Required Reagents Analysis
```python
constraints = {
    "required_reagents": ["copper catalyst", "CuI"]
}
# → Detects "copper" or "cu" in text
# → Routes to C_N_Coupling_Cu (Ullmann)
```

### 3. Structural Analysis
If no explicit metal is specified, the router analyzes:
- Reactant functional groups (aryl halide + amine → C-N coupling)
- Agent/catalyst components in reaction SMILES
- Historical patterns in reaction database

## Implementation Details

### Function: `determine_final_reaction_type()`

Located in `app/cli_recommend.py`, this function:

1. Calls `chemtools.router.detect_family_from_reaction()` for structural analysis
2. Extracts metal preference from constraints
3. Applies routing logic for C-N coupling reactions
4. Normalizes to canonical family name via `chemtools.recommend.utils.canonical_family()`

### Returns

```python
{
    "final_type": "C_N_Coupling_Cu",           # Final determined type
    "initial_type": "C_N_Coupling",            # User-provided type
    "detected_family": "Ullmann_CN",           # Structurally detected
    "confidence": 0.90,                        # Detection confidence
    "detected_metal": "Cu",                    # Extracted metal catalyst
    "routing": "metal-specific: Ullmann...",   # Routing explanation
    "detection_info": {...}                    # Full detection details
}
```

## Testing

### Unit Test for Routing Logic

Run the test script:
```bash
python test_cli_routing.py
```

This demonstrates routing for:
- ✓ C-N Coupling with Copper → Ullmann (rule-based)
- ✓ C-N Coupling with Palladium → Buchwald (ML)
- ✓ C-N Coupling auto-detect → Based on structure
- ✓ Suzuki Coupling → Auto-detected
- ✓ User-specified with metal in reagents → Properly routed

### Interactive Testing

```bash
# Start in test mode
python app/cli_recommend.py --test

# Try these test cases:
# 1. Ullmann: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
#    Requirements: "use copper catalyst"

# 2. Buchwald: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
#    Requirements: "prefer palladium catalyst"

# 3. Suzuki: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
#    Requirements: (leave empty)
```

## Benefits

1. **Accurate Routing**: Ensures recommendations use the correct engine (rule-based vs ML)
2. **Better Results**: Ullmann reactions get deterministic constraint matching; Buchwald gets similarity-based ML
3. **User Transparency**: Clear indication of which system is being used
4. **Safe Testing**: Test mode validates input without actual API calls
5. **Visual Feedback**: Color-coded output for easy interpretation

## Integration with Existing System

This enhancement integrates seamlessly with:

- `chemtools/router.py`: Structural reaction family detection
- `chemtools/recommend/core.py`: Recommendation engine
- `chemtools/recommend/utils.py`: Family name normalization
- `llmtools/`: LLM-powered natural language parsing

## Future Enhancements

- [ ] Support for more metal catalysts (Ir, Rh, Au, etc.)
- [ ] Enhanced confidence scoring for routing decisions
- [ ] User feedback loop for improving detection accuracy
- [ ] Integration with reaction database statistics

## Related Documentation

- `app/CLI_IMPLEMENTATION_SUMMARY.md`: General CLI documentation
- `app/CLI_QUICKSTART.md`: Quick start guide
- `AGENTS.md`: Repository structure and guidelines
- `chemtools/router.py`: Reaction family detection logic
