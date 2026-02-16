# Interactive CLI - Quick Start Guide

## Installation

Ensure you have the required dependencies:
```bash
pip install rdkit openai rxnmapper
```

Set your API key:
```bash
export OPENAI_API_KEY='sk-...'
```

## Usage Modes

### 1. Interactive Mode (Recommended for Exploration)

Start the interactive shell:
```bash
python reaction_agent/cli.py
```

Then enter reactions at the prompt:
```
reaction> CCBr>>CCN
reaction> Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
reaction> quit
```

**Commands in interactive mode:**
- `quit` or `exit` - Exit the program
- `config` - Show current configuration
- `help` - Show help message
- `batch <file>` - Analyze reactions from file

### 2. Single Reaction Mode

Analyze one reaction directly:
```bash
python reaction_agent/cli.py --reaction "CCBr>>CCN"
```

Save output to file:
```bash
python reaction_agent/cli.py --reaction "CCBr>>CCN" --output result.json
```

### 3. Batch Mode

Analyze multiple reactions from a file:
```bash
python reaction_agent/cli.py --batch reaction_agent/examples/sample_reactions.txt
```

Results are saved to `results/<filename>/` by default. Specify custom directory:
```bash
python reaction_agent/cli.py --batch reactions.txt --output-dir my_results/
```

### 4. Deterministic Only Mode (No LLM)

Run without requiring an API key:
```bash
python reaction_agent/cli.py --reaction "CCBr>>CCN" --no-llm
```

This performs:
- Cleaning and canonicalization
- Spectator detection
- Atom mapping (if rxnmapper available)
- Bond change extraction

## Configuration Options

### Model Selection

```bash
# Fast and cheap (default)
python reaction_agent/cli.py --model gpt-4o-mini

# Better quality
python reaction_agent/cli.py --model gpt-4o

# Latest GPT-5
python reaction_agent/cli.py --model gpt-5.2

# Reasoning model
python reaction_agent/cli.py --model o3-mini

# DeepSeek alternative
python reaction_agent/cli.py --provider aliyun --model deepseek-r1
```

### Analysis Options

```bash
# Skip atom mapping (faster)
python reaction_agent/cli.py --reaction "CCBr>>CCN" --skip-mapping

# Keep spectators in analysis
python reaction_agent/cli.py --reaction "CCBr>>CCN" --keep-spectators

# Adjust temperature (ignored for GPT-5/o-series)
python reaction_agent/cli.py --reaction "CCBr>>CCN" --temperature 0.7

# Adjust max tokens
python reaction_agent/cli.py --reaction "CCBr>>CCN" --max-tokens 3000
```

### Output Options

```bash
# Disable colors
python reaction_agent/cli.py --no-color

# Minimal output
python reaction_agent/cli.py --reaction "CCBr>>CCN" --quiet

# Save to custom location
python reaction_agent/cli.py --reaction "CCBr>>CCN" --output my_result.json
```

## Examples

### Example 1: Quick Analysis

```bash
python reaction_agent/cli.py --reaction "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
```

**Output:**
```
================================================================================
  INPUT
================================================================================
Raw:   Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
Clean: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1

================================================================================
  DETERMINISTIC ANALYSIS
================================================================================
Mapping QC: ✓ OK
  Confidence: 0.95

Bond Changes (4):
  BC1: broken bond between :7 and :13 (single)
  BC2: formed bond between :7 and :14 (single)
  BC3: broken bond between :14 and :15 (single)
  BC4: order_change bond between :8 and :9 (aromatic->single)

Reaction Center: 7, 8, 9, 13, 14, 15

================================================================================
  LLM INTERPRETATION
================================================================================
Reaction Class: cross_coupling
Tags: Buchwald-Hartwig, C-N coupling

Confidence: 0.85

Roles:
  electrophile: bromobenzene
  nucleophile: aniline
  leaving_group: bromide

Mechanistic Events (1):

  E1: substitution
    Bond changes: BC1, BC2
    Palladium-catalyzed C-N bond formation via oxidative addition, amine coordination, and reductive elimination
    Confidence: 0.88

Mechanism Summary:
  1. Oxidative addition of bromobenzene to Pd(0) catalyst
  2. Coordination and deprotonation of aniline
  3. Reductive elimination forming C-N bond and regenerating Pd(0)

================================================================================
  METADATA
================================================================================
Model: gpt-4o-mini (openai)
Tokens: 1245 (prompt: 567, completion: 678)
Latency: 2340 ms
Speed: 532.1 tokens/sec
```

### Example 2: Batch Processing

Create `my_reactions.txt`:
```
CCBr>>CCN
Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
CN(C)C.Clc1cc(Cl)nc(Cl)c1>>CN(C)c1cc(Cl)nc(Cl)c1
```

Run batch analysis:
```bash
python reaction_agent/cli.py --batch my_reactions.txt
```

**Output:**
```
================================================================================
  BATCH ANALYSIS (3 reactions)
================================================================================

[1/3] CCBr>>CCN...
  → nucleophilic_substitution (confidence: 0.82)

[2/3] Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1...
  → cross_coupling (confidence: 0.85)

[3/3] CN(C)C.Clc1cc(Cl)nc(Cl)c1>>CN(C)c1cc(Cl)nc(Cl)c1...
  → nucleophilic_substitution (confidence: 0.79)

================================================================================
  BATCH SUMMARY
================================================================================
Total: 3
Successful: 3
Failed: 0
✓ Batch results saved to results/my_reactions
```

### Example 3: Interactive Session

```bash
$ python reaction_agent/cli.py

================================================================================
  Reaction SMILES Analysis Agent - Interactive CLI
================================================================================

Model: gpt-4o-mini (openai)
ℹ Initializing LLM client...
✓ Initialized: gpt-4o-mini

================================================================================
  INTERACTIVE MODE
================================================================================

Enter reaction SMILES to analyze (or 'quit' to exit)
Commands:
  'quit' or 'exit' - Exit the program
  'config' - Show current configuration
  'help' - Show this help message
  'batch <file>' - Analyze reactions from file

reaction> CCBr>>CCN

Analyzing: CCBr>>CCN
--------------------------------------------------------------------------------
ℹ Step 1/2: Running deterministic analysis...
✓ Deterministic analysis complete
ℹ Step 2/2: LLM interpretation complete

[... results displayed ...]

reaction> config

Current Configuration:
  Model: gpt-4o-mini
  Provider: openai
  Temperature: 0.0
  Max Tokens: 2000
  Drop Spectators: True

reaction> quit

Goodbye!
```

### Example 4: Using Different Models

```bash
# Compare results across models
for model in gpt-4o-mini gpt-4o o3-mini gpt-5.2; do
  echo "Testing $model..."
  python reaction_agent/cli.py \
    --reaction "CCBr>>CCN" \
    --model $model \
    --quiet \
    --output "results_${model}.json"
done
```

### Example 5: Deterministic Analysis Only

For testing without API costs:
```bash
python reaction_agent/cli.py \
  --reaction "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1" \
  --no-llm
```

**Output:**
```json
{
  "input": {
    "rxn_smiles_raw": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "rxn_smiles_clean": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "reactants_clean": ["Brc1ccccc1", "Nc1ccccc1"],
    "products_clean": ["c1ccc(Nc2ccccc2)cc1"],
    "spectators": [],
    "parse_warnings": [],
    "standardization_actions": []
  },
  "tool_facts": {
    "mapped_rxn_smiles": "...",
    "mapping_qc": {"ok": true, "confidence": 0.95},
    "bond_changes": [...],
    "reaction_center_atoms": [7, 8, 9, 13, 14, 15]
  }
}
```

## File Format for Batch Processing

Create a text file with one reaction per line:

```txt
# Comments start with #
# Empty lines are ignored

CCBr>>CCN
Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1
CN(C)C.Clc1cc(Cl)nc(Cl)c1>>CN(C)c1cc(Cl)nc(Cl)c1

# More reactions...
```

## Troubleshooting

### Issue: API Key Not Found

```bash
$ python reaction_agent/cli.py
✗ OPENAI_API_KEY environment variable not set
```

**Solution:**
```bash
export OPENAI_API_KEY='sk-...'
```

### Issue: Colors Not Displaying

On Windows or when piping output, colors may not work. Disable them:
```bash
python reaction_agent/cli.py --no-color
```

### Issue: rxnmapper Not Installed

If you see warnings about missing rxnmapper:
```bash
pip install rxnmapper
```

Or skip mapping:
```bash
python reaction_agent/cli.py --reaction "CCBr>>CCN" --skip-mapping
```

### Issue: Import Errors

Make sure you're running from the project root:
```bash
cd /path/to/Condition-agent
python reaction_agent/cli.py
```

## Advanced Usage

### Piping Output

Process output with other tools:
```bash
python reaction_agent/cli.py --reaction "CCBr>>CCN" --quiet --no-color | jq .
```

### Scripting

Use in shell scripts:
```bash
#!/bin/bash

REACTIONS=(
  "CCBr>>CCN"
  "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
)

for rxn in "${REACTIONS[@]}"; do
  python reaction_agent/cli.py \
    --reaction "$rxn" \
    --output "results/$(echo $rxn | md5sum | cut -d' ' -f1).json"
done
```

### Custom Configuration

Set environment variables for defaults:
```bash
export OPENAI_API_KEY='sk-...'
export DEFAULT_MODEL='gpt-5.2'

python reaction_agent/cli.py  # Uses gpt-5.2 by default
```

## Tips

1. **Start with interactive mode** to explore and test reactions
2. **Use `--no-llm` mode** for quick validation without API costs
3. **Save outputs** with `--output` for later analysis
4. **Use batch mode** for processing multiple reactions efficiently
5. **Try different models** to compare quality and speed
6. **Use `--quiet`** when piping output or in scripts

## Getting Help

```bash
# Show all options
python reaction_agent/cli.py --help

# In interactive mode
reaction> help
```

## Next Steps

After using the CLI:
- Analyze your own reactions
- Process reaction databases
- Compare models and configurations
- Integrate into your workflow
- Build custom scripts using the Python API

For programmatic usage, see the main README and API documentation.
