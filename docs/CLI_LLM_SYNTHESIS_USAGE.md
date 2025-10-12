# CLI Tool - LLM Multi-Source Synthesis Testing

## Overview

The `local_recommendation_cli.py` tool has been updated to support testing the new LLM-enhanced multi-source synthesis recommendation system. This allows you to test the integration of ML precedents, rule-based matching, and protocol-based recommendations through LLM synthesis.

## Installation

Ensure all dependencies are installed:

```powershell
pip install -r requirements.txt
```

## Basic Usage

### Test LLM Synthesis Mode

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
```

### With Custom Constraints

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --constraints '{\"scale\": \"multigram\", \"cost\": \"low\"}'
```

### Using V1 Prompt (Original)

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --llm-prompt-version v1
```

### Custom LLM Provider/Model

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --llm-provider openai `
  --llm-model gpt-4
```

## Command-Line Arguments

### Required

- `--rxn SMILES`: Reaction SMILES string (or use `-i` for interactive mode)
- `--strategy STRATEGY`: Recommendation strategy

### Strategy Options

- `all` - Run all recommendation modes (rule, ML, fusion, protocol, LLM)
- `rule` - Rule-based matching only
- `ml` - ML-based precedent search only
- `fusion` - Fusion of rule + ML
- `protocol` - Protocol-based recommendations
- `llm` - **NEW** LLM multi-source synthesis (combines rule + ML + protocol)

### LLM-Specific Arguments

- `--llm-provider {aliyun,openai}` - LLM provider (default: aliyun)
- `--llm-model MODEL` - LLM model name (default: deepseek-v3.2-exp)
- `--llm-prompt-version {v1,v2}` - Prompt version (default: v2)
- `--constraints JSON` - User constraints as JSON string

### Other Arguments

- `--rxn-type {suzuki,ullmann,buchwald}` - Reaction type (auto-detected if omitted)
- `--label LABEL` - Custom label for output files
- `--scdb PATH` - Path to SCDB database (for rule-based)
- `-k K` - Number of ML neighbors (default: 3)
- `--limit LIMIT` - Max protocol matches (default: 5)
- `-i, --interactive` - Interactive mode (prompts for reaction)

## Output

### Console Output

The tool displays a summary of the LLM synthesis result:

```
LLM Multi-Source Synthesis (deepseek-v3.2-exp):
  Status: success
  Total processing time: 45823.4 ms
  LLM latency: 43.6s
  Tokens used: 1480
  Sources: ML=1, Rule=1, Protocol=1

  Confidence: HIGH
  Reasoning: All three sources agree on Pd(PPh3)4 catalyst...

  RECOMMENDED CONDITIONS:
    Catalyst: Pd(PPh3)4
    Ligand: PPh3 (pre-complexed)
    Solvent: THF
    Temperature: 80°C
    Base: K2CO3
    Rationale: High consensus across all sources. ML precedent shows 0.95 similarity...

  BACKUP CONDITIONS (2):
    1. Pd(dppf)Cl2
       When: If Pd(PPh3)4 is unavailable or cost-sensitive...
    2. PdCl2(dppf)
       When: Alternative to Pd(dppf)Cl2 with similar performance...

  WARNINGS (1):
    • Electron-poor aryl halides may require elevated temperature...

  Consensus: Catalyst=high, Solvent=high
```

### JSON Output

Results are saved to `results/` directory:

- `{timestamp}_{label}_llm_synthesis_local.json` - Full LLM synthesis result

JSON structure:

```json
{
  "status": "success",
  "synthesis": {
    "consensus_analysis": {
      "catalyst": {
        "agreement": "high",
        "consensus_value": "Pd(PPh3)4",
        "notes": "All sources recommend Pd(0) catalysts..."
      },
      "solvent": {...},
      "temperature": {...},
      "base": {...}
    },
    "confidence_level": "high",
    "confidence_reasoning": "...",
    "recommended_condition": {
      "catalyst": "Pd(PPh3)4",
      "ligand": "PPh3 (pre-complexed)",
      "solvent": "THF",
      "temperature": "80°C",
      "base": "K2CO3",
      "additive": "None",
      "rationale": "..."
    },
    "backup_conditions": [...],
    "warnings": [...],
    "source_comparison": {
      "ml_contribution": "...",
      "rule_contribution": "...",
      "protocol_contribution": "..."
    }
  },
  "sources_used": {
    "ml_precedents": 1,
    "rule_matches": 1,
    "protocol_procedures": 1
  },
  "llm_metadata": {
    "model": "deepseek-v3.2-exp",
    "provider": "aliyun",
    "tokens": 1480,
    "latency_ms": 43600,
    "processing_time_ms": 45823
  }
}
```

## Examples

### Example 1: Suzuki Coupling (High Confidence)

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --rxn-type suzuki
```

Expected: High confidence, Pd(PPh3)4 recommendation, THF/K2CO3

### Example 2: Buchwald-Hartwig with Cost Constraint

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1" `
  --rxn-type buchwald `
  --constraints '{\"cost\": \"low\"}'
```

Expected: Medium confidence, cost-effective catalyst suggestion

### Example 3: Ullmann with Scale Constraint

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1" `
  --rxn-type ullmann `
  --constraints '{\"scale\": \"multigram\", \"cost\": \"low\"}'
```

Expected: Recommendations optimized for larger scale and cost

### Example 4: Compare V1 vs V2 Prompts

Run with V1:
```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --llm-prompt-version v1
```

Run with V2:
```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --llm-prompt-version v2
```

Compare latency and quality (V2 should be ~28% faster with same quality)

### Example 5: Run All Strategies for Comparison

```powershell
python scripts/local_recommendation_cli.py `
  --strategy all `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
```

Generates 5 JSON files:
- `*_rule_local.json`
- `*_ml_local.json`
- `*_fusion_local.json`
- `*_protocol_local.json`
- `*_llm_synthesis_local.json`

## Troubleshooting

### Issue: LLM imports not available

**Error**: `ModuleNotFoundError: No module named 'llmtools'`

**Solution**: Install LLM dependencies:
```powershell
pip install openai python-dotenv
```

### Issue: API key not configured

**Error**: `LLMClient initialization failed`

**Solution**: Set environment variables:
```powershell
$env:OPENAI_API_KEY = "sk-..."  # For OpenAI
$env:ALIYUN_API_KEY = "sk-..."  # For Aliyun/DeepSeek
```

### Issue: Constraints JSON parsing error

**Error**: `Warning: Failed to parse constraints JSON`

**Solution**: Ensure JSON is properly escaped:
```powershell
# Correct (PowerShell)
--constraints '{\"cost\": \"low\"}'

# Wrong
--constraints '{"cost": "low"}'  # Missing escapes
```

### Issue: Protocol database not built

**Warning**: Protocol recommendations may be empty

**Solution**: Build protocol index:
```powershell
python -m chemtools.protocol.cli build
```

## Performance Notes

- **V2 Prompt**: ~28% faster than V1 (43.6s vs 60.4s average)
- **Token Usage**: V2 uses ~1480 tokens vs V1's ~1632 tokens
- **Quality**: Both versions maintain 5.0/5.0 quality score
- **Recommendation**: Use V2 by default for better performance

## Integration with Step 5

This CLI tool enables rapid testing before production API integration (Step 5). Once validated:

1. Test various reaction types and constraints
2. Verify chemistry guidelines are working correctly
3. Validate confidence thresholds
4. Ensure backup conditions are appropriate
5. Check warnings are triggered correctly

Then proceed to integrate into FastAPI (`app/main.py`) with confidence.

## Related Documentation

- [Step 4: Prompt Refinement](./STEP4_PROMPT_REFINEMENT_COMPLETE.md)
- [LLM Tools README](../llmtools/README.md)
- [API Documentation](./API_DOCUMENTATION.md)
