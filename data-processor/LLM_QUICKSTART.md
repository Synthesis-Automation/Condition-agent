# Quick Start: LLM-Enhanced Reagent Generator

## Setup (One-time)

### 1. Set API Key

**PowerShell**:
```powershell
[System.Environment]::SetEnvironmentVariable('OPENAI_API_KEY', 'sk-your-key-here', 'User')
```

Or for Aliyun/DeepSeek:
```powershell
[System.Environment]::SetEnvironmentVariable('ALIYUN_API_KEY', 'sk-your-key-here', 'User')
```

**Restart** your terminal after setting keys.

### 2. Verify Installation

```powershell
python -c "from llmtools.clients import LLMClient; print('✓ LLM tools ready')"
```

## Basic Usage

### Launch the Generator

```powershell
cd c:\Git-softwares\Condition-agent
python data-processor\reagent_taxonomy_qt.py
```

### Enable LLM

1. **LLM assistance** dropdown → Select **"Use LLM"**
2. **Provider** → Select **"Openai"** or **"Aliyun"**
3. **Model** → Select recommended model (marked with "(recommended)")

### Generate Entry

1. **CAS number**: Enter reagent CAS (e.g., `121-44-8`)
2. **Reagent role**: Select role (or "Other Reagent" to test auto-upgrade)
3. Click **Generate**
4. Review JSON output in the editor
5. Click **Save** to persist to registry

## Three Key Features

### ✨ Feature 1: Auto-upgrade "other_reagent"

**Try this**: CAS `7664-38-2` (phosphoric acid) with role "Other Reagent"

**Expected**: LLM upgrades to → `acid` with family `mineral_acids`

### ✨ Feature 2: Smart field assignment

**Try this**: CAS `121-44-8` (triethylamine) with role "Base"

**Expected**: LLM adds:
- `basicity: "strong"`
- `nucleophilicity: "moderate"`
- `sterics: "unhindered"`

### ✨ Feature 3: Review-ready output

**Output shows**:
- `review_summary`: Original vs. LLM suggestions
- `entry_original`: Deterministic assignment
- `entry_revised`: LLM-enhanced entry with fields
- `changes_applied`: What changed and why

## Example Session

```text
1. Launch generator
2. Enable LLM (OpenAI, gpt-4o-mini)
3. Enter CAS: 121-44-8
4. Select role: Other Reagent
5. Click Generate

Output:
{
  "review_summary": {
    "original_role": "other_reagent",
    "llm_status": "needs_review",
    "confidence": 0.95,
    "justification": "Tertiary aliphatic amine with strong basicity",
    "auto_upgrade": {
      "from": "other_reagent",
      "to": "base"
    }
  },
  "entry_revised": {
    "name": "Triethylamine",
    "cas": "121-44-8",
    "roles": {
      "base": {
        "families": ["tertiary_amines_aliphatic"],
        "basicity": "strong",
        "nucleophilicity": "moderate",
        "sterics": "unhindered"
      }
    }
  }
}

6. Review output (edit if needed)
7. Click Save
8. Done! Entry saved to data/reagents/base.json
```

## Test Cases

Try these CAS numbers to see LLM enhancements:

| CAS | Name | Role | LLM Enhancement |
|-----|------|------|-----------------|
| `7664-38-2` | Phosphoric acid | other_reagent | → acid |
| `121-44-8` | Triethylamine | base | + basicity/nucleophilicity |
| `14024-61-4` | Pd(OAc)₂ | metal_precursor | + metal/oxidation_states |
| `67-56-1` | Methanol | solvent | + proticity/polarity |
| `1310-73-2` | NaOH | base | + basicity (superbase) |

## Recommended Settings

**For accuracy**:
- Model: `gpt-4o` (OpenAI) or `deepseek-r1` (Aliyun)
- Temperature: 0.0 (deterministic)

**For speed/cost**:
- Model: `gpt-4o-mini` (OpenAI) or `deepseek-r1-distill-qwen-7b` (Aliyun)
- Temperature: 0.0

**For experimentation**:
- Model: `gpt-5-mini` (OpenAI) or `deepseek-v3.2-exp` (Aliyun)
- Temperature: 0.2 (slightly creative)

## Common Issues

### "LLM support unavailable"

**Fix**: 
```powershell
pip install openai  # or your chosen provider SDK
```

### "API key not found"

**Fix**: Set environment variable and **restart terminal**

### "Parse error" in LLM response

**Fix**: Try different model or check `llm_review.raw_response` for debugging

## Full Documentation

See `data-processor/LLM_ENHANCEMENTS.md` for complete guide.

## Support

- LLM module docs: `llmtools/README.md`
- API key setup: `docs/WINDOWS_ENV_SETUP.md`
- Project structure: `AGENTS.md`
