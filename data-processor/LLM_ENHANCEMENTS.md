# LLM-Enhanced Reagent Taxonomy Generator

## Overview

The reagent taxonomy generator (`reagent_taxonomy_qt.py`) has been enhanced with advanced LLM capabilities to improve accuracy and reliability of reagent classification and field assignments.

## Three Key Enhancements

### 1. Automatic Role Upgrade for "other_reagent"

**Problem**: When reagents are classified as `other_reagent`, it indicates the deterministic system couldn't identify a more specific role.

**LLM Solution**: When LLM is enabled and the initial role is `other_reagent`, the LLM automatically analyzes the reagent chemistry and suggests a more appropriate role.

**How it works**:
- LLM receives full reagent context (name, CAS, synonyms, structure, family)
- Evaluates against all available role types in the taxonomy
- Provides confidence score and justification for suggested role
- System automatically applies the upgrade if valid

**Example Output**:
```json
{
  "llm_auto_upgrade": {
    "from": "other_reagent",
    "to": "base",
    "reason": "This compound is a tertiary amine with strong basicity",
    "confidence": 0.95
  }
}
```

### 2. Intelligent Field Assignment

**Problem**: Deterministic rules may miss critical chemical properties like `basicity`, `oxidation_states`, `polarity`, etc.

**LLM Solution**: LLM analyzes reagent structure and chemistry to suggest reliable values for role-specific fields.

**Supported Fields by Role**:

**Base**:
- `basicity`: "weak" | "moderate" | "strong" | "superbase"
- `nucleophilicity`: "weak" | "moderate" | "strong"
- `sterics`: "unhindered" | "moderate" | "hindered"

**Metal Precursor/Preformed Catalyst**:
- `metal`: Element symbol (e.g., "Pd", "Ni", "Cu")
- `oxidation_states`: List of oxidation states (e.g., [0], [2], [0, 2])
- `ligand_type`: Type of ligand coordination

**Solvent**:
- `proticity`: "aprotic" | "protic"
- `polarity`: "polar" | "nonpolar" | "intermediate"
- `coordination`: "coordinating" | "weakly_coordinating" | "non_coordinating"

**Oxidant/Reductant/Condensation Agent**:
- `strength_band`: "weak" | "moderate" | "strong"

**Ligand**:
- `donors`: List of donor atoms (e.g., ["P"], ["N", "N"], ["P", "P"])
- `denticity`: Number (1, 2, 3, etc.)

**Organo-catalyst**:
- `activation_mode`: Type of catalysis (e.g., "lewis_base", "phase_transfer")
- `chirality`: "achiral" | "chiral"

**Enzyme**:
- `source`: Origin (e.g., "bacterial", "fungal", "mammalian")
- `cofactor_requirement`: Required cofactors

**Example Output**:
```json
{
  "changes_applied": {
    "field_suggestions_applied": {
      "basicity": "strong",
      "nucleophilicity": "moderate",
      "sterics": "hindered"
    }
  }
}
```

### 3. Enhanced Review Output Format

**Problem**: Need clear visibility into what LLM changed and why.

**Solution**: Structured output showing original vs. revised entries with comprehensive review metadata.

**Output Structure**:

```json
{
  "review_summary": {
    "original_role": "other_reagent",
    "original_family": "misc_general",
    "llm_status": "needs_review",
    "confidence": 0.92,
    "justification": "Compound is clearly a tertiary aliphatic amine with strong basicity",
    "auto_upgrade": {
      "from": "other_reagent",
      "to": "base",
      "reason": "Strong tertiary amine characteristics",
      "confidence": 0.92
    }
  },
  "entry_original": {
    "id": "cas-121-44-8",
    "name": "Triethylamine",
    "cas": "121-44-8",
    "roles": {
      "other_reagent": {
        "families": ["misc_general"]
      }
    }
  },
  "entry_revised": {
    "id": "cas-121-44-8",
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
  },
  "llm_review_details": {
    "model": "gpt-4o-mini",
    "provider": "openai",
    "tokens_used": 456,
    "latency_ms": 1823
  },
  "llm_alerts": [
    "Consider adding 'TEA' as common abbreviation"
  ],
  "changes_applied": {
    "role": {
      "from": "other_reagent",
      "to": "base"
    },
    "family": {
      "from": "misc_general",
      "to": "tertiary_amines_aliphatic"
    },
    "field_suggestions_applied": {
      "basicity": "strong",
      "nucleophilicity": "moderate",
      "sterics": "unhindered"
    },
    "synonyms_added": ["TEA"]
  }
}
```

## Using the Enhanced Generator

### 1. Enable LLM in the UI

1. Launch the generator: `python data-processor/reagent_taxonomy_qt.py`
2. In the "LLM assistance" dropdown, select **"Use LLM"**
3. Choose your provider (OpenAI or Aliyun/DeepSeek)
4. Select a model (recommended models are marked)

### 2. Configure API Keys

**Windows (PowerShell)**:
```powershell
[System.Environment]::SetEnvironmentVariable('OPENAI_API_KEY', 'sk-...', 'User')
[System.Environment]::SetEnvironmentVariable('ALIYUN_API_KEY', 'sk-...', 'User')
```

See `docs/WINDOWS_ENV_SETUP.md` for detailed instructions.

### 3. Generate Entry with LLM Review

1. Enter CAS number (e.g., `121-44-8` for triethylamine)
2. Select role (try `other_reagent` to test auto-upgrade)
3. Click **Generate**
4. Review the comprehensive output showing:
   - Original deterministic assignment
   - LLM suggestions and confidence
   - Revised entry with field assignments
   - Change summary

### 4. Review and Edit

The output is editable JSON. You can:
- Modify field values
- Adjust role/family assignments
- Add/remove synonyms
- Change any property

### 5. Save to Registry

Click **Save** to persist the entry (uses `entry_revised` if available, otherwise `entry`).

## Recommended Models

### OpenAI
- **Fast**: `gpt-4o-mini` (cost-effective, good accuracy)
- **Balanced**: `gpt-4o` (better reasoning)
- **Advanced**: `gpt-5-mini` (latest features)

### Aliyun/DeepSeek
- **Fast**: `deepseek-r1-distill-qwen-7b` (cost-effective)
- **Balanced**: `deepseek-v3` (production-ready)
- **Reasoning**: `deepseek-r1` (strong chemistry knowledge)

## Example Workflows

### Workflow 1: Auto-upgrade "other_reagent"

```
Input:
- CAS: 7664-38-2
- Role: other_reagent
- LLM: Enabled

Output:
✓ LLM identifies phosphoric acid as an acid
✓ Auto-upgrade: other_reagent → acid
✓ Family: mineral_acids
✓ Fields: acidity = "strong"
```

### Workflow 2: Enhance base classification

```
Input:
- CAS: 1310-73-2
- Role: base
- LLM: Enabled

Output:
✓ Confirms role: base
✓ Refined family: metal_hydroxides
✓ Fields: basicity = "superbase", nucleophilicity = "moderate"
✓ Added synonym: "caustic soda"
```

### Workflow 3: Metal precursor with oxidation states

```
Input:
- CAS: 14024-61-4
- Role: metal_precursor
- LLM: Enabled

Output:
✓ Confirms role: metal_precursor
✓ Fields: metal = "Pd", oxidation_states = [2]
✓ Family: pd_ii_salts
```

## Prompt Engineering Details

The LLM receives a structured prompt with:

1. **Reagent Identity**: Name, CAS, synonyms, abbreviations
2. **Deterministic Resolution**: Source, SMILES, token matches
3. **Proposed Assignment**: Role, family, keywords, examples
4. **Existing Conflicts**: Registry duplicates

The prompt explicitly instructs the LLM to:
- **Strongly recommend** specific roles for `other_reagent`
- Suggest **most specific family** that fits the chemistry
- Identify **critical missing fields** for the role
- Add **well-known synonyms** or trade names

See `llmtools/prompts.py` → `REAGENT_REGISTRY_REVIEW` for full template.

## LLM Response Schema

```json
{
  "status": "confirm" | "needs_review" | "reject",
  "proposed_role": "string",
  "proposed_family": "string",
  "confidence": 0.0-1.0,
  "justification": "short explanation",
  "alerts": ["list of warnings"],
  "suggested_synonyms": ["additional names"],
  "field_suggestions": {
    "field_name": "suggested_value"
  }
}
```

## Integration Architecture

```
reagent_taxonomy_qt.py (UI)
    ↓
generate_taxonomy_entry() (Core logic)
    ↓
review_taxonomy_proposal() (llmtools/reagent_review.py)
    ↓
LLMClient.chat() (llmtools/clients.py)
    ↓
OpenAI/Aliyun API
```

## Performance Notes

**Typical LLM call**:
- Tokens: 300-600 (prompt) + 200-400 (completion)
- Latency: 1-5 seconds
- Cost: ~$0.001-0.005 per entry (gpt-4o-mini)

**Batch processing**: Consider caching or batch API calls for large datasets.

## Troubleshooting

### LLM suggests invalid role

**Issue**: `proposed_role` not in `ROLE_CONFIG`

**Solution**: System automatically rejects invalid roles with error message in `adjustment_errors`

### Field suggestions don't apply

**Issue**: Field name not in role's expected fields

**Solution**: Check `ROLE_PAYLOAD_FIELDS` in `reagent_taxonomy_qt.py` for valid fields per role

### LLM parse error

**Issue**: LLM didn't return valid JSON

**Solution**: 
- Check `llm_review.status == "parse_error"`
- Review `llm_review.raw_response` for debugging
- Try different model or increase `max_tokens`

## Testing

Test the enhancements with these CAS numbers:

1. **Auto-upgrade test**: `7664-38-2` (phosphoric acid as "other_reagent")
2. **Field assignment test**: `121-44-8` (triethylamine)
3. **Metal precursor test**: `14024-61-4` (Pd(OAc)₂)
4. **Solvent test**: `67-56-1` (methanol)

## Related Documentation

- `llmtools/README.md` - LLM client and agent documentation
- `llmtools/reagent_review.py` - Review function implementation
- `docs/WINDOWS_ENV_SETUP.md` - API key configuration
- `AGENTS.md` - Project structure and conventions
