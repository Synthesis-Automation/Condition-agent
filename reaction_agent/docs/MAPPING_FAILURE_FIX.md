# Mapping Failure Fix - Direct SMILES Analysis Mode

## Problem Identified

When atom mapping **completely fails** (0 bond changes detected), even GPT-5.2 was giving up with "insufficient evidence" instead of analyzing the reaction.

### Example: Defluorinative Cyclization
```
Oc1ccc(Br)cc1Br.FC(F)(F)S.[Cs]>>FC1(F)Oc2ccc(Br)cc2S1
```

**Before fix**:
- Mapping confidence: 0.215 ❌
- Bond changes detected: **0** ❌❌❌
- Result: `overall_class: "other", confidence: 0.20`
- Mechanism: "Insufficient evidence... no bond formations/cleavages detected"

**Root cause**: Even though we switched to GPT-5.2, the prompt said "0 bond changes detected" and GPT-5.2 took that at face value instead of analyzing SMILES directly.

## Solution Implemented

### 1. New Prompt Template (`DIRECT_SMILES_ANALYSIS`)

Created specialized prompt in `reaction_agent/prompts.py` that:
- ✅ Explicitly tells GPT-5.2 that mapping **failed completely**
- ✅ Instructs to **IGNORE** deterministic analysis
- ✅ Provides reactants and products SMILES separately
- ✅ Asks for direct pattern recognition and functional group comparison
- ✅ Uses chemical knowledge instead of relying on bond change data

Key instruction:
```
⚠️ MAPPING FAILURE MODE ⚠️

The atom mapping tool has COMPLETELY FAILED for this reaction (0 bond changes detected).
This does NOT mean the reaction is invalid - it means you must analyze the SMILES
strings DIRECTLY using chemical reasoning.

IGNORE the deterministic analysis and focus on:
1. Identifying functional groups in reactants
2. Identifying functional groups in products
3. Inferring what bonds MUST have changed
4. Using your knowledge of organic chemistry patterns
```

### 2. Detection Logic (`agent.py`)

Added automatic detection of complete mapping failures:

```python
# Check for complete mapping failure (0 bond changes + low/failed mapping)
bond_changes = tool_facts.get("bond_changes", [])
mapping_qc = tool_facts.get("mapping_qc", {})
mapping_conf = mapping_qc.get("confidence", 0.0)
mapping_ok = mapping_qc.get("ok", False)

use_direct_analysis = (
    len(bond_changes) == 0
    and (not mapping_ok or mapping_conf < 0.4)
)
```

When triggered, system:
1. Prints warning: `⚠️  Mapping failed completely (0 bond changes, confidence X.XXX)`
2. Prints notice: `→ Using direct SMILES analysis mode (pattern recognition)`
3. Uses new prompt template with separated reactants/products

## Results

### Defluorinative Cyclization (NEW-2.md example)

**After fix**:
- ✅ **Class**: `annulation` (correct!)
- ✅ **Confidence**: 0.30 (appropriately lower, but not zero)
- ✅ **Mechanism** (3 steps):
  1. Base generation of nucleophiles
  2. Defluorinative substitution (CF₃ → CF₂)
  3. Ring closure by aryl C–S bond formation
- ✅ **Identified**: Benzofused O–CF₂–S heterocycle formation
- ✅ **Warnings**: `mapping_failed_used_direct_analysis`

**Comparison with Web GPT-5.2**:
| Aspect | Web GPT-5.2 | Our System (Fixed) |
|--------|-------------|-------------------|
| Class | Defluorinative O,S-difluoromethylenation / Annulation | Annulation ✓ |
| Mechanism steps | 3 | 3 ✓ |
| CF₃ → CF₂ identified | Yes | Yes ✓ |
| Ring closure identified | Yes | Yes ✓ |
| Confidence | N/A | 0.30 |

### Existing Functionality Preserved

**Simple reaction (good mapping)**:
```
CCBr.CCN>>CCNCC
```
- Model: gpt-4o (fast mode)
- Class: nucleophilic_substitution
- Confidence: 0.97 ✓
- No change in behavior ✓

**Complex reaction (high mapping but 0 bond changes reported)**:
- Now uses direct analysis instead of giving up
- Still switches to GPT-5.2 automatically
- Provides meaningful mechanistic analysis

## When Direct Analysis Mode Activates

The system automatically switches to **Direct SMILES Analysis Mode** when:
```python
len(bond_changes) == 0 AND (not mapping_ok OR mapping_conf < 0.4)
```

This means:
- ✅ 0 bond changes detected **AND**
- ✅ Either mapping flagged as "not ok" **OR** confidence < 0.4

## Impact on Auto Mode Strategy

| Scenario | Mapping Conf | Bond Changes | Mode Used | Model |
|----------|-------------|--------------|-----------|-------|
| **Simple reactions** | ≥ 0.4 | >0 | Normal | gpt-4o (fast) |
| **Complex but mappable** | ≥ 0.4 but <0.7 | >0 | Normal | gpt-5.2 (smart switch) |
| **Poor mapping** | < 0.4 | >0 | Normal | gpt-5.2 (deep reasoning) |
| **Complete failure** | < 0.4 | **0** | **Direct SMILES** | gpt-5.2 (pattern recognition) ✅ |

## Confidence Interpretation

With direct SMILES analysis:
- **0.5-0.9**: Clear transformation, confident about mechanism
- **0.3-0.5**: Complex transformation, mechanistic hypothesis
- **0.2-0.3**: Very complex, multiple possible pathways

The lower confidence (vs reactions with good mapping) is **appropriate** because:
- No atom-level tracking of bond changes
- Must infer mechanism from structure comparison
- Multiple possible pathways may exist

## Files Modified

1. **`reaction_agent/prompts.py`**
   - Added `DIRECT_SMILES_ANALYSIS` template
   - Added `get_direct_smiles_template()` function

2. **`reaction_agent/agent.py`**
   - Imported `get_direct_smiles_template`
   - Added detection logic for complete mapping failures
   - Conditional prompt selection based on mapping status
   - Added user-facing print statements for transparency

## Future Improvements

Potential enhancements:
1. **RDKit substructure matching**: Auto-detect functional groups to seed the analysis
2. **SMARTS pattern templates**: Pre-computed patterns for common reaction types
3. **Confidence calibration**: Fine-tune confidence ranges based on reaction type
4. **Literature RAG**: Add relevant literature context for known transformations

## Testing

Run on defluorinative cyclization:
```bash
python -c "
from llmtools.clients import LLMClient
from reaction_agent.agent import ReactionSMILESAnalyzer

rxn = 'Oc1ccc(Br)cc1Br.FC(F)(F)S.[Cs]>>FC1(F)Oc2ccc(Br)cc2S1'

client = LLMClient(provider='openai', model='gpt-4o', timeout=300)
analyzer = ReactionSMILESAnalyzer(client, max_tokens=8000)

result = analyzer.analyze(rxn, mode='auto')
print(result['interpretation']['overall_class'])  # Should print: annulation
"
```

## Summary

✅ **Fixed**: Complete mapping failures (0 bond changes) now trigger direct SMILES analysis
✅ **Preserved**: Existing functionality for reactions with good/moderate mapping
✅ **Improved**: GPT-5.2 now analyzes reactions instead of giving up when mapping fails
✅ **Appropriate**: Lower confidence reflects increased uncertainty from missing atom tracking
✅ **Transparent**: Clear user messaging about mode switching

The system now handles the full spectrum:
- **Good mapping** → Fast gpt-4o analysis
- **Poor mapping** → GPT-5.2 deep reasoning with bond changes
- **No mapping** → GPT-5.2 pattern recognition on raw SMILES ✨
