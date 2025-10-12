# Proposed LLM Workflow - Executive Summary

## The Problem

Current LLM integration is **messy**:
- Output has 7+ JSON keys (entry_original, llm_adjusted_entry, changes_applied, etc.)
- Unclear which entry is final
- Deterministic inference runs first, then LLM tries to fix it (redundant)
- No auto-detect role option

## The Solution

**Pure LLM workflow** with 4 clean steps:

1. **Resolve Identity** (PubChem) → Get name, SMILES, synonyms
2. **Classify Role** (LLM) → Auto-detect role from chemistry
3. **Assign Fields** (LLM) → Pick family + populate all required fields
4. **Verify Entry** (LLM) → Check for errors before saving

## Output Comparison

### Current (Messy)
```json
{
  "entry_original": {...},
  "llm_review": {
    "raw_response": "...",
    "analysis": {...}
  },
  "llm_adjusted_entry": {...},
  "llm_auto_upgrade": {...},
  "llm_applied_changes": {...},
  "changes_applied": {...},
  "llm_adjustment_errors": [...]
}
```
❌ Which entry do I save?

### Proposed (Clean)
```json
{
  "status": "ready_to_save",
  "workflow": {
    "step1_identity": {...},
    "step2_role": {confidence: 0.95},
    "step3_fields": {confidence: 0.92},
    "step4_verification": {approved: true}
  },
  "entry": {
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
```
✅ Clear: Save `entry`, review `workflow` if needed

## Key Benefits

1. ✅ **Auto-detect role** - User can select "Auto (LLM)" instead of guessing
2. ✅ **Clean output** - Single entry + optional workflow audit trail
3. ✅ **Better accuracy** - LLM builds from scratch, not fixing deterministic errors
4. ✅ **Full validation** - Step 4 catches mistakes (wrong oxidation states, etc.)
5. ✅ **Transparent** - Confidence score at each step
6. ✅ **Matches DB format** - Direct mapping to existing registry JSON

## UI Changes

**Before**:
```
CAS: [____]
Role: [Base ▼]  ← User must know
LLM: [x] Use LLM
```

**After**:
```
Mode: [●] LLM-powered  [ ] Deterministic
CAS: [____]
Role: [Auto-detect (LLM) ▼]  ← LLM picks!
      └ Or: Base, Ligand, Metal precursor...
```

## Example: Triethylamine

**Input**: CAS `121-44-8`, Role "Auto-detect"

**Step 1**: Resolve → Name: "Triethylamine", SMILES: `CCN(CC)CC`

**Step 2**: LLM classifies → Role: "base" (95% confidence)

**Step 3**: LLM assigns → Family: "tertiary_amines_aliphatic"
- basicity: "strong"
- nucleophilicity: "moderate"  
- sterics: "unhindered"

**Step 4**: LLM verifies → Approved ✓ (no issues)

**Output**: Clean entry ready to save!

## Cost Trade-off

- **Current**: 1 LLM call, ~$0.002/entry
- **Proposed**: 3 LLM calls, ~$0.003/entry (+$0.001)

**Worth it?** YES - Much better UX and accuracy for pennies.

## Implementation

**Effort**: 2-3 days
- Create 3 new prompts (role classification, field assignment, verification)
- Create `generate_taxonomy_entry_llm()` function
- Update UI with mode toggle
- Test with 5 examples

## Decision

**Recommend**: Implement proposed workflow

**Why**: Current mixed approach is confusing. Pure LLM workflow is cleaner, more accurate, and enables auto-detection.

---

## Files to Review

1. **PROPOSED_LLM_WORKFLOW.md** - Detailed specification
2. **LLM_WORKFLOW_VISUAL.md** - Visual diagrams and examples

## Next Steps

If approved:
1. Create prompt templates
2. Implement LLM classifier functions
3. Update UI
4. Test and iterate

Ready to proceed? 🚀
