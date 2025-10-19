# Role Override Enhancement - Dual Role Support

**Feature**: Keep both user-selected role AND LLM-suggested role  
**Date**: October 19, 2025  
**Status**: ✅ IMPLEMENTED

## Overview

When a user manually selects a role (e.g., "Metal precursor"), the system now:
1. ✅ Uses the user's selection as the **PRIMARY role**
2. ✅ Asks LLM for its suggestion in the background
3. ✅ If LLM suggests a different role, adds it as a **SECONDARY role**
4. ✅ Provides transparency: user sees both their choice AND the LLM's opinion

## Example Output

### User Selects "Metal precursor" for Copper

```json
{
  "name": "Copper",
  "cas": "7440-50-8",
  "roles": {
    "metal_precursor": {
      "metal": "Cu",
      "oxidation_states": [0, 1, 2]
    },
    "reductant": {
      "families": ["metal_powders"],
      "strength_band": "moderate",
      "_note": "LLM-suggested alternative role"
    }
  }
}
```

**Key points:**
- `metal_precursor` is listed FIRST (user's choice, primary role)
- `reductant` is listed SECOND (LLM suggestion, alternative role)
- `_note` field identifies it as LLM-suggested
- User can see both perspectives

## Benefits

### 1. User Control + AI Insight 🎯
- User's expertise is respected (primary role)
- AI provides additional perspective (secondary role)
- Best of both worlds!

### 2. Transparency 🔍
- User sees if their selection differs from LLM opinion
- Builds trust in the system
- Helps users learn when they might reconsider

### 3. Flexibility 📋
- Entry supports multiple roles (schema-compliant)
- Downstream systems can use primary, secondary, or both
- No information loss

### 4. Error Detection 🚨
- If roles differ significantly, may indicate:
  - User knows better (domain expertise)
  - LLM misunderstood (needs better prompts)
  - Reagent truly has dual roles (both valid)

## Workflow Details

### When User Selects a Role

```
[1] User selects "metal_precursor"
      ↓
[2] System uses "metal_precursor" as PRIMARY
      ↓
[3] System asks LLM: "What role do you think this is?"
      ↓
[4] LLM responds: "reductant"
      ↓
[5] System creates entry with BOTH roles:
    - metal_precursor (user's choice) - PRIMARY
    - reductant (LLM suggestion) - SECONDARY
      ↓
[6] User sees full picture in output
```

### When LLM Agrees with User

```
[1] User selects "catalyst"
      ↓
[2] System uses "catalyst" as PRIMARY
      ↓
[3] System asks LLM: "What role do you think this is?"
      ↓
[4] LLM responds: "catalyst" (same!)
      ↓
[5] System creates entry with ONLY ONE role:
    - catalyst (both user and LLM agree)
      ↓
[6] No secondary role added (no conflict)
```

### When User Uses Auto-detect

```
[1] User selects "🤖 Auto-detect (LLM)"
      ↓
[2] System asks LLM: "What role is this?"
      ↓
[3] LLM responds: "reductant"
      ↓
[4] System creates entry with ONE role:
    - reductant (LLM-determined)
      ↓
[5] No secondary role (nothing to compare)
```

## Implementation Details

### Changes to `generate_taxonomy_entry_llm()`

**Step 2: Role Classification Enhanced**

```python
if role_override:
    # Use user's role as PRIMARY
    role = role_override
    
    # ALSO get LLM suggestion (background)
    llm_result = classify_role(resolved_identity, llm_client)
    llm_suggested_role = llm_result.get("role")
    
    # Store both in workflow
    workflow["step2_role"] = {
        "role": role,
        "source": "user_override",
        "llm_suggestion": llm_suggested_role,
        "llm_confidence": llm_result.get("confidence"),
    }
```

**Step 3: Build Entry with Dual Roles**

```python
# Always add user's role first (PRIMARY)
entry["roles"] = {
    role: {
        "families": [fields_result["family"]],
        **fields_result.get("fields", {}),
    }
}

# If LLM suggested something different, add as SECONDARY
if llm_suggested_role and llm_suggested_role != role:
    llm_fields = assign_fields(resolved_identity, llm_suggested_role, ...)
    entry["roles"][llm_suggested_role] = {
        "families": [llm_fields["family"]],
        **llm_fields.get("fields", {}),
        "_note": "LLM-suggested alternative role",
    }
```

## Use Cases

### Case 1: User Knows Better (Domain Expert)
**Scenario**: Chemist knows Copper should be "metal_precursor" for their C-N coupling reactions

**Action**: Select "Metal precursor"

**Result**:
```json
{
  "roles": {
    "metal_precursor": { /* user's expert choice */ },
    "reductant": { /* LLM thought this, but user knows better */ }
  }
}
```

**Value**: User's expertise preserved, LLM provides alternative view

### Case 2: Reagent Has Dual Roles
**Scenario**: Some reagents can play multiple roles depending on reaction

**Action**: Select primary role, see LLM's alternative

**Result**: Both roles captured, reflects chemical reality

**Value**: No information loss, flexible for different uses

### Case 3: Learning Opportunity
**Scenario**: User unsure, wants to see what LLM thinks

**Action**: Select best guess, compare with LLM

**Result**: 
- If LLM agrees → confidence boost
- If LLM differs → learn something new

**Value**: Educational, builds understanding

### Case 4: Quality Control
**Scenario**: Systematic review of reagent database

**Action**: Review entries where user role ≠ LLM suggestion

**Result**: Identify misclassifications, improve prompts

**Value**: Continuous improvement

## Configuration

No configuration needed! This is automatic behavior when:
- ✅ User selects a specific role (not auto-detect)
- ✅ Pure LLM workflow is enabled
- ✅ LLM client is available

## Performance Impact

**Minimal**: 
- LLM call happens anyway for field assignment
- Role classification is fast (~1-2 seconds)
- Added latency: ~2 seconds per entry
- Benefit: Transparency and dual perspectives

**Optimization**: LLM calls could be parallelized in future

## Error Handling

If LLM suggestion fails:
- ✅ User's role is still used (no blocking)
- ✅ Error is logged in workflow
- ✅ No secondary role added
- ✅ Entry generation continues normally

## Backward Compatibility

- ✅ Existing entries unaffected
- ✅ Single-role entries still valid
- ✅ Schema supports multiple roles (always has)
- ✅ No breaking changes

## Future Enhancements

1. **Visual Indicators in UI**
   - Highlight when user role ≠ LLM suggestion
   - Color-code confidence levels
   - Add explanation tooltips

2. **Role Comparison Mode**
   - Show side-by-side comparison
   - Explain differences
   - Recommend which to use when

3. **Batch Review Tool**
   - Filter entries with role conflicts
   - Bulk accept/reject LLM suggestions
   - Export conflict reports

4. **Confidence Scoring**
   - Show LLM confidence for suggestion
   - Warn if user role has low confidence
   - Suggest review if conflict + low confidence

5. **Learning from Corrections**
   - Track when users override LLM
   - Improve LLM prompts based on patterns
   - Fine-tune role classification model

## Testing

Test with Copper (CAS: 7440-50-8):

1. Select "Metal precursor" → Should show both roles
2. Select "Reductant" → Might show only one (if LLM agrees)
3. Select "Auto-detect" → Should show only LLM's choice

Expected behavior verified ✅

## Summary

This enhancement provides:
- ✅ **Transparency**: User sees both their choice and LLM opinion
- ✅ **Control**: User's selection is always primary
- ✅ **Insight**: LLM provides alternative perspective
- ✅ **Flexibility**: Entry supports multiple valid roles
- ✅ **Trust**: System shows its reasoning, not just results

**Result**: Better decisions through human expertise + AI assistance! 🎯
