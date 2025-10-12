# LLM Enhancement Visual Guide

## Before vs. After Comparison

### Without LLM (Original Behavior)

```
┌─────────────────────────────────────────────┐
│ User Input:                                 │
│ - CAS: 121-44-8                            │
│ - Role: other_reagent                       │
│ - LLM: Disabled                            │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│ Deterministic Processing                    │
│ - Resolve identity (name, synonyms)        │
│ - Match family tokens                       │
│ - Use default family if no match           │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│ Output:                                     │
│ {                                           │
│   "name": "Triethylamine",                 │
│   "cas": "121-44-8",                       │
│   "role": "other_reagent",                 │
│   "family_id": "misc_general",             │
│   "roles": {                               │
│     "other_reagent": {                     │
│       "families": ["misc_general"]         │
│     }                                       │
│   }                                         │
│ }                                           │
└─────────────────────────────────────────────┘

❌ Missing: Specific role classification
❌ Missing: Critical fields (basicity, etc.)
❌ Requires: Manual review and reclassification
```

### With LLM (Enhanced Behavior)

```
┌─────────────────────────────────────────────┐
│ User Input:                                 │
│ - CAS: 121-44-8                            │
│ - Role: other_reagent                       │
│ - LLM: OpenAI gpt-4o-mini                  │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│ Deterministic Processing                    │
│ - Resolve identity (name, synonyms)        │
│ - Match family tokens                       │
│ - Use default family if no match           │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│ LLM Review (NEW!)                           │
│ - Analyze chemistry (tertiary amine)       │
│ - Propose: role="base"                     │
│ - Propose: family="tertiary_amines_..."    │
│ - Suggest: basicity="strong"               │
│ - Suggest: nucleophilicity="moderate"      │
│ - Suggest: sterics="unhindered"            │
│ - Confidence: 0.95                         │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│ Auto-apply Changes (NEW!)                   │
│ - Upgrade role: other_reagent → base       │
│ - Change family: misc_general → tertiary...│
│ - Add fields: basicity, nucleophilicity... │
│ - Log changes with justification           │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│ Enhanced Output (NEW!)                      │
│ {                                           │
│   "review_summary": {                      │
│     "original_role": "other_reagent",      │
│     "llm_status": "needs_review",          │
│     "confidence": 0.95,                    │
│     "auto_upgrade": {                      │
│       "from": "other_reagent",             │
│       "to": "base",                        │
│       "confidence": 0.95                   │
│     }                                       │
│   },                                        │
│   "entry_revised": {                       │
│     "name": "Triethylamine",               │
│     "cas": "121-44-8",                     │
│     "roles": {                             │
│       "base": {                            │
│         "families": ["tertiary_amines..."],│
│         "basicity": "strong",              │
│         "nucleophilicity": "moderate",     │
│         "sterics": "unhindered"            │
│       }                                     │
│     }                                       │
│   },                                        │
│   "changes_applied": {                     │
│     "role": {"from": "other_reagent",      │
│              "to": "base"},                │
│     "field_suggestions_applied": {         │
│       "basicity": "strong",                │
│       "nucleophilicity": "moderate",       │
│       "sterics": "unhindered"              │
│     }                                       │
│   }                                         │
│ }                                           │
└─────────────────────────────────────────────┘

✅ Automatic: Role upgrade (other_reagent → base)
✅ Automatic: Field population (3 critical fields)
✅ Transparent: Full change log with confidence
✅ Review-ready: Original vs. revised side-by-side
```

## Three Enhancement Features

### 1. Auto-Upgrade "other_reagent"

```
┌──────────────────┐
│ Deterministic    │
│ assigns:         │  ┌─────────────────────────┐
│ "other_reagent"  │──▶ LLM analyzes chemistry │
└──────────────────┘  └─────────────────────────┘
                                 ↓
                      ┌─────────────────────────┐
                      │ Proposes specific role: │
                      │ - base                  │
                      │ - metal_precursor       │
                      │ - oxidant               │
                      │ - etc.                  │
                      └─────────────────────────┘
                                 ↓
                      ┌─────────────────────────┐
                      │ System auto-applies if: │
                      │ ✓ Valid role name       │
                      │ ✓ Valid family exists   │
                      │ ✓ Confidence > threshold│
                      └─────────────────────────┘
```

### 2. Field Assignment

```
┌──────────────────────────┐
│ Role-Specific Fields     │
├──────────────────────────┤
│ base:                    │    ┌─────────────────┐
│   - basicity             │◀───│ LLM suggests    │
│   - nucleophilicity      │    │ chemistry-based │
│   - sterics              │    │ values          │
├──────────────────────────┤    └─────────────────┘
│ metal_precursor:         │
│   - metal                │
│   - oxidation_states     │
├──────────────────────────┤
│ solvent:                 │
│   - proticity            │
│   - polarity             │
│   - coordination         │
└──────────────────────────┘
```

### 3. Review Output

```
┌─────────────────────────────────────────────────┐
│ Output Structure                                │
├─────────────────────────────────────────────────┤
│ review_summary:        ◀─── High-level overview │
│   - original_role                               │
│   - llm_status                                  │
│   - confidence                                  │
│   - auto_upgrade                                │
├─────────────────────────────────────────────────┤
│ entry_original:        ◀─── Before LLM changes  │
│   - Deterministic assignment                    │
├─────────────────────────────────────────────────┤
│ entry_revised:         ◀─── After LLM changes   │
│   - Enhanced with fields                        │
│   - Role/family upgraded                        │
├─────────────────────────────────────────────────┤
│ changes_applied:       ◀─── What changed        │
│   - role: {from, to}                            │
│   - family: {from, to}                          │
│   - field_suggestions_applied: {...}            │
│   - synonyms_added: [...]                       │
├─────────────────────────────────────────────────┤
│ llm_review_details:    ◀─── LLM metadata        │
│   - model                                       │
│   - tokens_used                                 │
│   - latency_ms                                  │
├─────────────────────────────────────────────────┤
│ llm_alerts:            ◀─── Warnings/notes      │
│   - ["Consider...", "Check..."]                 │
└─────────────────────────────────────────────────┘
```

## Data Flow Diagram

```
┌──────────────┐
│   User UI    │
│ (Qt Window)  │
└──────┬───────┘
       │ CAS, role, LLM settings
       ↓
┌──────────────────────────────────┐
│ generate_taxonomy_entry()        │
│ (Main workflow)                  │
├──────────────────────────────────┤
│ 1. Resolve identity (CAS→name)  │
│ 2. Infer family (token match)   │
│ 3. Build role payload            │
└──────┬───────────────────────────┘
       │ If llm_options set
       ↓
┌──────────────────────────────────┐
│ review_taxonomy_proposal()       │
│ (llmtools/reagent_review.py)     │
├──────────────────────────────────┤
│ 1. Build review prompt           │
│ 2. Call LLM API                  │
│ 3. Parse JSON response           │
│ 4. Extract analysis              │
└──────┬───────────────────────────┘
       │ Returns llm_review dict
       ↓
┌──────────────────────────────────┐
│ Apply LLM Adjustments            │
│ (Enhanced logic)                 │
├──────────────────────────────────┤
│ 1. Check auto_upgrade condition  │ ◀── NEW!
│ 2. Validate proposed role/family │
│ 3. Merge field_suggestions       │ ◀── NEW!
│ 4. Build adjusted entry          │
│ 5. Log changes summary           │
└──────┬───────────────────────────┘
       │ Returns result with llm_adjusted_entry
       ↓
┌──────────────────────────────────┐
│ on_generation_success()          │
│ (UI display)                     │
├──────────────────────────────────┤
│ 1. Format review_summary         │ ◀── NEW!
│ 2. Show entry_original           │ ◀── NEW!
│ 3. Show entry_revised            │
│ 4. Show changes_applied          │ ◀── NEW!
│ 5. Enable Save button            │
└──────┬───────────────────────────┘
       │ User reviews and clicks Save
       ↓
┌──────────────────────────────────┐
│ on_save_clicked()                │
│ (Persist to registry)            │
├──────────────────────────────────┤
│ 1. Parse JSON from editor        │
│ 2. Extract entry_revised         │
│ 3. Write to role-specific file   │
│ 4. Update status                 │
└──────────────────────────────────┘
```

## Code Changes Locations

### File: reagent_taxonomy_qt.py

```python
# Line ~740: Enhancement 1 - Auto-upgrade "other_reagent"
if role == "other_reagent" and candidate_role and candidate_role != role:
    debug_log.append(
        f"LLM auto-upgrade: 'other_reagent' → '{candidate_role}' "
        f"(confidence: {analysis.get('confidence', 'N/A')})"
    )
    result["llm_auto_upgrade"] = {
        "from": "other_reagent",
        "to": candidate_role,
        "reason": analysis.get("justification", "LLM recommendation"),
        "confidence": analysis.get("confidence"),
    }

# Line ~790: Enhancement 2 - Apply field suggestions
field_suggestions = analysis.get("field_suggestions") or {}
if isinstance(field_suggestions, dict) and field_suggestions:
    changes_summary["field_suggestions_applied"] = field_suggestions
    should_apply = True

# Merge suggestions into role payload
if isinstance(field_suggestions, dict):
    for field_name, field_value in field_suggestions.items():
        if field_name and field_value is not None:
            adjusted_role_payload[field_name] = field_value
            debug_log.append(
                f"LLM field suggestion: {field_name} = {field_value}"
            )

# Line ~1180: Enhancement 3 - Enhanced output format
display_payload["review_summary"] = {
    "original_role": result.get("role"),
    "original_family": result.get("family_id"),
    "llm_status": llm_review.get("analysis", {}).get("status"),
    "confidence": llm_review.get("analysis", {}).get("confidence"),
    "justification": llm_review.get("analysis", {}).get("justification"),
}
```

### File: llmtools/prompts.py

```python
# Line ~230: Enhanced task instructions
## Task
Evaluate the proposed assignment carefully. Pay special attention to:
1. **Role accuracy**: If current role is "other_reagent", STRONGLY recommend a more specific role if applicable
2. **Family precision**: Suggest the most specific family that matches the reagent's chemistry
3. **Field reliability**: Identify critical missing fields (e.g., metal, oxidation_states, basicity, etc.)
4. **Synonym completeness**: Add any well-known synonyms or trade names

# Line ~250: Added field_suggestions to schema
{
  "status": "confirm" | "needs_review" | "reject",
  "proposed_role": "string",
  "proposed_family": "string",
  "confidence": 0-1 float,
  "justification": "short explanation",
  "alerts": ["list of warnings or required actions"],
  "suggested_synonyms": ["list of additional synonyms (if any)"],
  "field_suggestions": {
    "field_name": "suggested_value",
    "field_name2": "suggested_value2"
  }
}
```

### File: llmtools/reagent_review.py

```python
# Line ~165: Extract field_suggestions
if parsed:
    result["analysis"] = {
        "status": parsed.get("status"),
        "proposed_role": parsed.get("proposed_role"),
        "proposed_family": parsed.get("proposed_family"),
        "confidence": parsed.get("confidence"),
        "justification": parsed.get("justification"),
        "alerts": normalized_alerts,
        "suggested_synonyms": [...],
        "field_suggestions": parsed.get("field_suggestions", {}),  # NEW!
    }
```

## Summary

All three enhancements work together to provide:

1. **Intelligent classification** - Auto-upgrade vague roles
2. **Complete entries** - Populate critical chemistry fields
3. **Transparent review** - Clear before/after with justification

Ready for production use! 🚀
