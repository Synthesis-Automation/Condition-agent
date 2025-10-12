# LLM Workflow Visual Guide

## Current Problem (Mixed Mode)

```
User Input (CAS + Role)
        ↓
┌───────────────────────────────────────┐
│  DETERMINISTIC PROCESSING             │
│  - Resolve identity                   │
│  - Token matching                     │
│  - Infer family                       │
│  - Build role payload                 │
└───────────────┬───────────────────────┘
                ↓
        [entry_original]
                ↓
┌───────────────────────────────────────┐
│  LLM REVIEW & ADJUSTMENT              │
│  - Review deterministic assignment    │
│  - Suggest corrections                │
│  - Auto-upgrade if "other_reagent"    │
└───────────────┬───────────────────────┘
                ↓
        [llm_adjusted_entry]
                ↓
┌───────────────────────────────────────┐
│  MESSY OUTPUT                         │
│  - entry_original                     │
│  - llm_review (raw_response, etc.)   │
│  - llm_adjusted_entry                 │
│  - llm_auto_upgrade                   │
│  - llm_applied_changes                │
│  - changes_applied                    │
│  - llm_adjustment_errors              │
└───────────────────────────────────────┘

❌ Confusing: Which entry is final?
❌ Redundant: Why run deterministic if LLM fixes it?
❌ Complex: Multiple adjustment logic paths
```

---

## Proposed Streamlined Workflow (LLM Mode)

```
User Input (CAS only or CAS + Role hint)
        ↓
┌──────────────────────────────────────────────┐
│ STEP 1: RESOLVE CHEMICAL IDENTITY           │
│ Source: PubChem / CAS Registry              │
│ Output: name, SMILES, synonyms, InChI key   │
└──────────────┬───────────────────────────────┘
               ↓
        {identity bundle}
               ↓
┌──────────────────────────────────────────────┐
│ STEP 2: LLM ROLE CLASSIFICATION              │
│ Trigger: User selects "Auto" or "other"      │
│ LLM Task: Pick correct role from 13 options  │
│ Output: {role, confidence, reasoning}        │
└──────────────┬───────────────────────────────┘
               ↓
        {role: "base", confidence: 0.95}
               ↓
┌──────────────────────────────────────────────┐
│ STEP 3: LLM FAMILY & FIELD ASSIGNMENT        │
│ Input: identity + role                       │
│ LLM Task: Pick family + assign all fields    │
│ Schema: Role-specific field definitions      │
│ Output: {family, fields, abbreviations}      │
└──────────────┬───────────────────────────────┘
               ↓
        {family: "tertiary_amines_aliphatic",
         fields: {basicity: "strong", ...}}
               ↓
┌──────────────────────────────────────────────┐
│ STEP 4: LLM VERIFICATION                     │
│ Input: Complete entry JSON                   │
│ LLM Task: Check for errors/inconsistencies   │
│ Output: {approved: true/false, issues: [...]}│
└──────────────┬───────────────────────────────┘
               ↓
┌──────────────────────────────────────────────┐
│ CLEAN OUTPUT                                 │
│ {                                            │
│   "status": "ready_to_save",                 │
│   "workflow": {                              │
│     "step1_identity": {...},                 │
│     "step2_role": {...},                     │
│     "step3_fields": {...},                   │
│     "step4_verification": {...}              │
│   },                                          │
│   "entry": {                                 │
│     "id": "...",                             │
│     "name": "...",                           │
│     "roles": { "base": {...} }               │
│   }                                           │
│ }                                             │
└───────────────────────────────────────────────┘

✅ Clear: Single entry + workflow audit trail
✅ Efficient: Direct LLM path, no fixing
✅ Transparent: Confidence at each step
```

---

## Step-by-Step Example: Triethylamine

### Input
```
CAS: 121-44-8
User selection: "Auto-detect role (LLM)"
```

### Step 1: Identity Resolution
```json
{
  "cas": "121-44-8",
  "name": "Triethylamine",
  "smiles": "CCN(CC)CC",
  "synonyms": ["TEA", "N,N-Diethylethanamine"],
  "inchi_key": "ZMANZCXQSJIPKH-UHFFFAOYSA-N",
  "molecular_formula": "C6H15N",
  "source": "pubchem"
}
```

### Step 2: Role Classification (LLM)

**LLM sees**:
- Name: Triethylamine
- SMILES: CCN(CC)CC (tertiary amine structure)
- Available roles: [ligand, base, metal_precursor, ...]

**LLM returns**:
```json
{
  "role": "base",
  "confidence": 0.95,
  "reasoning": "Tertiary aliphatic amine (R3N) with basic nitrogen"
}
```

### Step 3: Field Assignment (LLM)

**LLM sees**:
- Role: base
- Available families: [tertiary_amines_aliphatic, amidine_guanidine_bases, ...]
- Required fields: {basicity, nucleophilicity, sterics}
- Field options: {basicity: [weak, moderate, strong, superbase], ...}

**LLM returns**:
```json
{
  "family": "tertiary_amines_aliphatic",
  "fields": {
    "basicity": "strong",
    "nucleophilicity": "moderate",
    "sterics": "unhindered"
  },
  "abbreviations": ["TEA"],
  "additional_synonyms": [],
  "confidence": 0.92,
  "reasoning": "Unhindered tertiary amine, typical pKa ~10.7"
}
```

### Step 4: Verification (LLM)

**LLM sees**:
```json
{
  "name": "Triethylamine",
  "smiles": "CCN(CC)CC",
  "roles": {
    "base": {
      "families": ["tertiary_amines_aliphatic"],
      "basicity": "strong",
      "nucleophilicity": "moderate",
      "sterics": "unhindered"
    }
  }
}
```

**LLM checks**:
- ✅ SMILES matches tertiary amine (3 alkyl groups on N)
- ✅ "strong" basicity reasonable for R3N
- ✅ "unhindered" correct (ethyl groups are small)
- ✅ No oxidation states (not a metal)
- ✅ No contradictions

**LLM returns**:
```json
{
  "approved": true,
  "issues": [],
  "suggestions": ["Consider adding pKa value in notes"]
}
```

### Final Output
```json
{
  "status": "ready_to_save",
  "cas": "121-44-8",
  "workflow": {
    "step1_identity": {
      "source": "pubchem",
      "name": "Triethylamine",
      "smiles": "CCN(CC)CC"
    },
    "step2_role": {
      "role": "base",
      "confidence": 0.95,
      "reasoning": "Tertiary aliphatic amine"
    },
    "step3_fields": {
      "family": "tertiary_amines_aliphatic",
      "confidence": 0.92,
      "fields_assigned": ["basicity", "nucleophilicity", "sterics"]
    },
    "step4_verification": {
      "approved": true,
      "issues": [],
      "suggestions": ["Add pKa value"]
    }
  },
  "entry": {
    "id": "ZMANZCXQSJIPKH-UHFFFAOYSA-N",
    "name": "Triethylamine",
    "abbreviation": ["TEA"],
    "aliases": ["N,N-Diethylethanamine"],
    "cas": "121-44-8",
    "inchi_key": "ZMANZCXQSJIPKH-UHFFFAOYSA-N",
    "smiles": "CCN(CC)CC",
    "roles": {
      "base": {
        "families": ["tertiary_amines_aliphatic"],
        "basicity": "strong",
        "nucleophilicity": "moderate",
        "sterics": "unhindered"
      }
    },
    "embedding_text": "type: BASE | family: Tertiary amines..."
  }
}
```

**User sees**:
- ✅ Clean entry ready to save
- ✅ Workflow audit trail (optional review)
- ✅ Confidence scores at each step
- ✅ LLM suggestions

---

## UI Changes

### Current UI (Mixed Mode)
```
┌─────────────────────────────────┐
│ CAS: [121-44-8            ]     │
│ Role: [Base ▼]                  │  ← User must know role
│ LLM: [x] Use LLM                │
└─────────────────────────────────┘
                ↓
        [Generate]
                ↓
┌─────────────────────────────────┐
│ Output (messy JSON):            │
│ - entry_original                │
│ - llm_review                    │
│ - llm_adjusted_entry            │
│ - changes_applied               │
│ (User confused which to save)   │
└─────────────────────────────────┘
```

### Proposed UI (Streamlined)
```
┌─────────────────────────────────────────┐
│ Mode: [●] LLM-powered  [ ] Deterministic│
└─────────────────────────────────────────┘
┌─────────────────────────────────────────┐
│ CAS: [121-44-8                    ]     │
│ Role: [Auto-detect (LLM) ▼]             │  ← Auto-detect option!
│       └ Or: Base, Ligand, etc.          │
└─────────────────────────────────────────┘
                ↓
        [Generate with LLM]
                ↓
┌─────────────────────────────────────────┐
│ ✅ Workflow Complete                    │
│                                          │
│ Step 1: Identity ✓ (pubchem)            │
│ Step 2: Role ✓ (base, 95% confidence)   │
│ Step 3: Fields ✓ (3 assigned)           │
│ Step 4: Verified ✓ (no issues)          │
│                                          │
│ [Show Entry] [Show Workflow] [Edit]     │
└─────────────────────────────────────────┘
                ↓
        [Save to Registry]
```

**Entry Tab**:
```json
{
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
```

**Workflow Tab** (optional):
```
Step 1: Identity Resolution
  Source: pubchem
  Name: Triethylamine
  SMILES: CCN(CC)CC

Step 2: Role Classification
  Role: base
  Confidence: 95%
  Reasoning: Tertiary aliphatic amine

Step 3: Field Assignment
  Family: tertiary_amines_aliphatic
  Confidence: 92%
  Fields: basicity=strong, nucleophilicity=moderate, sterics=unhindered

Step 4: Verification
  Status: Approved ✓
  Issues: None
  Suggestions: Consider adding pKa value
```

---

## Comparison Table

| Feature | Current (Mixed) | Proposed (Pure LLM) |
|---------|----------------|---------------------|
| **Output complexity** | 7+ JSON keys | 3 keys (status, workflow, entry) |
| **User confusion** | High (which entry is final?) | Low (entry is final) |
| **LLM calls** | 1 (review only) | 3 (classify, assign, verify) |
| **Accuracy** | LLM fixes deterministic errors | LLM builds correctly from start |
| **Transparency** | Partial (changes log) | Full (workflow audit trail) |
| **Auto-detect role** | ❌ No | ✅ Yes |
| **Validation** | ⚠️ Partial | ✅ Full LLM verification |
| **Edit before save** | ⚠️ Edit messy JSON | ✅ Edit clean entry |
| **Token usage** | ~1000 tokens | ~1500 tokens (3 calls) |
| **Cost** | ~$0.002/entry | ~$0.003/entry (+50%) |

**Trade-off**: Slightly higher cost (~$0.001 more) for **much better** UX and accuracy.

---

## Summary

### Problems Solved
1. ✅ **Messy output** → Clean entry + optional workflow view
2. ✅ **No auto-detect** → LLM can classify role from SMILES
3. ✅ **Weak validation** → Full LLM verification step
4. ✅ **Confusing workflow** → Linear 4-step pipeline

### Key Advantages
1. **Single source of truth**: Entry is built by LLM, not adjusted
2. **Confidence tracking**: Each step has confidence score
3. **Better UX**: Clear status (ready_to_save vs. needs_review)
4. **Flexible**: User can override any step
5. **Auditable**: Workflow shows all LLM decisions

### Implementation Effort
- **New code**: ~400 lines (3 LLM functions + UI updates)
- **Reuse**: Identity resolver, registry store, save logic
- **Testing**: 5 test cases (one per role type)
- **Timeline**: 2-3 days

Ready to implement? 🚀
