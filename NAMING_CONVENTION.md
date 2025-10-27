# Naming Convention: Dataset ↔ Taxonomy Mapping

## Overview
This document defines the unified naming convention for reaction types across the **dataset files** and **taxonomy system** to ensure consistent mapping and future scalability.

---

## Current State Analysis

### Dataset Labels (from `data/reaction_dataset/*.jsonl`)
- **Format**: PascalCase or space-separated words
- **Examples**: 
  - `"Suzuki"` (in Suzuki.jsonl)
  - `"Amide formation"` (in Amide_formation.jsonl)
  - `"original_dataset"` (placeholder in C_N_Coupling.jsonl, C_O_Coupling.jsonl, C_S_Coupling.jsonl)

### Taxonomy IDs (from `chemtools/taxonomy/data/reaction_types.json`)
- **Format**: `snake_case`
- **Examples**: 
  - `"suzuki_miyaura"`
  - `"buchwald_hartwig_c_n"`
  - `"ullmann_cn"`
  - `"snar_cn"` (newly created)

### Taxonomy Source IDs (mapping datasets → taxonomy)
- **Format**: PascalCase with hyphens (kebab-PascalCase)
- **Examples**:
  - `"Suzuki-Miyaura"`
  - `"Buchwald-Hartwig-C-N"`
  - `"Ullmann-CN"`
  - `"SNAr-CN"` (newly created)

---

## Unified Naming Convention (RECOMMENDED)

### 1. **Dataset `reaction_type` field** (what goes in JSONL files)

**Format**: **PascalCase with hyphens** for multi-word names

**Rules**:
- Use the **most common literature name** for the reaction
- For bond-specific variants, append the bond type: `-C-N`, `-C-O`, `-C-S`
- Match the `source_ids` field in taxonomy exactly

**Examples**:
```json
// C-C Coupling datasets
"reaction_type": "Suzuki-Miyaura"
"reaction_type": "Negishi"
"reaction_type": "Heck"
"reaction_type": "Sonogashira"

// C-N Coupling datasets  
"reaction_type": "Buchwald-Hartwig-C-N"
"reaction_type": "Ullmann-CN"
"reaction_type": "Chan-Lam"
"reaction_type": "SNAr-CN"

// C-O Coupling datasets
"reaction_type": "Ullmann-Ether"
"reaction_type": "SNAr-CO"

// C-S Coupling datasets
"reaction_type": "CS-Coupling"
"reaction_type": "SNAr-CS"

// Amide/Ester formation
"reaction_type": "Amide-Formation"
"reaction_type": "Steglich-Esterification"
```

---

### 2. **Taxonomy `id` field** (internal identifier)

**Format**: **snake_case**

**Rules**:
- Descriptive, lowercase with underscores
- Must be unique across the entire taxonomy
- Should match the common name but in snake_case
- For variants, append bond type: `_c_n`, `_c_o`, `_c_s`

**Examples**:
```json
"id": "suzuki_miyaura"
"id": "buchwald_hartwig_c_n"
"id": "ullmann_cn"
"id": "snar_cn"
"id": "snar_co"
"id": "snar_cs"
"id": "ullmann_ether"
"id": "cs_coupling"
```

---

### 3. **Taxonomy `source_ids` field** (dataset mapping)

**Format**: **Array of strings matching dataset labels**

**Rules**:
- Must match exactly what appears in dataset `reaction_type` fields
- One source_id = one dataset label pattern
- For new synthetic categories (like split SNAr), use descriptive hyphenated names
- This is the **bridge** between datasets and taxonomy

**Examples**:
```json
"source_ids": ["Suzuki-Miyaura"]
"source_ids": ["Buchwald-Hartwig-C-N", "Buchwald"]  // Can have aliases
"source_ids": ["SNAr-CN"]
"source_ids": ["SNAr-CO"]
"source_ids": ["SNAr-CS"]
```

---

### 4. **Taxonomy `name` field** (human-readable display)

**Format**: **Proper Case with full descriptive text**

**Rules**:
- Full, human-readable name for UI display
- Can include chemical notation and bond types in parentheses
- More verbose than source_ids if needed

**Examples**:
```json
"name": "Suzuki-Miyaura Coupling"
"name": "Buchwald-Hartwig C-N Amination"
"name": "Nucleophilic Aromatic Substitution (C-N)"
"name": "Nucleophilic Aromatic Substitution (C-O)"
"name": "Nucleophilic Aromatic Substitution (C-S)"
```

---

## Migration Plan for SNAr

### Current State
- No standalone SNAr dataset yet
- Three SNAr variants created in taxonomy: `snar_cn`, `snar_co`, `snar_cs`

### When Creating SNAr Dataset(s)

#### Option A: **Separate dataset files** (RECOMMENDED)
```
data/reaction_dataset/SNAr-CN.jsonl     → "reaction_type": "SNAr-CN"
data/reaction_dataset/SNAr-CO.jsonl     → "reaction_type": "SNAr-CO"
data/reaction_dataset/SNAr-CS.jsonl     → "reaction_type": "SNAr-CS"
```

**Pros**: 
- Clean separation by bond type
- Easy to manage and extend
- Matches existing pattern (separate files for C_N, C_O, C_S coupling)

#### Option B: **Single dataset file with variants**
```
data/reaction_dataset/SNAr.jsonl
```
Each reaction labeled based on actual nucleophile:
```json
{"reaction_type": "SNAr-CN", ...}  // When nucleophile is amine
{"reaction_type": "SNAr-CO", ...}  // When nucleophile is alcohol
{"reaction_type": "SNAr-CS", ...}  // When nucleophile is thiol
```

**Pros**:
- All SNAr reactions in one place
- Easier initial collection

---

## Consistency Checks

### Files to Update for Full Compliance

1. **`data/reaction_dataset/Suzuki.jsonl`**
   - ❌ Current: `"reaction_type": "Suzuki"`
   - ✅ Should be: `"reaction_type": "Suzuki-Miyaura"`
   - **Action**: Update all records OR add `"Suzuki"` to `source_ids` in taxonomy

2. **`data/reaction_dataset/Amide_formation.jsonl`**
   - ❌ Current: `"reaction_type": "Amide formation"` (with space)
   - ✅ Should be: `"reaction_type": "Amide-Formation"`
   - **Action**: Update all records to use hyphen

3. **`data/reaction_dataset/C_N_Coupling.jsonl`** (and C_O, C_S)
   - ❌ Current: `"reaction_type": "original_dataset"` (placeholder)
   - ✅ Should be: Actual reaction types like `"Buchwald-Hartwig-C-N"`, `"Ullmann-CN"`, etc.
   - **Action**: Re-process with proper reaction type detection

---

## Implementation Guidelines

### For New Datasets
1. Use the naming convention from the start
2. Label `reaction_type` with PascalCase-Hyphenated format
3. Ensure it matches a `source_ids` entry in taxonomy

### For New Taxonomy Entries
1. Choose descriptive `snake_case` for `id`
2. Set `source_ids` to match expected dataset labels
3. Add bond-type suffix for variants: `-CN`, `-CO`, `-CS`
4. Write clear `name` for UI display

### For Mapping Dataset → Taxonomy
```python
# Pseudo-code for matching
def find_taxonomy_entry(dataset_reaction_type: str, taxonomy: List[ReactionType]):
    for rxn in taxonomy:
        if dataset_reaction_type in rxn.source_ids:
            return rxn
    return None
```

---

## Summary Table

| Component            | Format                     | Example                                    | Purpose                        |
|----------------------|----------------------------|--------------------------------------------|--------------------------------|
| Dataset `reaction_type` | PascalCase-Hyphenated   | `"SNAr-CN"`                                | Label in JSONL files           |
| Taxonomy `id`        | snake_case                 | `"snar_cn"`                                | Internal unique identifier     |
| Taxonomy `source_ids`| Array of PascalCase-Hyphen | `["SNAr-CN"]`                              | Maps datasets to taxonomy      |
| Taxonomy `name`      | Proper Case (descriptive)  | `"Nucleophilic Aromatic Substitution (C-N)"` | Human-readable display name |

---

## Benefits of This Convention

1. **Consistency**: Clear mapping between datasets and taxonomy
2. **Extensibility**: Easy to add new reaction variants (e.g., SNAr-CN, SNAr-CO)
3. **Readability**: Each format serves a specific purpose
4. **Compatibility**: Hyphenated format works in filenames, JSON, and code
5. **Future-proof**: Can add aliases to `source_ids` without breaking existing data

---

## Next Steps

1. ✅ **Already Done**: Created `snar_cn`, `snar_co`, `snar_cs` in taxonomy with source_ids `["SNAr-CN"]`, `["SNAr-CO"]`, `["SNAr-CS"]`

2. **TODO**: When creating SNAr datasets, use matching reaction_type labels:
   - `"reaction_type": "SNAr-CN"` for C-N bond formation
   - `"reaction_type": "SNAr-CO"` for C-O bond formation  
   - `"reaction_type": "SNAr-CS"` for C-S bond formation

3. **Consider**: Update existing datasets to follow convention (or add backwards-compatible source_ids)
