# Taxonomy Simplification Summary

## Changes Made (January 21, 2026)

### 1. Simplified Group IDs

**Before:**

- Scaffolds: `Ar`, `Alkyl`, `RCH2`, etc. (no prefix)
- Substituents: `Cl`, `Br`, `NH2`, `OH`, etc. (no prefix)
- Special variants: `Ar_Subst`, `AromN_Subst` (ugly suffixes)

**After:**

- Scaffolds: `Ar`, `Alkyl`, `RCH2`, etc. (unchanged)
- Substituents: `-Cl`, `-Br`, `-NH2`, `-OH`, etc. (all have `-` prefix)
- Clean variants: `-Ar`, `-AromN`, `-Alkenyl`, `-Alkyl` (no `_Subst` suffix)

### 2. Removed "name" Field

**Before:** Each group had both `id` and `name` fields (redundant)

```json
{
  "id": "Cl",
  "name": "-Cl",
  "kind": "substituent",
  ...
}
```

**After:** Only `id` field needed (serves as display name)

```json
{
  "id": "-Cl",
  "kind": "substituent",
  ...
}
```

### 3. Simplified Compound ID Generation

**Before:** Complex logic with `-` insertion and `_Subst` stripping

```python
display_a = group_a.replace("_Subst", "")
display_b = group_b.replace("_Subst", "")
entry["id"] = f"{display_a}-{display_b}"
```

**After:** Simple concatenation (dash already in substituent ID)

```python
entry["id"] = f"{group_a}{group_b}"  # Just A+B!
```

### 4. Example Compound Structures

**Biaryl (your original question):**

```json
{
  "A": "Ar",
  "B": "-Ar",
  "description": "Biaryl (Ar-Ar) bond."
}
```

→ Auto-generated ID: **`Ar-Ar`** ✓

**Aryl chloride:**

```json
{
  "A": "Ar",
  "B": "-Cl"
}
```

→ Auto-generated ID: **`Ar-Cl`** ✓

**Alkyl alcohol:**

```json
{
  "A": "Alkyl",
  "B": "-OH"
}
```

→ Auto-generated ID: **`Alkyl-OH`** ✓

## Benefits

1. **Cleaner IDs**: No more `_Subst` suffixes
2. **Less redundancy**: Removed `name` field from 92 groups, removed `id` field from 364 compounds
3. **Simpler logic**: ID generation is now just `A+B` concatenation
4. **Self-documenting**: Substituent IDs start with `-`, making them instantly recognizable
5. **Maintainable**: Fewer fields = fewer opportunities for inconsistency

## Statistics

- **Groups**: 92 total (14 scaffolds, 78 substituents)
- **Renamed**: 78 substituent IDs (added `-` prefix, cleaned up `_Subst` variants)
- **Fields removed**: 92 `name` fields (groups) + 364 `id` fields (compounds) = **456 lines eliminated**
- **Compounds**: All 364 compounds work with simple A+B concatenation

## Validation Results

```
✓ All 364 compounds load correctly
✓ All 92 groups validated
✓ 85/92 groups actively used in compounds
✓ Zero errors, zero warnings
✓ Featurizer builds 364 compiled patterns successfully
```

## Key Design Principles

1. **Scaffolds have no prefix** (e.g., `Ar`, `Alkyl`) - they're the "left side"
2. **Substituents have `-` prefix** (e.g., `-Cl`, `-OH`) - they're the "right side"
3. **ID is the display name** - no separate `name` field needed
4. **Compound IDs are A+B** - the dash comes naturally from the substituent prefix
5. **No special cases** - consistent rules apply to all groups

## Migration Notes

If you have external code referencing group IDs:

- Update substituent references: `"Cl"` → `"-Cl"`, `"Ar_Subst"` → `"-Ar"`
- Scaffold IDs unchanged: `"Ar"`, `"Alkyl"` stay the same
- Compound IDs remain the same (e.g., `"Ar-Cl"` is still `"Ar-Cl"`)
