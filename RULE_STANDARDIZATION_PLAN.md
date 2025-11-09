# Rule File Standardization Plan - FUTURE-PROOF DESIGN

## Executive Summary

**Status:** Early in project lifecycle - perfect timing for standardization  
**Current State:** 9 rule files with 39 unique condition fields  
**Goal:** Establish robust standards now to avoid technical debt as you scale to 50+ rule files

## Key Findings from Analysis

### Issues Identified

1. **Mixed Data Types** (4 issues):
   - `catalyst_loading_molpct`: string + float
   - `base_equiv`: string + float  
   - `solvent_volume_mL_per_mmol`: string + float
   - Root cause: Sometimes single values, sometimes ranges

2. **Inconsistent Presence**:
   - `catalyst`: Only 5/62 occurrences (other reactions use `pd_precatalyst`, `ru_catalyst`, etc.)
   - Different reaction families need different fields

3. **Field Proliferation**: 39 unique fields across 9 files
   - Some overlap (e.g., `notes` vs `note`, `order_of_addition` exists once)
   - Risk of synonym explosion as you scale

## Recommended Standardization Strategy

### Phase 1: Standardize Existing Fields ✅ **DO THIS NOW**

#### 1.1 Data Type Consistency

**Rule:** All numeric fields that can be ranges MUST use **string type with standardized format**

```json
{
  "catalyst_loading_molpct": "0.5-2.0",     // ✅ Range as string
  "base_equiv": "2.0-3.0",                   // ✅ Range as string
  "temperature_C": "40-80",                  // ✅ Range as string
  "time_h": "1-8",                           // ✅ Range as string
  
  // When single value, still use string for consistency
  "catalyst_loading_molpct": "1.5",         // ✅ Single value as string
  "base_equiv": "2.0"                        // ✅ Single value as string
}
```

**Rationale:**
- Simpler parsing: Always expect string, check for "-" separator
- Forward compatible: Can add more complex expressions later (e.g., "1.0-2.0 (optimal: 1.5)")
- Avoids JSON schema validation headaches

#### 1.2 Field Naming Conventions

**Rule:** Use descriptive snake_case with units in suffix

```json
{
  // ✅ STANDARDIZED - Use these
  "catalyst_loading_molpct": "1.0-3.0",     // mol% loading
  "ligand_loading_molpct": "2.0-5.0",       // mol% loading
  "base_equiv": "2.0-3.0",                  // equivalents
  "temperature_C": "60-100",                // degrees Celsius
  "time_h": "4-16",                         // hours
  "pressure_bar": "1-3",                    // bar
  "concentration_M": "0.2-0.5",             // molarity
  
  // ❌ AVOID - Don't use these variants
  "catalyst_mol_pct": "...",                // Missing "loading"
  "base_equivalents": "...",                // Too verbose
  "temp_C": "...",                          // Abbreviated
  "reaction_time_h": "...",                 // Redundant "reaction"
}
```

#### 1.3 Option Formatting

**Rule:** Multiple options separated by " or " (space-or-space)

```json
{
  "catalyst": "PdCl2(PPh3)2 or Pd(dppf)Cl2·DCM",  // ✅ Correct
  "solvent": "THF, toluene, or DMF",              // ✅ Correct (Oxford comma OK)
  
  "catalyst": "PdCl2(PPh3)2/Pd(dppf)Cl2·DCM",     // ❌ Wrong separator
  "solvent": "THF or toluene or DMF"              // ⚠️  Acceptable but prefer comma for 3+
}
```

#### 1.4 Array vs String for Lists

**Rule:** Use arrays for truly independent items; use comma-separated strings for alternative choices

```json
{
  // ✅ Array for independent additives (all can be used together)
  "additives": [
    "CuI (2-5 mol%)",
    "3 Å MS optional",
    "H2O (5-10 equiv)"
  ],
  
  // ✅ String for alternatives (pick one)
  "catalyst": "PdCl2(PPh3)2 or Pd(dppf)Cl2·DCM",
  
  // ✅ String for optional modifiers
  "atmosphere": "N2 or Ar; thoroughly degassed"
}
```

### Phase 2: Define Core Schema with Required/Optional Fields

Create `data/rule_db_v2/SCHEMA_CONDITIONS.md`:

```markdown
# Conditions Schema for Rule Files

## Core Fields (Present in most rules)

### Required
- `solvent`: String, options separated by " or "
- `temperature_C`: String, range "min-max" or single "value"
- `time_h`: String, range or single value

### Highly Recommended
- `base`: String, options (if reaction needs base)
- `catalyst` OR family-specific catalyst field (e.g., `pd_precatalyst`, `ru_catalyst`)
- `catalyst_loading_molpct`: String, range or single
- `atmosphere`: String (e.g., "N2 or Ar")

### Optional but Standardized
- `ligand`: String, options
- `ligand_loading_molpct`: String
- `base_equiv`: String
- `additives`: Array of strings
- `notes`: String (general guidance)

## Family-Specific Fields

### Palladium-Catalyzed
- `pd_precatalyst`: String (preferred over generic `catalyst`)
- `pd_source`: String
- `pd_loading_molpct`: String

### Copper-Catalyzed
- `cu_source`: String
- `cu_loading_molpct`: String

### Ruthenium-Catalyzed (Metathesis)
- `ru_catalyst`: String
- `substrate_concentration_mM`: String
- `ethylene_removal`: String

### Amide Formation
- `coupling_system`: String
- `stoichiometry_equiv`: String
- `order_of_addition`: String (optional)

### Reductive Amination
- `reducing_agent`: String
- `acid_or_buffer`: String

## Field Naming Rules

1. **Units in suffix**: `_molpct`, `_equiv`, `_C`, `_h`, `_mM`, `_M`, `_bar`
2. **Full words**: `catalyst_loading` not `cat_load`
3. **Consistent order**: `{component}_{property}_{unit}`
   - ✅ `catalyst_loading_molpct`
   - ✅ `base_equiv`
   - ✅ `temperature_C`

## Value Format Rules

1. **Ranges**: "min-max" with hyphen, no spaces
   - ✅ `"0.5-2.0"`
   - ❌ `"0.5 - 2.0"`, `"0.5–2.0"` (en-dash)

2. **Single values**: String for consistency
   - ✅ `"1.5"`
   - ⚠️  `1.5` (numeric acceptable but string preferred)

3. **Options**: " or " separator
   - ✅ `"THF or toluene"`
   - ✅ `"THF, toluene, or DMF"` (Oxford comma for 3+)

4. **Arrays**: For independent items only
   - ✅ `["CuI (2-5 mol%)", "3 Å MS optional"]`
```

### Phase 3: Create Validation & Migration Tools

#### 3.1 Schema Validator Enhancement

Update `chemtools/schema/validator.py` to check:

```python
# Add to RuleSourceValidator class

STANDARD_NUMERIC_FIELDS = {
    "catalyst_loading_molpct",
    "ligand_loading_molpct", 
    "pd_loading_molpct",
    "cu_loading_molpct",
    "base_equiv",
    "temperature_C",
    "time_h",
    "pressure_bar",
    "concentration_M",
    "substrate_concentration_mM",
}

def validate_field_types(self, conditions: Dict[str, Any], path: str):
    """Validate that numeric fields are consistently typed as strings."""
    for field in STANDARD_NUMERIC_FIELDS:
        if field in conditions:
            value = conditions[field]
            if not isinstance(value, str):
                self.warnings.append({
                    "path": f"{path}.{field}",
                    "message": f"Numeric field '{field}' should be string type, got {type(value).__name__}",
                    "suggestion": f"Change to string: \"{value}\""
                })

def validate_range_format(self, value: str, field: str, path: str):
    """Validate range format (e.g., '0.5-2.0')."""
    if "-" in value:
        parts = value.split("-")
        if len(parts) != 2:
            self.errors.append({
                "path": f"{path}.{field}",
                "message": f"Invalid range format: '{value}'",
                "expected": "Use 'min-max' format (e.g., '0.5-2.0')"
            })
        else:
            # Try to parse as numbers
            try:
                low = float(parts[0].strip())
                high = float(parts[1].strip())
                if low >= high:
                    self.warnings.append({
                        "path": f"{path}.{field}",
                        "message": f"Range min >= max: {low} >= {high}"
                    })
            except ValueError:
                self.errors.append({
                    "path": f"{path}.{field}",
                    "message": f"Non-numeric range: '{value}'"
                })
```

#### 3.2 Auto-Standardization Script

Create `scripts/standardize_rules.py`:

```python
"""
Auto-standardize rule files to conform to conventions.

Actions:
- Convert numeric fields to string type
- Standardize range format (remove spaces, use hyphen)
- Fix field naming inconsistencies
- Add missing required fields with placeholders
"""

import json
from pathlib import Path
from typing import Any, Dict

def standardize_numeric_field(value: Any) -> str:
    """Convert numeric field to standard string format."""
    if isinstance(value, str):
        # Already string, just clean up
        return value.replace(" - ", "-").replace("–", "-")  # en-dash to hyphen
    elif isinstance(value, (int, float)):
        return str(value)
    return str(value)

def standardize_conditions(conditions: Dict[str, Any]) -> Dict[str, Any]:
    """Standardize a conditions dictionary."""
    result = {}
    
    NUMERIC_FIELDS = {
        "catalyst_loading_molpct",
        "ligand_loading_molpct",
        "base_equiv",
        "temperature_C",
        "time_h",
        # ... add all standard numeric fields
    }
    
    for key, value in conditions.items():
        if key in NUMERIC_FIELDS:
            result[key] = standardize_numeric_field(value)
        else:
            result[key] = value
    
    return result

def standardize_rule_file(file_path: Path, dry_run: bool = True):
    """Standardize a single rule file."""
    with open(file_path, encoding='utf-8') as f:
        data = json.load(f)
    
    changes = []
    
    # Standardize default_rule conditions
    if "default_rule" in data and "conditions" in data["default_rule"]:
        old = data["default_rule"]["conditions"]
        new = standardize_conditions(old)
        if old != new:
            changes.append("default_rule.conditions")
            data["default_rule"]["conditions"] = new
    
    # Standardize base_rules conditions
    for idx, rule in enumerate(data.get("base_rules", [])):
        if "conditions" in rule:
            old = rule["conditions"]
            new = standardize_conditions(old)
            if old != new:
                changes.append(f"base_rules[{idx}].conditions")
                rule["conditions"] = new
    
    if changes:
        print(f"📝 {file_path.name}: {len(changes)} changes")
        for change in changes:
            print(f"   - {change}")
        
        if not dry_run:
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            print(f"   ✅ Saved")
    else:
        print(f"✅ {file_path.name}: Already standardized")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="Actually modify files")
    args = parser.parse_args()
    
    rule_dir = Path("data/rule_db_v2")
    
    print("=" * 80)
    print("RULE FILE STANDARDIZATION")
    print("=" * 80)
    print(f"Mode: {'EXECUTE' if args.execute else 'DRY RUN'}")
    print()
    
    for rule_file in sorted(rule_dir.glob("*.json")):
        standardize_rule_file(rule_file, dry_run=not args.execute)
    
    print()
    if not args.execute:
        print("💡 Run with --execute to apply changes")
```

### Phase 4: Addition Sequence Integration

With standardized conditions, the addition sequence generator (from previous proposal) becomes much simpler:

```python
# chemtools/formatters/addition_sequence.py

ROLE_FIELD_MAPPING = {
    # Direct field mappings
    "solvent": "solvent",
    "base": "base",
    "ligand": "ligand",
    
    # Metal catalyst - family specific
    "metal_catalyst": [
        "catalyst",
        "pd_precatalyst",
        "pd_source",
        "ru_catalyst",
        "cu_source",
        "metal_catalyst"
    ],
    
    # Additives
    "additive": "additives",
}

def extract_chemicals_from_conditions(
    conditions: Dict[str, Any],
    family: Optional[str] = None
) -> List[Dict[str, Any]]:
    """Extract chemicals from standardized conditions dict."""
    chemicals = []
    
    # Solvent (order 1)
    if "solvent" in conditions:
        chemicals.append({
            "name": _pick_first_option(conditions["solvent"]),
            "role": "solvent",
            "amount": {"volume_ml": "to_be_scaled"},
            "addition_order": 1
        })
    
    # Base (order 2)
    if "base" in conditions:
        base_equiv = _parse_range_midpoint(conditions.get("base_equiv", "2.0"))
        chemicals.append({
            "name": _pick_first_option(conditions["base"]),
            "role": "base",
            "amount": {"equivalents": base_equiv},
            "addition_order": 2
        })
    
    # ... etc (simplified because fields are now standardized)
```

## Implementation Timeline

### Week 1: Standards & Validation
- [ ] Create `SCHEMA_CONDITIONS.md` documentation
- [ ] Update `validator.py` with new checks
- [ ] Run validation on all 9 files, document issues

### Week 2: Standardization
- [ ] Create `standardize_rules.py` script
- [ ] Run dry-run on all files
- [ ] Review proposed changes
- [ ] Execute standardization

### Week 3: Addition Sequence
- [ ] Implement `addition_sequence.py` with standardized field mapping
- [ ] Integrate into `UnifiedRecommender`
- [ ] Test on all 9 rule families

### Week 4: Documentation & Training
- [ ] Update `AGENTS.md` with new standards
- [ ] Create authoring guide for new rule files
- [ ] Add examples to README

## Benefits of Standardization NOW

1. **Avoid Technical Debt**
   - 9 files is manageable
   - 50+ files would be painful to retrofit

2. **Cleaner Addition Sequence Generator**
   - Predictable field names
   - Consistent types simplify parsing
   - Less error handling needed

3. **Better Validation**
   - Catch errors at authoring time
   - Automated checks in CI/CD

4. **Easier Onboarding**
   - Clear documentation
   - Consistent patterns across families

5. **LLM-Friendly**
   - More predictable for LLM-based curation
   - Easier to generate synthetic rules

## Future-Proofing Strategies

### 1. Semantic Versioning for Conditions

Add version field to track schema evolution:

```json
{
  "schema_version": "2.0",
  "conditions_schema_version": "1.0",  // NEW
  "default_rule": {
    "conditions": {...}
  }
}
```

### 2. Extensibility Points

Reserve namespace for custom fields:

```json
{
  "conditions": {
    // Standard fields
    "catalyst": "...",
    
    // Custom/experimental fields (prefix with _)
    "_experimental_activation_method": "microwave"
  }
}
```

### 3. Deprecation Pathway

When fields need to change:

```json
{
  "conditions": {
    "catalyst_loading_molpct": "1.0-3.0",
    "_deprecated_cat_loading": "1.0-3.0"  // Mark old field
  },
  "_migration_notes": [
    "cat_loading renamed to catalyst_loading_molpct in conditions_schema v1.0"
  ]
}
```

## Recommendation

**✅ START WITH PHASE 1 NOW** (Week 1-2, ~8 hours total)

1. Document standards (2 hours)
2. Run analysis script (done ✅)
3. Create standardization script (3 hours)
4. Apply to all 9 files (1 hour)
5. Validate & test (2 hours)

This establishes foundation before you scale to more rule files.

**✅ THEN DO PHASE 3-4** (Week 3-4)

Integration with addition sequence generator becomes trivial with standardized fields.

## Questions to Resolve

1. **Numeric strings vs numbers**: Prefer strings for ranges, but allow numbers for single values?
   - **Recommendation**: All strings for consistency

2. **Oxford comma**: Required or optional in option lists?
   - **Recommendation**: Optional but encouraged for 3+ items

3. **Family-specific fields**: Document in SCHEMA_CONDITIONS.md or separate family guides?
   - **Recommendation**: Both - core schema + family supplements

4. **Validation strictness**: Warnings vs errors for non-standard fields?
   - **Recommendation**: Warnings for now, errors after migration period
