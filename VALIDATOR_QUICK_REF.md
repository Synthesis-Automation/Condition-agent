# Reagent Validator - Quick Reference

## Installation
Already installed - part of `chemtools.reagent` module.

## CLI Usage

```bash
# Validate entire database
python -m chemtools.reagent.validator data/reagents

# Or use the convenience script
python scripts/validate_reagent_db.py

# Validate specific role
python scripts/validate_reagent_db.py --role ligand

# Verbose output (show all issues)
python scripts/validate_reagent_db.py --verbose

# Strict type checking
python scripts/validate_reagent_db.py --strict
```

## Python API

### Import
```python
from chemtools.reagent import (
    validate_entry,
    validate_role_file,
    validate_database,
    print_validation_summary
)
```

### Validate Single Entry
```python
entry = {
    "id": "...",
    "name": "...",
    # ... other fields
}

issues = validate_entry(entry, role="ligand")

if not issues:
    print("✅ Valid!")
else:
    for issue in issues:
        print(f"{issue['severity']}: {issue['field']} - {issue['message']}")
```

### Validate Before Save (UI Integration)
```python
def on_save_clicked(self):
    entry = self._last_result.get("entry")
    role = self._last_result.get("role")
    
    # Validate before saving
    issues = validate_entry(entry, role=role)
    if issues:
        errors = [i for i in issues if i["severity"] == "error"]
        if errors:
            self.show_error(f"Cannot save: {errors[0]['message']}")
            return
    
    # Safe to save
    store.add_entry(role, entry)
```

### Validate Entire Database
```python
results = validate_database("data/reagents")

print(f"Total: {results['total_entries']}")
print(f"Valid: {results['valid_entries']}")
print(f"Invalid: {results['invalid_entries']}")
print(f"Errors: {results['error_count']}")
print(f"Warnings: {results['warning_count']}")

# Pretty print
print_validation_summary(results, verbose=True)
```

### Validate Specific Roles
```python
# Only validate ligands and bases
results = validate_database(
    "data/reagents",
    roles=["ligand", "base"]
)
```

## Return Values

### validate_entry()
Returns list of issues:
```python
[
    {
        "severity": "error",           # or "warning"
        "field": "roles.ligand.donors",
        "message": "Missing required field 'donors'"
    }
]
```

Empty list = valid entry.

### validate_database()
Returns comprehensive dict:
```python
{
    "registry_dir": "data/reagents",
    "roles_checked": ["ligand", "base", ...],
    "total_files": 8,
    "total_entries": 383,
    "valid_entries": 379,
    "invalid_entries": 4,
    "error_count": 6,
    "warning_count": 224,
    "by_role": {
        "ligand": {
            "total_entries": 155,
            "valid_entries": 154,
            "invalid_entries": 1,
            "errors": [...],
            "warnings": [...],
            "entry_issues": {
                "index_16": [{...}]
            }
        }
    }
}
```

## Required Fields

All entries must have these 9 fields:
1. `id` (string)
2. `name` (string)
3. `abbreviation` (array or null)
4. `aliases` (array or null)
5. `cas` (string or null)
6. `inchi_key` (string or null)
7. `smiles` (string or null)
8. `roles` (object)
9. `embedding_text` (string or null)

Plus role-specific fields in `roles.<role>`:
- **Ligand**: `families`, `donors`, `denticity`
- **Base**: `families`, `basicity`, `nucleophilicity`, `sterics`
- **Acid**: `families`, `acidity`
- **Oxidant/Reductant**: `families`, `strength_band`
- **Solvent**: `families`, `proticity`, `polarity`, `coordination`
- See full list in VALIDATOR_README.md

## Common Use Cases

### 1. Pre-Commit Hook
```bash
#!/bin/bash
python -m chemtools.reagent.validator data/reagents
exit $?  # Exit with validation status
```

### 2. CI/CD Pipeline
```yaml
# .github/workflows/validate.yml
- name: Validate Reagent Database
  run: |
    python -m chemtools.reagent.validator data/reagents
```

### 3. Data Quality Check
```python
# Check database health regularly
results = validate_database("data/reagents")
if results["invalid_entries"] > 0:
    print(f"⚠️ Warning: {results['invalid_entries']} invalid entries")
    # Send alert, create issue, etc.
```

### 4. Find Specific Issues
```python
# Find all entries missing InChIKey
results = validate_database("data/reagents")
for role, role_data in results["by_role"].items():
    for entry_id, issues in role_data.get("entry_issues", {}).items():
        for issue in issues:
            if "inchi_key" in issue["field"]:
                print(f"{role}/{entry_id}: {issue['message']}")
```

## Exit Codes (CLI)
- `0` - All entries valid
- `1` - Validation errors found or directory not found

## Performance
- ~1000 entries/second
- Current database (383 entries): < 1 second

## Current Database Status
As of last validation:
- **Total entries**: 383
- **Valid**: 379 (99.0%)
- **Invalid**: 4 (1.0%)
- **Errors**: 6
- **Warnings**: 224

## Files
- **Implementation**: `chemtools/reagent/validator.py`
- **Documentation**: `chemtools/reagent/VALIDATOR_README.md`
- **CLI Script**: `scripts/validate_reagent_db.py`
- **Summary**: `VALIDATOR_IMPLEMENTATION.md`
- **This guide**: `VALIDATOR_QUICK_REF.md`
