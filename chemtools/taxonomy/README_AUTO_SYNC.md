# Organic Taxonomy Auto-Synchronization System

This directory contains tools to automatically validate and synchronize the organic taxonomy system, ensuring that `organic_compounds.v1.3.json` stays consistent with `organic_groups.v1.3.json`.

## Why Auto-Sync?

When you update organic_groups (add new groups, change IDs, modify SMARTS), the compound definitions need to stay in sync. These tools automate:

1. **Validation** - Check that all compound references (A, B) exist in groups
2. **Detection** - Find new groups that aren't used in any compounds  
3. **Suggestion** - Auto-generate compound definitions for new groups
4. **Consistency** - Ensure naming follows the simplified ID-based system

## Tools

### 1. validate_and_sync.py

Validates the taxonomy and optionally fixes naming issues.

**Usage:**
```bash
# Validate only (show warnings/errors)
python chemtools/taxonomy/validate_and_sync.py

# Auto-fix naming mismatches (name != id)
python chemtools/taxonomy/validate_and_sync.py --fix

# Check mode for CI/CD (exit code 1 if issues found)
python chemtools/taxonomy/validate_and_sync.py --check-only
```

**What it checks:**
- ✅ All A/B group references exist in organic_groups
- ✅ Compound names match their IDs (simplified naming)
- ✅ Compound IDs follow A-B pattern
- ✅ Version dependencies are correct
- ✅ Which groups are used/unused

**Example output:**
```
✓ Loaded 92 organic groups
✓ Loaded 364 organic compounds
✓ All group references are valid
✓ All compound names match their IDs
⚠ 9 unused groups: SCN, NCS, SO2H, SO3H, ...
✓ Validation PASSED
```

### 2. suggest_compounds.py

Suggests new compound definitions for unused or new groups.

**Usage:**
```bash
# Show suggestions for all unused groups
python chemtools/taxonomy/suggest_compounds.py

# Generate JSON for suggestions
python chemtools/taxonomy/suggest_compounds.py --generate

# Suggest compounds for specific scaffold
python chemtools/taxonomy/suggest_compounds.py --scaffold Ar

# Suggest compounds for specific substituent  
python chemtools/taxonomy/suggest_compounds.py --substituent Cl

# Combined suggestions
python chemtools/taxonomy/suggest_compounds.py --scaffold Bn --substituent NCS
```

**Example output:**
```
✓ Loaded 92 groups (14 scaffolds, 78 substituents)
✓ Loaded 364 existing compounds

Unused Substituents (9):
  - SCN     (-SCN)     Thiocyanate core
  - NCS     (-NCS)     Isothiocyanate core
  - SO2H    (-SO2H)    Sulfinic acid core
  ...

Suggested Compounds (123 total):
1. Ar-SCN:
   Template: single_bond
   A: Ar (scaffold)
   B: SCN (substituent)
   Description: Ar- paired with -SCN substituent.

2. Ar-NCS:
   Template: single_bond
   A: Ar (scaffold)
   B: NCS (substituent)
   Description: Ar- paired with -NCS substituent.
...
```

### 3. fix_organometallic_refs.py

One-time fix for organometallic group references (Sn, Zn, Mg, Si → Sn*, Zn*, Mg*, Si*).

## Workflow: Adding New Groups

When you add new groups to `organic_groups.v1.3.json`:

### Step 1: Add groups
Edit `data/organic_groups.v1.3.json` and add your new groups:

```json
{
  "id": "MyNewGroup",
  "name": "-MyNewGroup",
  "kind": "substituent",
  "priority": 2,
  "smarts": "[...]",
  "description": "My new functional group",
  "tags": ["my_tag"]
}
```

### Step 2: Validate
Run validation to see what breaks:

```bash
python chemtools/taxonomy/validate_and_sync.py
```

If you see errors like:
```
✗ Compound 'Ar-MyNewGroup': Group B='MyNewGroup' not found
```

That means compounds reference the old group ID. Fix them manually or with a script.

### Step 3: Generate compound suggestions
See what compounds should be added for the new group:

```bash
python chemtools/taxonomy/suggest_compounds.py --substituent MyNewGroup
```

Or for all unused groups:

```bash
python chemtools/taxonomy/suggest_compounds.py
```

### Step 4: Add compounds
Copy the suggested JSON and add to `organic_compounds.v1.3.json`:

```json
{
  "id": "Ar-MyNewGroup",
  "name": "Ar-MyNewGroup",
  "template": "single_bond",
  "A": "Ar",
  "B": "MyNewGroup",
  "description": "Aryl MyNewGroup compound"
}
```

### Step 5: Re-validate
Confirm everything is consistent:

```bash
python chemtools/taxonomy/validate_and_sync.py
```

Should see:
```
✓ Validation PASSED - Taxonomy is consistent!
```

## CI/CD Integration

Add to your CI pipeline to prevent inconsistencies:

```yaml
# .github/workflows/taxonomy-validation.yml
- name: Validate Taxonomy
  run: python chemtools/taxonomy/validate_and_sync.py --check-only
```

This will fail the build if:
- Compounds reference non-existent groups
- Naming conventions are violated
- Dependencies are incorrect

## Naming Convention

After simplification, the system follows:
- **Compound ID** = `{A}-{B}` using exact group IDs
- **Compound name** = Same as ID (simplified)
- **Old names** = Preserved in `aliases` array

Examples:
- `Ar-Ar` (was "Biaryl", now in aliases)
- `Ar-Cl` (was already consistent)
- `Ar-OR` (was "Aryl-Ether", now in aliases)

## Files

```
chemtools/taxonomy/
├── data/
│   ├── organic_groups.v1.3.json      # Group definitions (80+ groups)
│   └── organic_compounds.v1.3.json   # Compound definitions (360+ compounds)
├── validate_and_sync.py              # Main validation tool
├── suggest_compounds.py              # Compound suggestion engine
├── fix_organometallic_refs.py        # One-time fix script
└── README_AUTO_SYNC.md               # This file
```

## Benefits

✅ **Automatic detection** of inconsistencies when groups change  
✅ **Reduced manual work** - auto-suggest compounds for new groups  
✅ **Prevents errors** - catch issues before they break the system  
✅ **Maintainability** - clear workflow for taxonomy updates  
✅ **CI/CD ready** - can be integrated into build pipelines  

## Example Workflow

```bash
# 1. Add new group "SCN" to organic_groups.v1.3.json
vim data/organic_groups.v1.3.json

# 2. Validate (will show SCN is unused)
python validate_and_sync.py

# 3. Get suggestions for SCN
python suggest_compounds.py --substituent SCN

# 4. Add suggested compounds to organic_compounds.v1.3.json  
vim data/organic_compounds.v1.3.json

# 5. Re-validate (should pass)
python validate_and_sync.py

# ✓ Done! Taxonomy is consistent and up-to-date
```

## Notes

- **Substituent variants**: Groups like `Ar_Subst` reference the base `Ar` group - validation handles this automatically
- **Organometallics**: Groups with asterisks (Sn*, Zn*, Mg*, Si*) should match in compound definitions
- **Templates**: Most compounds use `single_bond`, but some use `via_oxygen` for sulfonate esters
- **Priority**: Higher priority groups match first in detection (defined in organic_groups)
