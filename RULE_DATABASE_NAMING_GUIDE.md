# Rule Database File Naming Convention Guide

## 📋 Current Naming Pattern

Based on your codebase analysis, here's the **logical naming convention** your code expects:

### Pattern: `{ReactionFamily}_{Variant}_db.json`

---

## 📁 Expected File Names in `data/rule_db/`

### Currently Referenced by Code:

| Display Name | File Name Expected | Status |
|--------------|-------------------|---------|
| **Suzuki Coupling** | `Suzuki_db.json` | ❌ Missing (you have `suzuki.json`) |
| **C-N Coupling (Pd)** | `C_N_Coupling_Pd_db.json` | ❌ Missing (you have `buchwald_cn.json`) |
| **C-N Coupling (Cu)** | `C_N_Coupling_Cu_db.json` | ❌ Missing (you have `ullman_cn.json`) |
| **C-N Coupling (Ni)** | `C_N_Coupling_Ni_db.json` | ❌ Missing |
| **Amide Formation** | `Amide_formation_db.json` | ✅ Exists |

### Your Current Files:

| Current File Name | Should Be Renamed To | Reason |
|-------------------|---------------------|---------|
| `suzuki.json` | `Suzuki_db.json` | Match code expectations |
| `buchwald_cn.json` | `C_N_Coupling_Pd_db.json` | Buchwald-Hartwig is Pd-based C-N |
| `ullman_cn.json` | `C_N_Coupling_Cu_db.json` | Ullman is Cu-based C-N |
| `amide_formation_db.json` | ✅ Already correct | Keep as is |

---

## 🎯 Naming Rules

### 1. **Base Pattern**
```
{Family}_{Variant}_db.json
```

### 2. **Capitalization Rules**
- **First letter of each word**: UPPERCASE
- **Rest of letters**: lowercase  
- **Underscores**: Between words and before `_db.json`

### 3. **Examples**

✅ **Correct:**
```
Suzuki_db.json
C_N_Coupling_Pd_db.json
C_N_Coupling_Cu_db.json
Amide_formation_db.json
Heck_reaction_db.json
Sonogashira_db.json
```

❌ **Incorrect:**
```
suzuki.json                    # Missing _db suffix and wrong case
buchwald_cn.json               # Not descriptive enough
C-N-Coupling-Pd-db.json        # Hyphens instead of underscores
cn_coupling_pd.json            # All lowercase, missing _db
Suzuki_Miyaura_DB.json         # DB should be lowercase
```

---

## 🔧 Quick Rename Commands

Run these commands in PowerShell to fix your current files:

```powershell
cd c:\Git-softwares\Condition-agent\data\rule_db

# Rename existing files to match expected names
mv suzuki.json Suzuki_db.json
mv buchwald_cn.json C_N_Coupling_Pd_db.json
mv ullman_cn.json C_N_Coupling_Cu_db.json

# Verify
ls *.json
```

Expected output after rename:
```
Amide_formation_db.json
C_N_Coupling_Cu_db.json
C_N_Coupling_Pd_db.json
Suzuki_db.json
```

---

## 📊 Where These Names Are Used

### 1. **UI Simple** (`app/ui_simple.py` line 120)
```python
RULE_DATABASES = {
    "C-N Coupling (Cu)": str(SCDB_DIR / "C_N_Coupling_Cu_db.json"),
    "C-N Coupling (Pd)": str(SCDB_DIR / "C_N_Coupling_Pd_db.json"),
    "C-N Coupling (Ni)": str(SCDB_DIR / "C_N_Coupling_Ni_db.json"),
    "Amide Formation": str(SCDB_DIR / "Amide_formation_db.json"),
    "Suzuki Coupling": str(SCDB_DIR / "Suzuki_db.json"),
}
```

### 2. **Local Recommendation CLI** (`app/local_recommendation_cli.py` line 138)
```python
db_map = {
    "Suzuki": "data/rule_db/Suzuki_db.json",
    "Suzuki_CC": "data/rule_db/Suzuki_db.json",
    "suzuki_miyaura": "data/rule_db/Suzuki_db.json",
    "C_N_Coupling": "data/rule_db/C_N_Coupling_Cu_db.json",
    "cn_coupling": "data/rule_db/C_N_Coupling_Cu_db.json",
    "Amide_formation": "data/rule_db/Amide_formation_db.json",
    "amide_coupling": "data/rule_db/Amide_formation_db.json",
}
```

### 3. **Default Path** (`app/recommendation_cli_utils.py` line 15)
```python
DEFAULT_SCDB_PATH = "data/rule_db/Suzuki_db.json"
```

---

## 📝 Adding New Rule Databases

When adding a new reaction type, follow these steps:

### Step 1: Create the file with proper naming

**Template:**
```
{ReactionFamily}_{Variant}_db.json
```

**Examples:**
- Heck reaction: `Heck_db.json`
- Negishi coupling: `Negishi_db.json`
- Grignard addition: `Grignard_addition_db.json`
- Esterification: `Esterification_db.json`

### Step 2: Update the mappings in code

**In `app/ui_simple.py`:**
```python
RULE_DATABASES = {
    # ... existing entries ...
    "Heck Reaction": str(SCDB_DIR / "Heck_db.json"),
}
```

**In `app/local_recommendation_cli.py`:**
```python
db_map = {
    # ... existing entries ...
    "Heck": "data/rule_db/Heck_db.json",
    "heck_reaction": "data/rule_db/Heck_db.json",
}
```

### Step 3: Ensure JSON structure matches

Your JSON file must have these required fields:
```json
{
  "name": "Human Readable Name",
  "reaction_type": "machine_readable_type",
  "version": "2025-11-06",
  "applies_if": { ... },
  "default_rule": { ... },
  "base_rules": [ ... ],
  "modifiers": [ ... ]
}
```

---

## 🔍 Validation

After renaming, validate your setup:

```powershell
# Check files exist
cd c:\Git-softwares\Condition-agent
python -c "from pathlib import Path; import json; [print(f'✓ {f.name}' if f.exists() and json.loads(f.read_text()) else f'✗ {f.name}') for f in Path('data/rule_db').glob('*_db.json')]"

# Test loading in Python
python -c "from chemtools.rule import RuleEngine; engine = RuleEngine.from_file('data/rule_db/Suzuki_db.json'); print('✓ Suzuki loaded successfully')"
```

---

## 💡 Best Practices

### DO ✅
- Use `{Family}_{Variant}_db.json` pattern
- Use title case (First_Letter_Uppercase)
- Use underscores between words
- Always end with `_db.json`
- Match `reaction_type` in JSON to the family name
- Keep names descriptive but concise

### DON'T ❌
- Use hyphens instead of underscores
- Use all lowercase or all uppercase
- Omit the `_db` suffix
- Use spaces in filenames
- Use inconsistent capitalization
- Create duplicate mappings

---

## 📖 Summary Table

| Component | Pattern | Example |
|-----------|---------|---------|
| **File Name** | `{Family}_{Variant}_db.json` | `Suzuki_db.json` |
| **Capitalization** | Title_Case_With_Underscores | `C_N_Coupling_Pd_db.json` |
| **Suffix** | Always `_db.json` | `Amide_formation_db.json` |
| **reaction_type in JSON** | Lowercase with underscores | `"suzuki"` or `"c_n_coupling"` |
| **Display Name (UI)** | Human-friendly with spaces | `"Suzuki Coupling"` |

---

## ✅ Your To-Do Checklist

- [ ] Rename `suzuki.json` → `Suzuki_db.json`
- [ ] Rename `buchwald_cn.json` → `C_N_Coupling_Pd_db.json`
- [ ] Rename `ullman_cn.json` → `C_N_Coupling_Cu_db.json`
- [ ] Keep `amide_formation_db.json` as is ✓
- [ ] Test loading each database with RuleEngine
- [ ] Verify UI Simple can find all databases
- [ ] Update any custom scripts that reference old names

---

**Created:** 2025-11-06  
**Pattern Source:** `app/ui_simple.py`, `app/local_recommendation_cli.py`  
**Convention:** PascalCase_With_Underscores_db.json
