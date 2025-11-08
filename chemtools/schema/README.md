# Schema Validation & Build System

Complete system for validating and building unified indexes from **protocol** (specific reactions) and **rule** (reaction families) databases.

## 📁 Files

- `condition_source_v2.json` - Unified JSON schema for protocols and rules
- `validator.py` - Multi-level validation system with CLI
- `builder.py` - Build system for creating unified indexes
- `__init__.py` - Module exports

## 🎯 Core Concepts

### Protocols vs Rules

**Protocols** = Specific experimental procedures (narrow scope)
- Example: Single Org. Synth. procedure with exact amounts
- Required: `reaction_smiles`, `reaction_setup`, detailed `chemicals`
- Scope: 1 reaction or closely related variants

**Rules** = Generalized condition guidelines (broad scope)
- Example: "Suzuki coupling general conditions"
- Required: `reference_reactions` (for DRFP matching), `applies_if`, `base_rules`
- Scope: Entire reaction family (100s of variants)

## 🔍 Validation System

### Four Validation Levels

1. **Schema Validation** - Structure and types (JSON Schema)
2. **Chemical Validation** - SMILES/SMARTS validity (RDKit)
3. **Semantic Validation** - Business logic (protocols need specific reactions, rules need references)
4. **Quality Checks** - Best practices (tags, DOI, versioning)

### Usage

```bash
# Validate a single file
python -m chemtools.schema.validator data/protocol_db/my_protocol.json

# Validate with strict mode (warnings = errors)
python -m chemtools.schema.validator data/rule_db/sonogashira_v2.json --strict

# Batch validate directory
python -m chemtools.schema.validator data/protocol_db_v2 --batch

# Show info messages (suggestions)
python -m chemtools.schema.validator my_file.json --show-info
```

### Validation Output

```
================================================================================
Validation Report: alpha_arylation_dong_v100p0099_v2.json
================================================================================
Type: PROTOCOL
✅ VALID - No issues found
================================================================================
```

With errors/warnings:
```
❌ ERRORS (2):
  • reaction.reaction_smiles: Protocols must have reaction_smiles
    Fix: Add the specific reaction SMILES for this protocol
  • reaction.reference_reactions[1]: Invalid reactant SMILES: X#Y
    
⚠️  WARNINGS (1):
  • reaction.reference_reactions: Only 2 reference reactions - recommend 5-10
    Suggestion: Add more diverse reference reactions for better matching
```

## 🏗️ Build System

### Build Pipeline

5-step automated pipeline:
1. **Validate** all source files (optional: `--skip-validation`)
2. **Load** and parse protocols + rules
3. **Compute DRFP** fingerprints for similarity matching
4. **Build index** - create unified `index.json` + `fingerprints.npz`
5. **Generate statistics** - counts, families, coverage

### Usage

```bash
# Build unified index
python -m chemtools.schema.builder build \
    --protocols data/protocol_db_v2 \
    --rules data/rule_db_v2 \
    --output build/unified_index_v2.0 \
    --version 2.0

# Validate before building (batch validation)
python -m chemtools.schema.builder validate \
    --protocols data/protocol_db_v2 \
    --rules data/rule_db_v2 \
    --strict
```

### Build Output

```
================================================================================
Building Unified Index v2.0
================================================================================

[1/5] Validating source files...
  ✅ alpha_arylation_dong_v100p0099_v2.json
  ✅ sonogashira_v2.json

[2/5] Loading source files...
  Loaded 1 protocols
  Loaded 1 rules

[3/5] Computing DRFP fingerprints...
  Computed 1 protocol DRFPs (0 failed)
  Computed 6 rule reference DRFPs (0 failed)

[4/5] Building unified index...
  Index saved to build\unified_index_demo
  Index size: 0.00 MB

[5/5] Generating statistics...
  Statistics saved to build\unified_index_demo\stats.json

================================================================================
✅ BUILD SUCCESSFUL
================================================================================
```

### Output Files

- `index.json` - Metadata for all protocols and rules
- `fingerprints.npz` - DRFP vectors for similarity search
- `stats.json` - Build statistics (counts, families, DRFP success rate)

## 📝 Schema v2.0 Format

### Protocol Example

```json
{
  "schema_version": "2.0",
  "source_type": "protocol",
  "metadata": {
    "id": "alpha_arylation_dong_v100p0099",
    "name": "α-Arylation of Cyclopentanones",
    "version": "v1.0",
    "created_date": "2025-11-08",
    "tags": ["palladium", "enamine", "alpha-arylation"]
  },
  "source": {
    "title": "α-Arylation of Cyclopentanones...",
    "journal": "Organic Syntheses",
    "doi": "10.15227/orgsyn.100.0099"
  },
  "reaction": {
    "family": "Pd_Enamine_alpha_Arylation",
    "reaction_smiles": "O=C1CCCC1.Brc2ccc(C(C)=O)cc2>>...",
    "scope": {
      "scope_type": "narrow",
      "compatible_functional_groups": ["ketone", "aromatic"]
    }
  },
  "conditions": {
    "catalyst": "Pd(OAc)2",
    "catalyst_loading_molpct": 1.0,
    "ligand": "P(o-tol)3",
    "base": "NaOAc",
    "solvent": "1,4-dioxane",
    "temperature_C": 130,
    "time_h": 24.0
  },
  "reaction_setup": [
    {
      "chemicals": [
        {
          "name": "Cyclopentanone",
          "smiles": "O=C1CCCC1",
          "amount": {"weight_g": 1.68, "mmol": 20.0, "equivalents": 1.0},
          "role": "starting_material"
        }
      ]
    }
  ]
}
```

### Rule Example

```json
{
  "schema_version": "2.0",
  "source_type": "rule",
  "metadata": {
    "id": "sonogashira_v2",
    "name": "Sonogashira Coupling Guidelines",
    "version": "v2.0",
    "tags": ["sonogashira", "palladium", "alkyne"]
  },
  "reaction": {
    "family": "Sonogashira",
    "reference_reactions": [
      "Ic1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",
      "Brc1ccc(C(F)(F)F)cc1.C#CCCCCC>>..."
    ],
    "scope": {
      "scope_type": "broad",
      "substrate_limitations": "Requires terminal alkyne..."
    }
  },
  "applies_if": {
    "all": ["terminal_alkyne_present"],
    "any": ["aryl_halide_present", "vinyl_halide_present"]
  },
  "default_rule": {
    "id": "DEF_sonogashira_general",
    "conditions": {
      "catalyst": "PdCl2(PPh3)2",
      "catalyst_loading_molpct": "0.5-2.0",
      "base": "Et3N or DIPEA",
      "temperature_C": "40-80"
    }
  },
  "base_rules": [
    {
      "name": "Aryl iodides/bromides",
      "reactant_features": {"any": ["sp2_bromide_present"]},
      "conditions": {...}
    }
  ],
  "modifiers": [
    {
      "id": "MOD_glaser_suppression",
      "when": ["glaser_coupling_observed"],
      "suggest": "Use copper-free conditions..."
    }
  ]
}
```

## 🔧 Programmatic API

### Validation

```python
from chemtools.schema import ConditionSourceValidator

validator = ConditionSourceValidator()
report = validator.validate_file("my_protocol.json", strict=False)

if report.is_valid:
    print("✅ Valid!")
else:
    for error in report.errors:
        print(f"❌ {error.path}: {error.message}")
```

### Building

```python
from chemtools.schema import UnifiedIndexBuilder, BuildConfig
from pathlib import Path

config = BuildConfig(
    protocol_dir=Path("data/protocol_db_v2"),
    rule_dir=Path("data/rule_db_v2"),
    output_dir=Path("build/unified_index"),
    version="2.0"
)

builder = UnifiedIndexBuilder(config)
report = builder.build()

if report.success:
    print(f"✅ Built index with {report.num_protocols} protocols, {report.num_rules} rules")
    print(f"   DRFP computed: {report.drfp_computed}, failed: {report.drfp_failed}")
```

## 📊 Statistics

The build system generates comprehensive statistics:

```json
{
  "build_info": {
    "version": "2.0",
    "build_date": "2025-11-08T19:41:22",
    "total_sources": 2
  },
  "protocols": {
    "count": 1,
    "families": {"Pd_Enamine_alpha_Arylation": 1}
  },
  "rules": {
    "count": 1,
    "families": {"Sonogashira": 1}
  },
  "drfp": {
    "computed": 6,
    "failed": 0
  }
}
```

## 🎓 Migration Guide

### Converting v1 Protocol to v2

1. **Wrap in object** (v1 was array, v2 is object)
2. **Add metadata section**:
   ```json
   "metadata": {
     "id": "unique_id",
     "name": "Protocol name",
     "version": "v1.0",
     "created_date": "2025-11-08",
     "tags": ["tag1", "tag2"]
   }
   ```
3. **Add schema_version and source_type**:
   ```json
   "schema_version": "2.0",
   "source_type": "protocol"
   ```
4. **Add conditions summary** (extracted from reaction_setup)
5. **Add scope** to reaction

### Converting v1 Rule to v2

1. **Add metadata section** (same as protocols)
2. **Add reference_reactions**:
   ```json
   "reaction": {
     "reference_reactions": [
       "Ic1ccccc1.C#Cc1ccccc1>>...",
       "Brc1ccc(C)cc1.C#CCCCCC>>..."
     ]
   }
   ```
3. **Standardize conditions format** (use unified schema)
4. **Add schema_version and source_type**

## 🚀 Next Steps

1. **Migrate existing databases** - Convert all v1 files to v2 format
2. **Add CI/CD** - Automate validation on every commit
3. **Integrate with recommender** - Use unified index for similarity search
4. **Add more reference reactions** - Populate all rules with 5-10 examples

## 📚 References

- Full implementation plan: `SCHEMA_VALIDATION_BUILD_PLAN.md`
- JSON Schema spec: `condition_source_v2.json`
- Example protocol (v2): `data/protocol_db/alpha_arylation_dong_v100p0099_v2.json`
- Example rule (v2): `data/rule_db/sonogashira_v2.json`
