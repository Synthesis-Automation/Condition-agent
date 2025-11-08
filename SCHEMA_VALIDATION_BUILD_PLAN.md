# Implementation Plan: Unified Schema, Validation & Build System

## Vision Statement

**Goal:** Create a unified, validated database system where:
- **Protocols** = Specific experimental procedures for narrow reaction scopes (1 reaction or closely related variants)
- **Rules** = Generalized condition guidelines for broad reaction families (hundreds of variants)
- **Both** indexed and searchable using the same DRFP similarity infrastructure
- **Schema** standardized, validated, and version-controlled
- **Build tools** automate validation, indexing, and deployment

---

## Part 1: Unified Schema Design

### Core Principle: Semantic Distinction

```
┌─────────────────────────────────────────────────────────────┐
│                   Condition Source                          │
│  (Abstract base - all searchable entities)                  │
└─────────────────────────────────────────────────────────────┘
              ↓                           ↓
   ┌──────────────────────┐    ┌──────────────────────┐
   │    PROTOCOL          │    │       RULE           │
   │  (Specific Recipe)   │    │  (Family Guidelines) │
   ├──────────────────────┤    ├──────────────────────┤
   │ Scope: 1 reaction    │    │ Scope: 100s variants │
   │ Detail: Exact amounts│    │ Detail: Ranges       │
   │ Example: Org.Synth   │    │ Example: Buchwald CN │
   │ chemicals: [...]     │    │ base_rules: [...]    │
   │ workup: [...]        │    │ modifiers: [...]     │
   └──────────────────────┘    └──────────────────────┘
```

### Schema v2.0 - Unified Format

**File:** `chemtools/schema/condition_source_v2.json`

```json
{
  "$schema": "http://json-schema.org/draft-07/schema#",
  "$id": "https://condition-agent.io/schemas/condition-source-v2.json",
  "title": "Condition Source",
  "description": "Unified schema for protocols (specific reactions) and rules (reaction families)",
  "type": "object",
  
  "required": ["schema_version", "source_type", "metadata", "reaction"],
  
  "properties": {
    "schema_version": {
      "type": "string",
      "const": "2.0",
      "description": "Schema version for compatibility checking"
    },
    
    "source_type": {
      "type": "string",
      "enum": ["protocol", "rule"],
      "description": "Distinguishes specific protocols from generalized rules"
    },
    
    "metadata": {
      "type": "object",
      "required": ["id", "name", "version", "created_date"],
      "properties": {
        "id": {
          "type": "string",
          "pattern": "^[a-z0-9_-]+$",
          "description": "Unique identifier (e.g., 'sonogashira_v100p0099', 'suzuki_general_v1')"
        },
        "name": {
          "type": "string",
          "description": "Human-readable name"
        },
        "version": {
          "type": "string",
          "pattern": "^v?\\d+\\.\\d+(\\.\\d+)?$",
          "description": "Semantic version (e.g., 'v1.0', '2.1.3')"
        },
        "created_date": {
          "type": "string",
          "format": "date",
          "description": "ISO 8601 date (YYYY-MM-DD)"
        },
        "last_updated": {
          "type": "string",
          "format": "date"
        },
        "status": {
          "type": "string",
          "enum": ["active", "deprecated", "draft"],
          "default": "active"
        },
        "tags": {
          "type": "array",
          "items": {"type": "string"},
          "description": "Searchable tags (e.g., ['palladium', 'coupling', 'ligand-free'])"
        }
      }
    },
    
    "source": {
      "type": "object",
      "description": "Literature/source citation (primarily for protocols)",
      "properties": {
        "title": {"type": "string"},
        "journal": {"type": "string"},
        "volume": {"type": ["integer", "string"]},
        "pages": {"type": "string"},
        "year": {"type": "integer", "minimum": 1900, "maximum": 2100},
        "doi": {"type": "string", "pattern": "^10\\.\\d{4,}/.+$"},
        "url": {"type": "string", "format": "uri"},
        "authors": {
          "type": "array",
          "items": {"type": "string"}
        }
      }
    },
    
    "reaction": {
      "type": "object",
      "required": ["family"],
      "properties": {
        "family": {
          "type": "string",
          "description": "Reaction family (e.g., 'Suzuki', 'Sonogashira', 'C_N_Coupling')"
        },
        "reaction_smiles": {
          "type": "string",
          "description": "Specific reaction SMILES (required for protocols, optional for rules)"
        },
        "reaction_smarts": {
          "type": "array",
          "items": {"type": "string"},
          "description": "Reaction SMARTS patterns for structural filtering"
        },
        "reference_reactions": {
          "type": "array",
          "items": {"type": "string"},
          "description": "Representative SMILES for DRFP similarity matching (required for rules)",
          "minItems": 1
        },
        "scope": {
          "type": "object",
          "description": "Defines the applicability scope",
          "properties": {
            "scope_type": {
              "type": "string",
              "enum": ["specific", "narrow", "broad", "general"],
              "description": "specific=1 reaction, narrow=5-10 variants, broad=50-100, general=family-wide"
            },
            "compatible_functional_groups": {
              "type": "array",
              "items": {"type": "string"}
            },
            "incompatible_functional_groups": {
              "type": "array",
              "items": {"type": "string"}
            },
            "substrate_limitations": {
              "type": "string",
              "description": "Free-text description of limitations"
            }
          }
        },
        "notes": {"type": "string"}
      }
    },
    
    "conditions": {
      "type": "object",
      "description": "Common condition format (used by both protocols and rules)",
      "properties": {
        "catalyst": {"type": "string"},
        "catalyst_loading_molpct": {"type": ["string", "number"]},
        "ligand": {"type": "string"},
        "ligand_loading_molpct": {"type": ["string", "number"]},
        "base": {"type": "string"},
        "base_equiv": {"type": ["string", "number"]},
        "solvent": {"type": "string"},
        "solvent_volume_ml": {"type": "number"},
        "temperature_C": {"type": ["string", "number"]},
        "time_h": {"type": ["string", "number"]},
        "atmosphere": {"type": "string"},
        "concentration_M": {"type": ["string", "number"]},
        "additives": {
          "type": "array",
          "items": {"type": "string"}
        }
      }
    }
  },
  
  "allOf": [
    {
      "if": {
        "properties": {"source_type": {"const": "protocol"}}
      },
      "then": {
        "required": ["reaction_setup", "conditions"],
        "properties": {
          "reaction_setup": {
            "type": "array",
            "description": "Detailed experimental setup (protocol-specific)",
            "items": {
              "type": "object",
              "required": ["chemicals"],
              "properties": {
                "chemicals": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "required": ["name", "role"],
                    "properties": {
                      "name": {"type": "string"},
                      "abbreviation": {"type": ["string", "null"]},
                      "cas": {"type": ["string", "null"]},
                      "smiles": {"type": ["string", "null"]},
                      "amount": {
                        "type": "object",
                        "properties": {
                          "weight_g": {"type": ["number", "null"]},
                          "mmol": {"type": ["number", "null"]},
                          "volume_ml": {"type": ["number", "null"]},
                          "equivalents": {"type": ["number", "null"]}
                        }
                      },
                      "role": {
                        "type": "string",
                        "enum": ["starting_material", "reagent", "catalyst", "ligand", 
                                "base", "solvent", "additive", "workup"]
                      }
                    }
                  }
                },
                "conditions": {
                  "type": "array",
                  "items": {
                    "type": "object",
                    "properties": {
                      "temperature_C": {"type": "number"},
                      "time_h": {"type": "number"},
                      "atmosphere": {"type": "string"},
                      "stirring_rpm": {"type": "number"}
                    }
                  }
                }
              }
            }
          },
          "workup_and_purification": {
            "type": "array",
            "description": "Post-reaction workup steps",
            "items": {
              "type": "object",
              "properties": {
                "quench": {"type": "array", "items": {"type": "object"}},
                "workup": {"type": "array", "items": {"type": "object"}},
                "purification": {"type": "array", "items": {"type": "object"}},
                "notes": {"type": "array", "items": {"type": "string"}}
              }
            }
          }
        }
      }
    },
    {
      "if": {
        "properties": {"source_type": {"const": "rule"}}
      },
      "then": {
        "required": ["applies_if", "base_rules"],
        "properties": {
          "applies_if": {
            "type": "object",
            "description": "Feature-based applicability criteria",
            "properties": {
              "all": {
                "type": "array",
                "items": {"type": "string"},
                "description": "All features must be present"
              },
              "any": {
                "type": "array",
                "items": {"type": "string"},
                "description": "At least one feature must be present"
              },
              "none": {
                "type": "array",
                "items": {"type": "string"},
                "description": "None of these features should be present"
              }
            }
          },
          "default_rule": {
            "type": "object",
            "required": ["id", "description", "conditions"],
            "properties": {
              "id": {"type": "string"},
              "description": {"type": "string"},
              "conditions": {"$ref": "#/properties/conditions"}
            }
          },
          "base_rules": {
            "type": "array",
            "description": "Conditional rules based on substrate features",
            "items": {
              "type": "object",
              "required": ["name", "id", "conditions"],
              "properties": {
                "name": {"type": "string"},
                "id": {"type": "string"},
                "description": {"type": "string"},
                "reactant_features": {
                  "type": "object",
                  "properties": {
                    "all": {"type": "array", "items": {"type": "string"}},
                    "any": {"type": "array", "items": {"type": "string"}},
                    "none": {"type": "array", "items": {"type": "string"}}
                  }
                },
                "conditions": {"$ref": "#/properties/conditions"}
              }
            }
          },
          "modifiers": {
            "type": "array",
            "description": "Dynamic adjustments based on symptoms/features",
            "items": {
              "type": "object",
              "required": ["id", "when", "suggest"],
              "properties": {
                "id": {"type": "string"},
                "when": {
                  "type": "array",
                  "items": {"type": "string"},
                  "description": "Trigger conditions (features or symptoms)"
                },
                "suggest": {"type": "string"},
                "priority": {"type": "integer", "minimum": 1, "maximum": 10}
              }
            }
          }
        }
      }
    }
  ]
}
```

---

## Part 2: Validation System

### Architecture

```
┌────────────────────────────────────────────────────────────┐
│              Validation Pipeline                           │
├────────────────────────────────────────────────────────────┤
│                                                            │
│  Input: JSON file (protocol or rule)                      │
│     ↓                                                      │
│  1. Schema Validation (JSON Schema)                       │
│     - Structure correctness                               │
│     - Required fields present                             │
│     - Type checking                                       │
│     ↓                                                      │
│  2. Chemical Validation (RDKit)                           │
│     - SMILES parsing                                      │
│     - SMARTS compilation                                  │
│     - CAS number format                                   │
│     ↓                                                      │
│  3. Semantic Validation (Business Logic)                  │
│     - Protocols: Must have reaction_smiles               │
│     - Rules: Must have reference_reactions               │
│     - Conditions: Reasonable ranges                       │
│     - Cross-references valid                              │
│     ↓                                                      │
│  4. Quality Checks (Warnings)                             │
│     - Missing optional fields                             │
│     - Deprecated patterns                                 │
│     - Best practice violations                            │
│     ↓                                                      │
│  Output: Validation Report                                │
│     - Errors (must fix)                                   │
│     - Warnings (should fix)                               │
│     - Info (suggestions)                                  │
│                                                            │
└────────────────────────────────────────────────────────────┘
```

### Implementation: Validator Class

**File:** `chemtools/schema/validator.py`

```python
"""
Unified validator for protocol and rule database files.

Usage:
    python -m chemtools.schema.validator data/protocol_db/my_protocol.json
    python -m chemtools.schema.validator data/rule_db/sonogashira_db.json --strict
"""

from pathlib import Path
from typing import Dict, Any, List, Optional
import json
import jsonschema
from dataclasses import dataclass, field
from enum import Enum

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


class ValidationLevel(Enum):
    ERROR = "error"      # Must fix - breaks schema
    WARNING = "warning"  # Should fix - best practice
    INFO = "info"        # Optional - suggestion


@dataclass
class ValidationMessage:
    level: ValidationLevel
    message: str
    path: str = ""  # JSON path (e.g., "reaction.reaction_smiles")
    fix_suggestion: Optional[str] = None


@dataclass
class ValidationReport:
    file_path: Path
    source_type: Optional[str] = None
    is_valid: bool = False
    errors: List[ValidationMessage] = field(default_factory=list)
    warnings: List[ValidationMessage] = field(default_factory=list)
    info: List[ValidationMessage] = field(default_factory=list)
    
    @property
    def total_issues(self) -> int:
        return len(self.errors) + len(self.warnings)
    
    def add_error(self, message: str, path: str = "", fix: str = None):
        self.errors.append(ValidationMessage(
            ValidationLevel.ERROR, message, path, fix
        ))
    
    def add_warning(self, message: str, path: str = "", fix: str = None):
        self.warnings.append(ValidationMessage(
            ValidationLevel.WARNING, message, path, fix
        ))
    
    def add_info(self, message: str, path: str = "", fix: str = None):
        self.info.append(ValidationMessage(
            ValidationLevel.INFO, message, path, fix
        ))
    
    def print_report(self, show_info: bool = False):
        print(f"\n{'='*80}")
        print(f"Validation Report: {self.file_path.name}")
        print(f"{'='*80}")
        
        if self.source_type:
            print(f"Type: {self.source_type.upper()}")
        
        if self.is_valid and not self.warnings:
            print("✅ VALID - No issues found")
            return
        
        if self.errors:
            print(f"\n❌ ERRORS ({len(self.errors)}):")
            for err in self.errors:
                print(f"  • {err.path}: {err.message}")
                if err.fix_suggestion:
                    print(f"    Fix: {err.fix_suggestion}")
        
        if self.warnings:
            print(f"\n⚠️  WARNINGS ({len(self.warnings)}):")
            for warn in self.warnings:
                print(f"  • {warn.path}: {warn.message}")
                if warn.fix_suggestion:
                    print(f"    Suggestion: {warn.fix_suggestion}")
        
        if show_info and self.info:
            print(f"\nℹ️  INFO ({len(self.info)}):")
            for info_msg in self.info:
                print(f"  • {info_msg.path}: {info_msg.message}")
        
        print(f"\n{'='*80}")
        if self.errors:
            print("❌ VALIDATION FAILED")
        elif self.warnings:
            print("⚠️  VALIDATION PASSED WITH WARNINGS")
        else:
            print("✅ VALIDATION PASSED")
        print(f"{'='*80}\n")


class ConditionSourceValidator:
    """
    Validates protocol and rule database files against unified schema.
    """
    
    def __init__(self, schema_path: Optional[Path] = None):
        if schema_path is None:
            schema_path = Path(__file__).parent / "condition_source_v2.json"
        
        self.schema_path = Path(schema_path)
        self.schema = self._load_schema()
    
    def _load_schema(self) -> dict:
        """Load JSON schema"""
        if not self.schema_path.exists():
            raise FileNotFoundError(f"Schema not found: {self.schema_path}")
        
        with open(self.schema_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    def validate_file(
        self,
        file_path: Path,
        strict: bool = False
    ) -> ValidationReport:
        """
        Validate a protocol or rule database file.
        
        Args:
            file_path: Path to JSON file
            strict: If True, warnings become errors
        
        Returns:
            ValidationReport with errors, warnings, and info
        """
        report = ValidationReport(file_path=Path(file_path))
        
        # 1. Load and parse JSON
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            report.add_error(
                f"Invalid JSON: {e}",
                path="",
                fix="Fix JSON syntax errors"
            )
            return report
        except Exception as e:
            report.add_error(f"Failed to read file: {e}")
            return report
        
        # Detect source type
        source_type = data.get('source_type')
        report.source_type = source_type
        
        # 2. JSON Schema validation
        self._validate_schema(data, report)
        
        # 3. Chemical validation (SMILES, SMARTS)
        if RDKIT_AVAILABLE:
            self._validate_chemistry(data, report)
        else:
            report.add_warning(
                "RDKit not available - skipping chemical validation",
                fix="Install RDKit: pip install rdkit"
            )
        
        # 4. Semantic validation
        self._validate_semantics(data, report)
        
        # 5. Quality checks
        self._quality_checks(data, report)
        
        # Determine overall validity
        if strict:
            report.is_valid = len(report.errors) == 0 and len(report.warnings) == 0
        else:
            report.is_valid = len(report.errors) == 0
        
        return report
    
    def _validate_schema(self, data: dict, report: ValidationReport):
        """Validate against JSON schema"""
        try:
            jsonschema.validate(instance=data, schema=self.schema)
        except jsonschema.ValidationError as e:
            report.add_error(
                f"Schema validation failed: {e.message}",
                path='.'.join(str(p) for p in e.path),
                fix="Check schema requirements"
            )
        except jsonschema.SchemaError as e:
            report.add_error(f"Invalid schema: {e.message}")
    
    def _validate_chemistry(self, data: dict, report: ValidationReport):
        """Validate SMILES and SMARTS patterns"""
        reaction = data.get('reaction', {})
        
        # Validate reaction_smiles
        if 'reaction_smiles' in reaction:
            smiles = reaction['reaction_smiles']
            if '>>' not in smiles:
                report.add_error(
                    "reaction_smiles must contain '>>' separator",
                    path="reaction.reaction_smiles"
                )
            # TODO: Parse with RDKit
        
        # Validate reaction_smarts
        if 'reaction_smarts' in reaction:
            for idx, smarts in enumerate(reaction['reaction_smarts']):
                if '>>' not in smarts:
                    report.add_error(
                        f"SMARTS pattern must contain '>>' separator",
                        path=f"reaction.reaction_smarts[{idx}]"
                    )
                # TODO: Compile with RDKit
        
        # Validate reference_reactions
        if 'reference_reactions' in reaction:
            for idx, ref_smiles in enumerate(reaction['reference_reactions']):
                if '>>' not in ref_smiles:
                    report.add_error(
                        f"Reference reaction must contain '>>' separator",
                        path=f"reaction.reference_reactions[{idx}]"
                    )
        
        # Validate chemical SMILES (for protocols)
        if data.get('source_type') == 'protocol':
            for setup_idx, setup in enumerate(data.get('reaction_setup', [])):
                for chem_idx, chem in enumerate(setup.get('chemicals', [])):
                    smiles = chem.get('smiles')
                    if smiles:
                        # TODO: Validate with RDKit
                        pass
    
    def _validate_semantics(self, data: dict, report: ValidationReport):
        """Business logic validation"""
        source_type = data.get('source_type')
        reaction = data.get('reaction', {})
        
        if source_type == 'protocol':
            # Protocols must have specific reaction_smiles
            if not reaction.get('reaction_smiles'):
                report.add_error(
                    "Protocols must have reaction_smiles",
                    path="reaction.reaction_smiles",
                    fix="Add the specific reaction SMILES for this protocol"
                )
            
            # Protocols should have narrow scope
            scope = reaction.get('scope', {})
            scope_type = scope.get('scope_type')
            if scope_type in ('broad', 'general'):
                report.add_warning(
                    f"Protocol has '{scope_type}' scope - should be 'specific' or 'narrow'",
                    path="reaction.scope.scope_type",
                    fix="Protocols are for specific reactions. Consider creating a rule database instead."
                )
        
        elif source_type == 'rule':
            # Rules must have reference_reactions for DRFP matching
            if not reaction.get('reference_reactions'):
                report.add_error(
                    "Rules must have reference_reactions for similarity matching",
                    path="reaction.reference_reactions",
                    fix="Add 3-10 representative reaction SMILES examples"
                )
            elif len(reaction.get('reference_reactions', [])) < 3:
                report.add_warning(
                    f"Only {len(reaction['reference_reactions'])} reference reactions - recommend 5-10",
                    path="reaction.reference_reactions",
                    fix="Add more diverse reference reactions for better matching"
                )
            
            # Rules should have broader scope
            scope = reaction.get('scope', {})
            scope_type = scope.get('scope_type')
            if scope_type == 'specific':
                report.add_warning(
                    "Rule has 'specific' scope - consider protocol format instead",
                    path="reaction.scope.scope_type"
                )
    
    def _quality_checks(self, data: dict, report: ValidationReport):
        """Best practice and quality checks"""
        metadata = data.get('metadata', {})
        
        # Check for tags
        if not metadata.get('tags'):
            report.add_info(
                "No tags provided - tags improve searchability",
                path="metadata.tags",
                fix="Add descriptive tags like ['palladium', 'coupling', 'ligand']"
            )
        
        # Check for version
        version = metadata.get('version', '')
        if not version.startswith('v'):
            report.add_info(
                "Version should start with 'v' (e.g., 'v1.0')",
                path="metadata.version"
            )
        
        # Check for DOI (protocols)
        if data.get('source_type') == 'protocol':
            source = data.get('source', {})
            if not source.get('doi'):
                report.add_info(
                    "DOI not provided - recommended for literature protocols",
                    path="source.doi"
                )


def validate_batch(
    directory: Path,
    pattern: str = "*.json",
    strict: bool = False
) -> Dict[Path, ValidationReport]:
    """
    Validate all JSON files in a directory.
    
    Args:
        directory: Directory containing JSON files
        pattern: Glob pattern for files
        strict: Strict validation mode
    
    Returns:
        Dictionary mapping file paths to validation reports
    """
    validator = ConditionSourceValidator()
    results = {}
    
    for file_path in Path(directory).glob(pattern):
        report = validator.validate_file(file_path, strict=strict)
        results[file_path] = report
    
    return results


# CLI interface
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Validate protocol and rule database files"
    )
    parser.add_argument('path', type=Path, help="File or directory to validate")
    parser.add_argument('--strict', action='store_true', help="Treat warnings as errors")
    parser.add_argument('--show-info', action='store_true', help="Show info messages")
    parser.add_argument('--batch', action='store_true', help="Validate all JSON files in directory")
    
    args = parser.parse_args()
    
    if args.batch:
        results = validate_batch(args.path, strict=args.strict)
        
        print(f"\nBatch Validation: {len(results)} files\n")
        
        valid_count = sum(1 for r in results.values() if r.is_valid)
        error_count = sum(1 for r in results.values() if r.errors)
        warning_count = sum(1 for r in results.values() if r.warnings and not r.errors)
        
        for file_path, report in results.items():
            status = "✅" if report.is_valid else "❌"
            print(f"{status} {file_path.name} - {len(report.errors)} errors, {len(report.warnings)} warnings")
        
        print(f"\nSummary: {valid_count}/{len(results)} valid")
        
    else:
        validator = ConditionSourceValidator()
        report = validator.validate_file(args.path, strict=args.strict)
        report.print_report(show_info=args.show_info)
        
        exit(0 if report.is_valid else 1)
```

---

## Part 3: Build System

### Build Pipeline

```
┌─────────────────────────────────────────────────────────────┐
│                    Build Pipeline                           │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  Input: Source files (data/protocol_db/, data/rule_db/)    │
│     ↓                                                       │
│  1. Validate                                                │
│     - Run validator on all files                           │
│     - Fail build if errors found                           │
│     ↓                                                       │
│  2. Compute DRFP                                            │
│     - Protocols: Compute from reaction_smiles             │
│     - Rules: Compute from reference_reactions             │
│     ↓                                                       │
│  3. Build Unified Index                                     │
│     - Create single index.json (metadata)                  │
│     - Create single fingerprints.npz (DRFPs)              │
│     ↓                                                       │
│  4. Generate Statistics                                     │
│     - Count protocols/rules                                │
│     - Family distribution                                  │
│     - Coverage analysis                                     │
│     ↓                                                       │
│  5. Package                                                 │
│     - Version index                                        │
│     - Create manifest                                       │
│     - Generate changelog                                    │
│     ↓                                                       │
│  Output: unified_index_v{version}.tar.gz                   │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### Implementation: Builder Class

**File:** `chemtools/schema/builder.py`

```python
"""
Build system for creating unified index from protocols and rules.

Usage:
    # Build unified index
    python -m chemtools.schema.builder build \\
        --protocols data/protocol_db \\
        --rules data/rule_db \\
        --output data/unified_index_v1.0
    
    # Validate before building
    python -m chemtools.schema.builder validate \\
        --protocols data/protocol_db \\
        --rules data/rule_db
"""

from pathlib import Path
from typing import List, Dict, Any, Optional
import json
import numpy as np
from datetime import datetime
from dataclasses import dataclass, field

from .validator import ConditionSourceValidator, ValidationReport


@dataclass
class BuildConfig:
    """Configuration for build process"""
    protocol_dir: Path
    rule_dir: Path
    output_dir: Path
    version: str = "1.0"
    fail_on_warnings: bool = False
    skip_validation: bool = False


@dataclass
class BuildReport:
    """Report of build process"""
    success: bool = False
    version: str = ""
    timestamp: datetime = field(default_factory=datetime.now)
    
    # Counts
    num_protocols: int = 0
    num_rules: int = 0
    num_total: int = 0
    
    # Validation
    validation_errors: int = 0
    validation_warnings: int = 0
    
    # DRFP computation
    drfp_computed: int = 0
    drfp_failed: int = 0
    
    # Index stats
    index_size_mb: float = 0.0
    
    messages: List[str] = field(default_factory=list)
    
    def add_message(self, msg: str):
        self.messages.append(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")
        print(msg)


class UnifiedIndexBuilder:
    """
    Builds unified index from protocol and rule database files.
    """
    
    def __init__(self, config: BuildConfig):
        self.config = config
        self.validator = ConditionSourceValidator()
        self.report = BuildReport(version=config.version)
    
    def build(self) -> BuildReport:
        """Execute full build pipeline"""
        self.report.add_message("="*80)
        self.report.add_message(f"Building Unified Index v{self.config.version}")
        self.report.add_message("="*80)
        
        try:
            # 1. Validate
            if not self.config.skip_validation:
                self.report.add_message("\n[1/5] Validating source files...")
                if not self._validate_all():
                    self.report.add_message("❌ Validation failed - aborting build")
                    return self.report
            else:
                self.report.add_message("\n[1/5] Skipping validation (--skip-validation)")
            
            # 2. Load and parse
            self.report.add_message("\n[2/5] Loading source files...")
            protocols = self._load_protocols()
            rules = self._load_rules()
            
            self.report.num_protocols = len(protocols)
            self.report.num_rules = len(rules)
            self.report.num_total = len(protocols) + len(rules)
            
            self.report.add_message(f"  Loaded {len(protocols)} protocols")
            self.report.add_message(f"  Loaded {len(rules)} rules")
            
            # 3. Compute DRFP
            self.report.add_message("\n[3/5] Computing DRFP fingerprints...")
            protocol_drfps = self._compute_drfps_for_protocols(protocols)
            rule_drfps = self._compute_drfps_for_rules(rules)
            
            # 4. Build index
            self.report.add_message("\n[4/5] Building unified index...")
            self._build_index(protocols, rules, protocol_drfps, rule_drfps)
            
            # 5. Generate stats
            self.report.add_message("\n[5/5] Generating statistics...")
            self._generate_stats(protocols, rules)
            
            self.report.success = True
            self.report.add_message("\n" + "="*80)
            self.report.add_message("✅ BUILD SUCCESSFUL")
            self.report.add_message("="*80)
            
        except Exception as e:
            self.report.add_message(f"\n❌ BUILD FAILED: {e}")
            self.report.success = False
        
        return self.report
    
    def _validate_all(self) -> bool:
        """Validate all source files"""
        protocol_files = list(self.config.protocol_dir.glob("*.json"))
        rule_files = list(self.config.rule_dir.glob("*.json"))
        
        all_valid = True
        
        for file_path in protocol_files + rule_files:
            report = self.validator.validate_file(
                file_path,
                strict=self.config.fail_on_warnings
            )
            
            if not report.is_valid:
                all_valid = False
                self.report.add_message(f"  ❌ {file_path.name}: {len(report.errors)} errors")
            elif report.warnings:
                self.report.add_message(f"  ⚠️  {file_path.name}: {len(report.warnings)} warnings")
            else:
                self.report.add_message(f"  ✅ {file_path.name}")
            
            self.report.validation_errors += len(report.errors)
            self.report.validation_warnings += len(report.warnings)
        
        return all_valid
    
    def _load_protocols(self) -> List[Dict[str, Any]]:
        """Load all protocol files"""
        protocols = []
        for file_path in self.config.protocol_dir.glob("*.json"):
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                data['_source_file'] = str(file_path)
                protocols.append(data)
        return protocols
    
    def _load_rules(self) -> List[Dict[str, Any]]:
        """Load all rule files"""
        rules = []
        for file_path in self.config.rule_dir.glob("*.json"):
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                data['_source_file'] = str(file_path)
                rules.append(data)
        return rules
    
    def _compute_drfps_for_protocols(self, protocols: List[Dict]) -> Dict[str, np.ndarray]:
        """Compute DRFP for each protocol"""
        # TODO: Implement DRFP computation
        return {}
    
    def _compute_drfps_for_rules(self, rules: List[Dict]) -> Dict[str, List[np.ndarray]]:
        """Compute DRFP for each rule's reference reactions"""
        # TODO: Implement DRFP computation for reference reactions
        return {}
    
    def _build_index(
        self,
        protocols: List[Dict],
        rules: List[Dict],
        protocol_drfps: Dict,
        rule_drfps: Dict
    ):
        """Build unified index files"""
        output_dir = Path(self.config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save index.json (metadata)
        index = {
            'version': self.config.version,
            'build_date': datetime.now().isoformat(),
            'num_protocols': len(protocols),
            'num_rules': len(rules),
            'protocols': [self._extract_metadata(p) for p in protocols],
            'rules': [self._extract_metadata(r) for r in rules]
        }
        
        with open(output_dir / 'index.json', 'w', encoding='utf-8') as f:
            json.dump(index, f, indent=2)
        
        # Save fingerprints.npz (DRFPs)
        # TODO: Save all DRFPs to single NPZ file
        
        self.report.add_message(f"  Index saved to {output_dir}")
    
    def _extract_metadata(self, data: Dict) -> Dict:
        """Extract indexable metadata from source"""
        return {
            'id': data.get('metadata', {}).get('id'),
            'name': data.get('metadata', {}).get('name'),
            'source_type': data.get('source_type'),
            'family': data.get('reaction', {}).get('family'),
            'tags': data.get('metadata', {}).get('tags', []),
            'source_file': data.get('_source_file')
        }
    
    def _generate_stats(self, protocols: List[Dict], rules: List[Dict]):
        """Generate statistics report"""
        stats = {
            'protocols': {
                'count': len(protocols),
                'families': self._count_families(protocols)
            },
            'rules': {
                'count': len(rules),
                'families': self._count_families(rules)
            }
        }
        
        stats_path = Path(self.config.output_dir) / 'stats.json'
        with open(stats_path, 'w', encoding='utf-8') as f:
            json.dump(stats, f, indent=2)
        
        self.report.add_message(f"  Statistics saved to {stats_path}")
    
    def _count_families(self, items: List[Dict]) -> Dict[str, int]:
        """Count items per reaction family"""
        families = {}
        for item in items:
            family = item.get('reaction', {}).get('family', 'unknown')
            families[family] = families.get(family, 0) + 1
        return families


# CLI interface
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Build unified index from protocols and rules")
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Build command
    build_parser = subparsers.add_parser('build', help="Build unified index")
    build_parser.add_argument('--protocols', type=Path, required=True, help="Protocol directory")
    build_parser.add_argument('--rules', type=Path, required=True, help="Rule directory")
    build_parser.add_argument('--output', type=Path, required=True, help="Output directory")
    build_parser.add_argument('--version', default="1.0", help="Index version")
    build_parser.add_argument('--fail-on-warnings', action='store_true')
    build_parser.add_argument('--skip-validation', action='store_true')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help="Validate source files")
    validate_parser.add_argument('--protocols', type=Path, required=True)
    validate_parser.add_argument('--rules', type=Path, required=True)
    validate_parser.add_argument('--strict', action='store_true')
    
    args = parser.parse_args()
    
    if args.command == 'build':
        config = BuildConfig(
            protocol_dir=args.protocols,
            rule_dir=args.rules,
            output_dir=args.output,
            version=args.version,
            fail_on_warnings=args.fail_on_warnings,
            skip_validation=args.skip_validation
        )
        
        builder = UnifiedIndexBuilder(config)
        report = builder.build()
        
        exit(0 if report.success else 1)
    
    elif args.command == 'validate':
        # Batch validate both directories
        from .validator import validate_batch
        
        print("Validating protocols...")
        protocol_results = validate_batch(args.protocols, strict=args.strict)
        
        print("\nValidating rules...")
        rule_results = validate_batch(args.rules, strict=args.strict)
        
        all_results = {**protocol_results, **rule_results}
        valid_count = sum(1 for r in all_results.values() if r.is_valid)
        
        print(f"\nValidation complete: {valid_count}/{len(all_results)} valid")
        
        exit(0 if valid_count == len(all_results) else 1)
```

---

## Part 4: Migration Plan

### Week-by-Week Roadmap

#### Week 1: Schema Design & Validation Foundation
**Days 1-2:** Schema design
- ✅ Finalize unified schema v2.0
- ✅ Write JSON schema file
- ✅ Document schema with examples

**Days 3-4:** Validator implementation
- ✅ Implement `ConditionSourceValidator` class
- ✅ JSON schema validation
- ✅ Chemical validation (RDKit)
- ✅ Semantic validation

**Day 5:** Testing
- ✅ Test validator on existing files
- ✅ Fix validation errors
- ✅ Document validation rules

#### Week 2: Build System & Index Generation
**Days 1-2:** Builder foundation
- ✅ Implement `UnifiedIndexBuilder` class
- ✅ File loading and parsing
- ✅ DRFP computation integration

**Days 3-4:** Index generation
- ✅ Unified index creation
- ✅ NPZ file generation
- ✅ Statistics and reporting

**Day 5:** Testing
- ✅ Test build pipeline
- ✅ Verify index format
- ✅ Performance benchmarking

#### Week 3: Database Migration
**Days 1-2:** Protocol migration
- ✅ Add schema_version to all protocols
- ✅ Standardize condition format
- ✅ Validate all protocols

**Days 3-4:** Rule migration
- ✅ Add reference_reactions to all rules
- ✅ Add schema_version
- ✅ Validate all rules

**Day 5:** Quality assurance
- ✅ Run full validation suite
- ✅ Fix all errors
- ✅ Address warnings

#### Week 4: Integration & Testing
**Days 1-2:** Unified recommender
- ✅ Integrate with UnifiedRecommender
- ✅ Update agent tools
- ✅ Backward compatibility

**Days 3-4:** End-to-end testing
- ✅ Test with 100+ reactions
- ✅ Performance testing
- ✅ Edge case validation

**Day 5:** Documentation
- ✅ Update user guides
- ✅ API documentation
- ✅ Migration guide

---

## Part 5: CI/CD Integration

### Automated Validation Pipeline

```yaml
# .github/workflows/validate-databases.yml
name: Validate Database Files

on:
  push:
    paths:
      - 'data/protocol_db/**.json'
      - 'data/rule_db/**.json'
  pull_request:
    paths:
      - 'data/protocol_db/**.json'
      - 'data/rule_db/**.json'

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      
      - name: Install dependencies
        run: |
          pip install -r requirements.txt
          pip install jsonschema rdkit
      
      - name: Validate protocol databases
        run: |
          python -m chemtools.schema.validator data/protocol_db --batch --strict
      
      - name: Validate rule databases
        run: |
          python -m chemtools.schema.validator data/rule_db --batch --strict
      
      - name: Build unified index
        run: |
          python -m chemtools.schema.builder build \\
            --protocols data/protocol_db \\
            --rules data/rule_db \\
            --output build/unified_index \\
            --fail-on-warnings
      
      - name: Upload index artifact
        uses: actions/upload-artifact@v3
        with:
          name: unified-index
          path: build/unified_index/
```

---

## Summary & Next Steps

### What We're Building

1. **Unified Schema v2.0**
   - Single schema for protocols (specific) and rules (general)
   - Versioned, extensible, well-documented

2. **Validation System**
   - Multi-level validation (schema, chemistry, semantics, quality)
   - Clear error messages with fix suggestions
   - CLI + programmatic API

3. **Build System**
   - Automated index generation
   - DRFP computation
   - Statistics and reporting

4. **CI/CD Pipeline**
   - Automated validation on every commit
   - Build verification
   - Quality gates

### Immediate Next Steps

**This Week:**
1. Review and approve schema design
2. Implement basic validator (schema validation only)
3. Test on 5-10 existing files

**Next 2 Weeks:**
1. Complete validator with all validation levels
2. Implement builder skeleton
3. Start migrating 3 example files to new schema

**Month 1:**
1. Full validation system
2. Full build system
3. Migrate all databases

**Month 2:**
1. Integration with unified recommender
2. CI/CD setup
3. Documentation

### Success Criteria

✅ All database files pass validation  
✅ Unified index builds successfully  
✅ Recommender works with unified index  
✅ Build time <5 minutes  
✅ Zero manual steps (fully automated)

---

## Conclusion

This plan provides a complete, phased approach to:
- **Standardizing** schema (protocols vs rules)
- **Validating** data quality (errors, warnings, best practices)
- **Building** unified index (DRFP + metadata)
- **Automating** the pipeline (CI/CD)

The key insight: **Protocols and rules are the same type (condition sources) with different scopes**. By recognizing this, we can unify 90% of the infrastructure while preserving the semantic distinction where it matters.

**Recommended action:** Start with Week 1 (schema + validator) to establish the foundation.
