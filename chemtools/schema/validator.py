"""
Unified validator for protocol and rule database files.

Usage:
    # Validate a single file
    python -m chemtools.schema.validator data/protocol_db/my_protocol.json
    
    # Validate with strict mode (warnings = errors)
    python -m chemtools.schema.validator data/rule_db/sonogashira_db.json --strict
    
    # Batch validate directory
    python -m chemtools.schema.validator data/protocol_db --batch
    
    # Show info messages
    python -m chemtools.schema.validator my_file.json --show-info
"""

from pathlib import Path
from typing import Dict, Any, List, Optional
import json
import jsonschema
from dataclasses import dataclass, field
from enum import Enum

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
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
                if err.path:
                    print(f"  • {err.path}: {err.message}")
                else:
                    print(f"  • {err.message}")
                if err.fix_suggestion:
                    print(f"    Fix: {err.fix_suggestion}")
        
        if self.warnings:
            print(f"\n⚠️  WARNINGS ({len(self.warnings)}):")
            for warn in self.warnings:
                if warn.path:
                    print(f"  • {warn.path}: {warn.message}")
                else:
                    print(f"  • {warn.message}")
                if warn.fix_suggestion:
                    print(f"    Suggestion: {warn.fix_suggestion}")
        
        if show_info and self.info:
            print(f"\nℹ️  INFO ({len(self.info)}):")
            for info_msg in self.info:
                if info_msg.path:
                    print(f"  • {info_msg.path}: {info_msg.message}")
                else:
                    print(f"  • {info_msg.message}")
        
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
        
        # Handle array format (wrap in container)
        if isinstance(data, list):
            # Multiple protocols/rules in array - validate first entry for source_type
            if len(data) == 0:
                report.add_error("Empty array in file")
                return report
            
            # Validate each item in the array
            all_valid = True
            for idx, item in enumerate(data):
                item_report = ValidationReport(file_path, self.schema_path)
                
                # Detect source type
                source_type = item.get('source_type')
                item_report.source_type = source_type
                
                # Run all validation levels on this item
                self._validate_schema(item, item_report)
                
                if RDKIT_AVAILABLE:
                    self._validate_chemistry(item, item_report)
                
                self._validate_semantics(item, item_report)
                self._quality_checks(item, item_report)
                
                # Merge results into main report
                if item_report.errors:
                    for error in item_report.errors:
                        report.add_error(f"[Item {idx}] {error.message}", error.path, error.fix_suggestion)
                    all_valid = False
                
                if item_report.warnings:
                    for warning in item_report.warnings:
                        report.add_warning(f"[Item {idx}] {warning.message}", warning.path, warning.fix_suggestion)
                
                # Use first item's source_type for report
                if idx == 0:
                    report.source_type = source_type
            
            # Determine overall validity
            if strict:
                report.is_valid = all_valid and len(report.warnings) == 0
            else:
                report.is_valid = all_valid
            
            return report
        
        # Single item format
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
            path = '.'.join(str(p) for p in e.path) if e.path else ""
            report.add_error(
                f"Schema validation failed: {e.message}",
                path=path,
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
            else:
                # Try to parse reaction SMILES
                try:
                    parts = smiles.split('>>')
                    if len(parts) != 2:
                        report.add_error(
                            "reaction_smiles must have exactly one '>>' separator",
                            path="reaction.reaction_smiles"
                        )
                    else:
                        # Parse reactants and products
                        reactants, products = parts
                        for mol_smiles in reactants.split('.'):
                            if mol_smiles.strip():
                                mol = Chem.MolFromSmiles(mol_smiles.strip())
                                if mol is None:
                                    report.add_error(
                                        f"Invalid reactant SMILES: {mol_smiles}",
                                        path="reaction.reaction_smiles"
                                    )
                        for mol_smiles in products.split('.'):
                            if mol_smiles.strip():
                                mol = Chem.MolFromSmiles(mol_smiles.strip())
                                if mol is None:
                                    report.add_error(
                                        f"Invalid product SMILES: {mol_smiles}",
                                        path="reaction.reaction_smiles"
                                    )
                except Exception as e:
                    report.add_warning(
                        f"Could not parse reaction_smiles: {e}",
                        path="reaction.reaction_smiles"
                    )
        
        # Validate reaction_smarts
        if 'reaction_smarts' in reaction:
            for idx, smarts in enumerate(reaction['reaction_smarts']):
                if '>>' not in smarts:
                    report.add_error(
                        f"SMARTS pattern must contain '>>' separator",
                        path=f"reaction.reaction_smarts[{idx}]"
                    )
                else:
                    # Try to parse reaction SMARTS
                    try:
                        rxn = AllChem.ReactionFromSmarts(smarts)
                        if rxn is None:
                            report.add_error(
                                f"Invalid SMARTS pattern: {smarts}",
                                path=f"reaction.reaction_smarts[{idx}]"
                            )
                    except Exception as e:
                        report.add_warning(
                            f"Could not parse SMARTS: {e}",
                            path=f"reaction.reaction_smarts[{idx}]"
                        )
        
        # Validate reference_reactions
        if 'reference_reactions' in reaction:
            for idx, ref_smiles in enumerate(reaction['reference_reactions']):
                if '>>' not in ref_smiles:
                    report.add_error(
                        f"Reference reaction must contain '>>' separator",
                        path=f"reaction.reference_reactions[{idx}]"
                    )
                else:
                    # Basic SMILES validation
                    try:
                        parts = ref_smiles.split('>>')
                        if len(parts) != 2:
                            report.add_error(
                                "Reference reaction must have exactly one '>>'",
                                path=f"reaction.reference_reactions[{idx}]"
                            )
                    except Exception as e:
                        report.add_warning(
                            f"Could not validate reference reaction: {e}",
                            path=f"reaction.reference_reactions[{idx}]"
                        )
        
        # Validate chemical SMILES (for protocols)
        if data.get('source_type') == 'protocol':
            for setup_idx, setup in enumerate(data.get('reaction_setup', [])):
                for chem_idx, chem in enumerate(setup.get('chemicals', [])):
                    smiles = chem.get('smiles')
                    if smiles and smiles.strip():
                        mol = Chem.MolFromSmiles(smiles)
                        if mol is None:
                            report.add_error(
                                f"Invalid SMILES for {chem.get('name', 'chemical')}: {smiles}",
                                path=f"reaction_setup[{setup_idx}].chemicals[{chem_idx}].smiles"
                            )
    
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
        
        # Check for version format
        version = metadata.get('version', '')
        if version and not version.startswith('v'):
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
        
        # Check reaction_smarts format
        reaction = data.get('reaction', {})
        if 'reaction_smarts' in reaction:
            smarts_list = reaction['reaction_smarts']
            if not smarts_list:
                report.add_info(
                    "Empty reaction_smarts array - consider removing or populating",
                    path="reaction.reaction_smarts"
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
