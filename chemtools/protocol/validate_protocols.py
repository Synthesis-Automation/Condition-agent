#!/usr/bin/env python3
"""
Protocol Validation CLI

Validates protocol JSON files to ensure reaction SMILES match their declared SMARTS patterns.

Usage:
    # Validate all protocols
    python -m chemtools.protocol.validate_protocols
    
    # Validate specific file
    python -m chemtools.protocol.validate_protocols --file "Aryl mesylate_Suzuki.json"
    
    # Show only errors
    python -m chemtools.protocol.validate_protocols --errors-only
    
    # Verbose output with details
    python -m chemtools.protocol.validate_protocols --verbose
    
    # Export validation report to JSON
    python -m chemtools.protocol.validate_protocols --output validation_report.json
"""

import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, asdict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    """Result of validating a single protocol"""
    filename: str
    valid: bool
    reaction_smiles: Optional[str] = None
    reaction_smarts: Optional[List[str]] = None
    errors: List[str] = None
    warnings: List[str] = None
    
    def __post_init__(self):
        if self.errors is None:
            self.errors = []
        if self.warnings is None:
            self.warnings = []
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def get_protocol_dir() -> Path:
    """Get the default protocol database directory"""
    # Assume we're in chemtools/protocol/, so go up to project root
    current_file = Path(__file__)
    project_root = current_file.parent.parent.parent
    protocol_dir = project_root / "data" / "protocol_db"
    
    if not protocol_dir.exists():
        raise FileNotFoundError(f"Protocol directory not found: {protocol_dir}")
    
    return protocol_dir


def load_protocol(filepath: Path) -> Optional[Dict[str, Any]]:
    """Load a protocol JSON file"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in {filepath.name}: {e}")
        return None
    except Exception as e:
        logger.error(f"Error loading {filepath.name}: {e}")
        return None


def match_reaction_smarts(reaction_smiles: str, smarts_patterns: List[str]) -> Tuple[bool, List[str]]:
    """
    Check if a reaction SMILES matches any of the provided SMARTS patterns.
    
    Returns:
        Tuple of (matches: bool, error_messages: List[str])
    """
    errors = []
    
    if not smarts_patterns:
        return True, []  # No patterns to validate against
    
    if not reaction_smiles or '>>' not in reaction_smiles:
        errors.append(f"Invalid reaction SMILES format: {reaction_smiles}")
        return False, errors
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        errors.append("RDKit not available - cannot validate SMARTS patterns")
        return False, errors
    
    # Parse the input reaction
    try:
        rxn_mol = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
        if rxn_mol is None:
            errors.append(f"Could not parse reaction SMILES: {reaction_smiles}")
            return False, errors
    except Exception as e:
        errors.append(f"Error parsing reaction SMILES: {e}")
        return False, errors
    
    # Try to match against each SMARTS pattern
    matched = False
    pattern_errors = []
    
    for pattern in smarts_patterns:
        if not pattern or '>>' not in pattern:
            pattern_errors.append(f"Invalid SMARTS pattern format: {pattern}")
            continue
        
        try:
            # Parse the SMARTS pattern as a reaction
            pattern_rxn = AllChem.ReactionFromSmarts(pattern)
            if pattern_rxn is None:
                pattern_errors.append(f"Could not parse reaction SMARTS: {pattern}")
                continue
            
            # Get reactants and products from both
            rxn_reactants = rxn_mol.GetReactants()
            rxn_products = rxn_mol.GetProducts()
            pattern_reactants = pattern_rxn.GetReactants()
            pattern_products = pattern_rxn.GetProducts()
            
            # Check if all pattern reactants match some rxn reactants
            reactants_match = all(
                any(rxn_r.HasSubstructMatch(pat_r) for rxn_r in rxn_reactants)
                for pat_r in pattern_reactants
            )
            
            # Check if all pattern products match some rxn products
            products_match = all(
                any(rxn_p.HasSubstructMatch(pat_p) for rxn_p in rxn_products)
                for pat_p in pattern_products
            )
            
            if reactants_match and products_match:
                matched = True
                break  # Found a matching pattern
                
        except Exception as e:
            pattern_errors.append(f"Error matching SMARTS pattern '{pattern}': {e}")
            continue
    
    if not matched:
        if pattern_errors:
            errors.extend(pattern_errors)
        else:
            errors.append(
                f"Reaction SMILES does not match any of the {len(smarts_patterns)} SMARTS pattern(s)"
            )
    
    return matched, errors


def validate_protocol(filepath: Path) -> ValidationResult:
    """Validate a single protocol file"""
    result = ValidationResult(
        filename=filepath.name,
        valid=False
    )
    
    # Load the protocol
    protocol = load_protocol(filepath)
    if protocol is None:
        result.errors.append("Failed to load protocol JSON")
        return result
    
    # Check structure
    if 'reaction' not in protocol:
        result.errors.append("Missing 'reaction' field")
        return result
    
    reaction = protocol['reaction']
    
    # Extract reaction_smiles
    reaction_smiles = reaction.get('reaction_smiles')
    if not reaction_smiles:
        result.errors.append("Missing 'reaction.reaction_smiles' field")
        return result
    
    result.reaction_smiles = reaction_smiles
    
    # Extract reaction_SMARTS
    reaction_smarts = reaction.get('reaction_SMARTS', [])
    if not isinstance(reaction_smarts, list):
        result.warnings.append("'reaction.reaction_SMARTS' should be a list")
        reaction_smarts = [reaction_smarts] if reaction_smarts else []
    
    result.reaction_smarts = reaction_smarts
    
    if not reaction_smarts:
        result.warnings.append("No reaction_SMARTS patterns defined")
        result.valid = True  # Not an error, just a warning
        return result
    
    # Validate SMARTS matching
    matched, errors = match_reaction_smarts(reaction_smiles, reaction_smarts)
    
    if matched:
        result.valid = True
    else:
        result.valid = False
        result.errors.extend(errors)
    
    return result


def validate_all_protocols(
    protocol_dir: Path,
    specific_file: Optional[str] = None
) -> List[ValidationResult]:
    """Validate all protocol files in the directory"""
    results = []
    
    if specific_file:
        # Validate single file
        filepath = protocol_dir / specific_file
        if not filepath.exists():
            logger.error(f"File not found: {specific_file}")
            return []
        
        result = validate_protocol(filepath)
        results.append(result)
    else:
        # Validate all JSON files
        json_files = sorted(protocol_dir.glob("*.json"))
        
        if not json_files:
            logger.warning(f"No JSON files found in {protocol_dir}")
            return []
        
        logger.info(f"Validating {len(json_files)} protocol files...")
        print()
        
        for filepath in json_files:
            result = validate_protocol(filepath)
            results.append(result)
    
    return results


def print_summary(results: List[ValidationResult], errors_only: bool = False, verbose: bool = False):
    """Print validation summary"""
    total = len(results)
    valid = sum(1 for r in results if r.valid)
    invalid = total - valid
    
    print("=" * 70)
    print("Protocol Validation Summary")
    print("=" * 70)
    print()
    print(f"Total protocols: {total}")
    print(f"✅ Valid: {valid}")
    print(f"❌ Invalid: {invalid}")
    print()
    
    if invalid > 0:
        print("=" * 70)
        print("Invalid Protocols:")
        print("=" * 70)
        print()
        
        for result in results:
            if not result.valid:
                print(f"❌ {result.filename}")
                print(f"   Reaction: {result.reaction_smiles}")
                print(f"   SMARTS patterns: {len(result.reaction_smarts or [])}")
                
                for error in result.errors:
                    print(f"   ERROR: {error}")
                
                if verbose and result.reaction_smarts:
                    print(f"   Patterns:")
                    for pattern in result.reaction_smarts:
                        print(f"     - {pattern}")
                
                print()
    
    # Show warnings if verbose
    if verbose or not errors_only:
        warnings_count = sum(len(r.warnings) for r in results if r.warnings)
        
        if warnings_count > 0:
            print("=" * 70)
            print("Warnings:")
            print("=" * 70)
            print()
            
            for result in results:
                if result.warnings:
                    print(f"⚠️  {result.filename}")
                    for warning in result.warnings:
                        print(f"   WARNING: {warning}")
                    print()
    
    # Show valid protocols if not errors_only
    if not errors_only and verbose:
        print("=" * 70)
        print("Valid Protocols:")
        print("=" * 70)
        print()
        
        for result in results:
            if result.valid:
                print(f"✅ {result.filename}")
                if verbose:
                    print(f"   Reaction: {result.reaction_smiles}")
                    if result.reaction_smarts:
                        print(f"   Matched {len(result.reaction_smarts)} SMARTS pattern(s)")
                    if result.warnings:
                        for warning in result.warnings:
                            print(f"   WARNING: {warning}")
                    print()


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Validate protocol JSON files for SMARTS pattern matching",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--file', '-f',
        help="Validate specific protocol file (filename only, not full path)"
    )
    
    parser.add_argument(
        '--protocol-dir',
        type=Path,
        help="Protocol directory (default: data/protocol_db)"
    )
    
    parser.add_argument(
        '--errors-only', '-e',
        action='store_true',
        help="Show only invalid protocols"
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Show detailed validation information"
    )
    
    parser.add_argument(
        '--output', '-o',
        type=Path,
        help="Export validation report to JSON file"
    )
    
    parser.add_argument(
        '--fail-on-error',
        action='store_true',
        help="Exit with code 1 if any protocols are invalid"
    )
    
    args = parser.parse_args(argv)
    
    # Get protocol directory
    try:
        protocol_dir = args.protocol_dir or get_protocol_dir()
    except FileNotFoundError as e:
        print(f"❌ Error: {e}")
        return 1
    
    # Validate protocols
    try:
        results = validate_all_protocols(protocol_dir, args.file)
    except Exception as e:
        print(f"❌ Validation failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    if not results:
        print("❌ No protocols to validate")
        return 1
    
    # Print summary
    print_summary(results, args.errors_only, args.verbose)
    
    # Export to JSON if requested
    if args.output:
        try:
            report = {
                'total': len(results),
                'valid': sum(1 for r in results if r.valid),
                'invalid': sum(1 for r in results if not r.valid),
                'results': [r.to_dict() for r in results]
            }
            
            with open(args.output, 'w', encoding='utf-8') as f:
                json.dump(report, f, indent=2, ensure_ascii=False)
            
            print(f"📄 Validation report saved to: {args.output}")
            print()
        except Exception as e:
            print(f"❌ Failed to save report: {e}")
            return 1
    
    # Exit with error code if requested
    if args.fail_on_error:
        invalid_count = sum(1 for r in results if not r.valid)
        if invalid_count > 0:
            return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
