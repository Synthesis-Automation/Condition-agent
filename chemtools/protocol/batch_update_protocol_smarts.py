"""
Batch processor for protocol JSON files - auto-generate reaction_smarts_applicability.

This tool scans all protocol JSON files in data/protocol_db/, extracts reaction_smiles,
generates chemistry-aware SMARTS patterns, and updates the files with the new patterns.

Usage:
    python -m chemtools.protocol.batch_update_protocol_smarts
    python -m chemtools.protocol.batch_update_protocol_smarts --dry-run
    python -m chemtools.protocol.batch_update_protocol_smarts --protocol-dir path/to/protocols
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

# Add parent directory to path if running as script (not as module)
if __name__ == "__main__" and __package__ is None:
    _script_dir = Path(__file__).resolve().parent
    _project_root = _script_dir.parent.parent
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))

from chemtools.util.smarts_builders import build_smarts_with_guards
from chemtools.util.substrate_classifier import classify_substrate

# RDKit for SMARTS validation (optional)
try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


@dataclass
class ProcessingResult:
    """Result of processing a single protocol file"""
    filename: str
    success: bool
    message: str
    old_pattern: Optional[Dict[str, Any]] = None
    new_pattern: Optional[Dict[str, Any]] = None


class ProtocolSmartsUpdater:
    """Batch updater for protocol JSON files"""
    
    def __init__(self, protocol_dir: Path, dry_run: bool = False):
        """Initialize updater.
        
        Args:
            protocol_dir: Directory containing protocol JSON files
            dry_run: If True, don't write changes, just report what would happen
        """
        self.protocol_dir = protocol_dir
        self.dry_run = dry_run
        self.results: List[ProcessingResult] = []
    
    def process_all_protocols(self) -> List[ProcessingResult]:
        """Process all protocol JSON files in the directory.
        
        Returns:
            List of ProcessingResult objects
        """
        if not self.protocol_dir.exists():
            print(f"Error: Directory not found: {self.protocol_dir}")
            return []
        
        json_files = list(self.protocol_dir.glob("*.json"))
        
        if not json_files:
            print(f"Warning: No JSON files found in {self.protocol_dir}")
            return []
        
        print(f"\nFound {len(json_files)} protocol files in {self.protocol_dir}")
        print("=" * 80)
        
        for json_file in sorted(json_files):
            result = self.process_protocol_file(json_file)
            self.results.append(result)
        
        return self.results
    
    def process_protocol_file(self, filepath: Path) -> ProcessingResult:
        """Process a single protocol JSON file.
        
        Args:
            filepath: Path to protocol JSON file
            
        Returns:
            ProcessingResult object
        """
        filename = filepath.name
        
        try:
            # Load protocol JSON
            with open(filepath, 'r', encoding='utf-8') as f:
                protocol = json.load(f)
            
            # Extract reaction SMILES
            reaction_smiles = protocol.get("reaction", {}).get("reaction_smiles")
            
            if not reaction_smiles:
                return ProcessingResult(
                    filename=filename,
                    success=False,
                    message="No reaction_smiles found in protocol"
                )
            
            print(f"\n[Processing] {filename}")
            print(f"   Reaction: {reaction_smiles}")
            
            # Analyze reactants to show substrate information
            if ">>" in reaction_smiles:
                parts = reaction_smiles.split(">>")
            elif ">" in reaction_smiles:
                parts = reaction_smiles.split(">")
            else:
                parts = [reaction_smiles, ""]
            
            reactants = parts[0].strip()
            reactant_smiles_list = [s.strip() for s in reactants.split('.') if s.strip()]
            
            # Show substrate classification for each reactant
            for i, smiles in enumerate(reactant_smiles_list, 1):
                try:
                    substrate_info = classify_substrate(smiles)
                    substrate_class = substrate_info.substrate_class
                    family = substrate_info.substrate_family
                    special_positions = substrate_info.special_positions
                    
                    # Get special positions if any
                    special_list = []
                    if hasattr(special_positions, 'benzylic') and special_positions.benzylic:
                        special_list.append('benzylic')
                    if hasattr(special_positions, 'allylic') and special_positions.allylic:
                        special_list.append('allylic')
                    if hasattr(special_positions, 'propargylic') and special_positions.propargylic:
                        special_list.append('propargylic')
                    
                    special_str = f" [{', '.join(special_list)}]" if special_list else ""
                    print(f"   Reactant {i}: {substrate_class} ({family}){special_str}")
                except Exception as e:
                    print(f"   Reactant {i}: {smiles[:40]}... - {str(e)[:50]}")
            
            # Check if pattern already exists
            old_pattern = protocol.get("reaction", {}).get("reaction_smarts_applicability")
            if old_pattern:
                print(f"   Info: Existing pattern found (will be overwritten)")
            
            # Generate new SMARTS pattern
            try:
                new_pattern = self.generate_smarts_pattern(reaction_smiles)
                
                # Validate the generated pattern
                validation_result = self.validate_smarts_pattern(new_pattern)
                if not validation_result['valid']:
                    return ProcessingResult(
                        filename=filename,
                        success=False,
                        message=f"Generated pattern is invalid: {validation_result['error']}",
                        old_pattern=old_pattern,
                        new_pattern=new_pattern
                    )
                
                # Print pattern details
                print(f"   [OK] Generated pattern:")
                print(f"      Core: {new_pattern['core']}")
                if new_pattern.get('guards_forbid'):
                    print(f"      Guards: {len(new_pattern['guards_forbid'])} exclusion patterns")
                    for guard in new_pattern['guards_forbid'][:3]:  # Show first 3
                        print(f"        - {guard}")
                    if len(new_pattern['guards_forbid']) > 3:
                        print(f"        ... and {len(new_pattern['guards_forbid']) - 3} more")
                        
            except Exception as e:
                return ProcessingResult(
                    filename=filename,
                    success=False,
                    message=f"Pattern generation failed: {e}",
                    old_pattern=old_pattern
                )
            
            # Update protocol with new pattern
            if "reaction" not in protocol:
                protocol["reaction"] = {}
            
            protocol["reaction"]["reaction_smarts_applicability"] = new_pattern
            
            # Write updated protocol (unless dry-run)
            if not self.dry_run:
                with open(filepath, 'w', encoding='utf-8') as f:
                    json.dump(protocol, f, indent=2, ensure_ascii=False)
                print(f"   [SUCCESS] Updated successfully")
            else:
                print(f"   [DRY RUN] Would update file")
            
            return ProcessingResult(
                filename=filename,
                success=True,
                message="Pattern generated and applied",
                old_pattern=old_pattern,
                new_pattern=new_pattern
            )
            
        except json.JSONDecodeError as e:
            return ProcessingResult(
                filename=filename,
                success=False,
                message=f"JSON parsing error: {e}"
            )
        except Exception as e:
            return ProcessingResult(
                filename=filename,
                success=False,
                message=f"Unexpected error: {e}"
            )
    
    def generate_smarts_pattern(self, reaction_smiles: str) -> Dict[str, Any]:
        """Generate SMARTS pattern from reaction SMILES focusing on reaction centers.
        
        Args:
            reaction_smiles: Reaction SMILES string
            
        Returns:
            Dictionary with 'core' and 'guards_forbid' keys
        """
        # Parse reaction SMILES
        if ">>" in reaction_smiles:
            parts = reaction_smiles.split(">>")
        elif ">" in reaction_smiles:
            parts = reaction_smiles.split(">")
        else:
            raise ValueError(f"Invalid reaction SMILES format: {reaction_smiles}")
        
        reactants = parts[0].strip()
        products = parts[-1].strip()
        
        # Generate patterns for each reactant with focus on reaction centers
        reactant_patterns = []
        reactant_smiles_list = [s.strip() for s in reactants.split('.') if s.strip()]
        reactant_info_list = []
        
        for smiles in reactant_smiles_list:
            try:
                # Get substrate classification to identify reaction center
                substrate_info = classify_substrate(smiles)
                reactant_info_list.append(substrate_info)
                
                # Generate pattern focused on reaction center
                result = self.build_reaction_center_pattern(smiles, substrate_info)
                reactant_patterns.append(result)
            except Exception as e:
                print(f"      Warning: Could not generate pattern for {smiles}: {e}")
                continue
        
        # Generate pattern for product (focus on main product)
        product_smiles_list = [s.strip() for s in products.split('.') if s.strip()]
        product_pattern = None
        if product_smiles_list:
            try:
                # Classify product to focus on reaction center
                product_info = classify_substrate(product_smiles_list[0])
                product_pattern = self.build_product_pattern(product_smiles_list[0], product_info, reactant_info_list)
            except Exception as e:
                print(f"      Warning: Could not generate product pattern: {e}")
        
        # Combine into reaction SMARTS focusing on reaction centers
        if reactant_patterns:
            # Combine reactant patterns with dots
            reactant_cores = [p['core'] for p in reactant_patterns if p.get('core')]
            core_pattern = '.'.join(reactant_cores) if len(reactant_cores) > 1 else reactant_cores[0]
            
            # Add product pattern if available
            if product_pattern:
                core_pattern = f"{core_pattern}>>{product_pattern}"
            
            # Collect all guard patterns
            all_guards = []
            for rp in reactant_patterns:
                all_guards.extend(rp.get('guards_forbid', []))
            
            # Remove duplicates while preserving order
            unique_guards = []
            seen = set()
            for guard in all_guards:
                guard_clean = guard.split('#')[0].strip()  # Remove comments for comparison
                if guard_clean not in seen:
                    seen.add(guard_clean)
                    unique_guards.append(guard)
            
            return {
                "core": core_pattern,
                "guards_forbid": unique_guards
            }
        else:
            raise ValueError("No patterns could be generated from reactants")
    
    def build_reaction_center_pattern(self, smiles: str, substrate_info) -> Dict[str, Any]:
        """Build pattern focused on the reaction center, not spectator groups.
        
        Args:
            smiles: SMILES string
            substrate_info: SubstrateInfo from classification
            
        Returns:
            Dictionary with 'core' and 'guards_forbid'
        """
        substrate_class = substrate_info.substrate_class
        family = substrate_info.substrate_family
        
        # For halides: focus on C-X bond
        if family == 'halide':
            if 'aryl' in substrate_class or 'heteroaryl' in substrate_class:
                # Aryl halide: focus on aromatic-halogen bond
                if 'iodide' in substrate_class:
                    return {
                        'core': 'c-[I]',
                        'guards_forbid': ['[CX4]-[I]  # Exclude aliphatic halides']
                    }
                elif 'bromide' in substrate_class:
                    return {
                        'core': 'c-[Br]',
                        'guards_forbid': ['[CX4]-[Br]  # Exclude aliphatic halides']
                    }
                elif 'chloride' in substrate_class:
                    return {
                        'core': 'c-[Cl]',
                        'guards_forbid': ['[CX4]-[Cl]  # Exclude aliphatic halides']
                    }
            elif 'alkyl' in substrate_class:
                # Alkyl halide: use existing detailed pattern
                result = build_smarts_with_guards(smiles)
                return result
        
        # For alkynes: focus on terminal alkyne
        if 'alkyne' in substrate_class or any('alkyne' in fg for fg in substrate_info.functional_groups):
            return {
                'core': 'C#C',
                'guards_forbid': []
            }
        
        # For amines: focus on N, not neighboring groups
        if family == 'amine':
            if 'aniline' in substrate_class:
                return {
                    'core': 'c-[NX3;H2;!$(NC=O)]',
                    'guards_forbid': ['[CX4]-[NX3;H2;!$(NC=O)]  # Exclude aliphatic primary amines']
                }
            else:
                # Use existing detailed pattern for aliphatic amines
                result = build_smarts_with_guards(smiles)
                return result
        
        # For boronic acids/esters: focus on B
        if family == 'boron':
            if 'boronic_ester' in substrate_class:
                return {
                    'core': '[B]1OC(C)(C)C(C)(C)O1',  # Pinacol ester pattern
                    'guards_forbid': []
                }
            elif 'boronic_acid' in substrate_class:
                return {
                    'core': '[B](O)O',
                    'guards_forbid': []
                }
            else:
                # Generic boron pattern
                return {
                    'core': '[B]',
                    'guards_forbid': []
                }
        
        # For other cases: fall back to detailed pattern
        result = build_smarts_with_guards(smiles)
        return result
    
    def build_product_pattern(self, smiles: str, product_info, reactant_info_list) -> str:
        """Build product pattern focusing on newly formed bonds.
        
        Args:
            smiles: Product SMILES
            product_info: SubstrateInfo for product
            reactant_info_list: List of SubstrateInfo for reactants
            
        Returns:
            SMARTS pattern string for product
        """
        # Identify what reaction type based on reactants
        reactant_families = [info.substrate_family for info in reactant_info_list]
        reactant_classes = [info.substrate_class for info in reactant_info_list]
        
        # Sonogashira: aryl halide + alkyne -> aryl alkyne
        if 'halide' in reactant_families and any('alkyne' in str(r) or any('alkyne' in fg for fg in info.functional_groups) for r, info in zip(reactant_classes, reactant_info_list)):
            return 'c-C#C'  # Aryl-alkyne bond
        
        # Suzuki: aryl halide + boronic acid/ester -> biaryl
        if 'halide' in reactant_families and 'boron' in reactant_families:
            # Check if product has biaryl
            if product_info.has_aromatic:
                return 'c-c'  # Biaryl bond
            else:
                return 'c-C'  # Aryl-alkyl bond
        
        # Buchwald-Hartwig: aryl halide + amine -> aryl amine
        if 'halide' in reactant_families and 'amine' in reactant_families:
            if 'aniline' in any(str(r) for r in reactant_classes):
                return 'c-[NX3;H1]'  # Aryl-NH-aryl (secondary amine)
            else:
                return 'c-[NX3]'  # Aryl-amine bond
        
        # Borylation: alkyl halide + boron -> alkyl-B
        if 'halide' in reactant_families and 'boron' in reactant_families:
            if any('alkyl' in str(r) for r in reactant_classes):
                if 'boronic_ester_pinacol' in str(product_info.substrate_class):
                    return '[CX4]-[B]1OC(C)(C)C(C)(C)O1'
                else:
                    return '[CX4]-[B]'
        
        # Generic: use first reactant's pattern with modifications
        # This is a fallback - may need refinement
        result = build_smarts_with_guards(smiles)
        return result['core']
    
    def validate_smarts_pattern(self, pattern: Dict[str, Any]) -> Dict[str, Any]:
        """Validate a SMARTS pattern.
        
        Args:
            pattern: Dictionary with 'core' and 'guards_forbid' keys
            
        Returns:
            Dictionary with 'valid' (bool) and optional 'error' (str) keys
        """
        if not HAS_RDKIT:
            # Can't validate without RDKit, assume valid
            return {'valid': True, 'warning': 'RDKit not available, validation skipped'}
        
        try:
            # Validate core pattern
            core = pattern.get('core', '')
            if not core:
                return {'valid': False, 'error': 'Core pattern is empty'}
            
            # Split reaction SMARTS if present
            if '>>' in core:
                reactant_part, product_part = core.split('>>')
                # Validate reactant part
                if reactant_part.strip():
                    mol = Chem.MolFromSmarts(reactant_part.strip())
                    if mol is None:
                        return {'valid': False, 'error': f'Invalid reactant SMARTS: {reactant_part}'}
                # Validate product part
                if product_part.strip():
                    mol = Chem.MolFromSmarts(product_part.strip())
                    if mol is None:
                        return {'valid': False, 'error': f'Invalid product SMARTS: {product_part}'}
            else:
                # Single SMARTS pattern
                mol = Chem.MolFromSmarts(core)
                if mol is None:
                    return {'valid': False, 'error': f'Invalid SMARTS: {core}'}
            
            # Validate guard patterns
            guards = pattern.get('guards_forbid', [])
            for i, guard in enumerate(guards):
                # Remove comments
                guard_clean = guard.split('#')[0].strip()
                if guard_clean:
                    mol = Chem.MolFromSmarts(guard_clean)
                    if mol is None:
                        return {'valid': False, 'error': f'Invalid guard pattern {i+1}: {guard_clean}'}
            
            return {'valid': True}
            
        except Exception as e:
            return {'valid': False, 'error': f'Validation error: {str(e)}'}
    
    def print_summary(self):
        """Print summary of processing results"""
        print("\n" + "=" * 80)
        print("PROCESSING SUMMARY")
        print("=" * 80)
        
        successful = [r for r in self.results if r.success]
        failed = [r for r in self.results if not r.success]
        
        print(f"\n[OK] Successful: {len(successful)}/{len(self.results)}")
        print(f"[FAIL] Failed: {len(failed)}/{len(self.results)}")
        
        if failed:
            print(f"\nFailed files:")
            for result in failed:
                print(f"   - {result.filename}: {result.message}")
        
        if successful:
            print(f"\nSuccessfully processed files:")
            for result in successful:
                status = "would be updated" if self.dry_run else "updated"
                had_existing = "(replaced existing)" if result.old_pattern else "(new pattern)"
                print(f"   - {result.filename} - {status} {had_existing}")
        
        if self.dry_run:
            print("\n[DRY RUN] No files were modified")
            print("   Run without --dry-run to apply changes")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        prog="batch-update-protocol-smarts",
        description="Batch generate and update reaction_smarts_applicability in protocol JSON files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all protocols in default directory (data/protocol_db)
  python -m chemtools.protocol.batch_update_protocol_smarts
  
  # Dry run - see what would change without modifying files
  python -m chemtools.protocol.batch_update_protocol_smarts --dry-run
  
  # Process protocols in custom directory
  python -m chemtools.protocol.batch_update_protocol_smarts --protocol-dir path/to/protocols
  
  # Dry run with custom directory
  python -m chemtools.protocol.batch_update_protocol_smarts --protocol-dir path/to/protocols --dry-run

Output format:
  Each protocol JSON file will have reaction_smarts_applicability added/updated:
  
  "reaction": {
    "reaction_smiles": "...",
    "reaction_smarts_applicability": {
      "core": "[CX4;H2,H3]-[I]>>[B]1OC(C)(C)C(C)(C)O1",
      "guards_forbid": [
        "[CX4;H1]-[I]  # Exclude secondary",
        "[CX4;H0]-[I]  # Exclude tertiary",
        ...
      ]
    }
  }
"""
    )
    
    parser.add_argument(
        "--protocol-dir",
        type=Path,
        default=Path("data/protocol_db"),
        help="Directory containing protocol JSON files (default: data/protocol_db)"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Dry run mode - show what would change without modifying files"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output with detailed pattern information"
    )
    
    args = parser.parse_args()
    
    # Show header
    print("\n" + "=" * 80)
    print("Protocol SMARTS Batch Updater")
    print("=" * 80)
    
    if args.dry_run:
        print("\nRunning in DRY RUN mode - no files will be modified")
    
    # Create updater and process
    updater = ProtocolSmartsUpdater(
        protocol_dir=args.protocol_dir,
        dry_run=args.dry_run
    )
    
    try:
        results = updater.process_all_protocols()
        updater.print_summary()
        
        # Exit with error code if any failed
        failed_count = len([r for r in results if not r.success])
        return 0 if failed_count == 0 else 1
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted by user")
        return 1
    except Exception as e:
        print(f"\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
