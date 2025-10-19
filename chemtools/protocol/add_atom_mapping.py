"""
Batch tool to add atom-mapped reaction SMILES to protocol JSON files.

Uses RXNMapper (IBM's Transformer-based attention-guided atom mapper) to
automatically generate atom mappings for reactions in protocol files.

Usage:
    python -m chemtools.protocol.add_atom_mapping
    python -m chemtools.protocol.add_atom_mapping --dry-run
    python -m chemtools.protocol.add_atom_mapping --protocol-dir path/to/protocols
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

# Add parent directory to path if running as script
if __name__ == "__main__" and __package__ is None:
    _script_dir = Path(__file__).resolve().parent
    _project_root = _script_dir.parent.parent
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))

# Try to import RXNMapper
try:
    from rxnmapper import RXNMapper
    HAS_RXNMAPPER = True
except ImportError:
    HAS_RXNMAPPER = False


@dataclass
class MappingResult:
    """Result of mapping a single protocol file"""
    filename: str
    success: bool
    message: str
    unmapped_smiles: Optional[str] = None
    mapped_smiles: Optional[str] = None
    confidence: Optional[float] = None


class ProtocolAtomMapper:
    """Batch atom mapper for protocol JSON files"""
    
    def __init__(self, protocol_dir: Path, dry_run: bool = False, force: bool = False):
        """Initialize mapper.
        
        Args:
            protocol_dir: Directory containing protocol JSON files
            dry_run: If True, don't write changes, just report
            force: If True, regenerate mappings even if they exist
        """
        self.protocol_dir = protocol_dir
        self.dry_run = dry_run
        self.force = force
        self.results: List[MappingResult] = []
        self.mapper = None
        
        if HAS_RXNMAPPER:
            print("Initializing RXNMapper (this may take a moment)...")
            self.mapper = RXNMapper()
            print("[OK] RXNMapper ready")
        else:
            print("[ERROR] RXNMapper not installed!")
            print("        Install with: pip install rxnmapper")
            sys.exit(1)
    
    def process_all_protocols(self) -> List[MappingResult]:
        """Process all protocol JSON files in directory.
        
        Returns:
            List of MappingResult objects
        """
        if not self.protocol_dir.exists():
            print(f"[ERROR] Directory not found: {self.protocol_dir}")
            return []
        
        json_files = list(self.protocol_dir.glob("*.json"))
        
        if not json_files:
            print(f"[WARNING] No JSON files found in {self.protocol_dir}")
            return []
        
        print(f"\nFound {len(json_files)} protocol files in {self.protocol_dir}")
        print("=" * 80)
        
        for json_file in sorted(json_files):
            result = self.process_protocol_file(json_file)
            self.results.append(result)
        
        return self.results
    
    def process_protocol_file(self, filepath: Path) -> MappingResult:
        """Process a single protocol JSON file.
        
        Args:
            filepath: Path to protocol JSON file
            
        Returns:
            MappingResult object
        """
        filename = filepath.name
        
        try:
            # Load protocol JSON
            with open(filepath, 'r', encoding='utf-8') as f:
                protocol = json.load(f)
            
            # Extract reaction SMILES
            reaction_smiles = protocol.get("reaction", {}).get("reaction_smiles")
            
            if not reaction_smiles:
                return MappingResult(
                    filename=filename,
                    success=False,
                    message="No reaction_smiles found in protocol"
                )
            
            print(f"\n[Processing] {filename}")
            print(f"   Reaction: {reaction_smiles[:80]}...")
            
            # Check if mapping already exists
            existing_mapping = protocol.get("reaction", {}).get("reaction_smiles_mapped")
            if existing_mapping and not self.force:
                print(f"   [SKIP] Mapping already exists (use --force to regenerate)")
                return MappingResult(
                    filename=filename,
                    success=True,
                    message="Mapping already exists",
                    unmapped_smiles=reaction_smiles,
                    mapped_smiles=existing_mapping
                )
            
            # Generate atom mapping
            try:
                mapped_result = self.generate_atom_mapping(reaction_smiles)
                
                if not mapped_result['success']:
                    print(f"   [FAIL] {mapped_result['message']}")
                    return MappingResult(
                        filename=filename,
                        success=False,
                        message=mapped_result['message'],
                        unmapped_smiles=reaction_smiles
                    )
                
                mapped_smiles = mapped_result['mapped_rxn']
                confidence = mapped_result.get('confidence', 0.0)
                
                print(f"   [OK] Mapped with confidence: {confidence:.1%}")
                print(f"   Mapped: {mapped_smiles[:80]}...")
                
                # Update protocol with mapping
                if "reaction" not in protocol:
                    protocol["reaction"] = {}
                
                protocol["reaction"]["reaction_smiles_mapped"] = mapped_smiles
                protocol["reaction"]["mapping_confidence"] = confidence
                protocol["reaction"]["mapping_tool"] = "RXNMapper-0.4.2"
                
                # Write updated protocol (unless dry-run)
                if not self.dry_run:
                    with open(filepath, 'w', encoding='utf-8') as f:
                        json.dump(protocol, f, indent=2, ensure_ascii=False)
                    print(f"   [SUCCESS] Updated file")
                else:
                    print(f"   [DRY RUN] Would update file")
                
                return MappingResult(
                    filename=filename,
                    success=True,
                    message="Mapping generated successfully",
                    unmapped_smiles=reaction_smiles,
                    mapped_smiles=mapped_smiles,
                    confidence=confidence
                )
                
            except Exception as e:
                print(f"   [ERROR] Mapping generation failed: {e}")
                return MappingResult(
                    filename=filename,
                    success=False,
                    message=f"Mapping generation error: {e}",
                    unmapped_smiles=reaction_smiles
                )
            
        except json.JSONDecodeError as e:
            return MappingResult(
                filename=filename,
                success=False,
                message=f"JSON parsing error: {e}"
            )
        except Exception as e:
            return MappingResult(
                filename=filename,
                success=False,
                message=f"Unexpected error: {e}"
            )
    
    def generate_atom_mapping(self, reaction_smiles: str) -> Dict[str, Any]:
        """Generate atom mapping for a reaction.
        
        Args:
            reaction_smiles: Unmapped reaction SMILES
            
        Returns:
            Dictionary with 'success', 'mapped_rxn', and 'confidence'
        """
        try:
            # RXNMapper expects a list of reactions
            results = self.mapper.get_attention_guided_atom_maps([reaction_smiles])
            
            if not results or len(results) == 0:
                return {
                    'success': False,
                    'message': 'RXNMapper returned no results'
                }
            
            result = results[0]
            mapped_rxn = result.get('mapped_rxn', '')
            confidence = result.get('confidence', 0.0)
            
            if not mapped_rxn:
                return {
                    'success': False,
                    'message': 'RXNMapper returned empty mapping'
                }
            
            return {
                'success': True,
                'mapped_rxn': mapped_rxn,
                'confidence': confidence
            }
            
        except Exception as e:
            return {
                'success': False,
                'message': f'RXNMapper error: {str(e)}'
            }
    
    def print_summary(self):
        """Print summary of processing results"""
        print("\n" + "=" * 80)
        print("ATOM MAPPING SUMMARY")
        print("=" * 80)
        
        successful = [r for r in self.results if r.success and r.mapped_smiles]
        skipped = [r for r in self.results if r.success and not r.mapped_smiles]
        failed = [r for r in self.results if not r.success]
        
        print(f"\n[OK] Successfully mapped: {len(successful)}/{len(self.results)}")
        print(f"[SKIP] Already mapped: {len(skipped)}/{len(self.results)}")
        print(f"[FAIL] Failed: {len(failed)}/{len(self.results)}")
        
        if successful:
            # Show confidence statistics
            confidences = [r.confidence for r in successful if r.confidence]
            if confidences:
                avg_conf = sum(confidences) / len(confidences)
                min_conf = min(confidences)
                max_conf = max(confidences)
                print(f"\nMapping Confidence:")
                print(f"   Average: {avg_conf:.1%}")
                print(f"   Range: {min_conf:.1%} - {max_conf:.1%}")
                
                # Warn about low confidence mappings
                low_conf = [r for r in successful if r.confidence and r.confidence < 0.7]
                if low_conf:
                    print(f"\n[WARNING] {len(low_conf)} mappings with confidence < 70%:")
                    for r in low_conf[:5]:  # Show first 5
                        print(f"   - {r.filename}: {r.confidence:.1%}")
                    if len(low_conf) > 5:
                        print(f"   ... and {len(low_conf) - 5} more")
        
        if failed:
            print(f"\n[FAIL] Failed files:")
            for result in failed:
                print(f"   - {result.filename}: {result.message}")
        
        if self.dry_run:
            print("\n[DRY RUN] No files were modified")
            print("   Run without --dry-run to apply changes")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        prog="add-atom-mapping",
        description="Batch add atom-mapped reaction SMILES to protocol JSON files using RXNMapper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all protocols in default directory (data/protocol_db)
  python -m chemtools.protocol.add_atom_mapping
  
  # Dry run - see what would change without modifying files
  python -m chemtools.protocol.add_atom_mapping --dry-run
  
  # Process protocols in custom directory
  python -m chemtools.protocol.add_atom_mapping --protocol-dir path/to/protocols
  
  # Force regenerate mappings even if they exist
  python -m chemtools.protocol.add_atom_mapping --force
  
Output:
  Each protocol JSON file will have atom-mapped SMILES added:
  
  "reaction": {
    "reaction_smiles": "NC(=O)c1ccccc1I.C#CCCCC>>NC(=O)c1ccccc1C#CCCCC",
    "reaction_smiles_mapped": "[NH2:7][C:8](=[O:9])[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[I:10].[C:11]#[C:12][CH2:13][CH2:14][CH2:15][CH3:16]>>[NH2:7][C:8](=[O:9])[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[C:11]#[C:12][CH2:13][CH2:14][CH2:15][CH3:16]",
    "mapping_confidence": 0.95,
    "mapping_tool": "RXNMapper-0.4.2"
  }

Note:
  - RXNMapper uses a Transformer model (first run downloads ~500MB model)
  - Mappings typically have >90% confidence for common reactions
  - Low confidence (<70%) may indicate unusual reactions or errors
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
        "--force",
        action="store_true",
        help="Force regenerate mappings even if they already exist"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output with detailed mapping information"
    )
    
    args = parser.parse_args()
    
    # Show header
    print("\n" + "=" * 80)
    print("Protocol Atom Mapping Tool (RXNMapper)")
    print("=" * 80)
    
    if not HAS_RXNMAPPER:
        print("\n[ERROR] RXNMapper not installed!")
        print("        Install with: pip install rxnmapper")
        return 1
    
    if args.dry_run:
        print("\nRunning in DRY RUN mode - no files will be modified")
    
    if args.force:
        print("Force mode: Will regenerate existing mappings")
    
    # Create mapper and process
    mapper = ProtocolAtomMapper(
        protocol_dir=args.protocol_dir,
        dry_run=args.dry_run,
        force=args.force
    )
    
    try:
        results = mapper.process_all_protocols()
        mapper.print_summary()
        
        # Exit with error code if any failed
        failed_count = len([r for r in results if not r.success])
        return 0 if failed_count == 0 else 1
        
    except KeyboardInterrupt:
        print("\n\n[INTERRUPTED] Stopped by user")
        return 1
    except Exception as e:
        print(f"\n[FATAL ERROR] {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
