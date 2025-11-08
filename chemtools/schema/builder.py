"""
Build system for creating unified index from protocols and rules.

Usage:
    # Build unified index
    python -m chemtools.schema.builder build \
        --protocols data/protocol_db \
        --rules data/rule_db \
        --output data/unified_index_v1.0
    
    # Validate before building
    python -m chemtools.schema.builder validate \
        --protocols data/protocol_db \
        --rules data/rule_db
"""

from pathlib import Path
from typing import List, Dict, Any, Optional
import json
import numpy as np
from datetime import datetime
from dataclasses import dataclass, field

from .validator import ConditionSourceValidator, ValidationReport

try:
    from drfp import DrfpEncoder
    DRFP_AVAILABLE = True
except ImportError:
    DRFP_AVAILABLE = False


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
        timestamp = datetime.now().strftime('%H:%M:%S')
        self.messages.append(f"[{timestamp}] {msg}")
        print(msg)


class UnifiedIndexBuilder:
    """
    Builds unified index from protocol and rule database files.
    """
    
    def __init__(self, config: BuildConfig):
        self.config = config
        self.validator = ConditionSourceValidator()
        self.report = BuildReport(version=config.version)
        
        if DRFP_AVAILABLE:
            self.drfp_encoder = DrfpEncoder()
        else:
            self.drfp_encoder = None
    
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
            if not DRFP_AVAILABLE:
                self.report.add_message("  ⚠️  DRFP library not available - skipping fingerprint computation")
                protocol_drfps = {}
                rule_drfps = {}
            else:
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
            import traceback
            traceback.print_exc()
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
                # Handle both v1 (array) and v2 (object) formats
                if isinstance(data, list):
                    # v1 format - take first item
                    if data:
                        data = data[0]
                        data['_source_file'] = str(file_path)
                        protocols.append(data)
                else:
                    # v2 format
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
        drfps = {}
        
        for protocol in protocols:
            protocol_id = protocol.get('metadata', {}).get('id', 'unknown')
            reaction = protocol.get('reaction', {})
            reaction_smiles = reaction.get('reaction_smiles')
            
            if not reaction_smiles:
                self.report.drfp_failed += 1
                continue
            
            try:
                # Encode reaction SMILES to DRFP
                fp = self.drfp_encoder.encode([reaction_smiles])[0]
                drfps[protocol_id] = fp
                self.report.drfp_computed += 1
            except Exception as e:
                self.report.add_message(f"  ⚠️  Failed to compute DRFP for {protocol_id}: {e}")
                self.report.drfp_failed += 1
        
        self.report.add_message(f"  Computed {self.report.drfp_computed} protocol DRFPs ({self.report.drfp_failed} failed)")
        return drfps
    
    def _compute_drfps_for_rules(self, rules: List[Dict]) -> Dict[str, List[np.ndarray]]:
        """Compute DRFP for each rule's reference reactions"""
        drfps = {}
        
        for rule in rules:
            rule_id = rule.get('metadata', {}).get('id', 'unknown')
            reaction = rule.get('reaction', {})
            reference_reactions = reaction.get('reference_reactions', [])
            
            if not reference_reactions:
                self.report.drfp_failed += 1
                continue
            
            rule_fps = []
            for ref_smiles in reference_reactions:
                try:
                    fp = self.drfp_encoder.encode([ref_smiles])[0]
                    rule_fps.append(fp)
                    self.report.drfp_computed += 1
                except Exception as e:
                    self.report.add_message(f"  ⚠️  Failed to compute DRFP for {rule_id} reference: {e}")
                    self.report.drfp_failed += 1
            
            if rule_fps:
                drfps[rule_id] = rule_fps
        
        self.report.add_message(f"  Computed {self.report.drfp_computed} rule reference DRFPs ({self.report.drfp_failed} failed)")
        return drfps
    
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
        
        index_path = output_dir / 'index.json'
        with open(index_path, 'w', encoding='utf-8') as f:
            json.dump(index, f, indent=2)
        
        # Save fingerprints.npz (DRFPs)
        if protocol_drfps or rule_drfps:
            fp_data = {
                'protocol_ids': list(protocol_drfps.keys()),
                'protocol_fps': np.array([protocol_drfps[k] for k in protocol_drfps.keys()]) if protocol_drfps else np.array([]),
                'rule_ids': list(rule_drfps.keys()),
                # For rules, we store all reference fingerprints flattened
                # with metadata to track which belong to which rule
            }
            
            # Flatten rule fingerprints
            rule_fps_flat = []
            rule_fp_indices = []
            for rule_id in rule_drfps.keys():
                start_idx = len(rule_fps_flat)
                rule_fps_flat.extend(rule_drfps[rule_id])
                end_idx = len(rule_fps_flat)
                rule_fp_indices.append((rule_id, start_idx, end_idx))
            
            fp_data['rule_fps'] = np.array(rule_fps_flat) if rule_fps_flat else np.array([])
            fp_data['rule_fp_indices'] = rule_fp_indices
            
            fp_path = output_dir / 'fingerprints.npz'
            np.savez_compressed(fp_path, **fp_data)
            
            # Calculate size
            self.report.index_size_mb = (index_path.stat().st_size + fp_path.stat().st_size) / (1024 * 1024)
        else:
            self.report.index_size_mb = index_path.stat().st_size / (1024 * 1024)
        
        self.report.add_message(f"  Index saved to {output_dir}")
        self.report.add_message(f"  Index size: {self.report.index_size_mb:.2f} MB")
    
    def _extract_metadata(self, data: Dict) -> Dict:
        """Extract indexable metadata from source"""
        metadata = data.get('metadata', {})
        reaction = data.get('reaction', {})
        
        return {
            'id': metadata.get('id', 'unknown'),
            'name': metadata.get('name', ''),
            'source_type': data.get('source_type', 'unknown'),
            'family': reaction.get('family', 'unknown'),
            'tags': metadata.get('tags', []),
            'version': metadata.get('version', ''),
            'source_file': data.get('_source_file', '')
        }
    
    def _generate_stats(self, protocols: List[Dict], rules: List[Dict]):
        """Generate statistics report"""
        stats = {
            'build_info': {
                'version': self.config.version,
                'build_date': datetime.now().isoformat(),
                'total_sources': len(protocols) + len(rules)
            },
            'protocols': {
                'count': len(protocols),
                'families': self._count_families(protocols)
            },
            'rules': {
                'count': len(rules),
                'families': self._count_families(rules)
            },
            'drfp': {
                'computed': self.report.drfp_computed,
                'failed': self.report.drfp_failed
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
            reaction = item.get('reaction', {})
            family = reaction.get('family', 'unknown')
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
