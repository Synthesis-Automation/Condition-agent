"""
Demo: Schema Validation & Build System

This script demonstrates:
1. Validating protocol and rule files
2. Building unified index
3. Loading and using the index

Usage:
    python demo_schema_system.py
"""

from pathlib import Path
from chemtools.schema import ConditionSourceValidator, UnifiedIndexBuilder, BuildConfig
import json


def main():
    print("="*80)
    print("Schema Validation & Build System Demo")
    print("="*80)
    
    # 1. Validate protocol file
    print("\n[Step 1] Validating protocol file...")
    validator = ConditionSourceValidator()
    
    protocol_path = Path("data/protocol_db/alpha_arylation_dong_v100p0099_v2.json")
    protocol_report = validator.validate_file(protocol_path)
    
    print(f"  File: {protocol_path.name}")
    print(f"  Type: {protocol_report.source_type}")
    print(f"  Valid: {'✅ Yes' if protocol_report.is_valid else '❌ No'}")
    print(f"  Errors: {len(protocol_report.errors)}")
    print(f"  Warnings: {len(protocol_report.warnings)}")
    
    # 2. Validate rule file
    print("\n[Step 2] Validating rule file...")
    
    rule_path = Path("data/rule_db/sonogashira_v2.json")
    rule_report = validator.validate_file(rule_path)
    
    print(f"  File: {rule_path.name}")
    print(f"  Type: {rule_report.source_type}")
    print(f"  Valid: {'✅ Yes' if rule_report.is_valid else '❌ No'}")
    print(f"  Errors: {len(rule_report.errors)}")
    print(f"  Warnings: {len(rule_report.warnings)}")
    
    # 3. Build unified index
    print("\n[Step 3] Building unified index...")
    
    config = BuildConfig(
        protocol_dir=Path("data/protocol_db_v2"),
        rule_dir=Path("data/rule_db_v2"),
        output_dir=Path("build/unified_index_demo"),
        version="2.0"
    )
    
    builder = UnifiedIndexBuilder(config)
    build_report = builder.build()
    
    print(f"  Success: {'✅ Yes' if build_report.success else '❌ No'}")
    print(f"  Protocols: {build_report.num_protocols}")
    print(f"  Rules: {build_report.num_rules}")
    print(f"  DRFP computed: {build_report.drfp_computed}")
    print(f"  DRFP failed: {build_report.drfp_failed}")
    print(f"  Index size: {build_report.index_size_mb:.2f} MB")
    
    # 4. Load and inspect index
    print("\n[Step 4] Loading unified index...")
    
    index_path = Path("build/unified_index_demo/index.json")
    with open(index_path, 'r', encoding='utf-8') as f:
        index = json.load(f)
    
    print(f"  Version: {index['version']}")
    print(f"  Build date: {index['build_date']}")
    print(f"  Total sources: {index['num_protocols'] + index['num_rules']}")
    
    print("\n  Protocols:")
    for p in index['protocols']:
        print(f"    • {p['id']} ({p['family']})")
        print(f"      Tags: {', '.join(p['tags'][:3])}...")
    
    print("\n  Rules:")
    for r in index['rules']:
        print(f"    • {r['id']} ({r['family']})")
        print(f"      Tags: {', '.join(r['tags'][:3])}...")
    
    # 5. Load statistics
    print("\n[Step 5] Loading statistics...")
    
    stats_path = Path("build/unified_index_demo/stats.json")
    with open(stats_path, 'r', encoding='utf-8') as f:
        stats = json.load(f)
    
    print(f"  Protocol families: {stats['protocols']['families']}")
    print(f"  Rule families: {stats['rules']['families']}")
    print(f"  DRFP success rate: {stats['drfp']['computed']}/{stats['drfp']['computed'] + stats['drfp']['failed']}")
    
    print("\n" + "="*80)
    print("✅ Demo completed successfully!")
    print("="*80)
    print("\nNext steps:")
    print("  1. Migrate more databases to v2.0 format")
    print("  2. Integrate unified index with recommender system")
    print("  3. Set up CI/CD for automatic validation")
    print()


if __name__ == '__main__':
    main()
