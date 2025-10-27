#!/usr/bin/env python3
"""
Generate comprehensive summary of all reaction datasets.
Run with: python generate_dataset_summary.py
"""

import sys
from pathlib import Path

# Add parent directory to path so chemtools can be imported
ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from chemtools.dataset_analytics import (
    get_all_families,
    get_dataset_stats,
    get_common_catalysts,
    get_common_ligands,
    get_common_bases,
    get_common_solvents,
    get_common_reagents,
    get_condition_cores,
    get_common_catalytic_systems
)
import json


def generate_summary():
    """Generate and print comprehensive dataset summary."""
    
    families = get_all_families()
    
    print("=" * 80)
    print(" REACTION DATASET SUMMARY")
    print(" Updated: October 21, 2025")
    print("=" * 80)
    print()
    
    # Overview
    print("[DATASET OVERVIEW]")
    print()
    
    total_reactions = 0
    family_stats = {}
    
    for family in families:
        stats = get_dataset_stats(family)
        family_stats[family] = stats
        total_reactions += stats['total_reactions']
        
        print(f"📊 {family}")
        print(f"   Total reactions: {stats['total_reactions']:,}")
        print(f"   Unique condition cores: {stats['unique_condition_cores']}")
        print(f"   Unique catalysts: {stats['unique_catalysts']}")
        print(f"   Unique ligands: (extracted from catalytic systems)")
        print(f"   Unique bases: {stats['unique_bases']}")
        print(f"   Unique solvents: {stats['unique_solvents']}")
        
        if stats['yield_stats']:
            ys = stats['yield_stats']
            coverage = ys['count'] / stats['total_reactions'] * 100
            print(f"   Yield data: {ys['count']:,}/{stats['total_reactions']:,} ({coverage:.1f}%)")
            print(f"   Yield range: {ys['min']:.1f}% - {ys['max']:.1f}% (mean: {ys['mean']:.1f}%, median: {ys['median']:.1f}%)")
        
        if stats['temperature_stats']:
            ts = stats['temperature_stats']
            coverage = ts['count'] / stats['total_reactions'] * 100
            print(f"   Temperature data: {ts['count']:,}/{stats['total_reactions']:,} ({coverage:.1f}%)")
            print(f"   Temperature range: {ts['min']:.0f}°C - {ts['max']:.0f}°C (mean: {ts['mean']:.0f}°C)")
        
        if stats['time_stats']:
            times = stats['time_stats']
            coverage = times['count'] / stats['total_reactions'] * 100
            print(f"   Time data: {times['count']:,}/{stats['total_reactions']:,} ({coverage:.1f}%)")
            print(f"   Time range: {times['min']:.1f}h - {times['max']:.1f}h (mean: {times['mean']:.1f}h, median: {times['median']:.1f}h)")
        
        print()
    
    print(f"🎯 TOTAL REACTIONS ACROSS ALL FAMILIES: {total_reactions:,}")
    print()
    print("=" * 80)
    print()
    
    # Detailed analytics for each family
    for family in families:
        print("=" * 80)
        print(f" DETAILED ANALYTICS: {family}")
        print("=" * 80)
        print()
        
        # Top catalytic systems (complete catalyst + ligand combinations)
        print("[TOP 10 CATALYTIC SYSTEMS]")
        systems = get_common_catalytic_systems(family, top_n=10)
        if systems:
            for system, count, avg_yield in systems:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {system}")
        else:
            print("  No catalytic system data")
        print()
        
        # Top catalysts
        print("[TOP 10 CATALYSTS]")
        catalysts = get_common_catalysts(family, top_n=10)
        if catalysts:
            for name, count, avg_yield in catalysts:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
        else:
            print("  No catalyst data")
        print()
        
        # Top ligands
        print("[TOP 10 LIGANDS]")
        ligands = get_common_ligands(family, top_n=10)
        if ligands:
            for name, count, avg_yield in ligands:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
        else:
            print("  No ligand data")
        print()
        
        # Top bases
        print("[TOP 10 BASES]")
        bases = get_common_bases(family, top_n=10)
        if bases:
            for name, count, avg_yield in bases:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
        else:
            print("  No base data")
        print()
        
        # Top solvents
        print("[TOP 10 SOLVENTS]")
        solvents = get_common_solvents(family, top_n=10)
        if solvents:
            for name, count, avg_yield in solvents:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
        else:
            print("  No solvent data")
        print()
        
        # Top condensation agents (for amide formation)
        if family == "Amide_formation":
            print("[TOP 10 CONDENSATION AGENTS]")
            condensation_agents = get_common_reagents(family, role='CONDENSATION_AGENT', top_n=10)
            if condensation_agents:
                for name, role, count, avg_yield in condensation_agents:
                    yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                    print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
            else:
                print("  No condensation agent data")
            print()
        
        # Top condition cores
        print("[TOP 10 CONDITION CORES]")
        cores = get_condition_cores(family, top_n=10)
        if cores:
            for core, count, avg_yield in cores:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                core_display = core[:65] + "..." if len(core) > 65 else core
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {core_display}")
        else:
            print("  No condition core data")
        print()
    
    print("=" * 80)
    print(" SUMMARY COMPLETE")
    print("=" * 80)
    print()
    print("Dataset files located in: data/reaction_dataset/")
    print("- *.jsonl files contain full reaction data")
    print("- *_drfp.npz files contain DRFP fingerprints for similarity search")
    print()
    print("Supporting databases:")
    print("- data/rule_db/ - Curated condition databases by reaction type")
    print("- data/protocol_db/ - Protocol-based recommendation database")
    print("- data/reagent_db/ - Reagent taxonomy and classification")
    print()


if __name__ == "__main__":
    generate_summary()
