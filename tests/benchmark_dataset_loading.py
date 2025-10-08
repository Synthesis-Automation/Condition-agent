"""
Benchmark dataset loading AND search performance.

Tests:
1. How long it takes to load each of the 5 datasets into memory
2. How long it takes to perform precedent search operations (k=5)

This helps identify whether loading or search operations are the bottleneck.
"""

import time
import json
from pathlib import Path
import sys

# Add project root to path
project_root = Path(__file__).parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from chemtools.recommend.core import recommend_from_reaction


def count_reactions_in_file(file_path):
    """Count reactions in a JSONL file."""
    count = 0
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                count += 1
    return count


def load_dataset_simple(file_path):
    """Load dataset as simple list of dicts (minimal parsing)."""
    reactions = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                reactions.append(json.loads(line))
    return reactions


def benchmark_dataset(name, file_path, family_name, test_smiles):
    """Benchmark loading and searching a single dataset.
    
    Args:
        name: Display name
        file_path: Path to JSONL file
        family_name: Family name for precedent search
        test_smiles: Test reaction SMILES for search benchmark
    """
    print(f"\n{'=' * 70}")
    print(f"Dataset: {name}")
    print(f"{'=' * 70}")
    print(f"File: {file_path.name}")
    
    # Check if file exists
    if not file_path.exists():
        print(f"❌ File not found: {file_path}")
        return None
    
    # Get file size
    file_size_mb = file_path.stat().st_size / (1024 * 1024)
    print(f"File size: {file_size_mb:.2f} MB")
    
    # Count reactions (quick scan)
    print("\n⏱️  Step 1: Counting reactions...")
    start = time.perf_counter()
    count = count_reactions_in_file(file_path)
    count_time = time.perf_counter() - start
    print(f"   Reaction count: {count:,}")
    print(f"   Time: {count_time:.3f} seconds")
    
    # Load full dataset into memory
    print("\n⏱️  Step 2: Loading dataset into memory...")
    start = time.perf_counter()
    reactions = load_dataset_simple(file_path)
    load_time = time.perf_counter() - start
    print(f"   Loaded: {len(reactions):,} reactions")
    print(f"   Time: {load_time:.3f} seconds")
    print(f"   Rate: {len(reactions) / load_time:.0f} reactions/second")
    
    # Memory estimate (rough)
    if reactions:
        sample_size = sys.getsizeof(str(reactions[0]))
        estimated_mb = (sample_size * len(reactions)) / (1024 * 1024)
        print(f"   Estimated memory: ~{estimated_mb:.1f} MB")
    
    # Benchmark search operation (k=5)
    print("\n⏱️  Step 3: Precedent search (k=5)...")
    search_time = None
    search_success = False
    
    try:
        start = time.perf_counter()
        result = recommend_from_reaction(
            test_smiles,
            k=5,
            family_override=family_name
        )
        search_time = time.perf_counter() - start
        
        # Check if we got results
        if result:
            precedent_pack = result.get('precedent_pack', {})
            precedents = precedent_pack.get('precedents', []) if isinstance(precedent_pack, dict) else []
            precedent_count = len(precedents)
            
            print(f"   Found: {precedent_count} precedents")
            print(f"   Time: {search_time:.3f} seconds")
            if search_time > 0:
                print(f"   Rate: {count / search_time:.0f} reactions searched/second")
            search_success = True
        else:
            print(f"   ⚠️  No results returned")
            print(f"   Time: {search_time:.3f} seconds")
            search_success = False
            
    except Exception as e:
        print(f"   ❌ Search failed: {e}")
        search_time = 0
    
    return {
        'name': name,
        'file': file_path.name,
        'file_size_mb': file_size_mb,
        'reaction_count': count,
        'count_time': count_time,
        'load_time': load_time,
        'search_time': search_time if search_time else 0,
        'search_success': search_success,
        'rate_load': len(reactions) / load_time if load_time > 0 else 0,
        'rate_search': count / search_time if search_time and search_time > 0 else 0
    }


def main():
    """Run benchmarks for all 5 datasets."""
    print("=" * 70)
    print("DATASET LOADING & SEARCH BENCHMARK")
    print("=" * 70)
    print("\nMeasuring:")
    print("  1. Dataset loading time (JSON parsing)")
    print("  2. Precedent search time (k=5, no ML)")
    print()
    
    data_dir = Path(__file__).parent / "data" / "reaction_dataset"
    
    # Test SMILES for each reaction type
    # Using simple, typical reactions for each family
    datasets = [
        (
            "C-N Coupling (Cu)", 
            data_dir / "C_N_Coupling_Cu.jsonl",
            "C_N_Coupling_Cu",
            "Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"  # Aniline + chlorobenzene
        ),
        (
            "C-N Coupling (Pd)", 
            data_dir / "C_N_Coupling_Pd.jsonl",
            "C_N_Coupling_Pd",
            "Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"  # Same reaction, different catalyst
        ),
        (
            "C-N Coupling (Ni)", 
            data_dir / "C_N_Coupling_Ni.jsonl",
            "C_N_Coupling_Ni",
            "Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"  # Same reaction, Ni catalyst
        ),
        (
            "Suzuki", 
            data_dir / "Suzuki.jsonl",
            "Suzuki",
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"  # Bromobenzene + phenylboronic acid
        ),
        (
            "Amide Formation", 
            data_dir / "Amide_formation.jsonl",
            "Amide_formation",
            "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1"  # Acetic acid + benzylamine
        ),
    ]
    
    results = []
    total_start = time.perf_counter()
    
    for name, file_path, family_name, test_smiles in datasets:
        result = benchmark_dataset(name, file_path, family_name, test_smiles)
        if result:
            results.append(result)
        time.sleep(0.1)  # Brief pause between datasets
    
    total_time = time.perf_counter() - total_start
    
    # Summary table
    print("\n" + "=" * 70)
    print("LOADING PERFORMANCE SUMMARY")
    print("=" * 70)
    print()
    print(f"{'Dataset':<25} {'Count':>10} {'Size (MB)':>10} {'Load Time':>12} {'Rate':>12}")
    print("-" * 70)
    
    for r in results:
        print(f"{r['name']:<25} {r['reaction_count']:>10,} {r['file_size_mb']:>10.2f} "
              f"{r['load_time']:>11.3f}s {r['rate_load']:>11.0f}/s")
    
    print("-" * 70)
    total_reactions = sum(r['reaction_count'] for r in results)
    total_size = sum(r['file_size_mb'] for r in results)
    total_load_time = sum(r['load_time'] for r in results)
    
    if total_load_time > 0:
        print(f"{'TOTAL':<25} {total_reactions:>10,} {total_size:>10.2f} "
              f"{total_load_time:>11.3f}s {total_reactions/total_load_time:>11.0f}/s")
    else:
        print(f"{'TOTAL':<25} {total_reactions:>10,} {total_size:>10.2f} "
              f"{'0.000':>11}s {'N/A':>12}")
    
    # Search performance summary
    print("\n" + "=" * 70)
    print("SEARCH PERFORMANCE SUMMARY (k=5)")
    print("=" * 70)
    print()
    print(f"{'Dataset':<25} {'Count':>10} {'Search Time':>12} {'Rate':>12} {'Status':>10}")
    print("-" * 70)
    
    for r in results:
        status = "✅ OK" if r['search_success'] else "❌ FAIL"
        search_time_str = f"{r['search_time']:.3f}s" if r['search_time'] > 0 else "N/A"
        rate_str = f"{r['rate_search']:.0f}/s" if r['rate_search'] > 0 else "N/A"
        print(f"{r['name']:<25} {r['reaction_count']:>10,} {search_time_str:>12} {rate_str:>12} {status:>10}")
    
    print("-" * 70)
    total_search_time = sum(r['search_time'] for r in results if r['search_time'] > 0)
    successful_searches = sum(1 for r in results if r['search_success'])
    
    if total_search_time > 0:
        print(f"{'TOTAL':<25} {total_reactions:>10,} {total_search_time:>11.3f}s "
              f"{total_reactions/total_search_time:>11.0f}/s {successful_searches}/{len(results)} OK")
    
    print()
    print(f"Total benchmark time: {total_time:.2f} seconds")
    print(f"Total reactions: {total_reactions:,}")
    if total_load_time > 0:
        print(f"Average load rate: {total_reactions / total_load_time:.0f} reactions/second")
    else:
        print(f"Average load rate: N/A (no datasets loaded)")
    
    if total_search_time > 0:
        print(f"Average search rate: {total_reactions / total_search_time:.0f} reactions/second")
        print()
        print(f"⚡ Loading is {total_search_time / total_load_time:.1f}x faster than searching")
    
    # Identify slowest operations
    if results:
        slowest_load = max(results, key=lambda x: x['load_time'])
        fastest_load = min(results, key=lambda x: x['load_time'])
        
        search_results = [r for r in results if r['search_time'] > 0]
        if search_results:
            slowest_search = max(search_results, key=lambda x: x['search_time'])
            fastest_search = min(search_results, key=lambda x: x['search_time'])
            
            print()
            print("📊 Performance Analysis:")
            print(f"   Loading:")
            print(f"      ⚡ Fastest: {fastest_load['name']} ({fastest_load['load_time']:.3f}s)")
            print(f"      🐌 Slowest: {slowest_load['name']} ({slowest_load['load_time']:.3f}s)")
            print(f"      📈 Slowdown: {slowest_load['load_time'] / fastest_load['load_time']:.1f}x")
            print()
            print(f"   Searching (k=5):")
            print(f"      ⚡ Fastest: {fastest_search['name']} ({fastest_search['search_time']:.3f}s)")
            print(f"      🐌 Slowest: {slowest_search['name']} ({slowest_search['search_time']:.3f}s)")
            print(f"      📈 Slowdown: {slowest_search['search_time'] / fastest_search['search_time']:.1f}x")


if __name__ == "__main__":
    main()
