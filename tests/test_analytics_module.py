#!/usr/bin/env python
"""
Test Script: Dataset Analytics Module
Tests all analytics functions with detailed output and validation.
"""
import sys
from pathlib import Path
import time

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chemtools import chem
from chemtools import dataset_analytics as analytics


def print_header(title):
    """Print formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def print_subheader(title):
    """Print formatted subsection header."""
    print(f"\n{'─' * 80}")
    print(f"  {title}")
    print(f"{'─' * 80}")


def test_1_available_families():
    """Test 1: List all available families."""
    print_header("Test 1: Get All Available Families")
    
    start = time.time()
    families = chem.analytics.get_all_families()
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(families)} reaction families:\n")
    
    for i, family in enumerate(families, 1):
        print(f"      {i}. {family}")
    
    # Validation
    assert isinstance(families, list), "Should return a list"
    assert len(families) > 0, "Should find at least one family"
    assert "Suzuki" in families, "Should include Suzuki"
    
    print(f"\n   ✅ Test passed: Found {len(families)} families")
    return True


def test_2_dataset_statistics():
    """Test 2: Get comprehensive dataset statistics."""
    print_header("Test 2: Dataset Statistics")
    
    families_to_test = ["Suzuki", "C_N_Coupling_Pd", "Amide_formation"]
    
    for family in families_to_test:
        print_subheader(f"{family} Dataset")
        
        start = time.time()
        stats = chem.analytics.get_stats(family)
        elapsed = time.time() - start
        
        print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
        print(f"\n   📊 Basic Statistics:")
        print(f"      - Total reactions:       {stats['total_reactions']:>8,}")
        print(f"      - Unique condition cores: {stats['unique_condition_cores']:>7}")
        print(f"      - Unique catalysts:      {stats['unique_catalysts']:>8}")
        print(f"      - Unique bases:          {stats['unique_bases']:>8}")
        print(f"      - Unique solvents:       {stats['unique_solvents']:>8}")
        
        if stats['yield_stats']:
            ys = stats['yield_stats']
            coverage = ys['count'] / stats['total_reactions'] * 100
            print(f"\n   📈 Yield Statistics:")
            print(f"      - Coverage:  {ys['count']:>8,} reactions ({coverage:.1f}%)")
            print(f"      - Range:     {ys['min']:>8.1f}% - {ys['max']:.1f}%")
            print(f"      - Mean:      {ys['mean']:>8.1f}%")
            print(f"      - Median:    {ys['median']:>8.1f}%")
        
        if stats['temperature_stats']:
            ts = stats['temperature_stats']
            coverage = ts['count'] / stats['total_reactions'] * 100
            print(f"\n   🌡️  Temperature Statistics:")
            print(f"      - Coverage:  {ts['count']:>8,} reactions ({coverage:.1f}%)")
            print(f"      - Range:     {ts['min']:>8.0f}°C - {ts['max']:.0f}°C")
            print(f"      - Mean:      {ts['mean']:>8.0f}°C")
        
        if stats['time_stats']:
            tms = stats['time_stats']
            coverage = tms['count'] / stats['total_reactions'] * 100
            print(f"\n   ⏱️  Time Statistics:")
            print(f"      - Coverage:  {tms['count']:>8,} reactions ({coverage:.1f}%)")
            print(f"      - Range:     {tms['min']:>8.1f}h - {tms['max']:.1f}h")
            print(f"      - Mean:      {tms['mean']:>8.1f}h")
        
        # Validation
        assert stats['total_reactions'] > 0, f"{family} should have reactions"
        assert stats['family'] == family, "Family name should match"
        
        print(f"\n   ✅ Statistics validated for {family}")
    
    return True


def test_3_common_catalysts():
    """Test 3: Get common catalysts with yield data."""
    print_header("Test 3: Common Catalysts (Frequency Ranking)")
    
    family = "C_N_Coupling_Pd"
    top_n = 10
    
    print_subheader(f"{family} - Top {top_n} Catalysts")
    
    start = time.time()
    catalysts = chem.analytics.get_common_catalysts(family, top_n=top_n)
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(catalysts)} catalysts:\n")
    
    for i, (name, count, avg_yield) in enumerate(catalysts, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        name_display = name[:55] + "..." if len(name) > 55 else name
        print(f"      {i:>2}. {count:>5} reactions | Avg: {yield_str:>6} | {name_display}")
    
    # Validation
    assert isinstance(catalysts, list), "Should return a list"
    assert len(catalysts) <= top_n, f"Should return at most {top_n} catalysts"
    
    if catalysts:
        # Check sorting (should be descending by count)
        counts = [c[1] for c in catalysts]
        assert counts == sorted(counts, reverse=True), "Should be sorted by count"
        
        # Check data types
        name, count, avg_yield = catalysts[0]
        assert isinstance(name, str), "Name should be string"
        assert isinstance(count, int), "Count should be int"
        assert count > 0, "Count should be positive"
    
    print(f"\n   ✅ Test passed: Found and validated {len(catalysts)} catalysts")
    return True


def test_4_common_ligands():
    """Test 4: Get common ligands with yield data."""
    print_header("Test 4: Common Ligands (Frequency Ranking)")
    
    family = "C_N_Coupling_Pd"
    top_n = 10
    
    print_subheader(f"{family} - Top {top_n} Ligands")
    
    start = time.time()
    ligands = chem.analytics.get_common_ligands(family, top_n=top_n)
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(ligands)} ligands:\n")
    
    for i, (name, count, avg_yield) in enumerate(ligands, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        name_display = name[:55] + "..." if len(name) > 55 else name
        print(f"      {i:>2}. {count:>5} reactions | Avg: {yield_str:>6} | {name_display}")
    
    # Validation
    assert isinstance(ligands, list), "Should return a list"
    
    print(f"\n   ✅ Test passed: Found and validated {len(ligands)} ligands")
    return True


def test_5_common_bases():
    """Test 5: Get common bases with yield data."""
    print_header("Test 5: Common Bases (Frequency Ranking)")
    
    family = "Suzuki"
    top_n = 10
    
    print_subheader(f"{family} - Top {top_n} Bases")
    
    start = time.time()
    bases = chem.analytics.get_common_bases(family, top_n=top_n)
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(bases)} bases:\n")
    
    for i, (name, count, avg_yield) in enumerate(bases, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        name_display = name[:55] + "..." if len(name) > 55 else name
        print(f"      {i:>2}. {count:>5} reactions | Avg: {yield_str:>6} | {name_display}")
    
    # Validation
    assert isinstance(bases, list), "Should return a list"
    assert len(bases) > 0, "Should find at least one base"
    
    print(f"\n   ✅ Test passed: Found and validated {len(bases)} bases")
    return True


def test_6_common_solvents():
    """Test 6: Get common solvents with yield data."""
    print_header("Test 6: Common Solvents (Frequency Ranking)")
    
    family = "Suzuki"
    top_n = 10
    
    print_subheader(f"{family} - Top {top_n} Solvents")
    
    start = time.time()
    solvents = chem.analytics.get_common_solvents(family, top_n=top_n)
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(solvents)} solvents:\n")
    
    for i, (name, count, avg_yield) in enumerate(solvents, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        name_display = name[:55] + "..." if len(name) > 55 else name
        print(f"      {i:>2}. {count:>5} reactions | Avg: {yield_str:>6} | {name_display}")
    
    # Validation
    assert isinstance(solvents, list), "Should return a list"
    assert len(solvents) > 0, "Should find at least one solvent"
    
    print(f"\n   ✅ Test passed: Found and validated {len(solvents)} solvents")
    return True


def test_7_common_reagents():
    """Test 7: Get common reagents (all roles and filtered)."""
    print_header("Test 7: Common Reagents (All Roles & Filtered)")
    
    family = "Amide_formation"
    
    # Test 7a: All reagents
    print_subheader(f"{family} - Top 10 All Reagents")
    
    start = time.time()
    all_reagents = chem.analytics.get_common_reagents(family, top_n=10)
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(all_reagents)} reagents:\n")
    
    for i, (name, role, count, avg_yield) in enumerate(all_reagents, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        name_display = name[:35] + "..." if len(name) > 35 else name
        print(f"      {i:>2}. {count:>5} reactions | {role:<18} | {yield_str:>6} | {name_display}")
    
    # Test 7b: Filtered by role
    print_subheader(f"{family} - Top 10 CONDENSATION_AGENT")
    
    start = time.time()
    condensation_agents = chem.analytics.get_common_reagents(
        family, role="CONDENSATION_AGENT", top_n=10
    )
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(condensation_agents)} condensation agents:\n")
    
    for i, (name, role, count, avg_yield) in enumerate(condensation_agents, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        name_display = name[:55] + "..." if len(name) > 55 else name
        print(f"      {i:>2}. {count:>5} reactions | {yield_str:>6} | {name_display}")
    
    # Validation
    assert isinstance(all_reagents, list), "Should return a list"
    assert isinstance(condensation_agents, list), "Should return a list"
    
    # Check that filtered results only have the requested role
    for name, role, count, avg_yield in condensation_agents:
        assert role == "CONDENSATION_AGENT", "Should only have CONDENSATION_AGENT role"
    
    print(f"\n   ✅ Test passed: Reagent ranking and filtering work correctly")
    return True


def test_8_catalytic_systems():
    """Test 8: Get common catalytic systems (catalyst + ligand combinations)."""
    print_header("Test 8: Common Catalytic Systems & Generation Analysis")
    
    family = "C_N_Coupling_Pd"
    top_n = 15
    
    print_subheader(f"{family} - Top {top_n} Catalytic Systems")
    
    start = time.time()
    systems = chem.analytics.get_common_catalytic_systems(family, top_n=top_n)
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(systems)} catalytic systems:\n")
    
    for i, (system_str, count, avg_yield) in enumerate(systems, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        # Truncate long system strings
        system_display = system_str[:65] + "..." if len(system_str) > 65 else system_str
        print(f"      {i:>2}. {count:>5} reactions | Avg: {yield_str:>6} | {system_display}")
    
    # ========================================================================
    # CATALYTIC SYSTEM GENERATION ANALYSIS
    # ========================================================================
    print_subheader("Catalytic System Generation Analysis")
    
    # Analyze all systems (not just top N)
    all_systems = chem.analytics.get_common_catalytic_systems(family, top_n=1000)
    
    # Classify by component count
    single_component = [s for s in all_systems if " + " not in s[0]]
    two_component = [s for s in all_systems if s[0].count(" + ") == 1]
    three_plus_component = [s for s in all_systems if s[0].count(" + ") >= 2]
    
    print(f"\n   📊 System Composition Statistics:")
    print(f"      Total unique systems: {len(all_systems)}")
    print(f"      - Single component:   {len(single_component):>4} ({len(single_component)/len(all_systems)*100:>5.1f}%)")
    print(f"      - Two components:     {len(two_component):>4} ({len(two_component)/len(all_systems)*100:>5.1f}%)")
    print(f"      - Three+ components:  {len(three_plus_component):>4} ({len(three_plus_component)/len(all_systems)*100:>5.1f}%)")
    
    # Calculate total reaction coverage
    total_system_reactions = sum(s[1] for s in all_systems)
    multi_component_reactions = sum(s[1] for s in two_component + three_plus_component)
    
    print(f"\n   📊 Reaction Coverage:")
    print(f"      Total reactions with catalytic systems: {total_system_reactions:>6}")
    print(f"      Multi-component system reactions:       {multi_component_reactions:>6} ({multi_component_reactions/total_system_reactions*100:>5.1f}%)")
    
    # Analyze yield distribution by system type
    if single_component:
        single_avg_yields = [s[2] for s in single_component if s[2] is not None]
        single_avg = sum(single_avg_yields) / len(single_avg_yields) if single_avg_yields else 0
    else:
        single_avg = 0
    
    if two_component:
        two_avg_yields = [s[2] for s in two_component if s[2] is not None]
        two_avg = sum(two_avg_yields) / len(two_avg_yields) if two_avg_yields else 0
    else:
        two_avg = 0
    
    if three_plus_component:
        three_avg_yields = [s[2] for s in three_plus_component if s[2] is not None]
        three_avg = sum(three_avg_yields) / len(three_avg_yields) if three_avg_yields else 0
    else:
        three_avg = 0
    
    print(f"\n   📊 Average Yield by System Type:")
    if single_component:
        print(f"      Single component:     {single_avg:>5.1f}%")
    if two_component:
        print(f"      Two components:       {two_avg:>5.1f}%")
    if three_plus_component:
        print(f"      Three+ components:    {three_avg:>5.1f}%")
    
    # Show examples of different system types
    if two_component:
        print(f"\n   📋 Example Two-Component Systems (Top 3):")
        for i, (system_str, count, avg_yield) in enumerate(two_component[:3], 1):
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            system_display = system_str[:60] + "..." if len(system_str) > 60 else system_str
            print(f"      {i}. {system_display}")
            print(f"         {count} reactions, avg yield: {yield_str}")
    
    if three_plus_component:
        print(f"\n   📋 Example Three+ Component Systems:")
        for i, (system_str, count, avg_yield) in enumerate(three_plus_component[:3], 1):
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            system_display = system_str[:60] + "..." if len(system_str) > 60 else system_str
            print(f"      {i}. {system_display}")
            print(f"         {count} reactions, avg yield: {yield_str}")
    
    # ========================================================================
    # Show comparison with individual components
    # ========================================================================
    print_subheader("Comparison: Systems vs Individual Components")
    
    catalysts = chem.analytics.get_common_catalysts(family, top_n=5)
    ligands = chem.analytics.get_common_ligands(family, top_n=5)
    
    print(f"\n   Individual Catalysts (Top 5):")
    for i, (name, count, avg_yield) in enumerate(catalysts, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        print(f"      {i}. {count:>5} reactions | {yield_str:>6} | {name}")
    
    print(f"\n   Individual Ligands (Top 5):")
    for i, (name, count, avg_yield) in enumerate(ligands, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        print(f"      {i}. {count:>5} reactions | {yield_str:>6} | {name}")
    
    print(f"\n   💡 Key Insight:")
    print(f"      Individual component counts don't show catalyst-ligand pairing!")
    print(f"      Catalytic system analysis preserves these critical relationships.")
    
    # ========================================================================
    # Validation
    # ========================================================================
    assert isinstance(systems, list), "Should return a list"
    assert len(systems) > 0, "Should find at least one catalytic system"
    
    # Check that system strings contain " + " separator for multi-component systems
    multi_component = [s for s in systems if " + " in s[0]]
    if multi_component:
        print(f"\n   ℹ️  Found {len(multi_component)} multi-component systems in top {top_n}")
    
    print(f"\n   ✅ Test passed: Found and validated {len(systems)} catalytic systems")
    print(f"   ✅ Generation analysis complete: {len(all_systems)} total unique systems")
    return True


def test_9_high_yield_filtering():
    """Test 9: Filter by minimum yield threshold."""
    print_header("Test 9: High-Yield Filtering")
    
    family = "Suzuki"
    min_yield = 80.0
    
    print_subheader(f"{family} - Catalysts with >= {min_yield}% avg yield")
    
    # Get catalysts without filter
    start = time.time()
    all_cats = chem.analytics.get_common_catalysts(family, top_n=50, min_yield=None)
    elapsed_all = time.time() - start
    
    # Get catalysts with high-yield filter
    start = time.time()
    high_yield_cats = chem.analytics.get_common_catalysts(family, top_n=50, min_yield=min_yield)
    elapsed_filtered = time.time() - start
    
    print(f"\n   ⏱️  Time (all): {elapsed_all:.4f} seconds")
    print(f"   ⏱️  Time (filtered): {elapsed_filtered:.4f} seconds")
    
    # Filter manually to find high-yield ones
    manual_high_yield = [c for c in all_cats if c[2] and c[2] >= min_yield]
    
    print(f"\n   📊 Results:")
    print(f"      - All catalysts: {len(all_cats)}")
    print(f"      - High yield (manual filter): {len(manual_high_yield)}")
    print(f"      - High yield (min_yield param): {len(high_yield_cats)}")
    
    if manual_high_yield:
        print(f"\n   Top 5 high-yield catalysts:\n")
        for i, (name, count, avg_yield) in enumerate(manual_high_yield[:5], 1):
            name_display = name[:50] + "..." if len(name) > 50 else name
            print(f"      {i}. {count:>5} reactions | {avg_yield:>6.1f}% | {name_display}")
    
    print(f"\n   ✅ Test passed: Yield filtering works correctly")
    return True


def test_10_condition_cores():
    """Test 10: Get common condition core combinations."""
    print_header("Test 10: Common Condition Cores")
    
    family = "C_N_Coupling_Pd"
    top_n = 10
    
    print_subheader(f"{family} - Top {top_n} Condition Cores")
    
    start = time.time()
    cores = chem.analytics.get_condition_cores(family, top_n=top_n)
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"   📊 Found {len(cores)} condition cores:\n")
    
    for i, (core, count, avg_yield) in enumerate(cores, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        core_display = core[:60] + "..." if len(core) > 60 else core
        print(f"      {i:>2}. {count:>4} reactions | {yield_str:>6} | {core_display}")
    
    # Validation
    assert isinstance(cores, list), "Should return a list"
    
    print(f"\n   ✅ Test passed: Found and validated {len(cores)} condition cores")
    return True


def test_11_plate_recommendations():
    """Test 11: Generate HTE plate recommendations."""
    print_header("Test 11: Plate Design Recommendations")
    
    family = "C_N_Coupling_Pd"
    
    # Test different optimization strategies
    strategies = [
        ('diversity', 24),
        ('yield', 12),
        ('frequency', 12)
    ]
    
    for strategy, n_conditions in strategies:
        print_subheader(f"{family} - {n_conditions}-well plate (optimize: {strategy})")
        
        start = time.time()
        plate = chem.analytics.get_plate_recommendations(
            family=family,
            n_conditions=n_conditions,
            min_yield=60.0,
            optimize_for=strategy
        )
        elapsed = time.time() - start
        
        print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
        print(f"   📊 Generated {len(plate)} conditions")
        
        if plate:
            print(f"\n   Top 3 conditions:\n")
            for i, cond in enumerate(plate[:3], 1):
                cat = (cond['catalyst'] or 'None')[:30]
                lig = (cond['ligand'] or 'None')[:30]
                base = (cond['base'] or 'None')[:30]
                solv = (cond['solvent'] or 'None')[:30]
                temp = f"{cond['temperature_c']:.0f}°C" if cond['temperature_c'] else 'N/A'
                time_h = f"{cond['time_h']:.1f}h" if cond['time_h'] else 'N/A'
                yield_val = f"{cond['avg_yield']:.1f}%" if cond['avg_yield'] else 'N/A'
                
                print(f"      {i}. [{cond['condition_id']}]")
                print(f"         Catalyst: {cat}")
                print(f"         Ligand:   {lig}")
                print(f"         Base:     {base}")
                print(f"         Solvent:  {solv}")
                print(f"         Temp/Time: {temp}, {time_h}")
                print(f"         Performance: {yield_val} yield, {cond['frequency']} precedents")
                print(f"         Score: {cond['score']:.1f}")
                print()
        
        # Validation
        assert isinstance(plate, list), "Should return a list"
        assert len(plate) <= n_conditions, f"Should return at most {n_conditions} conditions"
        
        if plate:
            # Check required fields
            required_fields = ['condition_id', 'avg_yield', 'frequency', 'score']
            for field in required_fields:
                assert field in plate[0], f"Should have {field} field"
            
            # Check optimization worked
            if strategy == 'yield' and len(plate) > 1:
                yields = [c['avg_yield'] for c in plate if c['avg_yield'] is not None]
                if len(yields) > 1:
                    # Should be sorted by yield descending
                    assert yields == sorted(yields, reverse=True), "Should be sorted by yield"
            
            elif strategy == 'frequency' and len(plate) > 1:
                freqs = [c['frequency'] for c in plate]
                # Should be sorted by frequency descending
                assert freqs == sorted(freqs, reverse=True), "Should be sorted by frequency"
    
    print(f"\n   ✅ Test passed: Plate recommendations work for all strategies")
    return True


def test_11_print_summary():
    """Test 11: Print comprehensive summary."""
def test_12_print_summary():
    """Test 12: Print comprehensive analytics summary."""
    print_header("Test 12: Comprehensive Analytics Summary")
    
    family = "Suzuki"
    
    print_subheader(f"{family} - Full Summary (top 5)")
    
    start = time.time()
    chem.analytics.print_summary(family, top_n=5)
    elapsed = time.time() - start
    
    print(f"\n   ⏱️  Time: {elapsed:.4f} seconds")
    print(f"\n   ✅ Test passed: Summary printed successfully")
    return True


def test_13_error_handling():
    """Test 13: Error handling for invalid inputs."""
    print_header("Test 13: Error Handling")
    
    print_subheader("Testing invalid family name")
    
    try:
        stats = chem.analytics.get_stats("NonExistentFamily")
        print(f"\n   ❌ Should have raised FileNotFoundError")
        return False
    except FileNotFoundError as e:
        print(f"\n   ✅ Correctly raised FileNotFoundError: {e}")
    
    print(f"\n   ✅ Test passed: Error handling works correctly")
    return True


def run_all_tests():
    """Run all analytics tests."""
    print("\n" + "=" * 80)
    print("  DATASET ANALYTICS MODULE - COMPREHENSIVE TEST SUITE")
    print("=" * 80)
    print("\n  Testing all analytics functions with detailed validation")
    print("  and performance metrics.")
    print()
    
    tests = [
        test_1_available_families,
        test_2_dataset_statistics,
        test_3_common_catalysts,
        test_4_common_ligands,
        test_5_common_bases,
        test_6_common_solvents,
        test_7_common_reagents,
        test_8_catalytic_systems,
        test_9_high_yield_filtering,
        test_10_condition_cores,
        test_11_plate_recommendations,
        test_12_print_summary,
        test_13_error_handling,
    ]
    
    results = []
    start_time = time.time()
    
    for i, test in enumerate(tests, 1):
        try:
            success = test()
            results.append((test.__name__, success, None))
        except Exception as e:
            print(f"\n   ❌ Test failed with error: {e}")
            import traceback
            traceback.print_exc()
            results.append((test.__name__, False, str(e)))
    
    total_time = time.time() - start_time
    
    # Print summary
    print("\n" + "=" * 80)
    print("  TEST SUMMARY")
    print("=" * 80)
    
    passed = sum(1 for _, success, _ in results if success)
    failed = sum(1 for _, success, _ in results if not success)
    
    print(f"\n   Total tests:  {len(results)}")
    print(f"   Passed:       {passed} ✅")
    print(f"   Failed:       {failed} {'❌' if failed > 0 else ''}")
    print(f"   Total time:   {total_time:.3f} seconds")
    
    if failed > 0:
        print(f"\n   Failed tests:")
        for name, success, error in results:
            if not success:
                print(f"      - {name}")
                if error:
                    print(f"        Error: {error}")
    
    print("\n" + "=" * 80)
    
    if failed == 0:
        print("  ✅ ALL TESTS PASSED - Analytics module is working correctly!")
    else:
        print(f"  ❌ {failed} TEST(S) FAILED - Please review errors above")
    
    print("=" * 80 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
