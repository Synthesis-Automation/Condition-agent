"""
Dataset analytics for reaction condition analysis.

Provides statistical analysis of reaction datasets including:
- Common catalysts, ligands, bases, solvents, and reagents (with ranking)
- Temperature, time, and yield distributions
- Plate design recommendations based on frequency and success
- Condition combination analysis

This is useful for:
- Data-driven condition recommendations
- High-throughput experimentation (HTE) plate design
- Understanding chemical space coverage
- Identifying successful condition patterns
"""

import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from functools import lru_cache
from statistics import mean, median


def get_data_dir() -> Path:
    """Get the data directory path."""
    module_dir = Path(__file__).parent.parent
    return module_dir / "data" / "reaction_dataset"


def _load_dataset(family: str) -> List[Dict[str, Any]]:
    """Load a reaction dataset JSONL file.
    
    Args:
        family: Reaction family name (e.g., "Suzuki", "C_N_Coupling_Pd")
        
    Returns:
        List of reaction dictionaries
    """
    dataset_path = get_data_dir() / f"{family}.jsonl"
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset not found: {dataset_path}")
    
    reactions = []
    with open(dataset_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                reactions.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    
    return reactions


def get_dataset_stats(family: str) -> Dict[str, Any]:
    """Get basic statistics for a reaction dataset.
    
    Args:
        family: Reaction family name
        
    Returns:
        Dictionary with dataset statistics:
        - total_reactions: Total number of reactions
        - unique_condition_cores: Number of unique condition cores
        - unique_solvents: Number of unique solvents
        - unique_bases: Number of unique bases
        - unique_catalysts: Number of unique catalysts
        - yield_stats: {count, min, max, mean, median}
        - temperature_stats: {count, min, max, mean, median}
        - time_stats: {count, min, max, mean, median}
    """
    reactions = _load_dataset(family)
    
    cores = set()
    solvents = set()
    bases = set()
    catalysts = set()
    yields = []
    temps = []
    times = []
    
    for rxn in reactions:
        # Condition core
        core = rxn.get('condition_core', '')
        if core:
            cores.add(core)
        
        # Catalytic system
        for cat in rxn.get('catalytic_system', []):
            name = cat.get('name')
            if name:
                catalysts.add(name)
        
        # Reagents
        for reagent in rxn.get('reagents', []):
            role = reagent.get('role', '')
            name = reagent.get('name', '')
            if role == 'BASE' and name:
                bases.add(name)
        
        # Solvents
        for solv in rxn.get('solvents', []):
            name = solv.get('name', '')
            if name:
                solvents.add(name)
        
        # Conditions
        conditions = rxn.get('conditions', {})
        
        y = conditions.get('yield_pct')
        if y is not None:
            yields.append(y)
        
        t = conditions.get('temperature_c')
        if t is not None:
            temps.append(t)
        
        time = conditions.get('time_h')
        if time is not None:
            times.append(time)
    
    def _stats(values: List[float]) -> Optional[Dict[str, float]]:
        if not values:
            return None
        return {
            'count': len(values),
            'min': min(values),
            'max': max(values),
            'mean': mean(values),
            'median': median(values)
        }
    
    return {
        'family': family,
        'total_reactions': len(reactions),
        'unique_condition_cores': len(cores),
        'unique_solvents': len(solvents),
        'unique_bases': len(bases),
        'unique_catalysts': len(catalysts),
        'yield_stats': _stats(yields),
        'temperature_stats': _stats(temps),
        'time_stats': _stats(times),
    }


def get_common_catalysts(family: str, top_n: int = 10, min_yield: Optional[float] = None) -> List[Tuple[str, int, float]]:
    """Get most common catalysts with frequency and average yield.
    
    Args:
        family: Reaction family name
        top_n: Number of top catalysts to return
        min_yield: Optional minimum yield filter (0-100)
        
    Returns:
        List of tuples: (catalyst_name, count, avg_yield)
        Sorted by count (descending)
    """
    reactions = _load_dataset(family)
    
    catalyst_data = defaultdict(lambda: {'count': 0, 'yields': []})
    
    for rxn in reactions:
        # Filter by yield if requested
        if min_yield is not None:
            y = rxn.get('conditions', {}).get('yield_pct')
            if y is None or y < min_yield:
                continue
        
        # Extract catalysts
        for cat in rxn.get('catalytic_system', []):
            name = cat.get('name')
            if name:
                catalyst_data[name]['count'] += 1
                
                # Track yields
                y = rxn.get('conditions', {}).get('yield_pct')
                if y is not None:
                    catalyst_data[name]['yields'].append(y)
    
    # Compute averages and format
    results = []
    for name, data in catalyst_data.items():
        avg_yield = mean(data['yields']) if data['yields'] else None
        results.append((name, data['count'], avg_yield))
    
    # Sort by count descending
    results.sort(key=lambda x: x[1], reverse=True)
    
    return results[:top_n]


def get_common_solvents(family: str, top_n: int = 10, min_yield: Optional[float] = None) -> List[Tuple[str, int, float]]:
    """Get most common solvents with frequency and average yield.
    
    Args:
        family: Reaction family name
        top_n: Number of top solvents to return
        min_yield: Optional minimum yield filter (0-100)
        
    Returns:
        List of tuples: (solvent_name, count, avg_yield)
        Sorted by count (descending)
    """
    reactions = _load_dataset(family)
    
    solvent_data = defaultdict(lambda: {'count': 0, 'yields': []})
    
    for rxn in reactions:
        # Filter by yield if requested
        if min_yield is not None:
            y = rxn.get('conditions', {}).get('yield_pct')
            if y is None or y < min_yield:
                continue
        
        # Extract solvents
        for solv in rxn.get('solvents', []):
            name = solv.get('name')
            if name:
                solvent_data[name]['count'] += 1
                
                # Track yields
                y = rxn.get('conditions', {}).get('yield_pct')
                if y is not None:
                    solvent_data[name]['yields'].append(y)
    
    # Compute averages and format
    results = []
    for name, data in solvent_data.items():
        avg_yield = mean(data['yields']) if data['yields'] else None
        results.append((name, data['count'], avg_yield))
    
    # Sort by count descending
    results.sort(key=lambda x: x[1], reverse=True)
    
    return results[:top_n]


def get_common_catalytic_systems(family: str, top_n: int = 10, min_yield: Optional[float] = None) -> List[Tuple[str, int, float]]:
    """Get most common catalytic systems (complete catalyst + ligand combinations) with frequency and average yield.
    
    Analyzes catalytic_system as complete units (e.g., "Pd(OAc)2 + RuPhos") rather than individual components.
    This preserves the important relationship between catalysts and ligands used together.
    
    Args:
        family: Reaction family name
        top_n: Number of top catalytic systems to return
        min_yield: Optional minimum yield filter (0-100)
        
    Returns:
        List of tuples: (system_string, count, avg_yield)
        System string format: "Component1 + Component2 + ..."
        Sorted by count (descending)
    """
    reactions = _load_dataset(family)
    
    system_data = defaultdict(lambda: {'count': 0, 'yields': []})
    
    for rxn in reactions:
        # Filter by yield if requested
        if min_yield is not None:
            y = rxn.get('conditions', {}).get('yield_pct')
            if y is None or y < min_yield:
                continue
        
        # Extract catalytic system
        cat_sys = rxn.get('catalytic_system', [])
        if cat_sys:
            # Join component names to create system string
            component_names = [comp.get('name', '') for comp in cat_sys if comp.get('name')]
            if component_names:
                system_str = ' + '.join(component_names)
                system_data[system_str]['count'] += 1
                
                # Track yields
                y = rxn.get('conditions', {}).get('yield_pct')
                if y is not None:
                    system_data[system_str]['yields'].append(y)
    
    # Compute averages and format
    results = []
    for system_str, data in system_data.items():
        avg_yield = mean(data['yields']) if data['yields'] else None
        results.append((system_str, data['count'], avg_yield))
    
    # Sort by count descending
    results.sort(key=lambda x: x[1], reverse=True)
    
    return results[:top_n]


def get_common_bases(family: str, top_n: int = 10, min_yield: Optional[float] = None) -> List[Tuple[str, int, float]]:
    """Get most common bases with frequency and average yield.
    
    Args:
        family: Reaction family name
        top_n: Number of top bases to return
        min_yield: Optional minimum yield filter (0-100)
        
    Returns:
        List of tuples: (base_name, count, avg_yield)
        Sorted by count (descending)
    """
    reactions = _load_dataset(family)
    
    base_data = defaultdict(lambda: {'count': 0, 'yields': []})
    
    for rxn in reactions:
        # Filter by yield if requested
        if min_yield is not None:
            y = rxn.get('conditions', {}).get('yield_pct')
            if y is None or y < min_yield:
                continue
        
        # Extract bases
        for reagent in rxn.get('reagents', []):
            if reagent.get('role') == 'BASE':
                name = reagent.get('name')
                if name:
                    base_data[name]['count'] += 1
                    
                    # Track yields
                    y = rxn.get('conditions', {}).get('yield_pct')
                    if y is not None:
                        base_data[name]['yields'].append(y)
    
    # Compute averages and format
    results = []
    for name, data in base_data.items():
        avg_yield = mean(data['yields']) if data['yields'] else None
        results.append((name, data['count'], avg_yield))
    
    # Sort by count descending
    results.sort(key=lambda x: x[1], reverse=True)
    
    return results[:top_n]


def get_common_reagents(family: str, role: Optional[str] = None, top_n: int = 10, min_yield: Optional[float] = None) -> List[Tuple[str, str, int, float]]:
    """Get most common reagents with frequency and average yield.
    
    Args:
        family: Reaction family name
        role: Optional role filter (e.g., 'BASE', 'ACID', 'OXIDANT', 'COUPLING_REAGENT')
        top_n: Number of top reagents to return
        min_yield: Optional minimum yield filter (0-100)
        
    Returns:
        List of tuples: (reagent_name, role, count, avg_yield)
        Sorted by count (descending)
    """
    reactions = _load_dataset(family)
    
    reagent_data = defaultdict(lambda: {'role': '', 'count': 0, 'yields': []})
    
    for rxn in reactions:
        # Filter by yield if requested
        if min_yield is not None:
            y = rxn.get('conditions', {}).get('yield_pct')
            if y is None or y < min_yield:
                continue
        
        # Extract reagents
        for reagent in rxn.get('reagents', []):
            reagent_role = reagent.get('role', '')
            name = reagent.get('name', '')
            
            # Filter by role if requested
            if role and reagent_role != role:
                continue
            
            if name:
                key = (name, reagent_role)
                reagent_data[key]['role'] = reagent_role
                reagent_data[key]['count'] += 1
                
                # Track yields
                y = rxn.get('conditions', {}).get('yield_pct')
                if y is not None:
                    reagent_data[key]['yields'].append(y)
    
    # Compute averages and format
    results = []
    for (name, reagent_role), data in reagent_data.items():
        avg_yield = mean(data['yields']) if data['yields'] else None
        results.append((name, reagent_role, data['count'], avg_yield))
    
    # Sort by count descending
    results.sort(key=lambda x: x[2], reverse=True)
    
    return results[:top_n]


def get_common_ligands(family: str, top_n: int = 10, min_yield: Optional[float] = None) -> List[Tuple[str, int, float]]:
    """Get most common ligands with frequency and average yield.
    
    Args:
        family: Reaction family name
        top_n: Number of top ligands to return
        min_yield: Optional minimum yield filter (0-100)
        
    Returns:
        List of tuples: (ligand_name, count, avg_yield)
        Sorted by count (descending)
    """
    reactions = _load_dataset(family)
    
    ligand_data = defaultdict(lambda: {'count': 0, 'yields': []})
    
    for rxn in reactions:
        # Filter by yield if requested
        if min_yield is not None:
            y = rxn.get('conditions', {}).get('yield_pct')
            if y is None or y < min_yield:
                continue
        
        # Extract ligands from catalytic system
        for cat in rxn.get('catalytic_system', []):
            role = cat.get('role', '')
            name = cat.get('name', '')
            
            if role == 'LIGAND' and name:
                ligand_data[name]['count'] += 1
                
                # Track yields
                y = rxn.get('conditions', {}).get('yield_pct')
                if y is not None:
                    ligand_data[name]['yields'].append(y)
    
    # Compute averages and format
    results = []
    for name, data in ligand_data.items():
        avg_yield = mean(data['yields']) if data['yields'] else None
        results.append((name, data['count'], avg_yield))
    
    # Sort by count descending
    results.sort(key=lambda x: x[1], reverse=True)
    
    return results[:top_n]


def get_condition_cores(family: str, top_n: int = 10, min_yield: Optional[float] = None) -> List[Tuple[str, int, float]]:
    """Get most common condition cores with frequency and average yield.
    
    Args:
        family: Reaction family name
        top_n: Number of top cores to return
        min_yield: Optional minimum yield filter (0-100)
        
    Returns:
        List of tuples: (core_string, count, avg_yield)
        Sorted by count (descending)
    """
    reactions = _load_dataset(family)
    
    core_data = defaultdict(lambda: {'count': 0, 'yields': []})
    
    for rxn in reactions:
        # Filter by yield if requested
        if min_yield is not None:
            y = rxn.get('conditions', {}).get('yield_pct')
            if y is None or y < min_yield:
                continue
        
        # Extract core
        core = rxn.get('condition_core', '')
        if core:
            core_data[core]['count'] += 1
            
            # Track yields
            y = rxn.get('conditions', {}).get('yield_pct')
            if y is not None:
                core_data[core]['yields'].append(y)
    
    # Compute averages and format
    results = []
    for core, data in core_data.items():
        avg_yield = mean(data['yields']) if data['yields'] else None
        results.append((core, data['count'], avg_yield))
    
    # Sort by count descending
    results.sort(key=lambda x: x[1], reverse=True)
    
    return results[:top_n]


def get_plate_recommendations(
    family: str,
    n_conditions: int = 96,
    min_yield: float = 60.0,
    optimize_for: str = 'diversity'
) -> List[Dict[str, Any]]:
    """Generate high-throughput plate recommendations.
    
    Recommends N conditions for a reaction plate based on:
    - Most frequent successful conditions
    - Yield performance
    - Chemical diversity (optional)
    
    Args:
        family: Reaction family name
        n_conditions: Number of conditions to recommend (e.g., 96 for 96-well plate)
        min_yield: Minimum yield threshold (0-100)
        optimize_for: Strategy - 'diversity', 'frequency', or 'yield'
        
    Returns:
        List of recommended conditions, each with:
        - condition_id: Unique identifier
        - catalyst: Catalyst name
        - ligand: Ligand name (if applicable)
        - base: Base name
        - solvent: Solvent name (or mix)
        - temperature_c: Temperature
        - time_h: Reaction time
        - avg_yield: Average historical yield
        - frequency: Number of precedents
        - score: Recommendation score (0-100)
    """
    reactions = _load_dataset(family)
    
    # Filter by minimum yield
    filtered = [r for r in reactions 
                if r.get('conditions', {}).get('yield_pct') is not None 
                and r.get('conditions', {}).get('yield_pct') >= min_yield]
    
    if len(filtered) < n_conditions:
        print(f"Warning: Only {len(filtered)} reactions meet min_yield={min_yield}. Lowering threshold.")
        filtered = sorted(reactions, 
                         key=lambda r: r.get('conditions', {}).get('yield_pct', 0),
                         reverse=True)[:n_conditions * 2]
    
    # Build condition fingerprints
    condition_groups = defaultdict(lambda: {
        'reactions': [],
        'yields': [],
        'temps': [],
        'times': []
    })
    
    for rxn in filtered:
        # Extract components
        catalysts = [c.get('name') for c in rxn.get('catalytic_system', []) 
                    if c.get('role') == 'CATALYST']
        ligands = [c.get('name') for c in rxn.get('catalytic_system', []) 
                  if c.get('role') == 'LIGAND']
        bases = [r.get('name') for r in rxn.get('reagents', []) 
                if r.get('role') == 'BASE']
        solvents = [s.get('name') for s in rxn.get('solvents', [])]
        
        catalyst = catalysts[0] if catalysts else None
        ligand = ligands[0] if ligands else None
        base = bases[0] if bases else None
        solvent = solvents[0] if solvents else None
        
        # Create fingerprint
        fingerprint = (catalyst, ligand, base, solvent)
        
        # Group reactions
        condition_groups[fingerprint]['reactions'].append(rxn)
        
        y = rxn.get('conditions', {}).get('yield_pct')
        if y is not None:
            condition_groups[fingerprint]['yields'].append(y)
        
        t = rxn.get('conditions', {}).get('temperature_c')
        if t is not None:
            condition_groups[fingerprint]['temps'].append(t)
        
        time = rxn.get('conditions', {}).get('time_h')
        if time is not None:
            condition_groups[fingerprint]['times'].append(time)
    
    # Score and rank conditions
    recommendations = []
    for fingerprint, data in condition_groups.items():
        catalyst, ligand, base, solvent = fingerprint
        
        frequency = len(data['reactions'])
        avg_yield = mean(data['yields']) if data['yields'] else None
        avg_temp = mean(data['temps']) if data['temps'] else None
        avg_time = mean(data['times']) if data['times'] else None
        
        # Skip if missing critical data
        if avg_yield is None:
            continue
        
        # Compute score based on strategy
        if optimize_for == 'frequency':
            score = frequency
        elif optimize_for == 'yield':
            score = avg_yield
        else:  # diversity
            # Balanced score: yield * log(frequency)
            import math
            score = avg_yield * math.log(frequency + 1)
        
        recommendations.append({
            'condition_id': f"{family}_{len(recommendations)+1}",
            'catalyst': catalyst,
            'ligand': ligand,
            'base': base,
            'solvent': solvent,
            'temperature_c': avg_temp,
            'time_h': avg_time,
            'avg_yield': avg_yield,
            'frequency': frequency,
            'score': score
        })
    
    # Sort by score and take top N
    recommendations.sort(key=lambda x: x['score'], reverse=True)
    
    return recommendations[:n_conditions]


def get_all_families() -> List[str]:
    """Get list of all available reaction families.
    
    Returns:
        List of family names (JSONL file basenames)
    """
    data_dir = get_data_dir()
    if not data_dir.exists():
        return []
    
    families = []
    for file in data_dir.glob("*.jsonl"):
        families.append(file.stem)
    
    return sorted(families)


def print_analytics_summary(family: str, top_n: int = 10):
    """Print a comprehensive analytics summary for a reaction family.
    
    Args:
        family: Reaction family name
        top_n: Number of top items to display
    """
    print("=" * 80)
    print(f" DATASET ANALYTICS: {family}")
    print("=" * 80)
    print()
    
    # Basic stats
    stats = get_dataset_stats(family)
    print("[DATASET STATISTICS]")
    print(f"  Total reactions: {stats['total_reactions']:,}")
    print(f"  Unique condition cores: {stats['unique_condition_cores']}")
    print(f"  Unique solvents: {stats['unique_solvents']}")
    print(f"  Unique bases: {stats['unique_bases']}")
    print(f"  Unique catalysts: {stats['unique_catalysts']}")
    print()
    
    if stats['yield_stats']:
        ys = stats['yield_stats']
        print(f"  Yield data: {ys['count']}/{stats['total_reactions']} reactions ({ys['count']/stats['total_reactions']*100:.1f}%)")
        print(f"    Range: {ys['min']:.1f}% - {ys['max']:.1f}%")
        print(f"    Mean: {ys['mean']:.1f}%")
        print(f"    Median: {ys['median']:.1f}%")
        print()
    
    if stats['temperature_stats']:
        ts = stats['temperature_stats']
        print(f"  Temperature data: {ts['count']}/{stats['total_reactions']} reactions ({ts['count']/stats['total_reactions']*100:.1f}%)")
        print(f"    Range: {ts['min']:.0f}°C - {ts['max']:.0f}°C")
        print(f"    Mean: {ts['mean']:.0f}°C")
        print()
    
    # Top catalysts
    print(f"[TOP {top_n} CATALYSTS]")
    catalysts = get_common_catalysts(family, top_n=top_n)
    if catalysts:
        for name, count, avg_yield in catalysts:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No catalyst data found")
    print()
    
    # Top ligands
    print(f"[TOP {top_n} LIGANDS]")
    ligands = get_common_ligands(family, top_n=top_n)
    if ligands:
        for name, count, avg_yield in ligands:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No ligand data found")
    print()
    
    # Top bases
    print(f"[TOP {top_n} BASES]")
    bases = get_common_bases(family, top_n=top_n)
    if bases:
        for name, count, avg_yield in bases:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No base data found")
    print()
    
    # Top solvents
    print(f"[TOP {top_n} SOLVENTS]")
    solvents = get_common_solvents(family, top_n=top_n)
    if solvents:
        for name, count, avg_yield in solvents:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No solvent data found")
    print()
    
    # Top condition cores
    print(f"[TOP {top_n} CONDITION CORES]")
    cores = get_condition_cores(family, top_n=top_n)
    if cores:
        for core, count, avg_yield in cores:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            core_display = core[:60] + "..." if len(core) > 60 else core
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {core_display}")
    else:
        print("  No condition core data found")
    print()
    
    print("=" * 80)
