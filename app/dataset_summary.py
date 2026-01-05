#!/usr/bin/env python3
"""
Generate comprehensive summary of all reaction datasets.
Run with: python app/dataset_summary.py
"""

import sys
import json
from pathlib import Path
from statistics import mean, median
from typing import Any, Dict, Iterable, List, Optional, Tuple

# Add parent directory to path so chemtools can be imported
ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

CATALYST_ROLES = {"CATALYST", "METAL_CATALYST", "METAL_PRECURSOR", "ORGANO_CATALYST"}
LIGAND_ROLES = {"LIGAND"}
BASE_ROLES = {"BASE"}
SOLVENT_ROLES = {"SOLVENT"}
CONDENSATION_ROLES = {"CONDENSATION_AGENT", "COUPLING_REAGENT"}


def _get_data_dir() -> Path:
    """Return the reaction dataset directory."""
    return ROOT_DIR / "data" / "reaction_dataset"


def _get_all_families() -> List[str]:
    """Return dataset family names based on JSONL filenames."""
    data_dir = _get_data_dir()
    if not data_dir.exists():
        return []
    return sorted(path.stem for path in data_dir.glob("*.jsonl"))


def _iter_dataset_records(family: str) -> Iterable[Dict[str, Any]]:
    """Yield parsed records for a reaction family."""
    dataset_path = _get_data_dir() / f"{family}.jsonl"
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset not found: {dataset_path}")

    with open(dataset_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(rec, dict):
                yield rec


def _coerce_float(value: Any) -> Optional[float]:
    """Convert value to float when possible."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None
    return None


def _normalize_name(value: Any) -> Optional[str]:
    """Normalize a display name into a non-empty string."""
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _extract_condition_list(conditions: Dict[str, Any], keys: Iterable[str]) -> List[str]:
    """Extract name lists from conditions for supported keys."""
    for key in keys:
        if key not in conditions:
            continue
        raw = conditions.get(key)
        if raw is None:
            return []
        if isinstance(raw, list):
            values = raw
        else:
            values = [raw]
        names = []
        for item in values:
            name = _normalize_name(item)
            if name:
                names.append(name)
        return names
    return []


def _extract_reagent_names(record: Dict[str, Any], roles: Iterable[str]) -> List[str]:
    """Extract reagent names that match any of the provided roles."""
    role_set = set(roles)
    names: List[str] = []
    reagents = record.get("reagents")
    if not isinstance(reagents, list):
        return names
    for reagent in reagents:
        if not isinstance(reagent, dict):
            continue
        role = reagent.get("role")
        if role not in role_set:
            continue
        name = _normalize_name(reagent.get("name"))
        if name:
            names.append(name)
    return names


def _extract_solvent_names(record: Dict[str, Any], conditions: Dict[str, Any]) -> List[str]:
    """Extract solvent names from record and conditions."""
    names: List[str] = []
    solvents = record.get("solvents")
    if isinstance(solvents, list):
        for solvent in solvents:
            if isinstance(solvent, dict):
                name = _normalize_name(solvent.get("name"))
            else:
                name = _normalize_name(solvent)
            if name:
                names.append(name)
    names.extend(_extract_reagent_names(record, SOLVENT_ROLES))
    names.extend(_extract_condition_list(conditions, ["solvent", "solvents"]))
    return names


def _extract_condition_values(record: Dict[str, Any]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Return yield, temperature, and time values from supported fields."""
    conditions = record.get("conditions")
    if not isinstance(conditions, dict):
        conditions = {}

    yield_value = _coerce_float(conditions.get("yield_pct"))
    if yield_value is None:
        yield_value = _coerce_float(conditions.get("yield"))
    if yield_value is None:
        yield_value = _coerce_float(record.get("yield"))

    temp_value = None
    for key in ("temperature_c", "temperature_C", "temp_c", "temp_C"):
        temp_value = _coerce_float(conditions.get(key))
        if temp_value is not None:
            break
    if temp_value is None:
        temp_value = _coerce_float(record.get("temperature_c"))

    time_value = None
    for key in ("time_h", "time_hr", "time_hours"):
        time_value = _coerce_float(conditions.get(key))
        if time_value is not None:
            break
    if time_value is None:
        time_value = _coerce_float(record.get("time_h"))

    return yield_value, temp_value, time_value


def _build_catalytic_system(
    record: Dict[str, Any],
    catalysts: Iterable[str],
    ligands: Iterable[str],
) -> Optional[str]:
    """Build a catalytic system label."""
    cat_sys = record.get("catalytic_system")
    if isinstance(cat_sys, list):
        components = []
        for comp in cat_sys:
            if isinstance(comp, dict):
                name = _normalize_name(comp.get("name"))
            else:
                name = _normalize_name(comp)
            if name:
                components.append(name)
        if components:
            return " + ".join(components)
    elif isinstance(cat_sys, str):
        text = cat_sys.strip()
        if text:
            return text

    combined = list(dict.fromkeys([*sorted(catalysts), *sorted(ligands)]))
    if combined:
        return " + ".join(combined)
    return None


def _build_condition_core(
    record: Dict[str, Any],
    catalysts: Iterable[str],
    ligands: Iterable[str],
    bases: Iterable[str],
    solvents: Iterable[str],
    additives: Iterable[str],
) -> Optional[str]:
    """Return condition core or a derived label when missing."""
    core = _normalize_name(record.get("condition_core"))
    if core:
        return core

    parts = []
    catalysts_list = sorted(catalysts)
    ligands_list = sorted(ligands)
    bases_list = sorted(bases)
    solvents_list = sorted(solvents)
    additives_list = sorted(additives)

    if catalysts_list:
        parts.append(f"Catalyst: {', '.join(catalysts_list)}")
    if ligands_list:
        parts.append(f"Ligand: {', '.join(ligands_list)}")
    if bases_list:
        parts.append(f"Base: {', '.join(bases_list)}")
    if solvents_list:
        parts.append(f"Solvent: {', '.join(solvents_list)}")
    if additives_list:
        parts.append(f"Additive: {', '.join(additives_list)}")

    return " | ".join(parts) if parts else None


def _update_item_stats(
    store: Dict[str, Dict[str, Any]],
    names: Iterable[str],
    yield_value: Optional[float],
) -> None:
    """Update per-item counts and yield lists."""
    for name in set(names):
        entry = store.setdefault(name, {"count": 0, "yields": []})
        entry["count"] += 1
        if yield_value is not None:
            entry["yields"].append(yield_value)


def _update_reagent_stats(
    store: Dict[Tuple[str, str], Dict[str, Any]],
    reagents: Iterable[Tuple[str, str]],
    yield_value: Optional[float],
) -> None:
    """Update reagent stats keyed by (name, role)."""
    for name, role in set(reagents):
        entry = store.setdefault((name, role), {"count": 0, "yields": []})
        entry["count"] += 1
        if yield_value is not None:
            entry["yields"].append(yield_value)


def _stats(values: List[float]) -> Optional[Dict[str, float]]:
    """Return summary stats for numeric values."""
    if not values:
        return None
    return {
        "count": len(values),
        "min": min(values),
        "max": max(values),
        "mean": mean(values),
        "median": median(values),
    }


def _finalize_item_stats(store: Dict[str, Dict[str, Any]]) -> List[Tuple[str, int, Optional[float]]]:
    """Return sorted item stats with average yields."""
    results: List[Tuple[str, int, Optional[float]]] = []
    for name, data in store.items():
        avg_yield = mean(data["yields"]) if data["yields"] else None
        results.append((name, data["count"], avg_yield))
    results.sort(key=lambda x: x[1], reverse=True)
    return results


def _finalize_reagent_stats(
    store: Dict[Tuple[str, str], Dict[str, Any]]
) -> List[Tuple[str, str, int, Optional[float]]]:
    """Return sorted reagent stats with average yields."""
    results: List[Tuple[str, str, int, Optional[float]]] = []
    for (name, role), data in store.items():
        avg_yield = mean(data["yields"]) if data["yields"] else None
        results.append((name, role, data["count"], avg_yield))
    results.sort(key=lambda x: x[2], reverse=True)
    return results


def _analyze_family(family: str, top_n: int = 20) -> Dict[str, Any]:
    """Analyze a reaction family and return summary stats plus rankings."""
    total_reactions = 0
    condition_cores: set[str] = set()
    catalysts_set: set[str] = set()
    ligands_set: set[str] = set()
    bases_set: set[str] = set()
    solvents_set: set[str] = set()
    yields: List[float] = []
    temps: List[float] = []
    times: List[float] = []

    system_stats: Dict[str, Dict[str, Any]] = {}
    catalyst_stats: Dict[str, Dict[str, Any]] = {}
    ligand_stats: Dict[str, Dict[str, Any]] = {}
    base_stats: Dict[str, Dict[str, Any]] = {}
    solvent_stats: Dict[str, Dict[str, Any]] = {}
    core_stats: Dict[str, Dict[str, Any]] = {}
    reagent_stats: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for record in _iter_dataset_records(family):
        total_reactions += 1
        conditions = record.get("conditions")
        if not isinstance(conditions, dict):
            conditions = {}

        yield_value, temp_value, time_value = _extract_condition_values(record)
        if yield_value is not None:
            yields.append(yield_value)
        if temp_value is not None:
            temps.append(temp_value)
        if time_value is not None:
            times.append(time_value)

        catalysts = set(_extract_reagent_names(record, CATALYST_ROLES))
        ligands = set(_extract_reagent_names(record, LIGAND_ROLES))
        bases = set(_extract_reagent_names(record, BASE_ROLES))
        solvents = set(_extract_solvent_names(record, conditions))
        additives = set(_extract_condition_list(conditions, ["additive", "additives"]))

        bases.update(_extract_condition_list(conditions, ["base", "bases"]))

        cat_sys = record.get("catalytic_system")
        if isinstance(cat_sys, list):
            for comp in cat_sys:
                if not isinstance(comp, dict):
                    continue
                name = _normalize_name(comp.get("name"))
                if not name:
                    continue
                role = comp.get("role")
                if role in CATALYST_ROLES:
                    catalysts.add(name)
                if role in LIGAND_ROLES:
                    ligands.add(name)

        if not catalysts and isinstance(cat_sys, str):
            cat_sys_name = _normalize_name(cat_sys)
            if cat_sys_name and cat_sys_name not in solvents:
                catalysts.add(cat_sys_name)

        catalysts_set.update(catalysts)
        ligands_set.update(ligands)
        bases_set.update(bases)
        solvents_set.update(solvents)

        _update_item_stats(catalyst_stats, catalysts, yield_value)
        _update_item_stats(ligand_stats, ligands, yield_value)
        _update_item_stats(base_stats, bases, yield_value)
        _update_item_stats(solvent_stats, solvents, yield_value)

        system_label = _build_catalytic_system(record, catalysts, ligands)
        if system_label:
            _update_item_stats(system_stats, [system_label], yield_value)

        core_label = _build_condition_core(record, catalysts, ligands, bases, solvents, additives)
        if core_label:
            condition_cores.add(core_label)
            _update_item_stats(core_stats, [core_label], yield_value)

        reagent_entries = []
        reagents = record.get("reagents")
        if isinstance(reagents, list):
            for reagent in reagents:
                if not isinstance(reagent, dict):
                    continue
                role = reagent.get("role")
                if role not in CONDENSATION_ROLES:
                    continue
                name = _normalize_name(reagent.get("name"))
                if name:
                    reagent_entries.append((name, role))
        _update_reagent_stats(reagent_stats, reagent_entries, yield_value)

    stats = {
        "family": family,
        "total_reactions": total_reactions,
        "unique_condition_cores": len(condition_cores),
        "unique_catalysts": len(catalysts_set),
        "unique_ligands": len(ligands_set),
        "unique_bases": len(bases_set),
        "unique_solvents": len(solvents_set),
        "yield_stats": _stats(yields),
        "temperature_stats": _stats(temps),
        "time_stats": _stats(times),
    }

    return {
        "stats": stats,
        "systems": _finalize_item_stats(system_stats)[:top_n],
        "catalysts": _finalize_item_stats(catalyst_stats)[:top_n],
        "ligands": _finalize_item_stats(ligand_stats)[:top_n],
        "bases": _finalize_item_stats(base_stats)[:top_n],
        "solvents": _finalize_item_stats(solvent_stats)[:top_n],
        "cores": _finalize_item_stats(core_stats)[:top_n],
        "condensation_agents": _finalize_reagent_stats(reagent_stats)[:top_n],
    }


def generate_summary() -> None:
    """Generate and print comprehensive dataset summary."""
    
    families = _get_all_families()
    family_analytics: Dict[str, Dict[str, Any]] = {}
    
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
        analytics = _analyze_family(family, top_n=20)
        family_analytics[family] = analytics
        stats = analytics["stats"]
        family_stats[family] = stats
        total_reactions += stats['total_reactions']
        
        print(f"📊 {family}")
        print(f"   Total reactions: {stats['total_reactions']:,}")
        print(f"   Unique condition cores: {stats['unique_condition_cores']}")
        print(f"   Unique catalysts: {stats['unique_catalysts']}")
        print(f"   Unique ligands: {stats['unique_ligands']}")
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
        print("[TOP 20 CATALYTIC SYSTEMS]")
        systems = family_analytics[family]["systems"]
        if systems:
            for system, count, avg_yield in systems:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {system}")
        else:
            print("  No catalytic system data")
        print()
        
        # Top catalysts
        print("[TOP 20 CATALYSTS]")
        catalysts = family_analytics[family]["catalysts"]
        if catalysts:
            for name, count, avg_yield in catalysts:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
        else:
            print("  No catalyst data")
        print()
        
        # Top ligands
        print("[TOP 20 LIGANDS]")
        ligands = family_analytics[family]["ligands"]
        if ligands:
            for name, count, avg_yield in ligands:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
        else:
            print("  No ligand data")
        print()
        
        # Top bases
        print("[TOP 20 BASES]")
        bases = family_analytics[family]["bases"]
        if bases:
            for name, count, avg_yield in bases:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
        else:
            print("  No base data")
        print()
        
        # Top solvents
        print("[TOP 20 SOLVENTS]")
        solvents = family_analytics[family]["solvents"]
        if solvents:
            for name, count, avg_yield in solvents:
                yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
        else:
            print("  No solvent data")
        print()
        
        # Top condensation agents (for amide formation)
        if family == "Amide_formation":
            print("[TOP 20 CONDENSATION AGENTS]")
            condensation_agents = family_analytics[family]["condensation_agents"]
            if condensation_agents:
                for name, role, count, avg_yield in condensation_agents:
                    yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
                    print(f"  {count:>5}× | Avg yield: {yield_str:>6} | {name}")
            else:
                print("  No condensation agent data")
            print()
        
        # Top condition cores
        print("[TOP 20 CONDITION CORES]")
        cores = family_analytics[family]["cores"]
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
