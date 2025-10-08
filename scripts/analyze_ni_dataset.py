"""
Analyze C_N_coupling_Ni dataset structure and vocabulary.

Quick analysis script to understand the Ni-catalyzed C-N coupling dataset.
"""

import json
from collections import Counter
from pathlib import Path

# Load dataset
dataset_path = Path("data/reaction_dataset/C_N_Coupling_Ni.jsonl")

print("=" * 80)
print("C-N COUPLING NI DATASET ANALYSIS")
print("=" * 80)
print()

reactions = []
with open(dataset_path, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
            reactions.append(rec)
        except:
            continue

print(f"Total reactions: {len(reactions)}")
print()

# Extract vocabulary
cores = []
bases = []
oxidants = []
solvents = []
yields = []
temps = []
times = []

for rxn in reactions:
    # Core
    core = rxn.get('condition_core', '')
    if core:
        cores.append(core)
    
    # Reagents
    for reagent in rxn.get('reagents', []):
        role = reagent.get('role', '')
        cas = reagent.get('cas', '')
        name = reagent.get('name', '')
        
        if role == 'BASE' and cas:
            bases.append(f"{name} ({cas})")
        elif role == 'OXIDANT' and cas:
            oxidants.append(f"{name} ({cas})")
    
    # Solvents
    for solv in rxn.get('solvents', []):
        cas = solv.get('cas', '')
        name = solv.get('name', '')
        if cas:
            solvents.append(f"{name} ({cas})")
    
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

print(f"Cores: {len(set(cores))} unique")
print(f"Bases: {len(set(bases))} unique")
print(f"Oxidants: {len(set(oxidants))} unique")
print(f"Solvents: {len(set(solvents))} unique")
print()

print(f"Yield data: {len(yields)}/{len(reactions)} ({len(yields)/len(reactions)*100:.1f}%)")
if yields:
    print(f"  Range: {min(yields):.1f}% - {max(yields):.1f}%")
    print(f"  Average: {sum(yields)/len(yields):.1f}%")
    print(f"  Median: {sorted(yields)[len(yields)//2]:.1f}%")
print()

print(f"Temperature data: {len(temps)}/{len(reactions)} ({len(temps)/len(reactions)*100:.1f}%)")
if temps:
    print(f"  Range: {min(temps):.0f}°C - {max(temps):.0f}°C")
    print(f"  Average: {sum(temps)/len(temps):.0f}°C")
print()

print(f"Time data: {len(times)}/{len(reactions)} ({len(times)/len(reactions)*100:.1f}%)")
if times:
    print(f"  Range: {min(times):.1f}h - {max(times):.1f}h")
    print(f"  Average: {sum(times)/len(times):.1f}h")
print()

# Top cores
core_counts = Counter(cores)
print("Top 10 cores:")
for core, count in core_counts.most_common(10):
    print(f"  {count:>4} - {core}")
print()

# Top bases
base_counts = Counter(bases)
print("Top 10 bases:")
for base, count in base_counts.most_common(10):
    print(f"  {count:>4} - {base}")
print()

# Top oxidants
oxidant_counts = Counter(oxidants)
print("Top 10 oxidants:")
for oxidant, count in oxidant_counts.most_common(10):
    print(f"  {count:>4} - {oxidant}")
print()

# Top solvents
solvent_counts = Counter(solvents)
print("Top 10 solvents:")
for solv, count in solvent_counts.most_common(10):
    print(f"  {count:>4} - {solv}")
print()

print("=" * 80)
print("COMPARISON: Ni vs Cu (Ullmann) vs Pd (Buchwald)")
print("=" * 80)
print()
print(f"{'Property':<25} {'Ni':<15} {'Cu (Ullmann)':<15} {'Pd (Buchwald)':<15}")
print("-" * 70)
print(f"{'Dataset size':<25} {len(reactions):<15} {'5,552':<15} {'1,343':<15}")
print(f"{'Unique cores':<25} {len(set(cores)):<15} {'27':<15} {'37':<15}")
print(f"{'Unique bases':<25} {len(set(bases)):<15} {'28':<15} {'20':<15}")
print(f"{'Unique solvents':<25} {len(set(solvents)):<15} {'49':<15} {'25':<15}")
if yields:
    print(f"{'Avg yield':<25} {f'{sum(yields)/len(yields):.1f}%':<15} {'74.2%':<15} {'~73%':<15}")
if temps:
    print(f"{'Avg temperature':<25} {f'{sum(temps)/len(temps):.0f}°C':<15} {'110°C':<15} {'100°C':<15}")

print()
print("Key Differences:")
print("  - Ni uses oxidants (Cu/Pd don't require them)")
print("  - Smaller dataset than Cu Ullmann, but larger than Pd Buchwald")
print("  - Ni appears to use more complex solvent systems (multiple solvents)")
print()

