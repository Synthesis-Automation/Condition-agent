"""Generate comprehensive Suzuki reaction conditions report."""

import sys
import os
sys.path.insert(0, '.')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from chemtools.scdb_matcher.matcher import match
from chemtools.scdb_matcher.loader import load_db

# Import test reactions
import importlib.util
spec = importlib.util.spec_from_file_location("sample_reactions", "tests/sample_reactions.py")
sample_reactions = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sample_reactions)
SUZUKI_DB_TEST_REACTIONS = sample_reactions.SUZUKI_DB_TEST_REACTIONS

# Load database
db = load_db("data/conditionDB/suzuki_db.json")

print("="*80)
print("SUZUKI COUPLING REACTION CONDITIONS REPORT")
print("="*80)
print(f"\nGenerated: October 5, 2025")
print(f"Database: suzuki_db.json ({len(db.entries)} entries)")
print(f"Test Cases: {len(SUZUKI_DB_TEST_REACTIONS)}")
print("\n" + "="*80)

# Process each test case
results = []
for rule_id, test_data in SUZUKI_DB_TEST_REACTIONS.items():
    smiles = test_data['smiles']
    description = test_data['description']
    
    # Match the reaction
    result = match(db, rxn_smiles=smiles)
    
    # Extract conditions
    matched_id = result.entry_id
    priority = result.priority
    conditions = result.conditions
    
    results.append({
        'test_id': rule_id,
        'description': description,
        'smiles': smiles,
        'matched_id': matched_id,
        'priority': priority,
        'conditions': conditions,
        'pass': matched_id == rule_id
    })

# Generate report
print("\n## RECOMMENDED CONDITIONS BY REACTION TYPE")
print("="*80)

for i, r in enumerate(results, 1):
    print(f"\n### {i}. {r['description']}")
    print(f"**Test Case ID:** `{r['test_id']}`")
    print(f"**SMILES:** `{r['smiles']}`")
    print(f"**Matched Rule:** `{r['matched_id']}` (Priority: {r['priority']})")
    
    if r['pass']:
        print(f"**Status:** [PASS] - Correct rule matched")
    else:
        print(f"**Status:** [FAIL] Expected `{r['test_id']}` but matched `{r['matched_id']}`")
    
    print(f"\n**Recommended Conditions:**")
    
    cond = r['conditions']
    
    # Pd source
    if 'pd_source' in cond and cond['pd_source']:
        pd_sources = cond['pd_source'] if isinstance(cond['pd_source'], list) else [cond['pd_source']]
        print(f"  - **Pd Source:** {', '.join(pd_sources)}")
    
    # Ligands
    if 'ligands' in cond and cond['ligands']:
        ligands = cond['ligands'] if isinstance(cond['ligands'], list) else [cond['ligands']]
        print(f"  - **Ligand:** {', '.join(ligands)}")
    
    # Boron partner
    if 'boron_partner' in cond and cond['boron_partner']:
        boron = cond['boron_partner'] if isinstance(cond['boron_partner'], list) else [cond['boron_partner']]
        print(f"  - **Boron Partner:** {', '.join(boron)}")
    
    # Base
    if 'base' in cond and cond['base']:
        bases = cond['base'] if isinstance(cond['base'], list) else [cond['base']]
        print(f"  - **Base:** {', '.join(bases)}")
    
    # Solvent
    if 'solvent' in cond and cond['solvent']:
        solvents = cond['solvent'] if isinstance(cond['solvent'], list) else [cond['solvent']]
        print(f"  - **Solvent:** {', '.join(solvents)}")
    
    # Temperature
    if 'temperature_C' in cond and cond['temperature_C']:
        temps = cond['temperature_C'] if isinstance(cond['temperature_C'], list) else [cond['temperature_C']]
        if len(temps) == 2:
            print(f"  - **Temperature:** {temps[0]}–{temps[1]}°C")
        elif len(temps) == 1:
            print(f"  - **Temperature:** {temps[0]}°C")
        else:
            print(f"  - **Temperature:** {', '.join(map(str, temps))}°C")
    
    # Time
    if 'time_h' in cond and cond['time_h']:
        times = cond['time_h'] if isinstance(cond['time_h'], list) else [cond['time_h']]
        if len(times) == 2:
            print(f"  - **Time:** {times[0]}–{times[1]} hours")
        elif len(times) == 1:
            print(f"  - **Time:** {times[0]} hours")
    
    # Loadings
    if 'loadings' in cond and cond['loadings']:
        loadings = cond['loadings']
        if 'Pd_mol%' in loadings:
            pd_mol = loadings['Pd_mol%']
            if isinstance(pd_mol, list) and len(pd_mol) == 2:
                print(f"  - **Pd Loading:** {pd_mol[0]}–{pd_mol[1]} mol%")
            else:
                print(f"  - **Pd Loading:** {pd_mol} mol%")
        
        if 'ligand_mol%' in loadings:
            lig_mol = loadings['ligand_mol%']
            if isinstance(lig_mol, list) and len(lig_mol) == 2:
                print(f"  - **Ligand Loading:** {lig_mol[0]}–{lig_mol[1]} mol%")
        
        if 'base_eq' in loadings:
            base_eq = loadings['base_eq']
            if isinstance(base_eq, list) and len(base_eq) == 2:
                print(f"  - **Base Equivalents:** {base_eq[0]}–{base_eq[1]} eq")

# Summary statistics
print("\n\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

passing = sum(1 for r in results if r['pass'])
failing = len(results) - passing

print(f"\nTotal Test Cases: {len(results)}")
print(f"Passing: {passing} ({passing/len(results)*100:.1f}%)")
print(f"Failing: {failing} ({failing/len(results)*100:.1f}%)")

print("\n### Rules Successfully Matched:")
for r in results:
    if r['pass']:
        print(f"  [PASS] {r['test_id']}")

print("\n### Rules Not Matched (Using Alternative Conditions):")
for r in results:
    if not r['pass']:
        print(f"  [FAIL] {r['test_id']} -> matched {r['matched_id']} instead")

print("\n" + "="*80)
print("END OF REPORT")
print("="*80)
