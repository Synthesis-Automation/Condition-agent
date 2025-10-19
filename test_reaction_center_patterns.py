"""Test reaction-center-focused pattern generation"""
from chemtools.protocol.batch_update_protocol_smarts import ProtocolSmartsUpdater
from pathlib import Path

test_files = [
    'Sonogashira-Coupling.json',
    'Pd_Buchwald_Arylsulfonate_Amination_CMPhos.json',
    'Alkyl_Iodide_Borylation.json',
    'Aryl mesylate_Suzuki.json',
]

updater = ProtocolSmartsUpdater(Path('data/protocol_db'), dry_run=True)

print("\n" + "="*80)
print("REACTION-CENTER-FOCUSED PATTERN GENERATION TEST")
print("="*80)

for filename in test_files:
    filepath = Path('data/protocol_db') / filename
    if not filepath.exists():
        print(f"\nSkipping {filename} (not found)")
        continue
    
    print(f"\n{filename}")
    print("-" * 80)
    
    result = updater.process_protocol_file(filepath)
    
    if result.success:
        print(f"  Core: {result.new_pattern['core']}")
        print(f"  Guards: {len(result.new_pattern.get('guards_forbid', []))} patterns")
    else:
        print(f"  FAILED: {result.message}")
