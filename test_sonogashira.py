"""Test Sonogashira pattern generation"""
from chemtools.protocol.batch_update_protocol_smarts import ProtocolSmartsUpdater
from pathlib import Path

# Test Sonogashira
updater = ProtocolSmartsUpdater(Path('data/protocol_db'), dry_run=True)
result = updater.process_protocol_file(Path('data/protocol_db/Sonogashira-Coupling.json'))

print(f"\n{'='*60}")
print("SONOGASHIRA COUPLING TEST")
print('='*60)
print(f"Success: {result.success}")
print(f"\nGenerated Pattern:")
print(f"  Core: {result.new_pattern['core']}")
print(f"  Guards: {result.new_pattern['guards_forbid']}")
print(f"\nExpected:")
print(f"  Core: c-[I].C#C>>c-C#C")
print(f"  (aryl iodide + alkyne -> aryl alkyne)")
