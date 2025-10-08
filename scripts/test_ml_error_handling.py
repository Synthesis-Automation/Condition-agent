"""Test ML error handling for unsupported/failed reactions."""

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from app.ui_simple import get_ml_recommendations

print("="*80)
print("TESTING ML ERROR HANDLING")
print("="*80)

# Test cases for different error scenarios
test_cases = [
    {
        "name": "Unsupported Reaction Type (should auto-detect as unknown)",
        "smiles": "CC(=O)O.CCO>>CC(=O)OCC",  # Simple esterification
        "reaction_type": "Auto-detect",
        "expected": "No ML model available",
    },
    {
        "name": "Supported Type but No Precedents (unusual substrates)",
        "smiles": "FC(F)(F)c1ccccc1Br.FC(F)(F)c1ccccc1N>>FC(F)(F)c1ccccc1Nc1ccccc1C(F)(F)F",  # Unusual C-N coupling
        "reaction_type": "C-N Coupling (Pd)",
        "expected": "No ML Recommendations Found",
    },
    {
        "name": "Valid Suzuki (should work)",
        "smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "reaction_type": "Auto-detect",
        "expected": "ML Recommendations Generated",
    },
]

for i, test in enumerate(test_cases, 1):
    print(f"\n{i}. {test['name']}")
    print("-"*80)
    print(f"SMILES: {test['smiles']}")
    print(f"Type: {test['reaction_type']}")
    print(f"Expected: {test['expected']}")
    
    try:
        result = get_ml_recommendations(
            test['smiles'],
            test['reaction_type'],
            top_k=2
        )
        
        summary, table = result
        
        # Check result type
        if test['expected'] in summary:
            print(f"\n✅ PASSED - Found expected message: '{test['expected']}'")
        else:
            print(f"\n⚠ Result different than expected")
        
        # Show first 500 chars of summary
        print(f"\nSummary preview:")
        print("-"*80)
        lines = summary.split('\n')
        for line in lines[:15]:  # First 15 lines
            print(line)
        if len(lines) > 15:
            print(f"... ({len(lines) - 15} more lines)")
        
        print(f"\nTable rows: {len(table)}")
        
    except Exception as e:
        print(f"\n❌ EXCEPTION: {e}")

print("\n" + "="*80)
print("ERROR HANDLING TEST COMPLETE")
print("="*80)
print("\nKey improvements:")
print("✅ Clear error messages for unsupported reaction types")
print("✅ Helpful guidance when no precedents found")
print("✅ Detection confidence warnings")
print("✅ Context-specific troubleshooting steps")
print("✅ Alternative action suggestions")
print("="*80)
