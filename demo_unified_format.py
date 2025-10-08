"""
Quick demo of the unified output formatter.

Shows how to convert ML-based recommendations to standardized JSON format
suitable for robotic execution systems.
"""

import json
import time
from chemtools.recommend.core import recommend_from_reaction
from chemtools.reaction_type_detector import detect_reaction_type
from chemtools.output_formatter import (
    format_ml_output,
    convert_raw_recommendation_to_standard,
    create_standard_output,
)


def demo_ullmann_coupling():
    """
    Demonstrate unified output format for Ullmann C-N coupling.
    """
    print("=" * 70)
    print("DEMO: Unified Output Format for Robotic Execution")
    print("=" * 70)
    
    # Test reaction
    smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    print(f"\n📋 Input Reaction: {smiles}")
    print(f"   Expected: C-N coupling (Ullmann or Buchwald-Hartwig)")
    
    # Run ML recommendation
    print("\n" + "-" * 70)
    print("Step 1: Running ML-based Recommendation")
    print("-" * 70)
    
    start = time.time()
    ml_result = recommend_from_reaction(smiles, k=50)
    elapsed_ms = (time.time() - start) * 1000
    
    print(f"\n✅ Complete in {elapsed_ms:.2f} ms")
    print(f"\n📊 Raw ML Output Structure:")
    print(f"   Keys: {list(ml_result.keys())}")
    print(f"   Family: {ml_result.get('family')}")
    print(f"   Recommendation core: {ml_result.get('recommendation', {}).get('core')}")
    
    # Check if formatted output exists
    if 'formatted' in ml_result:
        print(f"   Formatted output available: ✅")
        formatted = ml_result['formatted']
        if 'recommendations' in formatted:
            print(f"   Number of recommendations: {len(formatted['recommendations'])}")
    
    # Detect reaction type (use ML result if detector fails)
    detection = detect_reaction_type(smiles)
    if detection and detection.get("success"):
        detected_type = detection.get("mapped_family") or detection.get("rxn_class", "Unknown")
        confidence = detection.get("confidence") or 0.8
    else:
        # Fall back to ML family
        detected_type = ml_result.get('family', "Unknown")
        confidence = ml_result.get('recommendation', {}).get('confidence', 0.8)
    
    print(f"\n🔍 Detection: {detected_type} (confidence: {confidence:.3f})")
    
    # Try to convert to standardized format
    print("\n" + "-" * 70)
    print("Step 2: Examining Current Output Format")
    print("-" * 70)
    
    try:
        if 'formatted' in ml_result:
            formatted = ml_result['formatted']
            print("\n✅ Formatted output structure:")
            print(f"   Keys: {list(formatted.keys())}")
            
            # Save the current formatted output
            current_format_file = "output_ml_current_format.json"
            with open(current_format_file, "w") as f:
                json.dump(formatted, f, indent=2)
            print(f"\n💾 Saved current format to: {current_format_file}")
            
            # Check if there are recommended_conditions
            if 'recommended_conditions' in formatted:
                recs = formatted['recommended_conditions']
                print(f"\n📊 Recommended Conditions ({len(recs)} variants):")
                
                for i, rec in enumerate(recs[:3], 1):
                    print(f"\n  Variant {i}:")
                    print(f"    Keys: {list(rec.keys())}")
                    if 'chemicals' in rec:
                        print(f"    Chemicals: {len(rec['chemicals'])} reagents")
                        # Show first few chemicals
                        for chem in rec['chemicals'][:3]:
                            role = chem.get('role', 'unknown')
                            name = chem.get('name', 'N/A')
                            print(f"      - {role}: {name}")
                    if 'conditions' in rec:
                        cond = rec['conditions']
                        print(f"    Conditions: {list(cond.keys())}")
                
                # Print sample JSON
                print("\n" + "-" * 70)
                print("Sample Current Format (Variant #1)")
                print("-" * 70)
                rec1_json = json.dumps(recs[0], indent=2)
                if len(rec1_json) > 1500:
                    print(rec1_json[:1500] + "\n  ... (truncated)")
                else:
                    print(rec1_json)
            
            print("\n" + "=" * 70)
            print("✅ DEMO COMPLETE - Current format examined")
            print("=" * 70)
            print("\n📝 Next Steps:")
            print("   1. Current format saved to output_ml_current_format.json")
            print("   2. Review the structure")
            print("   3. Unified formatter can convert this to standard robot format")
            
        else:
            print("\n⚠️  No formatted output in ML result")
            print(f"   Available keys: {list(ml_result.keys())}")
            
    except Exception as e:
        print(f"\n❌ Examination failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    success = demo_ullmann_coupling()
    exit(0 if success else 1)
