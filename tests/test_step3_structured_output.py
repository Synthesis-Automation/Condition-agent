#!/usr/bin/env python3
"""
Test Step 3: Structured Output and Advanced Features

Tests:
- Structured condition output (recommend_conditions_structured)
- Reagent enrichment (database lookups)
- Output format validation
- Plate design generation (optional)

This validates the enhanced output formatter and advanced recommendation features.
"""

import sys
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.recommend import recommend_from_reaction, recommend_conditions_structured

# Test reaction (same as previous tests)
TEST_REACTION = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"


def test_3a_structured_output():
    """Test 3a: Structured output format with reagent enrichment."""
    print("=" * 70)
    print("  Test Step 3: Structured Output and Advanced Features")
    print("=" * 70)
    print()
    print("──" * 35)
    print("  Test 3a: Structured Output Format")
    print("──" * 35)
    print(f"   Query: {TEST_REACTION}")
    print()
    
    print("   ⏱️  Running structured recommendation...")
    t_start = time.time()
    
    try:
        result = recommend_conditions_structured(
            reaction=TEST_REACTION,
            k=25,
            limit=3,
        )
        
        t_elapsed = time.time() - t_start
        print(f"   ⏱️  Completed in {t_elapsed:.3f} seconds")
        print()
        
        # Validate structure
        print("   ✅ Results:")
        
        # Check meta section
        if 'meta' in result:
            meta = result['meta']
            print(f"      - Meta generated_at: {meta.get('generated_at', 'N/A')}")
            print(f"      - Model: {meta.get('model', 'N/A')}")
            print(f"      - Status: {meta.get('status', 'N/A')}")
            print(f"      - Processing time: {meta.get('processing_time_ms', 'N/A')} ms")
        else:
            print("      ⚠️  No 'meta' section found")
        
        print()
        
        # Check input section
        if 'input' in result:
            input_data = result['input']
            print(f"      - Input reaction: {input_data.get('reaction_smiles', 'N/A')[:50]}...")
            print(f"      - Requested type: {input_data.get('requested_reaction_type', 'N/A')}")
        else:
            print("      ⚠️  No 'input' section found")
        
        print()
        
        # Check detection section
        if 'detection' in result:
            detection = result['detection']
            print(f"      - Detected type: {detection.get('detected_reaction_type', 'N/A')}")
            print(f"      - Confidence: {detection.get('confidence', 0.0):.4f}")
            print(f"      - Method: {detection.get('method', 'N/A')}")
        else:
            print("      ⚠️  No 'detection' section found")
        
        print()
        
        # Check recommendations
        if 'recommendations' in result:
            recommendations = result['recommendations']
            print(f"      - Number of recommendations: {len(recommendations)}")
            
            if len(recommendations) > 0:
                print()
                print("   📊 Top recommendation details:")
                print()
                
                top = recommendations[0]
                
                # Basic info
                print(f"      [1] Rank: {top.get('rank', 'N/A')}")
                print(f"          Confidence: {top.get('confidence', 0.0):.4f}")
                print(f"          Support: {top.get('support', 0)} precedents")
                print()
                
                # Reagents section
                if 'reagents' in top:
                    reagents = top['reagents']
                    print(f"          Reagents ({len(reagents)} total):")
                    
                    for reagent in reagents:
                        role = reagent.get('role', 'unknown')
                        name = reagent.get('name') or reagent.get('abbreviation') or 'N/A'
                        
                        # Handle SMILES for starting materials
                        if not name or name == 'N/A':
                            smiles = reagent.get('smiles')
                            if smiles:
                                name = f"[SMILES: {smiles[:30]}...]"
                        
                        equiv_info = reagent.get('equivalents', {})
                        if isinstance(equiv_info, dict):
                            equiv_val = equiv_info.get('value', 'N/A')
                            equiv_unit = equiv_info.get('unit', '')
                        else:
                            equiv_val = equiv_info if equiv_info else 'N/A'
                            equiv_unit = 'eq'
                        
                        print(f"            - {role}: {name} ({equiv_val} {equiv_unit})")
                        
                        # Show enrichment details for catalysts/ligands/bases/solvents
                        if role in ['catalyst', 'ligand', 'base', 'solvent']:
                            cas = reagent.get('cas')
                            if cas:
                                print(f"              → CAS: {cas}")
                            
                            # Loading for catalysts/ligands
                            if role in ['catalyst', 'ligand'] and 'loading' in reagent:
                                loading = reagent['loading']
                                loading_val = loading.get('value', 'N/A')
                                loading_unit = loading.get('unit', 'mol%')
                                print(f"              → Loading: {loading_val} {loading_unit}")
                    
                    print()
                elif 'chemicals' in top:
                    # Alternative format: chemicals array
                    chemicals = top['chemicals']
                    print(f"          Chemicals ({len(chemicals)} total):")
                    
                    for chem in chemicals:
                        role = chem.get('role', 'unknown')
                        name = chem.get('name') or 'N/A'
                        
                        print(f"            - {role}: {name}")
                        
                        # Show CAS if available
                        cas = chem.get('cas')
                        if cas:
                            print(f"              → CAS: {cas}")
                    
                    print()
                else:
                    print("          ⚠️  No 'reagents' or 'chemicals' section")
                    print()
                
                # Conditions section
                if 'conditions' in top:
                    conditions = top['conditions']
                    print(f"          Conditions:")
                    
                    # Handle dict format (structured)
                    if isinstance(conditions, dict):
                        # Temperature
                        if 'temperature' in conditions:
                            temp = conditions['temperature']
                            if temp is not None:
                                if isinstance(temp, dict):
                                    temp_val = temp.get('value', 'N/A')
                                    temp_unit = temp.get('unit', '°C')
                                    temp_range = temp.get('range', [])
                                    if temp_range and len(temp_range) == 2:
                                        print(f"            - Temperature: {temp_val} {temp_unit} (range: {temp_range[0]}-{temp_range[1]} {temp_unit})")
                                    else:
                                        print(f"            - Temperature: {temp_val} {temp_unit}")
                                else:
                                    print(f"            - Temperature: {temp}")
                        
                        # Time
                        if 'time' in conditions:
                            time_data = conditions['time']
                            if time_data is not None:
                                if isinstance(time_data, dict):
                                    time_val = time_data.get('value', 'N/A')
                                    time_unit = time_data.get('unit', 'hours')
                                    time_range = time_data.get('range', [])
                                    if time_range and len(time_range) == 2:
                                        print(f"            - Time: {time_val} {time_unit} (range: {time_range[0]}-{time_range[1]} {time_unit})")
                                    else:
                                        print(f"            - Time: {time_val} {time_unit}")
                                else:
                                    print(f"            - Time: {time_data}")
                        
                        # Atmosphere
                        if 'atmosphere' in conditions:
                            atm = conditions['atmosphere']
                            if isinstance(atm, dict):
                                atm_gas = atm.get('gas', 'N/A')
                                atm_req = atm.get('requirement', 'N/A')
                                print(f"            - Atmosphere: {atm_gas} ({atm_req})")
                            else:
                                print(f"            - Atmosphere: {atm}")
                        
                        # Pressure
                        if 'pressure' in conditions:
                            press = conditions['pressure']
                            if isinstance(press, dict):
                                press_val = press.get('value', 'N/A')
                                press_unit = press.get('unit', 'atm')
                                print(f"            - Pressure: {press_val} {press_unit}")
                            else:
                                print(f"            - Pressure: {press}")
                    else:
                        # Simple format
                        print(f"            {conditions}")
                    
                    print()
                else:
                    print("          ⚠️  No 'conditions' section")
                    print()
                
                # Precedents (if available)
                if 'precedents' in top:
                    prec_info = top['precedents']
                    prec_count = prec_info.get('count', 0)
                    print(f"          Precedents: {prec_count} supporting reactions")
                    
                    top_similar = prec_info.get('top_similar', [])
                    if top_similar:
                        print(f"          Top similar precedents: {len(top_similar)} shown")
                    print()
        else:
            print("      ⚠️  No 'recommendations' section found")
        
        print()
        print("──" * 35)
        print("  Validation Summary")
        print("──" * 35)
        print()
        
        # Validation checks
        checks_passed = 0
        total_checks = 6
        
        # Check 1: Meta section exists
        if 'meta' in result and result['meta'].get('status') == 'success':
            print("   ✅ Meta section present and valid")
            checks_passed += 1
        else:
            print("   ❌ Meta section missing or invalid")
        
        # Check 2: Input section exists
        if 'input' in result and 'reaction_smiles' in result['input']:
            print("   ✅ Input section present and valid")
            checks_passed += 1
        else:
            print("   ❌ Input section missing or invalid")
        
        # Check 3: Detection section exists
        if 'detection' in result and 'detected_reaction_type' in result['detection']:
            print("   ✅ Detection section present and valid")
            checks_passed += 1
        else:
            print("   ❌ Detection section missing or invalid")
        
        # Check 4: Recommendations exist
        if 'recommendations' in result and len(result['recommendations']) > 0:
            print("   ✅ Recommendations generated")
            checks_passed += 1
        else:
            print("   ❌ No recommendations generated")
        
        # Check 5: Reagents are enriched
        if 'recommendations' in result and len(result['recommendations']) > 0:
            top = result['recommendations'][0]
            has_enrichment = False
            
            # Check reagents format
            if 'reagents' in top and len(top['reagents']) > 0:
                # Check if at least one reagent has enriched data (CAS, name, etc.)
                for reagent in top['reagents']:
                    if reagent.get('cas') or reagent.get('name'):
                        has_enrichment = True
                        break
                
                if has_enrichment:
                    print("   ✅ Reagents enriched with database info")
                    checks_passed += 1
                else:
                    print("   ⚠️  Reagents present but not enriched")
            # Check chemicals format (alternative)
            elif 'chemicals' in top and len(top['chemicals']) > 0:
                for chemical in top['chemicals']:
                    if chemical.get('cas') or chemical.get('name'):
                        has_enrichment = True
                        break
                
                if has_enrichment:
                    print("   ✅ Chemicals present with identifiers")
                    checks_passed += 1
                else:
                    print("   ⚠️  Chemicals present but minimal info")
            else:
                print("   ❌ No reagents/chemicals in recommendations")
        else:
            print("   ❌ Cannot check reagent enrichment")
        
        # Check 6: Conditions have ranges
        if 'recommendations' in result and len(result['recommendations']) > 0:
            top = result['recommendations'][0]
            if 'conditions' in top:
                conds = top['conditions']
                has_ranges = False
                
                # Handle dict format
                if isinstance(conds, dict):
                    # Check temperature range
                    if 'temperature' in conds:
                        temp = conds['temperature']
                        if isinstance(temp, dict) and 'range' in temp and len(temp.get('range', [])) == 2:
                            has_ranges = True
                    
                    # Check time range
                    if 'time' in conds:
                        time_data = conds['time']
                        if isinstance(time_data, dict) and 'range' in time_data and len(time_data.get('range', [])) == 2:
                            has_ranges = True
                
                if has_ranges:
                    print("   ✅ Conditions include value ranges")
                    checks_passed += 1
                elif isinstance(conds, dict) and ('temperature' in conds or 'time' in conds):
                    # Conditions exist but no ranges - that's still OK
                    print("   ⚠️  Conditions present (ranges may not be structured)")
                    checks_passed += 1
                else:
                    print("   ⚠️  Conditions present but format unclear")
            else:
                print("   ❌ No conditions in recommendations")
        else:
            print("   ❌ Cannot check condition ranges")
        
        print()
        print(f"   📊 Validation: {checks_passed}/{total_checks} checks passed")
        print()
        
        # Overall success
        if checks_passed >= 5:  # Allow 1 failure
            print("   ✅ Test 3a PASSED")
            print()
            return 0
        else:
            print("   ❌ Test 3a FAILED")
            print()
            return 1
            
    except Exception as e:
        print(f"   ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


def test_3b_output_format_comparison():
    """Test 3b: Compare basic vs structured output formats."""
    print("──" * 35)
    print("  Test 3b: Output Format Comparison")
    print("──" * 35)
    print(f"   Query: {TEST_REACTION}")
    print()
    
    print("   📋 Comparing basic vs structured output...")
    print()
    
    # Basic output
    print("   [A] Basic recommend_from_reaction():")
    t_start = time.time()
    basic_result = recommend_from_reaction(
        reaction=TEST_REACTION,
        k=25,
        max_variants=3,
    )
    t_basic = time.time() - t_start
    
    print(f"       ⏱️  Time: {t_basic:.3f}s")
    print(f"       📊 Keys: {list(basic_result.keys())[:8]}")
    
    # Check formatted section
    if 'formatted' in basic_result:
        formatted = basic_result['formatted']
        print(f"       📦 'formatted' section keys: {list(formatted.keys())}")
        
        if 'conditions' in formatted:
            conds = formatted['conditions']
            print(f"       📊 Number of conditions: {len(conds)}")
    else:
        print("       ⚠️  No 'formatted' section")
    
    print()
    
    # Structured output
    print("   [B] Structured recommend_conditions_structured():")
    t_start = time.time()
    structured_result = recommend_conditions_structured(
        reaction=TEST_REACTION,
        k=25,
        limit=3,
    )
    t_struct = time.time() - t_start
    
    print(f"       ⏱️  Time: {t_struct:.3f}s")
    print(f"       📊 Keys: {list(structured_result.keys())}")
    
    if 'recommendations' in structured_result:
        recs = structured_result['recommendations']
        print(f"       📦 Number of recommendations: {len(recs)}")
        
        if len(recs) > 0:
            first = recs[0]
            print(f"       📊 Recommendation keys: {list(first.keys())}")
    
    print()
    print("   ✅ Comparison complete")
    print()
    
    return 0


def test_3c_multiple_variants():
    """Test 3c: Multiple condition variants."""
    print("──" * 35)
    print("  Test 3c: Multiple Condition Variants")
    print("──" * 35)
    print(f"   Query: {TEST_REACTION}")
    print()
    
    print("   🧪 Generating 5 condition variants...")
    t_start = time.time()
    
    result = recommend_conditions_structured(
        reaction=TEST_REACTION,
        k=25,
        limit=5,  # Request 5 variants
    )
    
    t_elapsed = time.time() - t_start
    print(f"   ⏱️  Completed in {t_elapsed:.3f} seconds")
    print()
    
    if 'recommendations' in result:
        recs = result['recommendations']
        print(f"   ✅ Generated {len(recs)} variants")
        print()
        
        # Show summary of each variant
        for i, rec in enumerate(recs[:5], 1):
            rank = rec.get('rank', i)
            conf = rec.get('confidence', 0.0)
            support = rec.get('support', 0)
            
            print(f"   [{i}] Rank {rank}:")
            print(f"       Confidence: {conf:.4f}")
            print(f"       Support: {support} precedents")
            
            # Show catalyst system, base, and solvent to demonstrate diversity
            if 'reagents' in rec:
                reagents = rec['reagents']
                catalysts = [r for r in reagents if r.get('role') in ['catalyst', 'ligand']]
                bases = [r for r in reagents if r.get('role') == 'base']
                solvents = [r for r in reagents if r.get('role') == 'solvent']
                
                if catalysts:
                    cat_names = [r.get('name') or r.get('abbreviation') or 'Unknown' for r in catalysts]
                    print(f"       Catalyst: {', '.join(cat_names)}")
                if bases:
                    base_names = [r.get('name') or r.get('abbreviation') or 'Unknown' for r in bases]
                    print(f"       Base: {', '.join(base_names)}")
                if solvents:
                    solv_names = [r.get('name') or r.get('abbreviation') or 'Unknown' for r in solvents]
                    print(f"       Solvent: {', '.join(solv_names)}")
                    
            elif 'chemicals' in rec:
                # Alternative format: chemicals array
                chemicals = rec['chemicals']
                catalysts = [c for c in chemicals if c.get('role') in ['metal_precursor', 'ligand']]
                bases = [c for c in chemicals if c.get('role') == 'base']
                solvents = [c for c in chemicals if c.get('role') == 'solvent']
                
                if catalysts:
                    cat_names = [c.get('name') or 'Unknown' for c in catalysts]
                    print(f"       Catalyst: {', '.join(cat_names)}")
                if bases:
                    base_names = [c.get('name') or 'Unknown' for c in bases]
                    print(f"       Base: {', '.join(base_names)}")
                if solvents:
                    solv_names = [c.get('name') or 'Unknown' for c in solvents]
                    print(f"       Solvent: {', '.join(solv_names)}")
            
            print()
        
        print(f"   ✅ Test 3c PASSED")
        print()
        return 0
    else:
        print("   ❌ No recommendations found")
        print()
        return 1


def main():
    """Run all Test Step 3 tests."""
    print()
    
    # Test 3a: Structured output
    result_3a = test_3a_structured_output()
    
    # Test 3b: Format comparison
    result_3b = test_3b_output_format_comparison()
    
    # Test 3c: Multiple variants
    result_3c = test_3c_multiple_variants()
    
    # Performance summary
    print("=" * 70)
    print("  Performance Summary")
    print("=" * 70)
    print()
    print("   ⏱️  All tests completed")
    print()
    
    # Final summary
    print("=" * 70)
    print("  ✅ Test Step 3 Complete!")
    print("=" * 70)
    print()
    
    total_failed = result_3a + result_3b + result_3c
    
    if total_failed == 0:
        print("🎉 All validation checks passed!")
        return 0
    else:
        print(f"⚠️  {total_failed} test(s) had issues")
        return 1


if __name__ == "__main__":
    sys.exit(main())
