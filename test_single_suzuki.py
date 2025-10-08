"""
Quick Test: Single Suzuki Reaction Precedent Search
====================================================

Debug script to test precedent search on just one Suzuki reaction.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from chemtools import ChemTools

def test_single_suzuki():
    """Test precedent search for a single Suzuki reaction."""
    
    # Simple Suzuki coupling: PhBr + PhB(OH)2 -> biphenyl
    rxn_smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    
    print("=" * 80)
    print("SINGLE SUZUKI PRECEDENT SEARCH TEST")
    print("=" * 80)
    print(f"\nQuery Reaction: {rxn_smiles}")
    print("Description: Simple PhBr + PhB(OH)2 -> biphenyl\n")
    
    # Initialize ChemTools
    print("Initializing ChemTools...")
    chem = ChemTools()
    
    # Detect family
    print("\n1. Detecting reaction family...")
    detection = chem.router.detect_family(rxn_smiles)
    detected_family = detection.get("family") or detection.get("mapped_family", "Unknown")
    print(f"   Detected: {detected_family}")
    print(f"   Full detection result: {detection}")
    
    # Map to dataset family name
    family_map = {
        "Suzuki_CC": "Suzuki",
        "C_N_Coupling_Pd": "C_N_Coupling_Pd",
        "C_N_Coupling_Cu": "C_N_Coupling_Cu",
    }
    family = family_map.get(detected_family, detected_family)
    print(f"   Using family for search: {family}")
    
    # Configure search with DRFP enabled
    print("\n2. Configuring precedent search...")
    search_config = {
        "use_drfp": True,
        "drfp_weight": 0.7,  # High weight on reaction center similarity
        "precompute_drfp": True,
        "debug_timing": True,
        "reaction_smiles": rxn_smiles,
        "selective_loading": True,  # Only load Suzuki dataset
    }
    print(f"   DRFP enabled: {search_config['use_drfp']}")
    print(f"   DRFP weight: {search_config['drfp_weight']}")
    
    # Minimal features for DRFP-based search
    features = {"bin": "", "LG": "", "nuc_class": ""}
    
    # Search for precedents
    print("\n3. Searching for precedents...")
    try:
        result = chem.precedent.knn(
            family=family,
            features=features,
            k=5,
            relax=search_config
        )
        
        support = result.get("support", 0)
        precedents = result.get("precedents", [])
        timing = result.get("timing", {})
        drfp_strategy = result.get("drfp_load_strategy", {})
        
        print(f"\n   Results:")
        print(f"   - Total precedents found: {support}")
        print(f"   - Returning top: {len(precedents)}")
        
        if timing:
            print(f"\n   Timing:")
            for key, value in timing.items():
                print(f"   - {key}: {value:.3f}s")
        
        if drfp_strategy:
            print(f"\n   DRFP Loading Strategy:")
            print(f"   - Binary NPZ: {drfp_strategy.get('binary', 0)}")
            print(f"   - JSONL embedded: {drfp_strategy.get('jsonl', 0)}")
            print(f"   - On-demand compute: {drfp_strategy.get('computed', 0)}")
        
        # Display precedents
        if precedents:
            print("\n" + "=" * 80)
            print("TOP 5 PRECEDENTS")
            print("=" * 80)
            
            for i, prec in enumerate(precedents[:5], 1):
                print(f"\n{i}. Reaction: {prec.get('reaction_smiles', 'N/A')}")
                
                core = prec.get('condition_core') or prec.get('core')
                if core:
                    print(f"   Core: {core}")
                
                catalyst = prec.get('catalyst')
                if catalyst:
                    print(f"   Catalyst: {catalyst}")
                
                cat_sys = prec.get('catalytic_system')
                if cat_sys:
                    print(f"   Catalytic System: {cat_sys}")
                
                solvents = prec.get('solvents')
                if solvents:
                    solv_str = solvents if isinstance(solvents, str) else ', '.join(solvents) if isinstance(solvents, list) else str(solvents)
                    print(f"   Solvents: {solv_str}")
                
                reagents = prec.get('reagents')
                if reagents:
                    reag_str = reagents if isinstance(reagents, str) else ', '.join(reagents) if isinstance(reagents, list) else str(reagents)
                    print(f"   Reagents: {reag_str}")
                
                temp = prec.get('T_C')
                time_h = prec.get('time_h')
                if temp is not None or time_h is not None:
                    cond_parts = []
                    if temp is not None:
                        cond_parts.append(f"{temp}°C")
                    if time_h is not None:
                        cond_parts.append(f"{time_h}h")
                    print(f"   Conditions: {', '.join(cond_parts)}")
                
                yield_val = prec.get('yield')
                if yield_val is not None:
                    print(f"   Yield: {yield_val}%")
                
                reference = prec.get('reference')
                if reference:
                    print(f"   Reference: {reference}")
        else:
            print("\n[!] No precedents found!")
            print("\nPossible reasons:")
            print("- Dataset not loaded correctly")
            print("- Family name mismatch")
            print("- No matching precedents in database")
        
    except Exception as e:
        print(f"\n[ERROR] {str(e)}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("TEST COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    test_single_suzuki()
