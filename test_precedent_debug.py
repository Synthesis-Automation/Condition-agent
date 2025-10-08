"""
Diagnostic script to identify performance bottleneck in precedent search
"""

import sys
from pathlib import Path
import time

project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

def time_step(description):
    """Decorator to time each step"""
    def decorator(func):
        def wrapper(*args, **kwargs):
            print(f"\n[STEP] {description}...", flush=True)
            start = time.time()
            result = func(*args, **kwargs)
            elapsed = time.time() - start
            print(f"[TIME] {description}: {elapsed:.3f}s", flush=True)
            return result
        return wrapper
    return decorator

@time_step("Import ChemTools")
def import_chemtools():
    from chemtools import ChemTools
    return ChemTools

@time_step("Initialize ChemTools")
def init_chemtools(ChemTools):
    return ChemTools()

@time_step("Detect reaction family")
def detect_family(chem, rxn):
    return chem.router.detect_family(rxn)

@time_step("Load dataset directly")
def load_dataset():
    from chemtools.precedent.loader import _load_selective
    return _load_selective(families=["Suzuki"])

@time_step("Check dataset size")
def check_dataset(rows):
    print(f"  Dataset has {len(rows)} reactions")
    if rows:
        print(f"  First reaction: {rows[0].get('reaction_id', 'N/A')}")
    return rows

@time_step("Run precedent search (minimal)")
def search_precedents(chem, family, rxn_smiles):
    features = {"bin": "", "LG": "", "nuc_class": ""}
    search_config = {
        "use_drfp": False,  # Disable DRFP first to see baseline
        "selective_loading": True,
        "reaction_smiles": rxn_smiles,
    }
    return chem.precedent.knn(
        family=family,
        features=features,
        k=5,
        relax=search_config
    )

@time_step("Run precedent search (with DRFP)")
def search_with_drfp(chem, family, rxn_smiles):
    features = {"bin": "", "LG": "", "nuc_class": ""}
    search_config = {
        "use_drfp": True,
        "drfp_weight": 0.7,
        "precompute_drfp": True,
        "selective_loading": True,
        "reaction_smiles": rxn_smiles,
        "debug_timing": True,
    }
    return chem.precedent.knn(
        family=family,
        features=features,
        k=5,
        relax=search_config
    )

def main():
    print("=" * 80)
    print("PRECEDENT SEARCH PERFORMANCE DIAGNOSTIC")
    print("=" * 80)
    
    rxn_smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    print(f"\nQuery: {rxn_smiles}")
    
    # Step 1: Import
    ChemTools = import_chemtools()
    
    # Step 2: Initialize
    chem = init_chemtools(ChemTools)
    
    # Step 3: Detect family
    detection = detect_family(chem, rxn_smiles)
    detected = detection.get("family") or detection.get("mapped_family", "Suzuki_CC")
    family = "Suzuki"  # Map to dataset name
    print(f"  Detected: {detected} -> Using: {family}")
    
    # Step 4: Load dataset directly
    rows = load_dataset()
    rows = check_dataset(rows)
    
    # Step 5: Search without DRFP (baseline)
    print("\n" + "-" * 80)
    print("TEST 1: WITHOUT DRFP (baseline)")
    print("-" * 80)
    result1 = search_precedents(chem, family, rxn_smiles)
    print(f"  Found: {result1.get('support', 0)} precedents")
    print(f"  Timing: {result1.get('timing', {})}")
    
    # Step 6: Search with DRFP
    print("\n" + "-" * 80)
    print("TEST 2: WITH DRFP")
    print("-" * 80)
    result2 = search_with_drfp(chem, family, rxn_smiles)
    print(f"  Found: {result2.get('support', 0)} precedents")
    print(f"  Timing: {result2.get('timing', {})}")
    print(f"  DRFP strategy: {result2.get('drfp_load_strategy', {})}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("\nIf loading took >5s: Dataset loading is slow")
    print("If search without DRFP took >1s: Feature similarity is slow")
    print("If search with DRFP took >5s: DRFP computation is slow")
    print("\nCheck timing above to identify bottleneck.")
    
if __name__ == "__main__":
    main()
