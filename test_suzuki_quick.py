#!/usr/bin/env python
"""
Quick Suzuki Precedent Search Report
Tests precedent search for all Suzuki reactions and generates a report.
"""
import sys
from pathlib import Path
import time

# Add parent directory to path
parent_dir = Path(__file__).parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chemtools.precedent import knn
import importlib.util

# Load sample reactions
spec = importlib.util.spec_from_file_location("sample_reactions", parent_dir / "tests" / "sample_reactions.py")
sample_reactions_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sample_reactions_module)
SAMPLE_REACTIONS = sample_reactions_module.SAMPLE_REACTIONS


def extract_suzuki_reactions():
    """Extract all Suzuki coupling reactions from sample data."""
    suzuki_reactions = []
    for rxn in SAMPLE_REACTIONS:
        if "Suzuki" in rxn and ">>" in rxn:
            parts = rxn.split(">>")
            if len(parts) == 2:
                description = ""
                if "(" in parts[1]:
                    prod_parts = parts[1].rsplit("(", 1)
                    products = prod_parts[0].strip()
                    if len(prod_parts) == 2:
                        description = "(" + prod_parts[1]
                else:
                    products = parts[1].strip()
                
                reaction_smiles = parts[0] + ">>" + products
                suzuki_reactions.append({
                    "smiles": reaction_smiles,
                    "description": description or rxn,
                })
    return suzuki_reactions


def test_suzuki_precedents():
    """Test precedent search for Suzuki reactions."""
    
    print("\n" + "="*80)
    print("  SUZUKI PRECEDENT SEARCH QUALITY TEST")
    print("="*80)
    
    suzuki_reactions = extract_suzuki_reactions()
    print(f"\nFound {len(suzuki_reactions)} Suzuki reactions in sample data.")
    print("Testing precedent search with DRFP-based similarity...\n")
    
    # Configure DRFP search
    relax = {
        "use_drfp": True,
        "drfp_weight": 0.7,  # High weight on reaction center similarity
        "precompute_drfp": False,  # Binary NPZ files should exist
        "selective_loading": True,
        "debug_timing": False,  # Disable verbose timing for cleaner output
    }
    
    # Minimal features (DRFP-driven search)
    features = {
        "LG": "",
        "nuc_class": "",
        "bin": ""
    }
    
    results = []
    total_time = 0
    
    print("─" * 80)
    print("Running precedent search for each reaction...")
    print("─" * 80)
    
    for i, rxn_data in enumerate(suzuki_reactions, 1):  # Test ALL Suzuki reactions
        rxn_smiles = rxn_data["smiles"]
        description = rxn_data["description"]
        
        print(f"\n[{i}/{len(suzuki_reactions)}] {description[:60]}")
        print(f"    Reaction: {rxn_smiles[:70]}...")
        
        try:
            # Update relax with current reaction
            search_config = {**relax, "reaction_smiles": rxn_smiles}
            
            start_time = time.time()
            result = knn(family="Suzuki", features=features, k=5, relax=search_config)
            elapsed = time.time() - start_time
            total_time += elapsed
            
            precedents = result.get('precedents', [])
            support = result.get('support', 0)
            
            print(f"    Found: {support} total, showing top 5 | Time: {elapsed:.3f}s")
            
            if precedents:
                print(f"\n    Top 5 Precedents:")
                for j, prec in enumerate(precedents[:5], 1):
                    core = prec.get('condition_core') or prec.get('core', 'N/A')
                    yield_val = prec.get('yield', 'N/A')
                    print(f"      {j}. Core: {core:30s} | Yield: {yield_val}%")
            else:
                print(f"    [!] No precedents found")
            
            results.append({
                "index": i,
                "description": description,
                "smiles": rxn_smiles,
                "support": support,
                "precedents": precedents[:5],
                "time": elapsed
            })
            
        except Exception as e:
            print(f"    [ERROR] {str(e)}")
            results.append({
                "index": i,
                "description": description,
                "smiles": rxn_smiles,
                "error": str(e)
            })
    
    # Summary
    print("\n" + "=" * 80)
    print("  SUMMARY")
    print("=" * 80)
    
    successful = sum(1 for r in results if "error" not in r and r.get("support", 0) > 0)
    failed = sum(1 for r in results if "error" in r)
    no_precedents = sum(1 for r in results if "error" not in r and r.get("support", 0) == 0)
    
    print(f"\nTested: {len(results)} reactions")
    print(f"  [OK] Successful with precedents: {successful}")
    print(f"  [!]  No precedents found: {no_precedents}")
    print(f"  [ERROR] Errors: {failed}")
    print(f"\nTotal time: {total_time:.3f}s")
    print(f"Average per reaction: {total_time/len(results):.3f}s")
    
    # Save detailed markdown report
    report_file = parent_dir / "SUZUKI_PRECEDENT_REPORT.md"
    with open(report_file, "w", encoding="utf-8") as f:
        f.write("# Suzuki Precedent Search Quality Report\n\n")
        f.write(f"**Tested:** {len(results)} reactions\n")
        f.write(f"**Successful:** {successful}\n")
        f.write(f"**No Precedents:** {no_precedents}\n")
        f.write(f"**Errors:** {failed}\n")
        f.write(f"**Avg Time:** {total_time/len(results):.3f}s per reaction\n\n")
        f.write("---\n\n")
        
        for result in results:
            idx = result["index"]
            desc = result["description"]
            smiles = result["smiles"]
            
            f.write(f"## {idx}. {desc}\n\n")
            f.write(f"**Reaction:** `{smiles}`\n\n")
            
            if "error" in result:
                f.write(f"**Error:** {result['error']}\n\n")
                continue
            
            support = result.get("support", 0)
            prec_time = result.get("time", 0)
            f.write(f"**Found:** {support} total precedents | **Time:** {prec_time:.3f}s\n\n")
            
            precedents = result.get("precedents", [])
            if precedents:
                f.write("**Top 5 Precedents:**\n\n")
                for j, prec in enumerate(precedents, 1):
                    f.write(f"### Precedent {j}\n\n")
                    f.write(f"- **Reaction:** `{prec.get('reaction_smiles', 'N/A')}`\n")
                    
                    core = prec.get('condition_core') or prec.get('core')
                    if core:
                        f.write(f"- **Core:** {core}\n")
                    
                    catalyst = prec.get('catalyst')
                    if catalyst:
                        f.write(f"- **Catalyst:** {catalyst}\n")
                    
                    cat_sys = prec.get('catalytic_system')
                    if cat_sys:
                        f.write(f"- **Catalytic System:** {cat_sys}\n")
                    
                    solvents = prec.get('solvents')
                    if solvents:
                        solv_str = solvents if isinstance(solvents, str) else ', '.join([s.get('name', str(s)) if isinstance(s, dict) else str(s) for s in solvents]) if isinstance(solvents, list) else str(solvents)
                        f.write(f"- **Solvents:** {solv_str}\n")
                    
                    reagents = prec.get('reagents')
                    if reagents:
                        reag_str = reagents if isinstance(reagents, str) else ', '.join([r.get('name', str(r)) if isinstance(r, dict) else str(r) for r in reagents]) if isinstance(reagents, list) else str(reagents)
                        f.write(f"- **Reagents:** {reag_str}\n")
                    
                    temp = prec.get('T_C')
                    time_h = prec.get('time_h')
                    if temp is not None or time_h is not None:
                        cond_parts = []
                        if temp is not None:
                            cond_parts.append(f"{temp}C")
                        if time_h is not None:
                            cond_parts.append(f"{time_h}h")
                        f.write(f"- **Conditions:** {', '.join(cond_parts)}\n")
                    
                    yield_val = prec.get('yield')
                    if yield_val is not None:
                        f.write(f"- **Yield:** {yield_val}%\n")
                    
                    f.write("\n")
            else:
                f.write("**No precedents found**\n\n")
            
            f.write("---\n\n")
    
    print(f"\n[REPORT] Saved to: {report_file}")
    print("\n" + "=" * 80)
    print("  TEST COMPLETE")
    print("=" * 80 + "\n")
    
    return {
        "total": len(results),
        "successful": successful,
        "failed": failed,
        "avg_time": total_time/len(results) if results else 0
    }


if __name__ == "__main__":
    summary = test_suzuki_precedents()
    sys.exit(0 if summary["failed"] == 0 else 1)
