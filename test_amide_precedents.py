#!/usr/bin/env python
"""
Amide Formation Precedent Search Quality Test
Tests precedent search for all amide formation reactions and generates a report.
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


def extract_amide_reactions():
    """Extract all amide formation reactions from sample data."""
    amide_reactions = []
    for rxn in SAMPLE_REACTIONS:
        if "(Amide:" in rxn and ">>" in rxn:
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
                amide_reactions.append({
                    "smiles": reaction_smiles,
                    "description": description or rxn,
                })
    return amide_reactions


def test_amide_precedents():
    """Test precedent search for amide formation reactions."""
    
    print("\n" + "="*80)
    print("  AMIDE FORMATION PRECEDENT SEARCH QUALITY TEST")
    print("="*80)
    
    amide_reactions = extract_amide_reactions()
    print(f"\nFound {len(amide_reactions)} amide formation reactions in sample data.")
    print("Testing precedent search with DRFP-based similarity...\n")
    
    # Configure DRFP search
    relax = {
        "use_drfp": True,
        "drfp_weight": 0.7,  # High weight on reaction center similarity
        "precompute_drfp": False,  # Binary NPZ files should exist
        "selective_loading": True,
        "debug_timing": False,
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
    
    for i, rxn_data in enumerate(amide_reactions, 1):
        rxn_smiles = rxn_data["smiles"]
        description = rxn_data["description"]
        
        print(f"\n[{i}/{len(amide_reactions)}] {description[:60]}")
        print(f"    Reaction: {rxn_smiles[:70]}...")
        
        try:
            # Update relax with current reaction
            search_config = {**relax, "reaction_smiles": rxn_smiles}
            
            start_time = time.time()
            result = knn(family="Amide_formation", features=features, k=5, relax=search_config)
            elapsed = time.time() - start_time
            total_time += elapsed
            
            precedents = result.get('precedents', [])
            support = result.get('support', 0)
            
            print(f"    Found: {support} total, showing top 5 | Time: {elapsed:.3f}s")
            
            if precedents:
                print(f"\n    Top 5 Precedents:")
                for j, prec in enumerate(precedents[:5], 1):
                    # For amide formation, core is usually coupling reagent
                    core = prec.get('condition_core') or prec.get('core', 'N/A')
                    yield_val = prec.get('yield', 'N/A')
                    
                    # Try to extract coupling reagent info
                    reagents = prec.get('reagents', [])
                    coupling_reagent = "N/A"
                    if reagents:
                        if isinstance(reagents, list):
                            # Look for common coupling reagents
                            for r in reagents:
                                name = r.get('name', '') if isinstance(r, dict) else str(r)
                                if any(x in name.upper() for x in ['EDC', 'HATU', 'T3P', 'DCC', 'COMU', 'PyBOP']):
                                    coupling_reagent = name
                                    break
                        elif isinstance(reagents, str):
                            coupling_reagent = reagents[:30]
                    
                    print(f"      {j}. Core: {core:25s} | Reagent: {coupling_reagent:20s} | Yield: {yield_val}%")
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
    if results:
        print(f"\nTotal time: {total_time:.3f}s")
        print(f"Average per reaction: {total_time/len(results):.3f}s")
    
    # Save detailed markdown report
    report_file = parent_dir / "AMIDE_PRECEDENT_REPORT.md"
    with open(report_file, "w", encoding="utf-8") as f:
        f.write("# Amide Formation Precedent Search Quality Report\n\n")
        f.write(f"**Tested:** {len(results)} reactions\n")
        f.write(f"**Successful:** {successful}\n")
        f.write(f"**No Precedents:** {no_precedents}\n")
        f.write(f"**Errors:** {failed}\n")
        if results:
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
                    
                    # Coupling reagents
                    reagents = prec.get('reagents')
                    if reagents:
                        reag_list = []
                        if isinstance(reagents, list):
                            for r in reagents:
                                if isinstance(r, dict):
                                    name = r.get('name', 'N/A')
                                    role = r.get('role', '')
                                    if role:
                                        reag_list.append(f"{name} ({role})")
                                    else:
                                        reag_list.append(name)
                                else:
                                    reag_list.append(str(r))
                        elif isinstance(reagents, str):
                            reag_list.append(reagents)
                        
                        if reag_list:
                            f.write(f"- **Reagents:** {', '.join(reag_list[:5])}\n")  # Limit to 5
                    
                    # Solvents
                    solvents = prec.get('solvents')
                    if solvents:
                        solv_list = []
                        if isinstance(solvents, list):
                            for s in solvents:
                                if isinstance(s, dict):
                                    solv_list.append(s.get('name', str(s)))
                                else:
                                    solv_list.append(str(s))
                        elif isinstance(solvents, str):
                            solv_list.append(solvents)
                        
                        if solv_list:
                            f.write(f"- **Solvents:** {', '.join(solv_list)}\n")
                    
                    # Conditions
                    temp = prec.get('T_C')
                    time_h = prec.get('time_h')
                    if temp is not None or time_h is not None:
                        cond_parts = []
                        if temp is not None:
                            cond_parts.append(f"{temp}C")
                        if time_h is not None:
                            cond_parts.append(f"{time_h}h")
                        f.write(f"- **Conditions:** {', '.join(cond_parts)}\n")
                    
                    # Yield
                    yield_val = prec.get('yield')
                    if yield_val is not None:
                        f.write(f"- **Yield:** {yield_val}%\n")
                    
                    # Reference
                    reference = prec.get('reference')
                    if reference:
                        f.write(f"- **Reference:** {reference}\n")
                    
                    f.write("\n")
            else:
                f.write("**No precedents found**\n\n")
            
            f.write("---\n\n")
    
    print(f"\n[REPORT] Saved to: {report_file}")
    
    # Also save a summary analysis
    analysis_file = parent_dir / "AMIDE_PRECEDENT_ANALYSIS.md"
    with open(analysis_file, "w", encoding="utf-8") as f:
        f.write("# Amide Formation Precedent Search - Analysis\n\n")
        f.write("## Summary\n\n")
        f.write(f"**Performance:**\n\n")
        f.write(f"- Tested: {len(results)} amide formation reactions\n")
        f.write(f"- Successful: {successful}\n")
        f.write(f"- No precedents: {no_precedents}\n")
        f.write(f"- Errors: {failed}\n")
        if results:
            f.write(f"- Avg time: {total_time/len(results):.3f}s per reaction\n\n")
        
        f.write("## Key Findings\n\n")
        f.write("### Reaction Center Focus\n\n")
        f.write("The DRFP-based precedent search focuses on:\n")
        f.write("- Carboxylic acid + amine transformation\n")
        f.write("- Similar acid structures (aromatic, aliphatic, heteroaromatic)\n")
        f.write("- Similar amine types (primary, secondary, aromatic, aliphatic)\n")
        f.write("- Amide bond formation chemistry\n\n")
        
        f.write("### Common Coupling Reagents Found\n\n")
        
        # Extract coupling reagents from results
        all_reagents = []
        for r in results:
            if "precedents" in r:
                for prec in r["precedents"]:
                    reagents = prec.get('reagents', [])
                    if isinstance(reagents, list):
                        for reag in reagents:
                            if isinstance(reag, dict):
                                name = reag.get('name', '')
                                if name:
                                    all_reagents.append(name)
        
        # Count common coupling reagents
        from collections import Counter
        reagent_counter = Counter(all_reagents)
        common_reagents = reagent_counter.most_common(10)
        
        if common_reagents:
            f.write("Top 10 coupling reagents/additives in precedents:\n\n")
            for reag, count in common_reagents:
                f.write(f"- {reag}: {count} occurrences\n")
        
        f.write("\n---\n\n")
        f.write(f"**Test Date:** 2025-10-08\n")
        f.write(f"**Dataset:** Amide_formation.jsonl with precomputed DRFP fingerprints\n")
    
    print(f"[ANALYSIS] Saved to: {analysis_file}")
    
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
    summary = test_amide_precedents()
    sys.exit(0 if summary["failed"] == 0 else 1)
