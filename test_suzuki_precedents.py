"""
Test Precedent Search Quality for Suzuki Reactions
===================================================

This script evaluates how well the precedent search focuses on reaction centers
for all Suzuki reactions in the sample data. It generates a detailed report
showing the top 5 precedents for each query reaction.

The precedent search should prioritize:
1. Reaction center similarity (DRFP-based)
2. Similar functional groups near the reaction site
3. Ignore remote functional groups
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Import from local module
import importlib.util
spec = importlib.util.spec_from_file_location("sample_reactions", project_root / "tests" / "sample_reactions.py")
sample_reactions_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sample_reactions_module)
SAMPLE_REACTIONS = sample_reactions_module.SAMPLE_REACTIONS

from chemtools import ChemTools
import json
from datetime import datetime


def extract_suzuki_reactions():
    """Extract all Suzuki coupling reactions from sample data."""
    suzuki_reactions = []
    for rxn in SAMPLE_REACTIONS:
        if "Suzuki" in rxn and ">>" in rxn:
            # Parse reaction SMILES and description
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
                    "original": rxn
                })
    return suzuki_reactions


def format_precedent(prec, index):
    """Format a single precedent for display."""
    lines = []
    lines.append(f"\n  {index}. Reaction: {prec.get('reaction_smiles', 'N/A')}")
    
    # Core conditions
    core = prec.get('condition_core') or prec.get('core') or 'N/A'
    lines.append(f"     Core: {core}")
    
    # Catalyst details
    catalyst = prec.get('catalyst')
    if catalyst:
        lines.append(f"     Catalyst: {catalyst}")
    
    catalytic_system = prec.get('catalytic_system')
    if catalytic_system:
        lines.append(f"     Catalytic System: {catalytic_system}")
    
    # Solvents
    solvents = prec.get('solvents')
    if solvents:
        solvent_str = solvents if isinstance(solvents, str) else ', '.join(solvents) if isinstance(solvents, list) else str(solvents)
        lines.append(f"     Solvents: {solvent_str}")
    
    # Reagents
    reagents = prec.get('reagents')
    if reagents:
        reagent_str = reagents if isinstance(reagents, str) else ', '.join(reagents) if isinstance(reagents, list) else str(reagents)
        lines.append(f"     Reagents: {reagent_str}")
    
    # Temperature and time
    temp = prec.get('T_C')
    time_h = prec.get('time_h')
    if temp is not None or time_h is not None:
        cond_parts = []
        if temp is not None:
            cond_parts.append(f"{temp}°C")
        if time_h is not None:
            cond_parts.append(f"{time_h}h")
        lines.append(f"     Conditions: {', '.join(cond_parts)}")
    
    # Yield
    yield_val = prec.get('yield')
    if yield_val is not None:
        lines.append(f"     Yield: {yield_val}%")
    
    # Reference
    reference = prec.get('reference')
    if reference:
        lines.append(f"     Reference: {reference}")
    
    return '\n'.join(lines)


def test_all_suzuki_reactions():
    """Test precedent search for all Suzuki reactions and generate report."""
    print("=" * 80)
    print("SUZUKI PRECEDENT SEARCH QUALITY TEST")
    print("=" * 80)
    print("\nThis test evaluates how well the precedent search focuses on")
    print("reaction centers rather than remote functional groups.\n")
    
    # Extract Suzuki reactions
    suzuki_reactions = extract_suzuki_reactions()
    print(f"Found {len(suzuki_reactions)} Suzuki reactions in sample data.\n")
    
    # Initialize ChemTools
    print("Initializing ChemTools...")
    chem = ChemTools()
    
    # Test configuration with DRFP enabled (reaction center focused)
    search_config = {
        "use_drfp": True,
        "drfp_weight": 0.7,  # High weight on reaction fingerprint similarity
        "precompute_drfp": True,  # Use precomputed fingerprints if available
        "debug_timing": True,  # Show timing information
    }
    
    # Results storage
    results = []
    
    print("\nRunning precedent search for each Suzuki reaction...")
    print("-" * 80)
    
    for i, rxn_data in enumerate(suzuki_reactions, 1):
        rxn_smiles = rxn_data["smiles"]
        description = rxn_data["description"]
        
        print(f"\n[{i}/{len(suzuki_reactions)}] Testing: {description}")
        print(f"    Reaction: {rxn_smiles}")
        
        try:
            # Detect reaction family
            detection = chem.router.detect_family(rxn_smiles)
            detected_family = detection.get("family") or detection.get("mapped_family", "Suzuki_CC")
            
            # Map to dataset family name (loader expects "Suzuki" not "Suzuki_CC")
            family_map = {
                "Suzuki_CC": "Suzuki",
                "C_N_Coupling_Pd": "C_N_Coupling_Pd",
                "C_N_Coupling_Cu": "C_N_Coupling_Cu",
            }
            family = family_map.get(detected_family, detected_family)
            
            print(f"    Detected Family: {detected_family} -> Using: {family}")
            
            # Get precedents directly using precedent.knn
            # First we need features - for DRFP-only search, features can be minimal
            features = {"bin": "", "LG": "", "nuc_class": ""}
            
            # Build search config with reaction_smiles for DRFP
            full_search_config = {
                **search_config,
                "reaction_smiles": rxn_smiles
            }
            
            precedent_result = chem.precedent.knn(
                family=family,
                features=features,
                k=5,  # Top 5 precedents
                relax=full_search_config
            )
            
            precedents = precedent_result.get("precedents", [])
            support = precedent_result.get("support", 0)
            timing = precedent_result.get("timing", {})
            
            print(f"    Found {support} total precedents, showing top 5")
            if timing:
                total_time = timing.get('total', 0)
                print(f"    Search time: {total_time:.3f}s")
            
            # Store results
            result_entry = {
                "query_index": i,
                "query_smiles": rxn_smiles,
                "query_description": description,
                "detected_family": family,
                "support": support,
                "top_precedents": precedents[:5],
                "timing": timing
            }
            results.append(result_entry)
            
            # Display top precedents
            if precedents:
                print(f"\n    Top 5 Precedents:")
                for j, prec in enumerate(precedents[:5], 1):
                    formatted = format_precedent(prec, j)
                    # Indent the precedent details
                    indented = '\n'.join('    ' + line for line in formatted.split('\n'))
                    print(indented)
            else:
                print("    [!] No precedents found!")
            
        except Exception as e:
            print(f"    [ERROR] {str(e)}")
            import traceback
            traceback.print_exc()
            result_entry = {
                "query_index": i,
                "query_smiles": rxn_smiles,
                "query_description": description,
                "error": str(e)
            }
            results.append(result_entry)
    
    # Generate summary report
    print("\n" + "=" * 80)
    print("SUMMARY REPORT")
    print("=" * 80)
    
    successful = sum(1 for r in results if "error" not in r and r.get("support", 0) > 0)
    failed = sum(1 for r in results if "error" in r)
    no_precedents = sum(1 for r in results if "error" not in r and r.get("support", 0) == 0)
    
    print(f"\nTotal Reactions Tested: {len(results)}")
    print(f"  [OK] Successful: {successful}")
    print(f"  [!] No Precedents Found: {no_precedents}")
    print(f"  [ERROR] Errors: {failed}")
    
    # Average timing
    timings = [r.get("timing", {}).get("total", 0) for r in results if "timing" in r]
    if timings:
        avg_time = sum(timings) / len(timings)
        print(f"\nAverage Search Time: {avg_time:.3f}s")
    
    # Save detailed report to JSON
    report_file = project_root / "suzuki_precedent_report.json"
    with open(report_file, "w") as f:
        json.dump({
            "metadata": {
                "test_date": datetime.now().isoformat(),
                "total_reactions": len(results),
                "successful": successful,
                "no_precedents": no_precedents,
                "errors": failed,
                "search_config": search_config
            },
            "results": results
        }, f, indent=2)
    
    print(f"\n[REPORT] Detailed report saved to: {report_file}")
    
    # Save human-readable markdown report
    md_report_file = project_root / "SUZUKI_PRECEDENT_REPORT.md"
    with open(md_report_file, "w", encoding="utf-8") as f:
        f.write("# Suzuki Precedent Search Quality Report\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**Configuration:**\n")
        f.write(f"- DRFP Enabled: `{search_config['use_drfp']}`\n")
        f.write(f"- DRFP Weight: `{search_config['drfp_weight']}`\n")
        f.write(f"- Precompute DRFP: `{search_config['precompute_drfp']}`\n\n")
        
        f.write("## Summary\n\n")
        f.write(f"- **Total Reactions:** {len(results)}\n")
        f.write(f"- **Successful:** {successful}\n")
        f.write(f"- **No Precedents:** {no_precedents}\n")
        f.write(f"- **Errors:** {failed}\n")
        if timings:
            f.write(f"- **Avg Search Time:** {avg_time:.3f}s\n")
        f.write("\n---\n\n")
        
        f.write("## Detailed Results\n\n")
        for result in results:
            idx = result["query_index"]
            desc = result["query_description"]
            smiles = result["query_smiles"]
            
            f.write(f"### {idx}. {desc}\n\n")
            f.write(f"**Query Reaction:** `{smiles}`\n\n")
            
            if "error" in result:
                f.write(f"**❌ Error:** {result['error']}\n\n")
                continue
            
            family = result.get("detected_family", "Unknown")
            support = result.get("support", 0)
            f.write(f"**Detected Family:** {family}\n\n")
            f.write(f"**Total Precedents Found:** {support}\n\n")
            
            precedents = result.get("top_precedents", [])
            if precedents:
                f.write("**Top 5 Precedents:**\n\n")
                for j, prec in enumerate(precedents, 1):
                    f.write(f"#### Precedent {j}\n\n")
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
                        solv_str = solvents if isinstance(solvents, str) else ', '.join(solvents) if isinstance(solvents, list) else str(solvents)
                        f.write(f"- **Solvents:** {solv_str}\n")
                    
                    reagents = prec.get('reagents')
                    if reagents:
                        reag_str = reagents if isinstance(reagents, str) else ', '.join(reagents) if isinstance(reagents, list) else str(reagents)
                        f.write(f"- **Reagents:** {reag_str}\n")
                    
                    temp = prec.get('T_C')
                    time_h = prec.get('time_h')
                    if temp is not None or time_h is not None:
                        cond_parts = []
                        if temp is not None:
                            cond_parts.append(f"{temp}°C")
                        if time_h is not None:
                            cond_parts.append(f"{time_h}h")
                        f.write(f"- **Conditions:** {', '.join(cond_parts)}\n")
                    
                    yield_val = prec.get('yield')
                    if yield_val is not None:
                        f.write(f"- **Yield:** {yield_val}%\n")
                    
                    reference = prec.get('reference')
                    if reference:
                        f.write(f"- **Reference:** {reference}\n")
                    
                    f.write("\n")
            else:
                f.write("**⚠️ No precedents found**\n\n")
            
            f.write("---\n\n")
    
    print(f"[REPORT] Markdown report saved to: {md_report_file}")
    
    print("\n" + "=" * 80)
    print("TEST COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    test_all_suzuki_reactions()
