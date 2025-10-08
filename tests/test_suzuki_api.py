"""
Test Suzuki Reaction Recommendations via FastAPI Endpoints
===========================================================

This script tests the three main recommendation approaches:
1. Rule-Based: /match endpoint (SMARTS pattern matching)
2. ML-Based: /api/v1/recommend/conditions (DRFP k-NN precedent search)
3. Fusion: /api/v1/recommend/fusion (multi-source evidence with adaptive weights)

Usage:
    python test_suzuki_api.py

Requirements:
    - FastAPI server running on http://localhost:8000
    - Start server with: uvicorn app.main:app --reload --port 8000
"""

import sys
import io

# Set UTF-8 encoding for Windows console
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

import requests
import json
from typing import Dict, Any, List
from datetime import datetime


# API Configuration
BASE_URL = "http://localhost:8000"
SUZUKI_DB = "data/conditionDB/Suzuki_db.json"


def print_header(title: str):
    """Print formatted section header"""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def print_subheader(title: str):
    """Print formatted subsection header"""
    print("\n" + "-" * 80)
    print(f"  {title}")
    print("-" * 80)


def check_server_health() -> bool:
    """Check if FastAPI server is running"""
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=2)
        if response.status_code == 200:
            print("✓ Server is running")
            return True
        else:
            print(f"✗ Server returned status {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print(f"✗ Cannot connect to server at {BASE_URL}")
        print(f"  Error: {e}")
        print("\n  Please start the server with:")
        print("  uvicorn app.main:app --reload --port 8000")
        return False


def _extract_name(value: Any) -> str:
    """Extract name from various data structures (string, dict, list)"""
    if value is None:
        return "N/A"
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        return value.get("name", "N/A")
    if isinstance(value, list):
        if len(value) > 0:
            item = value[0]
            if isinstance(item, str):
                return item
            if isinstance(item, dict):
                return item.get("name", "N/A")
        return "N/A"
    return str(value)


def test_rule_based(reaction: str, description: str) -> Dict[str, Any]:
    """Test rule-based recommendation using /match endpoint"""
    print_subheader(f"Rule-Based: {description}")
    print(f"Reaction: {reaction}\n")
    
    try:
        payload = {
            "reaction": reaction,
            "db": SUZUKI_DB,
            "include_trace": True
        }
        
        response = requests.post(
            f"{BASE_URL}/match",
            json=payload,
            timeout=10
        )
        
        if response.status_code == 200:
            result = response.json()
            
            # Extract key information
            match_type = result.get("match_type", "N/A")
            entry_name = result.get("entry_name", "N/A")
            conditions = result.get("conditions", {})
            
            print(f"✓ Match Type: {match_type}")
            print(f"  Entry: {entry_name}")
            
            if conditions:
                print("\n  Recommended Conditions:")
                
                # Handle both dict and list structures
                if isinstance(conditions, list):
                    # If conditions is a list, use the first item
                    if len(conditions) > 0:
                        conditions = conditions[0]
                    else:
                        print("  (Empty conditions list)")
                        return result
                
                # Catalyst
                catalyst = conditions.get("catalyst") if isinstance(conditions, dict) else None
                if catalyst:
                    if isinstance(catalyst, dict):
                        core = catalyst.get("core", "N/A")
                        ligand = catalyst.get("ligand", "")
                        if ligand and ligand != "N/A" and ligand != "":
                            print(f"    Catalyst: {core} + {ligand}")
                        else:
                            print(f"    Catalyst: {core}")
                    elif isinstance(catalyst, list) and len(catalyst) > 0:
                        # Handle list of catalysts - show first one
                        cat_item = catalyst[0]
                        if isinstance(cat_item, dict):
                            core = cat_item.get("core", cat_item.get("name", "N/A"))
                            ligand = cat_item.get("ligand", "")
                            if ligand:
                                print(f"    Catalyst: {core} + {ligand}")
                            else:
                                print(f"    Catalyst: {core}")
                        else:
                            print(f"    Catalyst: {cat_item}")
                    else:
                        print(f"    Catalyst: {catalyst}")
                
                # Base
                base = conditions.get("base") if isinstance(conditions, dict) else None
                if base:
                    base_name = _extract_name(base)
                    if base_name != "N/A":
                        print(f"    Base: {base_name}")
                
                # Solvent
                solvent = conditions.get("solvent") if isinstance(conditions, dict) else None
                if solvent:
                    solv_name = _extract_name(solvent)
                    if solv_name != "N/A":
                        print(f"    Solvent: {solv_name}")
                
                # Temperature
                temp = conditions.get("temperature") if isinstance(conditions, dict) else None
                if temp:
                    if isinstance(temp, dict):
                        temp_val = temp.get("value", "N/A")
                    else:
                        temp_val = temp
                    if temp_val != "N/A":
                        print(f"    Temperature: {temp_val}°C")
                
                # Time
                time_val = conditions.get("time") if isinstance(conditions, dict) else None
                if time_val:
                    if isinstance(time_val, dict):
                        t_val = time_val.get("value", "N/A")
                    else:
                        t_val = time_val
                    if t_val != "N/A":
                        print(f"    Time: {t_val} hours")
            else:
                print("  No conditions returned")
            
            return result
        else:
            print(f"✗ Error: HTTP {response.status_code}")
            print(f"  {response.text}")
            return {}
            
    except Exception as e:
        print(f"✗ Exception: {e}")
        import traceback
        print("\nTraceback:")
        traceback.print_exc()
        return {}


def test_ml_based(reaction: str, description: str, k: int = 50) -> Dict[str, Any]:
    """Test ML-based recommendation using /api/v1/recommend/conditions endpoint"""
    print_subheader(f"ML-Based (k-NN): {description}")
    print(f"Reaction: {reaction}")
    print(f"k={k} precedents\n")
    
    try:
        payload = {
            "reaction": reaction,
            "reaction_type": None,  # Auto-detect
            "k": k,
            "limit": 5,
            "relax": {},
            "constraints": {}
        }
        
        response = requests.post(
            f"{BASE_URL}/api/v1/recommend/conditions",
            json=payload,
            timeout=30
        )
        
        if response.status_code == 200:
            result = response.json()
            
            # Extract metadata
            detection = result.get("detection", {})
            detected_type = detection.get("type", "Unknown")
            confidence = detection.get("confidence", 0.0)
            
            print(f"✓ Detected Type: {detected_type} (confidence: {confidence:.2%})")
            
            # Show recommendations
            recommendations = result.get("recommendations", [])
            print(f"  Found {len(recommendations)} recommendations\n")
            
            for i, rec in enumerate(recommendations[:3], 1):  # Show top 3
                reagents = rec.get("reagents", [])
                conditions_info = rec.get("conditions", {})
                precedent_count = rec.get("precedent_count", 0)
                rec_confidence = rec.get("confidence", 0.0)
                
                print(f"  {i}. Confidence: {rec_confidence:.2%} | Support: {precedent_count} precedents")
                
                # Extract catalyst, base, solvent
                catalyst = None
                base = None
                solvent = None
                
                for r in reagents:
                    role = r.get("role", "").upper()
                    name = r.get("name", "N/A")
                    
                    if role == "CATALYST":
                        catalyst = name
                    elif role == "BASE":
                        base = name
                    elif role == "SOLVENT":
                        solvent = name
                
                if catalyst:
                    print(f"     Catalyst: {catalyst}")
                if base:
                    print(f"     Base: {base}")
                if solvent:
                    print(f"     Solvent: {solvent}")
                
                # Temperature and time
                temp = conditions_info.get("temperature")
                if temp:
                    temp_val = temp.get("value") if isinstance(temp, dict) else temp
                    print(f"     Temperature: {temp_val}°C")
                
                time_val = conditions_info.get("time")
                if time_val:
                    t_val = time_val.get("value") if isinstance(time_val, dict) else time_val
                    print(f"     Time: {t_val} hours")
                
                print()
            
            return result
        else:
            print(f"✗ Error: HTTP {response.status_code}")
            if response.status_code == 500:
                print(f"  Server Error: {response.text}")
                print(f"\n  ⚠️  This is a server-side error. Check the uvicorn terminal for:")
                print(f"     - Python traceback")
                print(f"     - FileNotFoundError, KeyError, or AttributeError")
                print(f"     - Missing data files or import errors")
            else:
                print(f"  {response.text}")
            return {}
            
    except Exception as e:
        print(f"✗ Exception: {e}")
        return {}


def test_fusion(reaction: str, description: str, k: int = 50) -> Dict[str, Any]:
    """Test fusion recommendation using /api/v1/recommend/fusion endpoint"""
    print_subheader(f"Fusion (Multi-Source): {description}")
    print(f"Reaction: {reaction}")
    print(f"k={k} precedents\n")
    
    try:
        payload = {
            "reaction": reaction,
            "k": k,
            "max_variants": 5,
            "relax": {},
            "constraints": {}
        }
        
        response = requests.post(
            f"{BASE_URL}/api/v1/recommend/fusion",
            json=payload,
            timeout=30
        )
        
        if response.status_code == 200:
            result = response.json()
            
            # Check for fusion metadata
            fusion_meta = result.get("fusion_meta", {})
            if fusion_meta and "error" not in fusion_meta:
                print("✓ Fusion System Active")
                
                # Show adaptive weights
                weights = fusion_meta.get("adaptive_weights", {})
                if weights:
                    alpha = weights.get("alpha", 0.0)
                    beta = weights.get("beta", 0.0)
                    gamma = weights.get("gamma", 0.0)
                    delta = weights.get("delta", 0.0)
                    
                    print(f"  Adaptive Weights:")
                    print(f"    α (precedents) = {alpha:.3f}")
                    print(f"    β (analytics)  = {beta:.3f}")
                    print(f"    γ (rules)      = {gamma:.3f}")
                    print(f"    δ (ML)         = {delta:.3f}")
                
                # Evidence summary
                evidence = fusion_meta.get("evidence_summary", {})
                if evidence:
                    prec_count = evidence.get("precedent_count", 0)
                    diversity = evidence.get("diversity_score", 0.0)
                    
                    print(f"\n  Evidence Quality:")
                    print(f"    Precedent count: {prec_count}")
                    print(f"    Diversity score: {diversity:.3f}", end="")
                    if diversity < 0.3:
                        print(" (LOW - potential batch effect)")
                    else:
                        print(" (OK)")
                
                # Reasoning
                reasoning = fusion_meta.get("reasoning", [])
                if reasoning:
                    print(f"\n  Reasoning:")
                    for reason in reasoning[:3]:  # Show top 3 reasons
                        print(f"    - {reason}")
            else:
                print("⚠ Fusion metadata not available")
                if "error" in fusion_meta:
                    print(f"  Error: {fusion_meta.get('error')}")
            
            # Show recommendations
            recs = result.get("formatted", {}).get("recommended_conditions", [])
            if not recs:
                recs = result.get("recommended_conditions", [])
            
            print(f"\n  Recommendations: {len(recs)}")
            
            for i, rec in enumerate(recs[:3], 1):  # Show top 3
                summary = rec.get("summary", {})
                component_scores = rec.get("component_scores", {})
                
                core = summary.get("core", "N/A")
                base = summary.get("base", {})
                solvent = summary.get("solvent", {})
                confidence = summary.get("confidence", "UNKNOWN")
                
                base_name = base.get("name", "N/A") if isinstance(base, dict) else str(base)
                solv_name = solvent.get("name", "N/A") if isinstance(solvent, dict) else str(solvent)
                
                print(f"\n  {i}. Core: {core} | Confidence: {confidence}")
                print(f"     Base: {base_name}")
                print(f"     Solvent: {solv_name}")
                
                # Component scores
                if component_scores:
                    ps = component_scores.get("PS", 0.0)
                    as_score = component_scores.get("AS", 0.0)
                    rs = component_scores.get("RS", 0.0)
                    ms = component_scores.get("MS", 0.0)
                    
                    print(f"     Scores: PS={ps:.3f}, AS={as_score:.3f}, RS={rs:.3f}, MS={ms:.3f}")
            
            return result
        else:
            print(f"✗ Error: HTTP {response.status_code}")
            print(f"  {response.text}")
            return {}
            
    except Exception as e:
        print(f"✗ Exception: {e}")
        return {}


def test_reaction_comprehensive(
    reaction: str,
    description: str,
    run_rule: bool = True,
    run_ml: bool = True,
    run_fusion: bool = True
):
    """Test a single reaction with all three approaches"""
    print_header(f"Testing: {description}")
    print(f"SMILES: {reaction}")
    
    results = {}
    
    if run_rule:
        print_header("APPROACH 1: RULE-BASED (SMARTS Pattern Matching)")
        results["rule"] = test_rule_based(reaction, description)
    
    if run_ml:
        print_header("APPROACH 2: ML-BASED (DRFP k-NN Precedent Search)")
        results["ml"] = test_ml_based(reaction, description, k=50)
    
    if run_fusion:
        print_header("APPROACH 3: FUSION (Multi-Source Evidence)")
        results["fusion"] = test_fusion(reaction, description, k=50)
    
    return results


def main():
    """Run comprehensive Suzuki reaction tests"""
    print_header("SUZUKI REACTION API ENDPOINT TESTING")
    print(f"Test Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Target API: {BASE_URL}")
    
    # Check server health
    print_header("SERVER HEALTH CHECK")
    if not check_server_health():
        return
    
    # Define test reactions
    test_reactions = [
        {
            "smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
            "description": "Simple Suzuki: Ph-Br + Ph-B(OH)2",
            "note": "Classic benchmark reaction"
        },
        {
            "smiles": "Clc1ccc(C#N)cc1.c1ccc(B(O)O)cc1>>N#Cc1ccc(-c2ccccc2)cc1",
            "description": "Electron-poor ArCl with CN group",
            "note": "Challenging substrate (Cl + electron-withdrawing)"
        },
        {
            "smiles": "Brc1ccc(OC)cc1.c1ccc(B(O)O)cc1>>COc1ccc(-c2ccccc2)cc1",
            "description": "Electron-rich ArBr with OMe",
            "note": "Donor group effects"
        },
        {
            "smiles": "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1",
            "description": "Heteroaryl pyridine coupling",
            "note": "Potential catalyst coordination issues"
        },
        {
            "smiles": "Brc1ccccc1OCC.c1ccc(B(O)O)cc1C>>CCOc1ccccc1-c1ccccc1C",
            "description": "Ortho-substituted (sterically hindered)",
            "note": "Steric challenge"
        }
    ]
    
    # Test each reaction
    all_results = []
    for i, rxn_data in enumerate(test_reactions, 1):
        print(f"\n\n{'#' * 80}")
        print(f"# REACTION {i}/{len(test_reactions)}")
        print(f"{'#' * 80}")
        print(f"Note: {rxn_data['note']}")
        
        results = test_reaction_comprehensive(
            reaction=rxn_data["smiles"],
            description=rxn_data["description"]
        )
        
        all_results.append({
            "reaction": rxn_data,
            "results": results
        })
        
        # Add separator between reactions (removed interactive pause)
        if i < len(test_reactions):
            print("\n")
    
    # Summary
    print_header("TEST SUMMARY")
    print(f"\nTested {len(test_reactions)} Suzuki reactions")
    print(f"Test Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    success_count = {"rule": 0, "ml": 0, "fusion": 0}
    for result_set in all_results:
        results = result_set["results"]
        if results.get("rule"):
            success_count["rule"] += 1
        if results.get("ml"):
            success_count["ml"] += 1
        if results.get("fusion"):
            success_count["fusion"] += 1
    
    print("Success Rates:")
    print(f"  Rule-Based: {success_count['rule']}/{len(test_reactions)}")
    print(f"  ML-Based:   {success_count['ml']}/{len(test_reactions)}")
    print(f"  Fusion:     {success_count['fusion']}/{len(test_reactions)}")
    
    print("\nKey Observations:")
    print("  - Rule-Based: Fast, deterministic, limited to known patterns")
    print("  - ML-Based: Data-driven, handles novel substrates, needs precedents")
    print("  - Fusion: Combines all evidence, adaptive weights, most robust")
    
    print(f"\n{'=' * 80}\n")


if __name__ == "__main__":
    main()
