"""
Quick CLI Test for Pure LLM Workflow
=====================================

Tests the pure LLM workflow components with diverse CAS numbers.
This is a streamlined version that tests the core LLM functions directly.

Run: python quick_llm_test.py
"""

import json
import sys
import time
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from llmtools.clients import LLMClient
from llmtools.reagent_classifier import classify_role, assign_fields, verify_entry

# Test samples covering different reagent types
TEST_SAMPLES = [
    {"cas": "121-44-8", "name": "Triethylamine", "expected_role": "base", "category": "Organic base"},
    {"cas": "1633-05-2", "name": "Potassium tert-butoxide", "expected_role": "base", "category": "Strong base"},
    {"cas": "14221-01-3", "name": "Pd(PPh3)4", "expected_role": "catalyst", "category": "Pd catalyst"},
    {"cas": "67-56-1", "name": "Methanol", "expected_role": "solvent", "category": "Protic solvent"},
    {"cas": "109-99-9", "name": "Tetrahydrofuran", "expected_role": "solvent", "category": "Aprotic solvent"},
    {"cas": "603-35-0", "name": "Triphenylphosphine", "expected_role": "ligand", "category": "Phosphine ligand"},
    {"cas": "16940-66-2", "name": "Sodium borohydride", "expected_role": "reducing_agent", "category": "Reducing agent"},
    {"cas": "7722-84-1", "name": "Hydrogen peroxide", "expected_role": "oxidizing_agent", "category": "Oxidizer"},
]


def test_role_classification(llm_client: LLMClient):
    """Test Step 2: Role classification."""
    print("\n" + "="*80)
    print("STEP 2: ROLE CLASSIFICATION TEST")
    print("="*80)
    
    results = []
    total_time = 0
    
    for i, sample in enumerate(TEST_SAMPLES, 1):
        print(f"\n[{i}/{len(TEST_SAMPLES)}] Testing: {sample['name']} (CAS: {sample['cas']})")
        print(f"  Category: {sample['category']}")
        print(f"  Expected role: {sample['expected_role']}")
        
        start = time.time()
        try:
            # Build identity dict for the classifier
            identity = {
                "name": sample["name"],
                "cas": sample["cas"],
                "smiles": "",  # Empty for this test
                "molecular_formula": "",
                "synonyms": []
            }
            
            result = classify_role(
                identity=identity,
                llm_client=llm_client
            )
            elapsed = time.time() - start
            total_time += elapsed
            
            status = result.get("status")
            role = result.get("role")
            confidence = result.get("confidence", 0)
            reasoning = result.get("reasoning", "")
            
            match = role == sample["expected_role"]
            
            print(f"  [OK] Detected role: {role} ({confidence:.0%} confidence)")
            print(f"  {'[MATCH]' if match else '[MISMATCH]'}")
            print(f"  Reasoning: {reasoning[:100]}...")
            print(f"  Time: {elapsed:.2f}s")
            
            results.append({
                "sample": sample,
                "detected_role": role,
                "confidence": confidence,
                "match": match,
                "elapsed": elapsed
            })
            
            time.sleep(0.5)  # Rate limiting
            
        except Exception as e:
            elapsed = time.time() - start
            print(f"  ✗ ERROR: {e}")
            results.append({
                "sample": sample,
                "error": str(e),
                "match": False,
                "elapsed": elapsed
            })
    
    # Summary
    print(f"\n{'='*80}")
    print("ROLE CLASSIFICATION SUMMARY")
    print(f"{'='*80}")
    
    matches = sum(1 for r in results if r.get("match", False))
    total = len(results)
    avg_conf = sum(r.get("confidence", 0) for r in results if "confidence" in r) / total if total > 0 else 0
    avg_time = total_time / total if total > 0 else 0
    
    print(f"Accuracy: {matches}/{total} ({matches/total*100:.1f}%)")
    print(f"Avg confidence: {avg_conf:.0%}")
    print(f"Avg time: {avg_time:.2f}s")
    print(f"Total time: {total_time:.1f}s")
    
    return results


def test_field_assignment(llm_client: LLMClient, role_results):
    """Test Step 3: Field assignment."""
    print("\n" + "="*80)
    print("STEP 3: FIELD ASSIGNMENT TEST")
    print("="*80)
    
    results = []
    registry_dir = Path("../data/reagent_db")
    
    for i, role_result in enumerate(role_results, 1):
        if "error" in role_result:
            continue
            
        sample = role_result["sample"]
        role = role_result["detected_role"]
        
        print(f"\n[{i}/{len(role_results)}] Assigning fields for: {sample['name']}")
        print(f"  Role: {role}")
        
        try:
            # Build identity dict
            identity = {
                "name": sample["name"],
                "cas": sample["cas"],
                "smiles": "",
                "molecular_formula": "",
                "synonyms": []
            }
            
            result = assign_fields(
                identity=identity,
                role=role,
                registry_dir=registry_dir,
                llm_client=llm_client
            )
            
            status = result.get("status")
            family = result.get("family")
            fields = result.get("fields", {})
            
            print(f"  ✓ Family: {family}")
            print(f"  ✓ Fields assigned: {len(fields)} keys")
            if fields:
                print(f"    {', '.join(list(fields.keys())[:5])}...")
            
            results.append({
                "sample": sample,
                "role": role,
                "family": family,
                "fields_count": len(fields),
                "status": status
            })
            
            time.sleep(0.5)
            
        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            results.append({
                "sample": sample,
                "error": str(e)
            })
    
    successful = len([r for r in results if 'family' in r])
    print(f"\n{'='*80}")
    print(f"Field assignment completed: {successful}/{len(results)} successful")
    
    return results


def test_verification(llm_client: LLMClient, field_results):
    """Test Step 4: Verification."""
    print("\n" + "="*80)
    print("STEP 4: VERIFICATION TEST")
    print("="*80)
    
    results = []
    
    for i, field_result in enumerate(field_results, 1):
        if "error" in field_result:
            continue
            
        sample = field_result["sample"]
        
        # Build mock entry for verification
        entry = {
            "cas": sample["cas"],
            "name": sample["name"],
            "smiles": "",
            "roles": [field_result.get("role", "")],
            "family": field_result.get("family", ""),
        }
        
        print(f"\n[{i}/{len(field_results)}] Verifying: {sample['name']}")
        
        try:
            result = verify_entry(
                entry=entry,
                llm_client=llm_client
            )
            
            approved = result.get("approved", False)
            issues = result.get("issues", [])
            suggestions = result.get("suggestions", [])
            
            print(f"  {'✓ APPROVED' if approved else '✗ ISSUES FOUND'}")
            if issues:
                print(f"  Issues: {', '.join(issues[:3])}")
            if suggestions:
                print(f"  Suggestions: {', '.join(suggestions[:2])}")
            
            results.append({
                "sample": sample,
                "approved": approved,
                "issues_count": len(issues),
                "suggestions_count": len(suggestions)
            })
            
            time.sleep(0.5)
            
        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            results.append({
                "sample": sample,
                "error": str(e)
            })
    
    approved_count = sum(1 for r in results if r.get("approved", False))
    print(f"\n{'='*80}")
    print(f"Verification: {approved_count}/{len(results)} entries approved")
    
    return results


def main():
    """Run comprehensive workflow test."""
    print("="*80)
    print("PURE LLM WORKFLOW - CLI TEST SUITE")
    print("="*80)
    print(f"Testing {len(TEST_SAMPLES)} samples")
    print(f"Model: deepseek-v3.2-exp")
    print()
    
    # Create LLM client
    llm_client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
    
    # Test each step
    print("\n🧪 Starting 3-step workflow test...")
    
    # Step 2: Role classification
    role_results = test_role_classification(llm_client)
    
    # Step 3: Field assignment
    field_results = test_field_assignment(llm_client, role_results)
    
    # Step 4: Verification
    verify_results = test_verification(llm_client, field_results)
    
    # Final summary
    print("\n\n" + "="*80)
    print("FINAL SUMMARY")
    print("="*80)
    
    total_role = len(role_results)
    total_field = len(field_results) if field_results else 1
    total_verify = len(verify_results) if verify_results else 1
    
    role_accuracy = sum(1 for r in role_results if r.get("match", False)) / total_role * 100 if total_role > 0 else 0
    field_success = len([r for r in field_results if "family" in r]) / total_field * 100 if field_results else 0
    verify_approved = sum(1 for r in verify_results if r.get("approved", False)) / total_verify * 100 if verify_results else 0
    
    print(f"✓ Role Classification Accuracy: {role_accuracy:.1f}%")
    print(f"✓ Field Assignment Success: {field_success:.1f}%")
    print(f"✓ Verification Approval Rate: {verify_approved:.1f}%")
    
    print("\n✅ Test suite complete!")


if __name__ == "__main__":
    main()
