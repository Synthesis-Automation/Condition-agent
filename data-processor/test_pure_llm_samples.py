"""
Test Pure LLM Workflow with Diverse Sample CAS Numbers
=======================================================

Tests the pure LLM workflow with various reagent types to validate:
1. Role classification accuracy
2. Family assignment correctness
3. Field completeness
4. Verification quality
5. Overall workflow robustness

Run: python test_pure_llm_samples.py
"""

import json
import sys
import time
from pathlib import Path
from typing import Dict, Any, List

# Add parent directory to path to import from llmtools
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import the classifier functions directly
from llmtools.clients import LLMClient
from llmtools.reagent_classifier import classify_role, assign_fields, verify_entry

# Sample CAS numbers covering different reagent types
TEST_SAMPLES = [
    # Bases
    {
        "cas": "121-44-8",
        "name": "Triethylamine",
        "expected_role": "base",
        "expected_family": "tertiary_amines_aliphatic",
        "category": "Common organic base"
    },
    {
        "cas": "1633-05-2",
        "name": "Potassium tert-butoxide",
        "expected_role": "base",
        "expected_family": "alkoxides",
        "category": "Strong base"
    },
    {
        "cas": "7664-38-2",
        "name": "Phosphoric acid",
        "expected_role": "acid",
        "expected_family": "inorganic_acids",
        "category": "Inorganic acid"
    },
    
    # Catalysts
    {
        "cas": "14221-01-3",
        "name": "Tetrakis(triphenylphosphine)palladium(0)",
        "expected_role": "catalyst",
        "expected_family": "palladium_catalysts",
        "category": "Pd(0) catalyst for cross-coupling"
    },
    {
        "cas": "31274-51-8",
        "name": "[1,1'-Bis(diphenylphosphino)ferrocene]dichloropalladium(II)",
        "expected_role": "catalyst",
        "expected_family": "palladium_catalysts",
        "category": "Pd(II) catalyst with bidentate ligand"
    },
    {
        "cas": "14808-79-8",
        "name": "Grubbs Catalyst 1st Generation",
        "expected_role": "catalyst",
        "expected_family": "ruthenium_catalysts",
        "category": "Olefin metathesis catalyst"
    },
    
    # Solvents
    {
        "cas": "67-56-1",
        "name": "Methanol",
        "expected_role": "solvent",
        "expected_family": "alcohols",
        "category": "Protic polar solvent"
    },
    {
        "cas": "109-99-9",
        "name": "Tetrahydrofuran",
        "expected_role": "solvent",
        "expected_family": "ethers",
        "category": "Aprotic polar solvent"
    },
    {
        "cas": "127-19-5",
        "name": "N,N-Dimethylacetamide",
        "expected_role": "solvent",
        "expected_family": "amides",
        "category": "Polar aprotic solvent"
    },
    
    # Ligands
    {
        "cas": "603-35-0",
        "name": "Triphenylphosphine",
        "expected_role": "ligand",
        "expected_family": "phosphine_ligands",
        "category": "Monodentate phosphine"
    },
    {
        "cas": "154-23-4",
        "name": "BINAP",
        "expected_role": "ligand",
        "expected_family": "phosphine_ligands",
        "category": "Chiral bidentate phosphine"
    },
    
    # Reducing agents
    {
        "cas": "16940-66-2",
        "name": "Sodium borohydride",
        "expected_role": "reducing_agent",
        "expected_family": "metal_hydrides",
        "category": "Mild reducing agent"
    },
    {
        "cas": "7446-70-0",
        "name": "Aluminum chloride",
        "expected_role": "lewis_acid",
        "expected_family": "metal_halides",
        "category": "Lewis acid catalyst"
    },
    
    # Oxidizing agents
    {
        "cas": "7722-84-1",
        "name": "Hydrogen peroxide",
        "expected_role": "oxidizing_agent",
        "expected_family": "peroxides",
        "category": "Common oxidizer"
    },
    
    # Coupling reagents
    {
        "cas": "2524-03-0",
        "name": "Dicyclohexylcarbodiimide (DCC)",
        "expected_role": "coupling_reagent",
        "expected_family": "carbodiimides",
        "category": "Peptide coupling"
    },
    
    # Additives/Bases
    {
        "cas": "7757-93-9",
        "name": "Calcium phosphate",
        "expected_role": "additive",
        "expected_family": "inorganic_salts",
        "category": "Inorganic additive"
    },
]


def run_test_sample(sample: Dict[str, Any], llm_client: LLMClient, registry_dir: Path) -> Dict[str, Any]:
    """Run pure LLM workflow on a single test sample."""
    print(f"\n{'='*80}")
    print(f"Testing: {sample['name']} (CAS: {sample['cas']})")
    print(f"Category: {sample['category']}")
    print(f"Expected role: {sample['expected_role']}, family: {sample['expected_family']}")
    print(f"{'='*80}")
    
    start_time = time.time()
    
    try:
        result = generate_taxonomy_entry_llm(
            cas=sample["cas"],
            registry_dir=registry_dir,
            llm_client=llm_client,
            name_override=None
        )
        
        elapsed = time.time() - start_time
        
        # Extract key results
        status = result.get("status", "unknown")
        workflow = result.get("workflow", {})
        entry = result.get("entry", {})
        
        # Analyze workflow steps
        step1 = workflow.get("step1_identity", {})
        step2 = workflow.get("step2_role", {})
        step3 = workflow.get("step3_fields", {})
        step4 = workflow.get("step4_verification", {})
        
        # Extract detected values
        detected_name = step1.get("name", "N/A")
        detected_role = step2.get("role", "N/A")
        role_confidence = step2.get("confidence", 0.0)
        detected_family = step3.get("family", "N/A")
        verification_approved = step4.get("approved", False)
        verification_issues = step4.get("issues", [])
        
        # Check accuracy
        role_match = detected_role == sample["expected_role"]
        family_match = detected_family == sample["expected_family"]
        
        print(f"\n📊 RESULTS:")
        print(f"  Status: {status}")
        print(f"  Time: {elapsed:.2f}s")
        print(f"  Name: {detected_name}")
        print(f"  Role: {detected_role} (expected: {sample['expected_role']}) {'✓' if role_match else '✗'}")
        print(f"  Confidence: {role_confidence:.0%}")
        print(f"  Family: {detected_family} (expected: {sample['expected_family']}) {'✓' if family_match else '✗'}")
        print(f"  Verification: {'✓ Approved' if verification_approved else '✗ Issues found'}")
        
        if verification_issues:
            print(f"  Issues: {', '.join(verification_issues)}")
        
        # Workflow step status
        print(f"\n🔄 WORKFLOW STEPS:")
        for step_name, step_data in workflow.items():
            step_status = step_data.get("status", "unknown")
            print(f"  {step_name}: {step_status}")
        
        return {
            "sample": sample,
            "result": result,
            "elapsed": elapsed,
            "role_match": role_match,
            "family_match": family_match,
            "role_confidence": role_confidence,
            "verification_approved": verification_approved,
            "status": status
        }
        
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n❌ ERROR: {str(e)}")
        return {
            "sample": sample,
            "error": str(e),
            "elapsed": elapsed,
            "role_match": False,
            "family_match": False,
            "status": "error"
        }


def print_summary(results: List[Dict[str, Any]]) -> None:
    """Print summary statistics."""
    print(f"\n\n{'='*80}")
    print("SUMMARY STATISTICS")
    print(f"{'='*80}")
    
    total = len(results)
    successful = sum(1 for r in results if r["status"] in ["ready_to_save", "needs_review"])
    errors = sum(1 for r in results if r["status"] == "error")
    
    role_matches = sum(1 for r in results if r.get("role_match", False))
    family_matches = sum(1 for r in results if r.get("family_match", False))
    
    avg_confidence = sum(r.get("role_confidence", 0) for r in results if "role_confidence" in r) / total if total > 0 else 0
    avg_time = sum(r["elapsed"] for r in results) / total if total > 0 else 0
    
    approved = sum(1 for r in results if r.get("verification_approved", False))
    
    print(f"\n📈 Overall:")
    print(f"  Total samples: {total}")
    print(f"  Successful: {successful} ({successful/total*100:.1f}%)")
    print(f"  Errors: {errors} ({errors/total*100:.1f}%)")
    print(f"  Average time: {avg_time:.2f}s")
    
    print(f"\n🎯 Accuracy:")
    print(f"  Role matches: {role_matches}/{total} ({role_matches/total*100:.1f}%)")
    print(f"  Family matches: {family_matches}/{total} ({family_matches/total*100:.1f}%)")
    print(f"  Average confidence: {avg_confidence:.0%}")
    
    print(f"\n✅ Quality:")
    print(f"  Verification approved: {approved}/{total} ({approved/total*100:.1f}%)")
    
    # Group by category
    print(f"\n📊 By Category:")
    categories = {}
    for r in results:
        cat = r["sample"]["category"]
        if cat not in categories:
            categories[cat] = {"total": 0, "role_match": 0, "family_match": 0}
        categories[cat]["total"] += 1
        if r.get("role_match", False):
            categories[cat]["role_match"] += 1
        if r.get("family_match", False):
            categories[cat]["family_match"] += 1
    
    for cat, stats in sorted(categories.items()):
        role_pct = stats["role_match"] / stats["total"] * 100 if stats["total"] > 0 else 0
        family_pct = stats["family_match"] / stats["total"] * 100 if stats["total"] > 0 else 0
        print(f"  {cat}:")
        print(f"    Role: {stats['role_match']}/{stats['total']} ({role_pct:.0f}%)")
        print(f"    Family: {stats['family_match']}/{stats['total']} ({family_pct:.0f}%)")


def main():
    """Run all test samples and generate report."""
    print("="*80)
    print("PURE LLM WORKFLOW - COMPREHENSIVE TEST SUITE")
    print("="*80)
    print(f"Testing {len(TEST_SAMPLES)} diverse reagent samples")
    print(f"Model: deepseek-v3.2-exp (Aliyun)")
    print()
    
    # Setup
    registry_dir = Path("../data/reagents")
    llm_client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
    
    # Run all tests
    results = []
    for i, sample in enumerate(TEST_SAMPLES, 1):
        print(f"\n[{i}/{len(TEST_SAMPLES)}]", end=" ")
        result = run_test_sample(sample, llm_client, registry_dir)
        results.append(result)
        time.sleep(1)  # Rate limiting
    
    # Print summary
    print_summary(results)
    
    # Save detailed results
    output_file = Path("test_pure_llm_results.json")
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)
    
    print(f"\n💾 Detailed results saved to: {output_file}")
    print("\n✅ Test suite complete!")


if __name__ == "__main__":
    main()
