"""
Test LLM synthesis standard format conversion for robotic execution.
Ensures LLM output matches the same structure as ML/Rule/Protocol outputs.
"""

import sys
sys.path.insert(0, '.')

from llmtools.recommendation_llm import convert_llm_synthesis_to_standard_format

print("=" * 70)
print("Testing LLM Synthesis Standard Format Conversion")
print("=" * 70)

# Mock LLM synthesis result (from synthesize_recommendations_llm)
mock_synthesis_result = {
    "status": "success",
    "synthesis": {
        "consensus_analysis": {
            "catalyst": {
                "agreement": "high",
                "consensus_value": "Pd(PPh3)4",
                "notes": "All sources recommend Pd(0) catalyst"
            },
            "solvent": {
                "agreement": "high",
                "consensus_value": "THF",
                "notes": "THF preferred for solubility"
            },
            "temperature": {
                "agreement": "medium",
                "consensus_value": "80°C",
                "notes": "ML suggests 80-90°C, rules suggest 70-85°C"
            },
            "base": {
                "agreement": "high",
                "consensus_value": "K2CO3",
                "notes": "Standard base for Suzuki coupling"
            }
        },
        "confidence_level": "high",
        "confidence_reasoning": "All sources agree on Pd(PPh3)4 catalyst system. ML similarity 0.95, rule-based high-priority match, protocol confirms conditions.",
        "recommended_condition": {
            "catalyst": "Pd(PPh3)4",
            "ligand": "PPh3 (pre-complexed)",
            "solvent": "THF",
            "temperature": "80°C",
            "base": "K2CO3 (2 equiv)",
            "additive": "None",
            "rationale": "High consensus across all sources for standard Suzuki coupling with electron-neutral aryl bromide."
        },
        "backup_conditions": [
            {
                "catalyst": "Pd(dppf)Cl2",
                "ligand": "dppf (pre-complexed)",
                "solvent": "Dioxane/H2O (4:1)",
                "temperature": "90°C",
                "base": "K3PO4 (2 equiv)",
                "when_to_use": "If Pd(PPh3)4 gives <30% conversion after 12h at 80°C"
            },
            {
                "catalyst": "Pd(OAc)2",
                "ligand": "XPhos (5 mol%)",
                "solvent": "Toluene/EtOH/H2O (4:1:1)",
                "temperature": "100°C",
                "base": "Cs2CO3 (2 equiv)",
                "when_to_use": "If backup 1 fails (<20% after 12h) - more active catalyst system"
            }
        ],
        "warnings": [
            "Monitor for proto-debromination if reaction stalls",
            "Degassed solvent recommended for Pd(0) catalyst",
            "Boronic acid may oxidize if exposed to air for prolonged periods"
        ],
        "source_comparison": {
            "ml_contribution": "High - similarity 0.95, top precedent uses Pd(PPh3)4/THF/K2CO3",
            "rule_contribution": "Primary - high-priority SCDB scheme match for aryl bromide coupling",
            "protocol_contribution": "Supporting - literature protocols confirm 80°C as optimal"
        }
    },
    "sources_used": {
        "ml_precedents": 3,
        "rule_matches": 1,
        "protocol_procedures": 2
    },
    "llm_metadata": {
        "model": "deepseek-v3.2-exp",
        "tokens": 1480,
        "latency_ms": 43600,
        "processing_time_ms": 45823
    }
}

reaction_smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

print("\n1. Converting LLM synthesis to standard format...")
print("-" * 70)

standard_output = convert_llm_synthesis_to_standard_format(
    reaction_smiles=reaction_smiles,
    synthesis_result=mock_synthesis_result,
    requested_type="Suzuki_Miyaura",
    processing_time_ms=45823.0
)

print("✅ Conversion successful\n")

# Verify structure
print("2. Verifying standard structure...")
print("-" * 70)

required_keys = {"meta", "input", "detection", "recommended_conditions"}
if not required_keys.issubset(standard_output.keys()):
    print(f"❌ Missing required keys! Expected: {required_keys}, Got: {set(standard_output.keys())}")
    sys.exit(1)
print(f"✅ All required keys present: {required_keys}")

# Check meta section
meta = standard_output["meta"]
assert meta["model"] == "LLM-synthesis", f"Wrong model: {meta['model']}"
assert meta["status"] == "success", f"Wrong status: {meta['status']}"
assert "processing_time_ms" in meta, "Missing processing_time_ms"
print("✅ Meta section valid")

# Check input section
input_data = standard_output["input"]
assert input_data["reaction_smiles"] == reaction_smiles, "Wrong reaction SMILES"
assert input_data.get("requested_reaction_type") == "Suzuki_Miyaura", "Wrong requested type"
print("✅ Input section valid")

# Check detection section
detection = standard_output["detection"]
assert "confidence" in detection, "Missing confidence in detection"
assert detection["method"] == "llm-multi-source", f"Wrong method: {detection['method']}"
print(f"✅ Detection section valid (confidence: {detection['confidence']})")

# Check recommended_conditions
recommendations = standard_output["recommended_conditions"]
assert isinstance(recommendations, list), "Recommendations should be a list"
assert len(recommendations) > 0, "Should have at least 1 recommendation"
print(f"✅ Recommendations present: {len(recommendations)} recommendation(s)")

# Check first recommendation structure
rec1 = recommendations[0]
required_rec_keys = {"rank", "chemicals", "conditions"}
if not required_rec_keys.issubset(rec1.keys()):
    print(f"❌ Missing keys in recommendation! Expected: {required_rec_keys}, Got: {set(rec1.keys())}")
    sys.exit(1)
print(f"✅ Recommendation 1 has required keys: {required_rec_keys}")

# Check chemicals array
chemicals = rec1["chemicals"]
assert isinstance(chemicals, list), "Chemicals should be a list"
assert len(chemicals) > 0, "Should have chemicals"
print(f"✅ Chemicals array present: {len(chemicals)} chemical(s)")

# Verify chemical structure
catalyst_found = False
for chem in chemicals:
    if chem.get("role") == "metal_precursor":
        catalyst_found = True
        print(f"   - Catalyst: {chem.get('name', 'N/A')}")
        assert "name" in chem, "Chemical should have name"
        assert "role" in chem, "Chemical should have role"
        break

if not catalyst_found:
    print("⚠️  Warning: No catalyst found in chemicals")

# Check conditions
conditions = rec1["conditions"]
assert isinstance(conditions, dict) or isinstance(conditions, list), "Conditions should be dict or list"
print(f"✅ Conditions present: {conditions}")

# Check extras section for original synthesis data
extras = standard_output.get("extras", {})
assert "llm_synthesis" in extras, "Should preserve original LLM synthesis in extras"
assert "llm_metadata" in extras, "Should preserve LLM metadata in extras"
print("✅ Original synthesis data preserved in extras")

print("\n3. Checking robotic execution compatibility...")
print("-" * 70)

# Verify format matches output_formatter.py standard
from chemtools.output_formatter import ensure_standard_output

try:
    validated = ensure_standard_output(
        standard_output,
        default_model="LLM-synthesis",
        fallback_reaction_smiles=reaction_smiles,
    )
    print("✅ Output passes ensure_standard_output validation")
except Exception as e:
    print(f"❌ Validation failed: {e}")
    sys.exit(1)

# Check if output can be used for robotic execution
print("\n4. Robotic execution readiness:")
print("-" * 70)

print(f"   Reaction: {standard_output['input']['reaction_smiles']}")
print(f"   Family: {standard_output['detection']['family']}")
print(f"   Confidence: {standard_output['detection']['confidence']:.2%}")
print(f"   Recommendations: {len(recommendations)}")

for i, rec in enumerate(recommendations[:3], 1):
    print(f"\n   Recommendation {i}:")
    print(f"     Rank: {rec['rank']}")
    print(f"     Chemicals: {len(rec['chemicals'])}")
    
    # List key chemicals
    for chem in rec["chemicals"]:
        role = chem.get("role", "unknown")
        name = chem.get("name") or chem.get("abbreviation", "N/A")
        if role in ["metal_precursor", "ligand", "base", "solvent"]:
            print(f"       - {role.capitalize()}: {name}")
    
    # Show conditions
    cond = rec["conditions"]
    if isinstance(cond, dict):
        temp = cond.get("temperature")
        if temp:
            print(f"       - Temperature: {temp}")
    
    # Show reasoning if available
    summary = rec.get("summary", {})
    if summary.get("rationale"):
        rationale = summary["rationale"]
        short_rationale = rationale[:80] + "..." if len(rationale) > 80 else rationale
        print(f"       - Rationale: {short_rationale}")

print("\n" + "=" * 70)
print("✅ All tests passed! LLM synthesis standard format ready for robots.")
print("=" * 70)
print("\nThe LLM synthesis output now matches the same structure as ML/Rule/Protocol,")
print("making it fully compatible with robotic execution systems.")
print()
