"""
Test the new enhanced output format for both ML and rule-based recommendations.

This script demonstrates the improved JSON output structure with:
- Enhanced metadata (timestamps, processing time)
- Structured numerical data (temperature, time as objects)
- Enriched reagent information (CAS, InChI, properties)
- Confidence scores and ranking
- Expected yields and precedent counts
"""

import json
import os
import sys
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from chemtools import recommend, router, output_formatter
from chemtools.scdb_matcher import load_db, match


def test_ml_recommendations():
    """Test ML-based recommendations with new format."""
    
    print("="*80)
    print("ML-BASED RECOMMENDATIONS - NEW FORMAT")
    print("="*80)
    
    # Test reaction: C-N coupling
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    
    print(f"\nReaction: {reaction}")
    print(f"Expected Type: C-N Coupling (Pd)")
    
    # Get recommendations
    data = recommend.recommend_conditions_structured(
        reaction=reaction,
        reaction_type="C_N_Coupling_Pd",
        k=10,
        limit=6,  # Get 6 to filter to top 3
        relax={
            "use_drfp": False,
            "precompute_drfp": False,
            "use_rxn_insight": True,
        },
        constraints=None,
    )
    
    # Build enhanced recommendations
    recommendations_data = []
    conditions_list = data.get("conditions", [])[:3]
    
    for i, cond in enumerate(conditions_list, 1):
        # Extract reagents
        reagents = []
        
        # Add starting materials
        reactants = reaction.split(">>")[0].split(".")[:2]
        for j, sm_smiles in enumerate(reactants, 1):
            reagents.append({
                "id": f"SM{j}",
                "name": None,
                "abbreviation": None,
                "cas": None,
                "smiles": sm_smiles,
                "inchi_key": None,
                "role": "electrophile" if j == 1 else "nucleophile",
                "equivalents": {
                    "value": 1.0 if j == 1 else 1.2,
                    "range": [1.0, 1.0] if j == 1 else [1.1, 1.5],
                    "unit": "eq"
                }
            })
        
        # Enrich base
        base_name = cond.get("base", "")
        if base_name:
            reagents.append(output_formatter.enrich_reagent(
                name=base_name,
                reagent_type="base",
                role="base",
                equivalents=2.0,
                equiv_range=(1.5, 2.5),
                reagent_id="BASE1"
            ))
        
        # Enrich solvent
        solvent_name = cond.get("solvent", "")
        if solvent_name:
            reagents.append(output_formatter.enrich_reagent(
                name=solvent_name,
                reagent_type="solvent",
                role="solvent",
                reagent_id="SOLV1"
            ))
        
        # Extract similarity score
        similarity = cond.get("similarity", 0.85)
        
        # Format conditions
        temp = cond.get("temperature_C", 110)
        time_h = cond.get("time_h", 12)
        
        formatted_cond = output_formatter.format_conditions(
            temperature=temp,
            temp_range=(temp - 20, temp + 20),
            time_hours=time_h,
            time_range=(time_h / 2, time_h * 1.5),
            atmosphere="N2",
        )
        
        # Build recommendation
        rec = output_formatter.format_recommendation(
            rank=i,
            confidence=0.90 - (i - 1) * 0.05,
            support=cond.get("support", 10),
            reaction_smiles=reaction,
            reagents=reagents,
            conditions=formatted_cond,
            similarity_score=similarity,
            expected_yield=75.0,
            yield_range=(60.0, 85.0),
        )
        
        recommendations_data.append(rec)
    
    # Build full output
    detection = data.get("detection", {})
    output = output_formatter.format_ml_output(
        reaction_smiles=reaction,
        requested_type="C-N Coupling (Pd)",
        detected_type=detection.get("family", "Unknown"),
        detection_confidence=detection.get("confidence", 0.0),
        recommendations_data=recommendations_data,
        processing_time_ms=500.0,
    )
    
    # Print JSON
    print("\n" + "="*80)
    print("JSON OUTPUT:")
    print("="*80)
    print(json.dumps(output, indent=2))
    
    # Save to file
    output_path = ROOT / "output_ml_recommendations.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n�?Saved to: {output_path}")
    
    return output


def test_rule_recommendations():
    """Test rule-based recommendations with new format."""
    
    print("\n" + "="*80)
    print("RULE-BASED RECOMMENDATIONS - NEW FORMAT")
    print("="*80)
    
    # Test reaction: Suzuki coupling
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    
    print(f"\nReaction: {reaction}")
    print(f"Expected Type: Suzuki Coupling")
    
    # Load database and match
    db_path = ROOT / "data" / "conditionDB" / "suzuki_db.json"
    
    if not db_path.exists():
        print(f"�?Database not found: {db_path}")
        return None
    
    db = load_db(str(db_path))
    result = match(db, reaction)
    
    # Extract conditions
    conditions_dict = {}
    if hasattr(result, 'conditions'):
        conditions_dict = result.conditions
    elif hasattr(result, 'to_json_dict'):
        result_dict = result.to_json_dict()
        if 'conditions' in result_dict:
            conditions_dict = result_dict['conditions']
    
    if not conditions_dict:
        print("�?No conditions found")
        return None
    
    # Use formatter to expand conditions into 3 recommendations
    recommendations_data = output_formatter.expand_rule_conditions_to_recommendations(
        reaction_smiles=reaction,
        conditions_dict=conditions_dict,
        num_recommendations=3,
    )
    
    # Build full output
    output = output_formatter.format_rule_output(
        reaction_smiles=reaction,
        requested_type="Suzuki Coupling",
        detected_type="Suzuki_CC",
        recommendations_data=recommendations_data,
        database_name="Suzuki Coupling",
        processing_time_ms=50.0,
    )
    
    # Print JSON
    print("\n" + "="*80)
    print("JSON OUTPUT:")
    print("="*80)
    print(json.dumps(output, indent=2))
    
    # Save to file
    output_path = ROOT / "output_rule_recommendations.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n�?Saved to: {output_path}")
    
    return output


def main():
    """Run both tests."""
    
    print("\n" + "="*80)
    print("TESTING NEW ENHANCED OUTPUT FORMAT")
    print("="*80)
    
    # Test ML
    ml_output = test_ml_recommendations()
    
    # Test Rule-based
    rule_output = test_rule_recommendations()
    
    print("\n" + "="*80)
    print("TESTS COMPLETE")
    print("="*80)
    
    print("\nKey Improvements Demonstrated:")
    print("�?Structured metadata (timestamps, model info, processing time)")
    print("�?Consistent field naming (reagents, not chemicals)")
    print("�?Numerical data as objects (temperature, time with ranges)")
    print("�?Enriched reagent information (CAS, InChI, properties)")
    print("�?Confidence scores and ranking")
    print("�?Expected yields with ranges")
    print("�?Support/precedent counts")
    
    print("\nOutput files:")
    if ml_output:
        print(f"  �?output_ml_recommendations.json")
    if rule_output:
        print(f"  �?output_rule_recommendations.json")


if __name__ == "__main__":
    main()
