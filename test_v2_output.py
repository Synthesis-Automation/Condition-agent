
import json
from data_processor.v2_processor_core import V2Processor

def test_v2_processor():
    processor = V2Processor()
    
    # A simple Buchwald-Hartwig reaction: Bromobenzene + Aniline
    reaction_data = {
        "reaction_id": "test_001",
        "reaction_smiles": "c1ccc(Br)cc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "reagents": [
            {"name": "Pd(PPh3)4", "cas": "14221-01-3"},
            {"name": "tBuOK", "cas": "865-47-4"},
            {"name": "XPhos", "cas": "564483-18-7"}
        ],
        "solvents": [
            {"name": "Toluene", "cas": "108-88-3"}
        ],
        "conditions": {
            "temperature_c": 80,
            "yield_pct": 95
        }
    }
    
    result = processor.process_reaction(reaction_data)
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    test_v2_processor()
