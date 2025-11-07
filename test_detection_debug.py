"""Debug detection results"""
from chemtools import detect_reaction

reactions = {
    "Suzuki": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    "Amide": "CC(=O)Cl.NCC>>CC(=O)NCC",
    "SN2": "CBr.NCC>>CCNCC",
    "SNAr": "[O-][N+](=O)c1ccc([N+](=O)[O-])cc1F.NH2C2H4OH>>[O-][N+](=O)c1ccc([N+](=O)[O-])cc1NCCO",
    "Radical": "CC.CCl4>>CCl+CCCl3",
}

for name, smiles in reactions.items():
    try:
        result = detect_reaction(smiles)
        print(f"\n{name}:")
        print(f"  Family: {result.get('family')}")
        print(f"  Confidence: {result.get('confidence', 0):.2f}")
        print(f"  Method: {result.get('method')}")
        print(f"  Details:")
        details = result.get('details', {})
        if 'rule_prediction' in details:
            rule_pred = details['rule_prediction']
            print(f"    Raw family: {rule_pred.get('raw_family')}")
            print(f"    Rule confidence: {rule_pred.get('confidence')}")
    except Exception as e:
        print(f"{name} -> ERROR: {e}")
        import traceback
        traceback.print_exc()
