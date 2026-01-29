import json

try:
    with open("chemtools/taxonomy/data/organic_compounds.v1.3.json", "r") as f:
        data = json.load(f)
    
    compounds = data.get("compounds", [])
    print(f"✓ Valid JSON: {len(compounds)} compounds")
    
    # Find our new Thiocarbonyl-SH entry
    thio_compound = next((c for c in compounds if c.get("id") == "Thiocarbonyl-SH"), None)
    if thio_compound:
        print(f"\n✓ Found Thiocarbonyl-SH:")
        print(f"  ID: {thio_compound.get('id')}")
        print(f"  SMARTS: {thio_compound.get('smarts')}")
        print(f"  Description: {thio_compound.get('description')}")
    else:
        print("\n❌ Thiocarbonyl-SH not found in compounds")
        
except json.JSONDecodeError as e:
    print(f"❌ JSON Parse Error: {e}")
except Exception as e:
    print(f"❌ Error: {e}")
