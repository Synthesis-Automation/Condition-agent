import json

try:
    with open("chemtools/taxonomy/data/compound_logic.json", "r", encoding="utf-8") as f:
        data = json.load(f)
    
    motif_sets = data.get("motif_sets", {})
    print(f"✓ Valid JSON: {len(motif_sets)} motif sets")
    
    # Check @thiols_sh group
    thiols_sh = motif_sets.get("thiols_sh", {})
    if thiols_sh:
        members = thiols_sh.get("members", [])
        print(f"\n✓ @thiols_sh group found:")
        print(f"  Members ({len(members)}):")
        for m in members:
            print(f"    - {m}")
        
        if "Thiocarbonyl-SH" in members:
            print(f"\n✓ Thiocarbonyl-SH is in @thiols_sh members")
        else:
            print(f"\n❌ Thiocarbonyl-SH NOT in @thiols_sh members")
    else:
        print("\n❌ @thiols_sh group not found")
        
except json.JSONDecodeError as e:
    print(f"❌ JSON Parse Error: {e}")
except Exception as e:
    print(f"❌ Error: {e}")