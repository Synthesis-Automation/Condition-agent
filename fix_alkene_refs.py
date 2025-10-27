"""Fix alkene reactant_type_id references from lowercase to capitalized"""
import json
from pathlib import Path

# Load reaction_types.json
filepath = Path("chemtools/taxonomy/data/reaction_types.json")
with open(filepath, encoding="utf-8") as f:
    data = json.load(f)

# Count and fix
count = 0
for reaction in data:
    if "reactants" in reaction:
        for reactant in reaction["reactants"]:
            if reactant.get("reactant_type_id") == "alkene":
                reactant["reactant_type_id"] = "Alkene"
                count += 1
                print(f"Fixed in {reaction['id']}: {reaction['name']}")

# Save back
with open(filepath, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2, ensure_ascii=False)

print(f"\n✓ Total fixed: {count}")
print("✓ reaction_types.json updated successfully")
