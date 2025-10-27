"""Analyze reactants field definitions using the shared ChemTools taxonomy."""

from collections import defaultdict

from chemtools.reagent import (
    get_reaction_type_definitions,
    get_reactant_type_definitions,
)

def load_json_files():
    """Load both JSON taxonomies from chemtools.reagent."""
    reaction_types = get_reaction_type_definitions()
    reactant_types = get_reactant_type_definitions()
    return reaction_types, reactant_types

def extract_valid_reactant_ids(reactant_types):
    """Extract all valid reactant type IDs and names."""
    valid_ids = set()
    valid_names = set()
    
    for category_key, category_data in reactant_types.items():
        # Add category key
        valid_ids.add(category_key)
        
        # Add all member IDs
        for member in category_data.get("members", []):
            valid_ids.add(member["id"])
            valid_names.add(member["name"])
    
    return valid_ids, valid_names

def analyze_reactants_field(reaction_types, valid_ids, valid_names):
    """Analyze reactants field in reaction_types and suggest corrections."""
    
    print("=" * 80)
    print("REACTANTS FIELD ANALYSIS")
    print("=" * 80)
    print(f"\nValid reactant type IDs: {len(valid_ids)}")
    print(f"Valid reactant type names: {len(valid_names)}")
    
    issues = []
    suggestions = defaultdict(list)
    total_reactions = 0
    
    for category_key, category_data in reaction_types.items():
        for reaction in category_data["reactions"]:
            total_reactions += 1
            reaction_id = reaction["id"]
            reactants = reaction.get("reactants", [])
            
            for reactant_str in reactants:
                # Try to map the reactant string to valid IDs
                matched = False
                
                # Direct match with category key
                if reactant_str in valid_ids:
                    matched = True
                    continue
                
                # Check if it contains a valid ID
                for valid_id in valid_ids:
                    if valid_id in reactant_str:
                        matched = True
                        break
                
                # Check if it contains a valid name
                if not matched:
                    for valid_name in valid_names:
                        if valid_name.lower() in reactant_str.lower():
                            matched = True
                            break
                
                # If no match found, record as issue
                if not matched:
                    issues.append({
                        "reaction": reaction_id,
                        "reactant_string": reactant_str,
                        "category": category_data["category"]
                    })
                    
                    # Try to suggest a match
                    suggestion = suggest_reactant_type(reactant_str, valid_ids)
                    if suggestion:
                        suggestions[reactant_str].append(suggestion)
    
    print(f"\nTotal reactions analyzed: {total_reactions}")
    print(f"Issues found: {len(issues)}")
    
    if issues:
        print("\n" + "=" * 80)
        print("ISSUES FOUND (reactants not mapping to reactant_types.json)")
        print("=" * 80)
        
        for issue in issues:
            print(f"\nReaction: {issue['reaction']} ({issue['category']})")
            print(f"  Reactant string: '{issue['reactant_string']}'")
            
            if issue['reactant_string'] in suggestions:
                print(f"  Suggested mapping: {suggestions[issue['reactant_string']]}")
    
    # Show reactant type coverage
    print("\n" + "=" * 80)
    print("REACTANT TYPE USAGE IN REACTIONS")
    print("=" * 80)
    
    used_types = set()
    for category_data in reaction_types.values():
        for reaction in category_data["reactions"]:
            for reactant_str in reaction.get("reactants", []):
                # Extract IDs mentioned
                for valid_id in valid_ids:
                    if valid_id in reactant_str:
                        used_types.add(valid_id)
    
    unused_types = valid_ids - used_types
    
    print(f"\nReactant types used: {len(used_types)}/{len(valid_ids)}")
    print(f"Coverage: {len(used_types)/len(valid_ids)*100:.1f}%")
    
    if unused_types:
        print(f"\nUnused reactant types ({len(unused_types)}):")
        for utype in sorted(unused_types):
            print(f"  - {utype}")
    
    return issues, suggestions

def suggest_reactant_type(reactant_str, valid_ids):
    """Suggest a matching reactant type based on the string."""
    reactant_lower = reactant_str.lower()
    
    # Common mappings
    mappings = {
        "arx": "ArX*",
        "aryl halide": "ArX*",
        "aryl halides": "ArX*",
        "vinyl halide": "VinylX*",
        "alkyl halide": "Alkyl-X",
        "r-x": "Alkyl-X",
        "amine": "RNH2/R2NH",
        "r2nh": "RNH2/R2NH",
        "rnh2": "RNH2/R2NH",
        "arnh2": "ArNH2/Ar2NH",
        "aniline": "ArNH2/Ar2NH",
        "alcohol": "ROH",
        "roh": "ROH",
        "phenol": "ArOH",
        "aroh": "ArOH",
        "thiol": "RSH",
        "boronic acid": "ArB*",
        "boronate": "ArB*",
        "grignard": "RMgX",
        "organometallic": "R-M",
        "organozinc": "RZnX",
        "organolithium": "RLi",
        "alkene": "alkene",
        "alkyne": "Alkyne",
        "aldehyde": "Aldehyde",
        "ketone": "Ketone",
        "carbonyl": ["Aldehyde", "Ketone"],
        "amide": "Amide",
        "lactam": "Amide",
        "urea": "Amide-type",
        "carboxylic acid": "Acyl-source",
        "ester": "Acyl-source",
    }
    
    for key, suggestion in mappings.items():
        if key in reactant_lower:
            return suggestion
    
    return None

if __name__ == "__main__":
    reaction_types, reactant_types = load_json_files()
    valid_ids, valid_names = extract_valid_reactant_ids(reactant_types)
    
    print(f"\nValid reactant type categories: {list(reactant_types.keys())}")
    print(f"\nSample valid IDs: {list(valid_ids)[:20]}")
    
    issues, suggestions = analyze_reactants_field(reaction_types, valid_ids, valid_names)
    
    # Generate recommended reactants field format
    print("\n" + "=" * 80)
    print("RECOMMENDED REACTANTS FIELD FORMAT")
    print("=" * 80)
    print("""
For consistency, reactants should reference reactant_types.json using:
1. Category keys (e.g., "ArX*", "ArB*", "Aliphatic-amine")
2. Member IDs (e.g., "ArBr", "ArCl", "RNH2")

Examples:
  "reactants": ["ArX*", "ArB*"]  # Uses category keys
  "reactants": ["ArBr", "ArB(OH)2"]  # Uses specific member IDs
  "reactants": ["ArX*", "Aliphatic-amine"]  # Mixed approach

Avoid free-text descriptions like:
  "reactants": ["ArX", "ArB(OH)2 or ArB(OR)2"]  # Too verbose
  "reactants": ["Aryl halides", "Boronic acid"]  # Not using defined IDs
    """)
