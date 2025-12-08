"""
Simplify reaction_types.json by removing redundant fields and null values.

Changes:
1. Remove 'required_roles' (always empty)
2. Remove null 'notes' and 'stoichiometry' from reactants - only keep if they have values
3. Remove 'original_tokens' if same as 'reactant_type_id'
4. Remove empty 'patterns' and 'example_reactants' arrays
5. Simplify pattern objects by removing redundant fields
6. Remove 'source_ids' if it matches the id
"""

import json
from pathlib import Path

INPUT_FILE = Path(__file__).parent.parent / "chemtools/taxonomy/data/reaction_types.json"
OUTPUT_FILE = Path(__file__).parent.parent / "chemtools/taxonomy/data/reaction_types_simplified.json"


def simplify_reactant(reactant: dict) -> dict:
    """Simplify a reactant object by removing null/redundant fields."""
    simplified = {}
    
    # Always keep reactant_type_id
    simplified["type"] = reactant.get("reactant_type_id", "")
    
    # Only keep original_tokens if different from type
    original = reactant.get("original_tokens", [])
    if original and original != [simplified["type"]]:
        simplified["aliases"] = original
    
    # Only keep notes if not null
    if reactant.get("notes"):
        simplified["notes"] = reactant["notes"]
    
    # Only keep stoichiometry if not null
    if reactant.get("stoichiometry"):
        simplified["stoichiometry"] = reactant["stoichiometry"]
    
    return simplified


def simplify_pattern(pattern: dict) -> dict:
    """Simplify a SMARTS pattern object."""
    simplified = {}
    
    # Keep essential fields
    if pattern.get("pattern_id"):
        simplified["id"] = pattern["pattern_id"]
    
    # smirks is the actual pattern - rename to 'smarts' for clarity
    if pattern.get("smirks"):
        simplified["smarts"] = pattern["smirks"]
    
    # Keep source and scope if informative
    if pattern.get("source"):
        simplified["source"] = pattern["source"]
    
    # Scope describes what the pattern matches
    if pattern.get("scope"):
        simplified["scope"] = pattern["scope"]
    
    return simplified


def simplify_example(example: dict) -> dict:
    """Simplify an example reactants object."""
    simplified = {}
    
    if example.get("educt1_smiles"):
        simplified["reactant1"] = example["educt1_smiles"]
    if example.get("educt2_smiles"):
        simplified["reactant2"] = example["educt2_smiles"]
    if example.get("source"):
        simplified["source"] = example["source"]
    
    return simplified


def simplify_entry(entry: dict) -> dict:
    """Simplify a single reaction type entry."""
    simplified = {}
    
    # Essential identification fields
    simplified["id"] = entry.get("id", "")
    simplified["name"] = entry.get("name", "")
    simplified["category"] = entry.get("category_id", "")
    
    # Keep aliases if non-empty
    aliases = entry.get("aliases", [])
    if aliases:
        simplified["aliases"] = aliases
    
    # Keep description
    if entry.get("description"):
        simplified["description"] = entry["description"]
    
    # Simplify metadata - flatten typical_catalysts
    metadata = entry.get("metadata", {})
    if metadata.get("typical_catalysts"):
        simplified["catalysts"] = metadata["typical_catalysts"]
    if metadata.get("typical_conditions"):
        simplified["conditions"] = metadata["typical_conditions"]
    
    # Simplify reactants
    reactants = entry.get("reactants", [])
    if reactants:
        simplified["reactants"] = [simplify_reactant(r) for r in reactants]
    
    # Simplify patterns - only include if non-empty
    patterns = entry.get("patterns", [])
    if patterns:
        simplified["patterns"] = [simplify_pattern(p) for p in patterns]
    
    # Simplify examples - only include if non-empty
    examples = entry.get("example_reactants", [])
    if examples:
        simplified["examples"] = [simplify_example(e) for e in examples]
    
    return simplified


def main():
    # Load original
    with open(INPUT_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    print(f"Loaded {len(data)} entries from {INPUT_FILE}")
    
    # Simplify
    simplified = [simplify_entry(entry) for entry in data]
    
    # Write output
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(simplified, f, indent=2, ensure_ascii=False)
    
    # Calculate size reduction
    original_size = INPUT_FILE.stat().st_size
    output_size = OUTPUT_FILE.stat().st_size
    reduction = (1 - output_size / original_size) * 100
    
    print(f"Written {len(simplified)} entries to {OUTPUT_FILE}")
    print(f"Size: {original_size:,} -> {output_size:,} bytes ({reduction:.1f}% reduction)")
    
    # Validate structure
    print("\n=== Validation ===")
    with_patterns = sum(1 for e in simplified if e.get("patterns"))
    with_examples = sum(1 for e in simplified if e.get("examples"))
    with_reactants = sum(1 for e in simplified if e.get("reactants"))
    print(f"Entries with patterns: {with_patterns}/{len(simplified)}")
    print(f"Entries with examples: {with_examples}/{len(simplified)}")
    print(f"Entries with reactants: {with_reactants}/{len(simplified)}")


if __name__ == "__main__":
    main()
