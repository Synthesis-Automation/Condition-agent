#!/usr/bin/env python3
"""
Map extracted reagents from z-Score Peaks dataset to reagent registry system roles.
Converts dataset roles (like "Coupling Reagent") to registry roles (like "condensation_agent").
"""

import re
from pathlib import Path
from typing import Dict, List, Tuple

# Mapping from dataset roles to registry system roles
ROLE_MAPPING = {
    "Additive": "additive",
    "Base": "base",
    "Catalyst": "metal_precursor",  # Most catalysts in the dataset are metal precursors
    "Coupling Reagent": "condensation_agent",  # As specified by user
    "Ligand": "ligand",
    "Solvent": "solvent",
}

# Additional role classifications based on reagent names (for catalysts)
# Some catalysts should be classified as preformed_metal_catalyst
PREFORMED_CATALYST_PATTERNS = [
    r"Pd\(PPh3\)",
    r"Pd\(OAc\)2",
    r"Pd\(dba\)",
    r"Pd\(acac\)",
    r"Pd2\(dba\)3",
    r"PdCl2\(PPh3\)2",
    r"PdCl2\(dppf\)",
    r"PdCl2\(AmPhos\)",
    r"PdCl2\(XPhos\)",
    r"Pd\(amphos\)",
    r"Pd\(dtbpf\)",
    r"Pd\(QPhos\)",
    r"Pd\(XantPhos\)",
    r"Pd\(cataCXium\)",
    r"tBuBrettPhos Pd",
    r"tBuXPhos Pd",
    r"RuPhos Pd",
    r"XPhos Pd",
    r"SPhos Pd",
    r"JohnPhos Pd",
    r"CyJohnPhos Pd",
    r"PCy3 Pd",
    r"PEPPSI",
]

def is_preformed_catalyst(reagent_name: str) -> bool:
    """Check if a catalyst should be classified as preformed_metal_catalyst."""
    for pattern in PREFORMED_CATALYST_PATTERNS:
        if re.search(pattern, reagent_name, re.IGNORECASE):
            return True
    return False

def parse_reagent_section(content: str, start_line: str) -> Tuple[str, List[Dict]]:
    """Parse a reagent section from the markdown file."""
    lines = content.split('\n')
    
    # Find the section
    section_start = None
    for i, line in enumerate(lines):
        if line.startswith(start_line):
            section_start = i
            break
    
    if section_start is None:
        return "", []
    
    # Extract original role and count from header
    header_match = re.match(r'## (.+?) \((\d+) reagents\)', lines[section_start])
    if header_match:
        original_role = header_match.group(1)
        count = header_match.group(2)
    else:
        return "", []
    
    # Find next section or end of file
    section_end = len(lines)
    for i in range(section_start + 1, len(lines)):
        if lines[i].startswith('## ') and not lines[i].startswith('### '):
            section_end = i
            break
    
    # Parse reagents in this section
    reagents = []
    current_reagent = None
    
    for i in range(section_start + 1, section_end):
        line = lines[i].strip()
        
        if line.startswith('### '):
            # New reagent
            if current_reagent:
                reagents.append(current_reagent)
            
            reagent_name = line[4:].strip()
            current_reagent = {
                'name': reagent_name,
                'role': None,
                'occurrences': None,
                'reactions': []
            }
        elif line.startswith('- **Role:**') and current_reagent:
            current_reagent['role'] = line.split('**Role:**')[1].strip()
        elif line.startswith('- **Occurrences:**') and current_reagent:
            occ_match = re.search(r'\d+', line)
            if occ_match:
                current_reagent['occurrences'] = int(occ_match.group())
        elif line.startswith('- **Used in reactions:**') and current_reagent:
            reactions_str = line.split('**Used in reactions:**')[1].strip()
            current_reagent['reactions'] = [r.strip() for r in reactions_str.split(',')]
    
    # Add last reagent
    if current_reagent:
        reagents.append(current_reagent)
    
    return original_role, reagents

def main():
    input_file = Path("extracted_reagents_registry.md")
    output_file = Path("reagents_mapped_to_registry_roles.md")
    
    if not input_file.exists():
        print(f"Error: {input_file} not found")
        return
    
    content = input_file.read_text(encoding='utf-8')
    
    # Parse all sections
    all_reagents_by_registry_role: Dict[str, List[Dict]] = {
        role: [] for role in set(ROLE_MAPPING.values())
    }
    all_reagents_by_registry_role["preformed_metal_catalyst"] = []
    
    original_roles_parsed = []
    total_reagents = 0
    
    for dataset_role, registry_role in ROLE_MAPPING.items():
        print(f"Processing {dataset_role} -> {registry_role}...")
        
        # Find the section header pattern
        section_pattern = f"## {dataset_role}"
        original_role, reagents = parse_reagent_section(content, section_pattern)
        
        if not reagents:
            print(f"  Warning: No reagents found for {dataset_role}")
            continue
        
        original_roles_parsed.append((dataset_role, len(reagents)))
        total_reagents += len(reagents)
        
        # Classify reagents
        for reagent in reagents:
            # Special handling for catalysts
            if registry_role == "metal_precursor" and is_preformed_catalyst(reagent['name']):
                reagent['registry_role'] = "preformed_metal_catalyst"
                all_reagents_by_registry_role["preformed_metal_catalyst"].append(reagent)
            else:
                reagent['registry_role'] = registry_role
                all_reagents_by_registry_role[registry_role].append(reagent)
        
        print(f"  Processed {len(reagents)} reagents")
    
    # Generate output markdown
    output_lines = [
        "# Reagent Registry Role Mapping",
        "",
        "## Overview",
        "",
        f"- **Source Dataset:** z-Score Peaks with FG.csv",
        f"- **Total Reactions:** 66,308",
        f"- **Total Reagents Mapped:** {total_reagents}",
        f"- **Registry System Roles:** {len([r for r in all_reagents_by_registry_role.values() if r])}",
        "",
        "## Role Mapping",
        "",
        "| Dataset Role | Registry Role | Count |",
        "|--------------|---------------|-------|",
    ]
    
    for dataset_role, count in original_roles_parsed:
        registry_role = ROLE_MAPPING[dataset_role]
        output_lines.append(f"| {dataset_role} | `{registry_role}` | {count} |")
    
    # Add note about catalyst split
    metal_precursor_count = len(all_reagents_by_registry_role["metal_precursor"])
    preformed_count = len(all_reagents_by_registry_role["preformed_metal_catalyst"])
    if preformed_count > 0:
        output_lines.append(f"| Catalyst (metal precursors) | `metal_precursor` | {metal_precursor_count} |")
        output_lines.append(f"| Catalyst (preformed) | `preformed_metal_catalyst` | {preformed_count} |")
    
    output_lines.extend([
        "",
        "---",
        "",
    ])
    
    # Role descriptions
    role_descriptions = {
        "additive": "Phase-transfer agents, halide scavengers, fluoride sources, and related modifiers.",
        "base": "Bronsted/Lewis bases spanning amides, alkoxides, carbonates, superbases.",
        "metal_precursor": "Metal salts or complexes that generate the catalytically active species.",
        "preformed_metal_catalyst": "Precatalysts supplied with ligands; typically used as-is.",
        "condensation_agent": "Carbodiimides, uronium/phosphonium activators, and similar condensers.",
        "ligand": "Ligands including phosphines, NHCs, diimines, and ancillary donor sets.",
        "solvent": "Reaction media categorized by polarity, coordinating ability, and safety profile.",
    }
    
    # Generate sections for each registry role
    for registry_role in sorted(all_reagents_by_registry_role.keys()):
        reagents = all_reagents_by_registry_role[registry_role]
        if not reagents:
            continue
        
        output_lines.extend([
            f"## {registry_role.replace('_', ' ').title()} ({len(reagents)} reagents)",
            "",
            f"**Registry Role:** `{registry_role}`",
            "",
            f"**Description:** {role_descriptions.get(registry_role, 'N/A')}",
            "",
        ])
        
        # Sort reagents by occurrence count (descending)
        reagents_sorted = sorted(reagents, key=lambda x: x['occurrences'] or 0, reverse=True)
        
        # Add reagent entries
        for reagent in reagents_sorted:
            output_lines.extend([
                f"### {reagent['name']}",
                "",
                f"- **Registry Role:** `{registry_role}`",
                f"- **Original Role:** {reagent['role']}",
                f"- **Occurrences:** {reagent['occurrences']}",
                f"- **Used in {len(reagent['reactions'])} reaction types**",
                "",
            ])
            
            # Add reaction types (max 10 to keep readable)
            if reagent['reactions']:
                reactions_display = reagent['reactions'][:10]
                if len(reagent['reactions']) > 10:
                    reactions_display.append(f"... and {len(reagent['reactions']) - 10} more")
                output_lines.append("  **Reaction types:** " + ", ".join(reactions_display))
                output_lines.append("")
        
        output_lines.append("---")
        output_lines.append("")
    
    # Add summary statistics
    output_lines.extend([
        "## Summary Statistics",
        "",
        "### Reagents by Registry Role",
        "",
        "| Registry Role | Count | Top 3 Reagents |",
        "|---------------|-------|----------------|",
    ])
    
    for registry_role in sorted(all_reagents_by_registry_role.keys()):
        reagents = all_reagents_by_registry_role[registry_role]
        if not reagents:
            continue
        
        # Get top 3 by occurrence
        top_3 = sorted(reagents, key=lambda x: x['occurrences'] or 0, reverse=True)[:3]
        top_3_names = [f"{r['name']} ({r['occurrences']})" for r in top_3]
        top_3_str = "; ".join(top_3_names)
        
        output_lines.append(f"| `{registry_role}` | {len(reagents)} | {top_3_str} |")
    
    output_lines.extend([
        "",
        "---",
        "",
        "## Notes",
        "",
        "1. **Catalyst Classification:** Catalysts were split into two categories:",
        "   - `metal_precursor`: Simple metal salts (e.g., PdCl2, CuI)",
        "   - `preformed_metal_catalyst`: Pre-coordinated complexes (e.g., Pd(PPh3)4, PEPPSI-IPr)",
        "",
        "2. **Coupling Reagents:** Mapped to `condensation_agent` as per registry system convention.",
        "",
        "3. **Data Source:** Extracted from z-Score Peaks dataset (66,308 reactions, 42 reaction types).",
        "",
        "4. **Next Steps:**",
        "   - Cross-reference with existing reagent registry to identify new entries",
        "   - Add CAS numbers and SMILES structures",
        "   - Standardize naming variations",
        "   - Assign family classifications within each role",
        "",
    ])
    
    # Write output
    output_file.write_text('\n'.join(output_lines), encoding='utf-8')
    
    print(f"\n✓ Mapping complete!")
    print(f"  Total reagents processed: {total_reagents}")
    print(f"  Output saved to: {output_file}")
    print(f"\nReagent counts by registry role:")
    for registry_role in sorted(all_reagents_by_registry_role.keys()):
        count = len(all_reagents_by_registry_role[registry_role])
        if count > 0:
            print(f"  {registry_role}: {count}")

if __name__ == "__main__":
    main()
