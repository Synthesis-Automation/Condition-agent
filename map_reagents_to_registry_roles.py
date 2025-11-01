#!/usr/bin/env python3
"""
Map extracted reagents from z-Score Peaks dataset to reagent registry system roles.
This version reflects the unified `metal_catalyst` role that covers both precursor salts
and pre-ligated complexes.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

# Mapping from dataset roles to registry system roles
ROLE_MAPPING: Dict[str, str] = {
    "Additive": "additive",
    "Base": "base",
    "Catalyst": "metal_catalyst",
    "Coupling Reagent": "condensation_agent",
    "Ligand": "ligand",
    "Solvent": "solvent",
}


def parse_reagent_section(content: str, start_line: str) -> Tuple[str, List[Dict]]:
    """Parse a reagent section from the markdown file."""
    lines = content.split("\n")

    # Find the section
    section_start: int | None = None
    for idx, line in enumerate(lines):
        if line.startswith(start_line):
            section_start = idx
            break

    if section_start is None:
        return "", []

    # Extract original role and count from header
    header_match = re.match(r"## (.+?) \((\d+) reagents\)", lines[section_start])
    if header_match:
        original_role = header_match.group(1)
        count = header_match.group(2)
    else:
        return "", []

    # Find next section or end of file
    section_end = len(lines)
    for idx in range(section_start + 1, len(lines)):
        if lines[idx].startswith("## ") and not lines[idx].startswith("### "):
            section_end = idx
            break

    # Parse reagents in this section
    reagents: List[Dict] = []
    current_reagent: Dict | None = None

    for idx in range(section_start + 1, section_end):
        line = lines[idx].strip()

        if line.startswith("### "):
            # New reagent
            if current_reagent:
                reagents.append(current_reagent)

            reagent_name = line[4:].strip()
            current_reagent = {
                "name": reagent_name,
                "role": None,
                "occurrences": None,
                "reactions": [],
            }
        elif line.startswith("- **Role:**") and current_reagent:
            current_reagent["role"] = line.split("**Role:**")[1].strip()
        elif line.startswith("- **Occurrences:**") and current_reagent:
            occ_match = re.search(r"\d+", line)
            if occ_match:
                current_reagent["occurrences"] = int(occ_match.group())
        elif line.startswith("- **Used in reactions:**") and current_reagent:
            reactions_str = line.split("**Used in reactions:**")[1].strip()
            current_reagent["reactions"] = [r.strip() for r in reactions_str.split(",")]

    if current_reagent:
        reagents.append(current_reagent)

    return original_role, reagents


def main() -> None:
    input_file = Path("extracted_reagents_registry.md")
    output_file = Path("reagents_mapped_to_registry_roles.md")

    if not input_file.exists():
        print(f"Error: {input_file} not found")
        return

    content = input_file.read_text(encoding="utf-8")

    # Prepare containers for results
    all_reagents_by_registry_role: Dict[str, List[Dict]] = {
        role: [] for role in sorted(set(ROLE_MAPPING.values()))
    }
    original_roles_parsed: List[Tuple[str, int]] = []
    total_reagents = 0

    # Parse each dataset role
    for dataset_role, registry_role in ROLE_MAPPING.items():
        print(f"Processing {dataset_role} -> {registry_role}...")
        section_pattern = f"## {dataset_role}"
        original_role, reagents = parse_reagent_section(content, section_pattern)

        if not reagents:
            print(f"  Warning: No reagents found for {dataset_role}")
            continue

        original_roles_parsed.append((original_role or dataset_role, len(reagents)))
        total_reagents += len(reagents)

        for reagent in reagents:
            reagent["registry_role"] = registry_role
            all_reagents_by_registry_role[registry_role].append(reagent)

        print(f"  Processed {len(reagents)} reagents")

    # Build markdown report
    output_lines = [
        "# Reagent Registry Role Mapping",
        "",
        "## Overview",
        "",
        "- **Source Dataset:** z-Score Peaks with FG.csv",
        "- **Total Reactions:** 66,308",
        f"- **Total Reagents Mapped:** {total_reagents}",
        (
            "- **Registry System Roles:** "
            f"{len([reagents for reagents in all_reagents_by_registry_role.values() if reagents])}"
        ),
        "",
        "## Role Mapping",
        "",
        "| Dataset Role | Registry Role | Count |",
        "|--------------|---------------|-------|",
    ]

    for dataset_role, count in original_roles_parsed:
        registry_role = ROLE_MAPPING.get(dataset_role, "metal_catalyst")
        output_lines.append(f"| {dataset_role} | `{registry_role}` | {count} |")

    output_lines.extend(["", "---", "", "## Registry Role Details", ""])

    role_descriptions = {
        "additive": (
            "Phase-transfer agents, halide scavengers, fluoride sources, and related modifiers."
        ),
        "base": (
            "Brønsted/Lewis bases spanning amides, alkoxides, carbonates, and superbases."
        ),
        "condensation_agent": (
            "Carbodiimides, uronium/phosphonium activators, and related coupling partners."
        ),
        "ligand": (
            "Ligands including phosphines, NHCs, diimines, and ancillary donor sets."
        ),
        "metal_catalyst": (
            "Metal salts or pre-ligated complexes that provide the catalytic metal source."
        ),
        "solvent": (
            "Reaction media categorized by polarity, coordinating ability, and safety profile."
        ),
    }

    for registry_role in sorted(all_reagents_by_registry_role.keys()):
        reagents = all_reagents_by_registry_role[registry_role]
        if not reagents:
            continue

        output_lines.extend(
            [
                f"## {registry_role.replace('_', ' ').title()} ({len(reagents)} reagents)",
                "",
                f"**Registry Role:** `{registry_role}`",
                "",
            ]
        )

        description = role_descriptions.get(
            registry_role, "Description forthcoming."
        )
        output_lines.append(f"{description}")
        output_lines.append("")

        for reagent in reagents:
            output_lines.extend(
                [
                    f"### {reagent['name']}",
                    "",
                    f"- **Original Role:** {reagent.get('role') or 'Unknown'}",
                    f"- **Occurrences:** {reagent.get('occurrences') or 'Unknown'}",
                    (
                        "- **Used in reactions:** "
                        + ", ".join(reagent.get("reactions") or ["Not specified"])
                    ),
                    "",
                ]
            )

    output_file.write_text("\n".join(output_lines), encoding="utf-8")
    print(f"\nReport written to {output_file.resolve()}")


if __name__ == "__main__":
    main()

