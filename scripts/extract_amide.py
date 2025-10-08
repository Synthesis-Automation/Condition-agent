"""
Extract amide-specific logic from recommend.py.

This script reads recommend.py and extracts all amide-related functions
into a new families/amide.py module.
"""

import re

# Read the original recommend.py
with open('chemtools/recommend.py', 'r', encoding='utf-8') as f:
    content = f.read()
    lines = content.split('\n')

# Extract specific line ranges (based on analysis)
sections = {
    'analyze_carboxylic_acid': (297, 411),  # 115 lines
    'analyze_amine': (412, 461),  # 50 lines
    'has_free_alcohol': (462, 506),  # 45 lines (moved to substrate_analysis)
    'acid_classification': (507, 543),  # 37 lines
    'amine_classification': (544, 595),  # 52 lines
    'amide_substrate_class': (596, 601),  # 6 lines
    'infer_amide_category': (602, 626),  # 25 lines
    'water_management': (627, 637),  # 11 lines
    'looks_like_aromatic_aniline': (638, 650),  # 13 lines
    'amide_rule_feature_builder': (729, 794),  # 66 lines
}

print("Extracting amide-specific functions from recommend.py...")
print(f"Total sections: {len(sections)}")

extracted = []
for name, (start, end) in sections.items():
    section_lines = lines[start-1:end]  # -1 for 0-indexing
    extracted.append('\n'.join(section_lines))
    print(f"✓ Extracted {name}: lines {start}-{end} ({end-start+1} lines)")

# Combine all sections
full_content = '\n\n'.join(extracted)

print(f"\nTotal extracted lines: {len(full_content.split(chr(10)))}")
print("\nPreview (first 500 chars):")
print(full_content[:500])
print("\n...")
print("\nExtraction complete! Review the output and create amide.py")
