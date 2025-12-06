"""Validate the expanded calculable_features.json"""
import json
from collections import Counter
from pathlib import Path

json_file = Path(__file__).parent.parent / "chemtools" / "taxonomy" / "data" / "calculable_features.json"

with open(json_file, "r", encoding="utf-8") as f:
    data = json.load(f)

print("Valid JSON!")
print(f"\nVersion: {data['version']}")
print(f"Description: {data.get('description', 'N/A')[:80]}...")
print(f"\nFeatures: {len(data['features'])}")
derived_shortcuts = data.get('derived_shortcuts') or []
print(f"Derived shortcuts: {len(derived_shortcuts)}")
print(f"Categories: {len(data['schema_notes']['categories'])}")
functional_groups = data.get('functional_groups')
if functional_groups is None:
    print("Functional groups: (missing section, synthesizing from boolean features)")
else:
    print(f"Functional groups: {len(functional_groups)}")

print("\n=== Feature Categories ===")
for cat in data["schema_notes"]["categories"]:
    print(f"  - {cat}")

print("\n=== Sample Features (first 10) ===")
for i, feat in enumerate(data["features"][:10], 1):
    category = feat.get("category", "N/A")
    print(f"{i}. {feat['token']} ({feat['type']}) - {category}")
    print(f"   Why: {feat.get('why', 'N/A')[:60]}...")

print("\n=== Sample Derived Shortcuts (first 5) ===")
if derived_shortcuts:
    for i, derived in enumerate(derived_shortcuts[:5], 1):
        print(f"{i}. {derived['token']}")
        print(f"   Derive: {derived.get('derive', 'N/A')[:80]}...")
        print(f"   Why: {derived.get('why', 'N/A')[:60]}...")
else:
    print("  (none defined)")

print("\n=== Category Distribution ===")
category_counts = Counter(feat.get("category", "uncategorized") for feat in data["features"])
for cat, count in sorted(category_counts.items(), key=lambda x: -x[1])[:15]:
    print(f"  {cat}: {count} features")

# ---------------------------------------------------------------------------
# Functional group validation
# ---------------------------------------------------------------------------
print("\n=== Functional Group Integrity ===")
if functional_groups is None:
    functional_groups = []
    for feat in data.get("features", []):
        if feat.get("type") != "bool":
            continue
        detect = feat.get("detect")
        if not isinstance(detect, dict):
            continue
        token = feat.get("token", "")
        name = feat.get("name") or (token[:-8] if token.endswith("_present") else token)
        if not name:
            continue
        functional_groups.append({
            "name": name,
            "label": feat.get("label") or feat.get("why") or name.replace("_", " ").title(),
            "detect": detect,
            "text_patterns": feat.get("text_patterns", []),
            "category_tags": feat.get("category_tags", []) or ([feat.get("category")] if feat.get("category") else []),
        })
if not functional_groups:
    raise SystemExit("No functional group definitions available in calculable_features.json")

required_keys = {"name", "label", "detect"}
seen_names = set()
allowed_detect_keys = {"smarts_any"}

for entry in functional_groups:
    missing = required_keys - entry.keys()
    if missing:
        raise SystemExit(f"Functional group entry missing keys {missing}: {entry}")
    name = entry["name"]
    if name in seen_names:
        raise SystemExit(f"Duplicate functional group '{name}'")
    seen_names.add(name)

    detect = entry["detect"]
    if not isinstance(detect, dict) or not detect:
        raise SystemExit(f"Functional group '{name}' has invalid detect block")

    detect_keys = set(detect.keys())
    if not detect_keys <= allowed_detect_keys:
        raise SystemExit(f"Functional group '{name}' uses unsupported detect keys: {detect_keys}")

    smarts_any = detect.get("smarts_any") or []
    if not smarts_any or not all(isinstance(s, str) and s.strip() for s in smarts_any):
        raise SystemExit(f"Functional group '{name}' requires non-empty smarts_any list")

    text_patterns = entry.get("text_patterns", [])
    if text_patterns and not all(isinstance(p, str) for p in text_patterns):
        raise SystemExit(f"Functional group '{name}' has invalid text_patterns entries")

    category_tags = entry.get("category_tags", [])
    if category_tags and not all(isinstance(tag, str) and tag for tag in category_tags):
        raise SystemExit(f"Functional group '{name}' has invalid category_tags entries")

print(f"  {len(functional_groups)} functional group entries validated")

print("\nExpansion successful!")
print(f"  Total coverage: {len(data['features']) + len(derived_shortcuts)} tokens")
