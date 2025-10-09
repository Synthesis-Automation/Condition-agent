"""Debug precedent structure."""

from chemtools import recommend

rxn = "C/C=C/Br.C=CB(O)O>>C/C=C/C=C"

result = recommend.recommend_from_reaction(
    reaction=rxn,
    k=50,
    family_override="Suzuki",
    max_variants=1,
    rerank_strategy='rule',
)

pack = result.get('precedent_pack', {})
precs = pack.get('precedents', [])

print("Top 10 Precedents - Core & Catalytic System Info:")
print("=" * 120)

for i, p in enumerate(precs[:10], 1):
    print(f"\n{i}. Reaction ID: {p.get('reaction_id')}")
    print(f"   Core: {p.get('core')}")
    print(f"   Has catalytic_system: {bool(p.get('catalytic_system'))}")
    
    cat_sys = p.get('catalytic_system')
    if cat_sys:
        print(f"   Catalytic system ({len(cat_sys)} items):")
        for cat in cat_sys:
            print(f"     - {cat.get('name')} (CAS: {cat.get('cas')})")
    else:
        print(f"   Catalytic system: None")

# Check chosen core and group
recommendation = result.get('recommendation', {})
chosen_core = recommendation.get('core')

print(f"\n{'=' * 120}")
print(f"Chosen Core: {chosen_core}")

# Find group
group = [p for p in precs if str(p.get('core') or '') == chosen_core]
print(f"Group size: {len(group)}")

if group:
    print(f"\nPrecedents in chosen core group ({chosen_core}):")
    for i, p in enumerate(group[:5], 1):
        print(f"\n  {i}. {p.get('reaction_id')}")
        print(f"     Core: {p.get('core')}")
        print(f"     Has catalytic_system: {bool(p.get('catalytic_system'))}")
        if p.get('catalytic_system'):
            for cat in p.get('catalytic_system', [])[:2]:
                print(f"       - {cat.get('name')}")
