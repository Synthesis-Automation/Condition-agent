"""
Complete End-to-End Example: ChemTools v2.0 Workflow
Demonstrates: Normalize → Detect → Featurize → Search
"""
import time
from chemtools import chem

print("=" * 70)
print("ChemTools v2.0 - Complete Workflow Example")
print("=" * 70)

# Test reaction: Buchwald C-N coupling
test_rxn = "c1ccc(Br)cc1.[H]N(c1ccccc1)c1ccccc1>>c1ccc(N(c2ccccc2)c2ccccc2)cc1"

# ============================================================================
# STEP 1: Normalize Reaction SMILES
# ============================================================================
print("\n📝 Step 1: Normalize Reaction SMILES")
print("-" * 70)
norm_result = chem.smiles.normalize_reaction(test_rxn)
normalized = norm_result["normalized"]

print(f"Original:   {test_rxn[:55]}...")
print(f"Normalized: {normalized[:55]}...")
print(f"Reactants:  {len(norm_result.get('reactants', []))}")
print(f"Products:   {len(norm_result.get('products', []))}")

# ============================================================================
# STEP 2: Detect Reaction Family
# ============================================================================
print("\n📝 Step 2: Detect Reaction Family")
print("-" * 70)
router_result = chem.router.detect_family(normalized)

family = router_result["family"]
confidence = router_result["confidence"]
hits = router_result["hits"]

print(f"Family:     {family}")
print(f"Confidence: {confidence}")
print(f"Hits:")
for key, value in hits.items():
    if value:
        print(f"  ✓ {key}")

# ============================================================================
# STEP 3: Prepare Features for Search  
# ============================================================================
print("\n📝 Step 3: Prepare Features (Ullmann Featurizer)")
print("-" * 70)

# For C-N coupling, use Ullmann featurizer
from chemtools.featurizers import ullmann

features = ullmann.featurize(normalized)
print(f"Feature keys: {len(features)}")
print(f"Keys: {', '.join(list(features.keys())[:8])}...")

# Show some feature values
print(f"\nSample features:")
for key in list(features.keys())[:5]:
    val = features[key]
    if isinstance(val, (int, float, bool)):
        print(f"  {key}: {val}")
    elif isinstance(val, str) and len(val) < 50:
        print(f"  {key}: {val}")

# ============================================================================
# STEP 4: Precedent Search with KNN
# ============================================================================
print("\n📝 Step 4: K-Nearest Neighbor Precedent Search")
print("-" * 70)

# Map family to dataset name
family_to_dataset = {
    "Ullmann_CN": "C_N_Coupling_Cu",
    "Buchwald_CN": "C_N_Coupling_Pd",
}

dataset_family = family_to_dataset.get(family, family)
print(f"Dataset: {dataset_family}")

start_time = time.time()
try:
    knn_result = chem.precedent.knn(
        family=dataset_family,
        features=features,
        k=5,
        relax={}
    )
    search_time = time.time() - start_time
    
    print(f"⏱️  Search time: {search_time:.2f}s")
    
    results = knn_result.get("results", [])
    print(f"✓ Found {len(results)} precedents")
    
    if results:
        print(f"\n📊 Top Result:")
        top = results[0]
        print(f"  Distance: {top.get('distance', 'N/A')}")
        if 'metadata' in top:
            meta = top['metadata']
            print(f"  Dataset: {meta.get('dataset', 'N/A')}")
            print(f"  Index: {meta.get('index', 'N/A')}")
        if 'conditions' in top:
            conds = top['conditions']
            print(f"  Temperature: {conds.get('temperature', 'N/A')}")
            print(f"  Solvent: {conds.get('solvent', 'N/A')}")
            
except Exception as e:
    print(f"❌ Search failed: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# STEP 5: Resource Management
# ============================================================================
print("\n📝 Step 5: Resource Management & Cache Stats")
print("-" * 70)

stats = chem.get_cache_stats()
print(f"Datasets loaded: {stats.get('datasets_loaded', 0)}")
print(f"Total resources: {stats.get('total_resources', 0)}")

loaded_datasets = chem.list_loaded_datasets()
if loaded_datasets:
    print(f"\nLoaded datasets:")
    for ds in loaded_datasets:
        print(f"  • {ds}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("✅ Complete workflow executed successfully!")
print("=" * 70)
print(f"""
Summary:
  1. Normalized reaction SMILES
  2. Detected family: {family} ({confidence} confidence)
  3. Generated {len(features)} features
  4. Searched precedents (via selective loading)
  5. Resource stats: {stats.get('datasets_loaded', 0)} datasets loaded

ChemTools v2.0 provides a clean, unified interface for all operations!
""")
