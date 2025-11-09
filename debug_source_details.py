"""
Debug get_source_details
"""

from chemtools.recommend import UnifiedRecommender

recommender = UnifiedRecommender()

# Get a rule result
rxn = "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"
results = recommender.recommend(rxn, top_k=1, source_types=['rule'])

if results:
    result = results[0]
    print(f"Result ID: {result.id}")
    print(f"Source file: {result.source_file}")
    print()
    
    # Try to load details directly
    print("Attempting to load source details...")
    details = recommender.get_source_details(result.id)
    
    if details:
        print(f"✅ Got details!")
        print(f"Keys: {list(details.keys())}")
        if 'default_rule' in details:
            print(f"Has default_rule: YES")
            if 'conditions' in details['default_rule']:
                print(f"Has conditions: YES")
                conds = details['default_rule']['conditions']
                print(f"Conditions keys: {list(conds.keys())[:5]}...")
            else:
                print("❌ No conditions in default_rule")
        else:
            print("❌ No default_rule")
    else:
        print("❌ get_source_details returned None")
        print(f"Checking if file exists: {result.source_file}")
        from pathlib import Path
        if Path(result.source_file).exists():
            print("✅ File exists")
        else:
            print("❌ File does not exist")
