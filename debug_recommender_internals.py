
print("Script start")
from chemtools.recommend.recommender import HTERecommender, PROJECT_ROOT
from chemtools.featurizers.unified import featurize_reaction
import pandas as pd

try:
    db_path = PROJECT_ROOT / "data" / "protocol_db_v2"
    recommender = HTERecommender(db_path)
    
    # 1. Check if DF is loaded
    if recommender.df is None or recommender.df.empty:
        print("DF is empty!")
        exit()
    
    print(f"DF shape: {recommender.df.shape}")
    
    # 2. Check Source_Group values
    if 'Source_Group' in recommender.df.columns:
        print(f"Source Groups: {recommender.df['Source_Group'].unique()}")
    else:
        print("Source_Group column MISSING")

    # 3. Check for the specific reaction
    # Reactants: O=C1CCCC1 and Brc2ccc(C(C)=O)cc2
    # Standard key generation
    
    # Generate key for the reactants
    reactant_a = "O=C1CCCC1"
    reactant_b = "Brc2ccc(C(C)=O)cc2"
    product = "O=C(C)c3ccc(C4C(CCC4)=O)cc3"
    rxn_smiles = f"{reactant_a}.{reactant_b}>>{product}"
    
    print(f"Featurizing: {rxn_smiles}")
    features = featurize_reaction(rxn_smiles, options={"confirm_coupling_products": True})
    key = features.get("reaction", {}).get("reaction_key")
    print(f"Generated Key: {key}")
    
    # Look for this key in the index
    if key in recommender.indexed_data:
        print("Key found in index!")
        indices = recommender.indexed_data[key]
        print(f"Number of indices: {len(indices)}")
        
        # Check source group of these matches
        subset = recommender.df.iloc[indices]
        print("Subset sources:")
        print(subset['Source_Group'].value_counts())
        
        # Check filtering logic manually
        min_experiments = 2
        
        condition_cols = ['Catalyst', 'Ligand', 'Base', 'Solvent']
        grouped = subset.groupby(condition_cols, dropna=False)
        print(f"Number of groups: {len(grouped)}")
        
        for condition_tuple, group_df in grouped:
            is_proto = (group_df['Source_Group'] == 'protocols').any()
            print(f"Group: {condition_tuple}, Size: {len(group_df)}, IsProto: {is_proto}")

    else:
        print("Key NOT found in index.")
        # Print some keys to see what they look like
        print("First 5 keys in index:")
        print(list(recommender.indexed_data.keys())[:5])
        
        # Check if we can find a partial match
        print("Searching for partial match...")
        for k in recommender.indexed_data.keys():
            if "Br" in k and "acidic-H" in k:
                 print(f"Potential candidate: {k}")


except Exception as e:
    import traceback
    traceback.print_exc()

print("Script end")
