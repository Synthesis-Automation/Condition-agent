"""
Demo: Summarizing All Conditions for Copper-Catalyzed C-N Coupling
===================================================================

This script demonstrates how to use HTE tools to summarize all conditions
for a given reactant pair (e.g., copper-catalyzed C-N coupling).

Features:
1. Get all condition recommendations for a reactant pair
2. Export filtered experiments to CSV
3. Analyze condition distributions
4. Generate summary statistics
"""

import pandas as pd
from chemtools.HTE import HTERecommender
from chemtools.HTE.analytics import HTEAnalytics


def demo_condition_recommendations(reactant_a: str, reactant_b: str, catalyst: str = "Cu"):
    """
    Get condition recommendations for a specific reactant pair with catalyst filter.
    
    Args:
        reactant_a: SMILES of first reactant
        reactant_b: SMILES of second reactant
        catalyst: Catalyst metal filter (default: "Cu")
    """
    print("\n" + "="*80)
    print("1. CONDITION RECOMMENDATIONS")
    print("="*80)
    
    recommender = HTERecommender()
    result = recommender.recommend(
        reactant_a_smiles=reactant_a,
        reactant_b_smiles=reactant_b,
        catalyst_filter=catalyst,
        top_k=20,
        min_experiments=1  # Get all conditions, even single experiments
    )
    
    print(f"\nReactant A: {reactant_a} → Type: {result.reactant_a_type}")
    print(f"Reactant B: {reactant_b} → Type: {result.reactant_b_type}")
    print(f"Predicted Reaction: {result.predicted_reaction_type} "
          f"(confidence: {result.reaction_type_confidence*100:.1f}%)")
    print(f"\n📊 Found {result.total_matching_experiments} matching experiments")
    print(f"   ({result.database_coverage:.2f}% of database)")
    print(f"\n🏆 Top {len(result.recommendations)} conditions (ranked by Z-Score):\n")
    
    # Create summary table
    summary_data = []
    for i, rec in enumerate(result.recommendations, 1):
        summary_data.append({
            'Rank': i,
            'Z-Score': f"{rec.avg_z_score:.2f}",
            'Catalyst': rec.catalyst,
            'Ligand': rec.ligand,
            'Base': rec.base,
            'Solvent': rec.solvent,
            'Experiments': rec.num_experiments,
            'Avg_Yield_%': f"{rec.avg_yield:.1f}",
            'Confidence': f"{rec.confidence_score:.1f}"
        })
    
    df = pd.DataFrame(summary_data)
    print(df.to_string(index=False))
    
    return result


def demo_export_filtered_data(reactant_a_type: str, reactant_b_type: str, 
                               catalyst: str = "Cu", output_file: str = None):
    """
    Export all experiments matching the reactant types and catalyst filter.
    
    Args:
        reactant_a_type: Reactant A type (e.g., "ArBr")
        reactant_b_type: Reactant B type (e.g., "ArNH2")
        catalyst: Catalyst metal filter
        output_file: Output CSV filename
    """
    print("\n" + "="*80)
    print("2. EXPORT FILTERED EXPERIMENTS")
    print("="*80)
    
    analytics = HTEAnalytics()
    
    if output_file is None:
        output_file = f"{catalyst.lower()}_{reactant_a_type}_{reactant_b_type}_conditions.csv"
    
    num_exported = analytics.export_subset(
        output_path=output_file,
        catalyst_filter=catalyst,
        reactant_a_type=reactant_a_type,
        reactant_b_type=reactant_b_type
    )
    
    print(f"\n✅ Exported {num_exported} experiments to: {output_file}")
    
    # Show sample of exported data
    df = pd.read_csv(output_file)
    print(f"\nColumns in exported file:")
    print(f"  {', '.join(df.columns.tolist())}")
    
    return output_file


def demo_condition_analytics(reactant_a_type: str, reactant_b_type: str, catalyst: str = "Cu"):
    """
    Analyze condition distributions for a reactant pair.
    
    Args:
        reactant_a_type: Reactant A type
        reactant_b_type: Reactant B type
        catalyst: Catalyst metal filter
    """
    print("\n" + "="*80)
    print("3. CONDITION ANALYTICS")
    print("="*80)
    
    # Load data
    df = pd.read_csv("data/HTE_db/HTE_0.csv")
    
    # Filter by reactant types and catalyst
    subset = df[
        (df['Reactant_A_Type'] == reactant_a_type) & 
        (df['Reactant_B_Type'] == reactant_b_type)
    ]
    
    if catalyst:
        subset = subset[subset['Catalyst'].str.contains(catalyst, case=False, na=False)]
    
    print(f"\nFiltered to {len(subset)} experiments")
    print(f"  Reactant A Type: {reactant_a_type}")
    print(f"  Reactant B Type: {reactant_b_type}")
    print(f"  Catalyst Filter: {catalyst}")
    
    # Unique conditions
    print(f"\n📊 Condition Component Analysis:")
    print(f"  Unique Catalysts: {subset['Catalyst'].nunique()}")
    print(f"  Unique Ligands: {subset['Ligand'].nunique()}")
    print(f"  Unique Bases: {subset['Base'].nunique()}")
    print(f"  Unique Solvents: {subset['Solvent'].nunique()}")
    
    # Top components
    print(f"\n🔝 Top 5 Most Common Components:")
    
    print("\n  Catalysts:")
    for cat, count in subset['Catalyst'].value_counts().head(5).items():
        print(f"    {cat}: {count} experiments")
    
    print("\n  Ligands:")
    for lig, count in subset['Ligand'].value_counts().head(5).items():
        print(f"    {lig}: {count} experiments")
    
    print("\n  Bases:")
    for base, count in subset['Base'].value_counts().head(5).items():
        print(f"    {base}: {count} experiments")
    
    print("\n  Solvents:")
    for solv, count in subset['Solvent'].value_counts().head(5).items():
        print(f"    {solv}: {count} experiments")
    
    # Yield statistics
    print(f"\n📈 Yield Statistics:")
    print(f"  Average Yield: {subset['AREA_TOTAL_REDUCED'].mean():.1f}%")
    print(f"  Median Yield: {subset['AREA_TOTAL_REDUCED'].median():.1f}%")
    print(f"  Max Yield: {subset['AREA_TOTAL_REDUCED'].max():.1f}%")
    print(f"  Min Yield: {subset['AREA_TOTAL_REDUCED'].min():.1f}%")
    
    # Success rate (yield > 50%)
    success_rate = (subset['AREA_TOTAL_REDUCED'] > 50).mean() * 100
    print(f"  Success Rate (>50% yield): {success_rate:.1f}%")
    
    return subset


def demo_full_condition_matrix(export_file: str):
    """
    Create a comprehensive summary of all unique condition combinations.
    
    Args:
        export_file: CSV file from export step
    """
    print("\n" + "="*80)
    print("4. FULL CONDITION MATRIX")
    print("="*80)
    
    df = pd.read_csv(export_file)
    
    # Group by all condition components
    condition_groups = df.groupby(['Catalyst', 'Ligand', 'Base', 'Solvent']).agg({
        'AREA_TOTAL_REDUCED': ['count', 'mean', 'median', 'max', 'min'],
        'z-Score': ['mean', 'max', 'min']
    }).reset_index()
    
    # Flatten column names
    condition_groups.columns = [
        'Catalyst', 'Ligand', 'Base', 'Solvent',
        'Num_Experiments', 'Avg_Yield', 'Median_Yield', 'Max_Yield', 'Min_Yield',
        'Avg_ZScore', 'Max_ZScore', 'Min_ZScore'
    ]
    
    # Sort by average yield
    condition_groups = condition_groups.sort_values('Avg_Yield', ascending=False)
    
    print(f"\nFound {len(condition_groups)} unique condition combinations")
    print(f"\nTop 10 conditions by average yield:\n")
    
    # Display top 10
    display_cols = ['Catalyst', 'Ligand', 'Base', 'Solvent', 
                    'Num_Experiments', 'Avg_Yield', 'Max_Yield']
    print(condition_groups[display_cols].head(10).to_string(index=False))
    
    # Save full matrix
    matrix_file = export_file.replace('.csv', '_matrix.csv')
    condition_groups.to_csv(matrix_file, index=False)
    print(f"\n✅ Saved full condition matrix to: {matrix_file}")
    
    return condition_groups


def main():
    """Run all demos"""
    print("\n" + "="*80)
    print("HTE CONDITION SUMMARIZATION DEMO")
    print("Copper-Catalyzed C-N Coupling: ArBr + ArNH2")
    print("="*80)
    
    # Example reactants
    reactant_a = "c1ccc(Br)cc1"  # Aryl bromide
    reactant_b = "c1ccc(N)cc1"   # Aryl amine
    catalyst = "Cu"
    
    # 1. Get condition recommendations
    result = demo_condition_recommendations(reactant_a, reactant_b, catalyst)
    
    # 2. Export all matching experiments
    export_file = demo_export_filtered_data(
        reactant_a_type=result.reactant_a_type,
        reactant_b_type=result.reactant_b_type,
        catalyst=catalyst,
        output_file="demo_cu_cn_conditions.csv"
    )
    
    # 3. Analyze condition distributions
    demo_condition_analytics(
        reactant_a_type=result.reactant_a_type,
        reactant_b_type=result.reactant_b_type,
        catalyst=catalyst
    )
    
    # 4. Create full condition matrix
    demo_full_condition_matrix(export_file)
    
    print("\n" + "="*80)
    print("✅ DEMO COMPLETE")
    print("="*80)
    print("\nGenerated files:")
    print("  1. demo_cu_cn_conditions.csv - All matching experiments")
    print("  2. demo_cu_cn_conditions_matrix.csv - Condition matrix with statistics")
    print("\nYou can now:")
    print("  - Analyze the exported data in Excel/Python")
    print("  - Use the condition matrix to select optimal conditions")
    print("  - Filter by yield thresholds or specific components")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
