"""
Testing Framework for Simple Condition Recommender
Provides multiple ways to test and validate the recommendation system.

UPDATED: Now uses cleaned dataset with generic reactant naming (Reactant_A, Reactant_B)
"""

import pandas as pd
import json
from pathlib import Path
from simple_condition_recommender import SimpleConditionRecommender
from collections import defaultdict
import random


class RecommenderTester:
    """Test the condition recommender with various scenarios."""
    
    def __init__(self, csv_path: str):
        self.csv_path = csv_path
        self.recommender = SimpleConditionRecommender(csv_path)
        self.df = pd.read_csv(csv_path)
        self.high_performers = self.df[self.df['z-Score'] > 1.0]
    
    def test_known_combinations(self, n_samples: int = 10):
        """
        Test 1: Sample real successful combinations from dataset and see if we recommend them.
        This validates that we can reproduce known good results.
        """
        print("\n" + "=" * 80)
        print(f"TEST 1: KNOWN SUCCESSFUL COMBINATIONS (n={n_samples})")
        print("=" * 80)
        print("Testing if we can recommend conditions that actually worked in the dataset\n")
        
        # Sample random high-performing reactions
        samples = self.high_performers.sample(n=min(n_samples, len(self.high_performers)))
        
        results = {
            'exact_match': 0,
            'category_match': 0,
            'reaction_match': 0,
            'found_in_top3': 0,
            'found_in_top5': 0,
            'not_found': 0
        }
        
        for idx, row in samples.iterrows():
            rxn_type = row['Reaction_Type_Standardized']
            reactant_a = row['Reactant_A']  # Updated: use Reactant_A
            reactant_b = row['Reactant_B']  # Updated: use Reactant_B
            
            # Skip if any are NaN
            if pd.isna(rxn_type) or pd.isna(reactant_a) or pd.isna(reactant_b) or reactant_a == '' or reactant_b == '':
                continue
            
            # Get recommendation
            rec = self.recommender.recommend(rxn_type, reactant_a, reactant_b, top_n=5)
            
            # Check if actual conditions are in recommendations
            actual_cat = row['Catalyst']
            actual_lig = row['Ligand']
            actual_base = row['Base']
            actual_solv = row['Solvent']
            
            found_rank = None
            for i, r in enumerate(rec['recommendations']):
                if (str(r['catalyst']) == str(actual_cat) and
                    str(r['ligand']) == str(actual_lig) and
                    str(r['base']) == str(actual_base) and
                    str(r['solvent']) == str(actual_solv)):
                    found_rank = i + 1
                    break
            
            # Update results
            results[f"{rec['match_level']}_match"] += 1
            
            if found_rank:
                if found_rank <= 3:
                    results['found_in_top3'] += 1
                if found_rank <= 5:
                    results['found_in_top5'] += 1
                print(f"✓ {rxn_type} ({reactant_a} + {reactant_b})")
                print(f"  Actual conditions found at rank {found_rank}")
                print(f"  Match level: {rec['match_level']}")
            else:
                results['not_found'] += 1
                print(f"✗ {rxn_type} ({reactant_a} + {reactant_b})")
                print(f"  Actual conditions NOT in top 5 recommendations")
                print(f"  Match level: {rec['match_level']}, Total precedents: {rec['total_precedents']}")
        
        print(f"\nRESULTS:")
        print(f"  Exact matches: {results['exact_match']}")
        print(f"  Category matches: {results['category_match']}")
        print(f"  Reaction-only matches: {results['reaction_match']}")
        print(f"  Found in top 3: {results['found_in_top3']} ({results['found_in_top3']/n_samples*100:.1f}%)")
        print(f"  Found in top 5: {results['found_in_top5']} ({results['found_in_top5']/n_samples*100:.1f}%)")
        print(f"  Not found: {results['not_found']}")
    
    def test_coverage_by_reaction_type(self):
        """
        Test 2: Check coverage for all major reaction types.
        Shows which reaction types have good/poor coverage.
        """
        print("\n" + "=" * 80)
        print("TEST 2: COVERAGE BY REACTION TYPE")
        print("=" * 80)
        print("How many substrate combinations exist for each reaction type?\n")
        
        # Get unique combinations per reaction type
        coverage = {}
        
        for rxn_type in self.df['Reaction_Type_Standardized'].dropna().unique():
            subset = self.high_performers[self.high_performers['Reaction_Type_Standardized'] == rxn_type]
            
            # Count unique substrate combinations (using Reactant_A and Reactant_B)
            combinations = subset.groupby(['Reactant_A', 'Reactant_B']).size()
            
            coverage[rxn_type] = {
                'total_high_performers': len(subset),
                'unique_combinations': len(combinations),
                'avg_precedents_per_combo': len(subset) / len(combinations) if len(combinations) > 0 else 0,
                'top_combinations': combinations.nlargest(3).to_dict()
            }
        
        # Sort by total high performers
        sorted_coverage = sorted(coverage.items(), key=lambda x: x[1]['total_high_performers'], reverse=True)
        
        for rxn_type, stats in sorted_coverage[:15]:  # Top 15
            print(f"{rxn_type}:")
            print(f"  High performers: {stats['total_high_performers']:,}")
            print(f"  Unique substrate combos: {stats['unique_combinations']}")
            print(f"  Avg precedents per combo: {stats['avg_precedents_per_combo']:.1f}")
            print(f"  Top combinations:")
            for (reactant_a, reactant_b), count in stats['top_combinations'].items():
                print(f"    {reactant_a} + {reactant_b}: {count} precedents")
            print()
    
    def test_interactive(self):
        """
        Test 3: Interactive testing - user can input their own queries.
        """
        print("\n" + "=" * 80)
        print("TEST 3: INTERACTIVE TESTING")
        print("=" * 80)
        print("Enter your reaction details to get recommendations")
        print("Type 'quit' to exit\n")
        
        # Show available options
        print("Available Reaction Types (top 10):")
        rxn_types = self.df['Reaction_Type_Standardized'].value_counts().head(10)
        for i, (rxn, count) in enumerate(rxn_types.items(), 1):
            print(f"  {i}. {rxn} ({count:,} reactions)")
        
        print("\nAvailable Reactant A Types (top 10):")
        reactant_a_types = self.df['Reactant_A'].value_counts().head(10)
        for i, (reactant_a, count) in enumerate(reactant_a_types.items(), 1):
            if pd.notna(reactant_a) and reactant_a != '':
                print(f"  {i}. {reactant_a} ({count:,} uses)")
        
        print("\nAvailable Reactant B Types (top 10):")
        reactant_b_types = self.df['Reactant_B'].value_counts().head(10)
        for i, (reactant_b, count) in enumerate(reactant_b_types.items(), 1):
            if pd.notna(reactant_b) and reactant_b != '':
                print(f"  {i}. {reactant_b} ({count:,} uses)")
        
        while True:
            print("\n" + "-" * 80)
            rxn = input("\nReaction Type (or 'quit'): ").strip()
            if rxn.lower() == 'quit':
                break
            
            reactant_a = input("Reactant A: ").strip()
            reactant_b = input("Reactant B: ").strip()
            
            try:
                result = self.recommender.recommend(rxn, reactant_a, reactant_b, top_n=5)
                print("\n" + json.dumps(result, indent=2))
            except Exception as e:
                print(f"Error: {e}")
    
    def test_edge_cases(self):
        """
        Test 4: Test edge cases and error handling.
        """
        print("\n" + "=" * 80)
        print("TEST 4: EDGE CASES & ERROR HANDLING")
        print("=" * 80)
        
        test_cases = [
            {
                'name': 'Non-existent reaction type',
                'rxn': 'FakeReaction',
                'reactant_a': 'ArBr',
                'reactant_b': 'RNH2'
            },
            {
                'name': 'Non-existent substrate combination',
                'rxn': 'Buchwald-Hartwig-C-N',
                'reactant_a': 'FakeReactant',
                'reactant_b': 'FakeReactant2'
            },
            {
                'name': 'Rare combination (should fall back)',
                'rxn': 'Stille',
                'reactant_a': 'ArBr',
                'reactant_b': 'alkyne'
            },
            {
                'name': 'Empty strings',
                'rxn': '',
                'reactant_a': '',
                'reactant_b': ''
            },
            {
                'name': 'Amide coupling (no electrophile)',
                'rxn': 'Amide-coupling',
                'reactant_a': 'RCO2H or M',
                'reactant_b': 'RNH2'
            }
        ]
        
        for tc in test_cases:
            print(f"\n{tc['name']}:")
            print(f"  Input: {tc['rxn']} + {tc['reactant_a']} + {tc['reactant_b']}")
            try:
                result = self.recommender.recommend(tc['rxn'], tc['reactant_a'], tc['reactant_b'])
                print(f"  Result: {result['match_level']}")
                print(f"  Recommendations: {len(result['recommendations'])}")
                if result['recommendations']:
                    print(f"  Top recommendation: {result['recommendations'][0]['catalyst']}")
            except Exception as e:
                print(f"  Error: {type(e).__name__}: {e}")
    
    def test_recommendation_quality(self, n_samples: int = 20):
        """
        Test 5: Analyze quality of recommendations by z-score.
        Are we recommending high z-score conditions?
        """
        print("\n" + "=" * 80)
        print(f"TEST 5: RECOMMENDATION QUALITY ANALYSIS (n={n_samples})")
        print("=" * 80)
        print("Checking if recommended conditions have high z-scores\n")
        
        # Sample random substrate combinations
        unique_combos = self.high_performers.groupby([
            'Reaction_Type_Standardized',
            'Reactant_A', 
            'Reactant_B'
        ]).size().reset_index()
        
        samples = unique_combos.sample(n=min(n_samples, len(unique_combos)))
        
        quality_metrics = {
            'avg_top1_zscore': [],
            'avg_top3_zscore': [],
            'max_possible_zscore': []
        }
        
        for idx, row in samples.iterrows():
            rxn = row['Reaction_Type_Standardized']
            reactant_a = row['Reactant_A']
            reactant_b = row['Reactant_B']
            
            if pd.isna(rxn) or pd.isna(reactant_a) or pd.isna(reactant_b):
                continue
            
            # Get recommendations
            result = self.recommender.recommend(rxn, reactant_a, reactant_b, top_n=3)
            
            if result['recommendations']:
                # Get z-scores from recommendations
                top1_z = result['recommendations'][0]['evidence']['avg_zscore']
                top3_z = [r['evidence']['avg_zscore'] for r in result['recommendations']]
                
                # Get actual max z-score for this combination
                actual_data = self.high_performers[
                    (self.high_performers['Reaction_Type_Standardized'] == rxn) &
                    (self.high_performers['Reactant_A'] == reactant_a) &
                    (self.high_performers['Reactant_B'] == reactant_b)
                ]
                max_z = actual_data['z-Score'].max() if len(actual_data) > 0 else 0
                
                quality_metrics['avg_top1_zscore'].append(top1_z)
                quality_metrics['avg_top3_zscore'].extend(top3_z)
                quality_metrics['max_possible_zscore'].append(max_z)
                
                print(f"{rxn} ({reactant_a} + {reactant_b}):")
                print(f"  Top recommendation z-score: {top1_z:.2f}")
                print(f"  Max possible z-score: {max_z:.2f}")
                print(f"  Gap: {max_z - top1_z:.2f}")
        
        print(f"\n\nQUALITY METRICS:")
        print(f"  Average top-1 recommendation z-score: {sum(quality_metrics['avg_top1_zscore'])/len(quality_metrics['avg_top1_zscore']):.2f}")
        print(f"  Average top-3 recommendations z-score: {sum(quality_metrics['avg_top3_zscore'])/len(quality_metrics['avg_top3_zscore']):.2f}")
        print(f"  Average max possible z-score: {sum(quality_metrics['max_possible_zscore'])/len(quality_metrics['max_possible_zscore']):.2f}")
        print(f"  Average gap (max - top1): {sum(quality_metrics['max_possible_zscore'])/len(quality_metrics['max_possible_zscore']) - sum(quality_metrics['avg_top1_zscore'])/len(quality_metrics['avg_top1_zscore']):.2f}")
    
    def test_all_quick(self):
        """Run all non-interactive tests quickly."""
        self.test_known_combinations(n_samples=20)
        self.test_coverage_by_reaction_type()
        self.test_edge_cases()
        self.test_recommendation_quality(n_samples=20)
    
    def generate_test_report(self, output_file: str = "test_report.md"):
        """Generate a comprehensive test report."""
        print(f"\n\nGenerating test report to {output_file}...")
        
        with open(output_file, 'w') as f:
            f.write("# Simple Condition Recommender - Test Report\n\n")
            f.write(f"Dataset: {self.csv_path}\n")
            f.write(f"Total reactions: {len(self.df):,}\n")
            f.write(f"High performers (z > 1.0): {len(self.high_performers):,}\n\n")
            
            # Coverage summary
            f.write("## Coverage Summary\n\n")
            rxn_types = self.df['Reaction_Type_Standardized'].nunique()
            reactant_a_types = self.df['Reactant_A'].nunique()
            reactant_b_types = self.df['Reactant_B'].nunique()
            
            f.write(f"- Unique reaction types: {rxn_types}\n")
            f.write(f"- Unique reactant A types: {reactant_a_types}\n")
            f.write(f"- Unique reactant B types: {reactant_b_types}\n")
            f.write(f"- Possible combinations: {rxn_types * reactant_a_types * reactant_b_types:,}\n\n")
            
            # Actual combinations
            actual_combos = self.df.groupby([
                'Reaction_Type_Standardized',
                'Reactant_A',
                'Reactant_B'
            ]).size()
            
            f.write(f"- Actual combinations in dataset: {len(actual_combos):,}\n")
            f.write(f"- Coverage: {len(actual_combos)/(rxn_types*reactant_a_types*reactant_b_types)*100:.2f}%\n\n")
            
            f.write("## Top Reaction Types by Coverage\n\n")
            for rxn_type in self.df['Reaction_Type_Standardized'].value_counts().head(10).index:
                subset = self.df[self.df['Reaction_Type_Standardized'] == rxn_type]
                high_subset = self.high_performers[self.high_performers['Reaction_Type_Standardized'] == rxn_type]
                combos = subset.groupby(['Reactant_A', 'Reactant_B']).size()
                
                f.write(f"### {rxn_type}\n")
                f.write(f"- Total reactions: {len(subset):,}\n")
                f.write(f"- High performers: {len(high_subset):,} ({len(high_subset)/len(subset)*100:.1f}%)\n")
                f.write(f"- Unique substrate combinations: {len(combos)}\n\n")
        
        print(f"✓ Test report saved to {output_file}")


def main():
    """Main testing interface."""
    import sys
    
    csv_path = Path(__file__).parent / "z-Score Peaks CLEANED.csv"
    tester = RecommenderTester(csv_path)
    
    print("=" * 80)
    print("SIMPLE CONDITION RECOMMENDER - TESTING FRAMEWORK")
    print("=" * 80)
    print("\nAvailable tests:")
    print("  1. Test known combinations (validate we can reproduce results)")
    print("  2. Coverage by reaction type")
    print("  3. Interactive testing (manual queries)")
    print("  4. Edge cases & error handling")
    print("  5. Recommendation quality analysis")
    print("  6. Run all tests (non-interactive)")
    print("  7. Generate test report")
    print("  0. Exit")
    
    while True:
        choice = input("\nSelect test (0-7): ").strip()
        
        if choice == '0':
            break
        elif choice == '1':
            n = input("Number of samples (default 20): ").strip()
            n = int(n) if n else 20
            tester.test_known_combinations(n_samples=n)
        elif choice == '2':
            tester.test_coverage_by_reaction_type()
        elif choice == '3':
            tester.test_interactive()
        elif choice == '4':
            tester.test_edge_cases()
        elif choice == '5':
            n = input("Number of samples (default 20): ").strip()
            n = int(n) if n else 20
            tester.test_recommendation_quality(n_samples=n)
        elif choice == '6':
            tester.test_all_quick()
        elif choice == '7':
            output = input("Output file (default test_report.md): ").strip()
            output = output if output else "test_report.md"
            tester.generate_test_report(output)
        else:
            print("Invalid choice")


if __name__ == "__main__":
    main()
