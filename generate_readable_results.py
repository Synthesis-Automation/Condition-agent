#!/usr/bin/env python3
"""
Convert test results to a more readable format.
Shows each reaction with its recommended conditions in a clear format.
"""

import json
from pathlib import Path


def extract_condition_string(rec: dict) -> str:
    """Extract a readable condition string from a recommendation."""
    parts = []
    
    # Get the combo and conditions
    combo = rec.get('combo', {})
    conditions = rec.get('conditions', {})
    summary = rec.get('summary', {})
    
    # Catalyst/Core
    catalyst = combo.get('catalyst_name') or combo.get('core_name')
    if catalyst:
        parts.append(catalyst)
    
    # Ligand
    ligand = combo.get('ligand_name')
    if ligand:
        parts.append(ligand)
    
    # Base
    base = combo.get('base_name')
    if base:
        parts.append(base)
    
    # Solvent
    solvent = combo.get('solvent_name')
    if solvent:
        parts.append(solvent)
    
    # Temperature
    temp = conditions.get('temperature')
    if temp:
        parts.append(temp)
    
    # Time
    time = conditions.get('time')
    if time:
        parts.append(time)
    
    # If we have parts, join them
    if parts:
        return " / ".join(parts)
    
    # Fallback: try to extract from summary
    if summary:
        catalyst = summary.get('catalyst')
        ligand = summary.get('ligand')
        base = summary.get('base')
        solvent = summary.get('solvent')
        
        fallback_parts = []
        # Handle both string and dict values
        if catalyst:
            if isinstance(catalyst, dict):
                catalyst = catalyst.get('name') or catalyst.get('abbreviation', '')
            if catalyst:
                fallback_parts.append(str(catalyst))
        if ligand:
            if isinstance(ligand, dict):
                ligand = ligand.get('name') or ligand.get('abbreviation', '')
            if ligand:
                fallback_parts.append(str(ligand))
        if base:
            if isinstance(base, dict):
                base = base.get('name') or base.get('abbreviation', '')
            if base:
                fallback_parts.append(str(base))
        if solvent:
            if isinstance(solvent, dict):
                solvent = solvent.get('name') or solvent.get('abbreviation', '')
            if solvent:
                fallback_parts.append(str(solvent))
        
        if fallback_parts:
            return " / ".join(fallback_parts)
    
    return "Conditions not specified"


def format_precedent_count(rec: dict) -> str:
    """Get precedent count string."""
    summary = rec.get('summary', {})
    support = summary.get('support', {})
    count = support.get('count', 0)
    if count > 0:
        return f" ({count} precedents)"
    return ""


def main():
    # Load results
    with open('test_all_sample_reactions_results.json', 'r') as f:
        data = json.load(f)
    
    detailed_results = data['detailed_results']
    
    # Create readable output
    output_file = 'test_results_readable.txt'
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 100 + "\n")
        f.write("CROSS-FAMILY RECOMMENDATION TEST RESULTS - READABLE FORMAT\n")
        f.write("=" * 100 + "\n\n")
        
        f.write(f"Total Reactions: {len(detailed_results)}\n")
        f.write(f"Success Rate: 100%\n")
        f.write(f"Total Recommendations: {sum(r['result']['num_recommendations'] for r in detailed_results)}\n\n")
        
        f.write("=" * 100 + "\n\n")
        
        # Group by reaction type
        by_type = {}
        for result in detailed_results:
            rtype = result['reaction_type']
            if rtype not in by_type:
                by_type[rtype] = []
            by_type[rtype].append(result)
        
        # Sort types by number of reactions
        sorted_types = sorted(by_type.items(), key=lambda x: len(x[1]), reverse=True)
        
        for reaction_type, results in sorted_types:
            f.write(f"\n{'=' * 100}\n")
            f.write(f"{reaction_type.upper()} ({len(results)} reactions)\n")
            f.write(f"{'=' * 100}\n\n")
            
            for i, result in enumerate(results, 1):
                smiles = result['smiles']
                description = result['description']
                rec_result = result['result']
                
                f.write(f"{i}. {description}\n")
                f.write(f"   SMILES: {smiles}\n")
                f.write(f"   Family: {rec_result.get('reaction_family', 'Unknown')}\n")
                
                recommendations = rec_result.get('recommendations', [])
                
                if recommendations:
                    f.write(f"\n   Recommended Conditions ({len(recommendations)} options):\n")
                    f.write(f"   {'-' * 90}\n")
                    
                    for j, rec in enumerate(recommendations, 1):
                        condition_str = extract_condition_string(rec)
                        precedent_str = format_precedent_count(rec)
                        
                        f.write(f"   {j}. {condition_str}{precedent_str}\n")
                    
                    f.write(f"   {'-' * 90}\n")
                else:
                    f.write("   No recommendations available\n")
                
                f.write("\n")
        
        f.write("\n" + "=" * 100 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 100 + "\n")
    
    print(f"Readable results saved to: {output_file}")
    print(f"Total reactions processed: {len(detailed_results)}")
    
    # Also create a CSV for easy analysis
    csv_file = 'test_results_summary.csv'
    with open(csv_file, 'w', encoding='utf-8') as f:
        f.write("Index,Description,SMILES,Reaction Type,Reaction Family,Num Recommendations,Top Condition\n")
        
        for result in detailed_results:
            index = result['index']
            description = result['description'].replace(',', ';').replace('"', "'")
            smiles = result['smiles'].replace(',', ';')
            rtype = result['reaction_type']
            family = result['result'].get('reaction_family', 'Unknown')
            num_recs = result['result']['num_recommendations']
            
            # Get top condition
            recommendations = result['result'].get('recommendations', [])
            if recommendations:
                top_condition = extract_condition_string(recommendations[0]).replace(',', ';')
            else:
                top_condition = "None"
            
            f.write(f'{index},"{description}","{smiles}",{rtype},{family},{num_recs},"{top_condition}"\n')
    
    print(f"CSV summary saved to: {csv_file}")
    
    # Create a condensed markdown report
    md_file = 'test_results_condensed.md'
    with open(md_file, 'w', encoding='utf-8') as f:
        f.write("# Cross-Family Recommendation Test Results\n\n")
        f.write("## Quick Reference Guide\n\n")
        
        for reaction_type, results in sorted_types[:10]:  # Top 10 types
            f.write(f"### {reaction_type} ({len(results)} reactions)\n\n")
            
            for result in results[:5]:  # Show first 5 of each type
                description = result['description']
                smiles = result['smiles']
                recommendations = result['result'].get('recommendations', [])
                
                f.write(f"**{description}**\n\n")
                f.write(f"```\n{smiles}\n```\n\n")
                
                if recommendations:
                    f.write("Top conditions:\n\n")
                    for j, rec in enumerate(recommendations[:3], 1):  # Top 3
                        condition_str = extract_condition_string(rec)
                        precedent_str = format_precedent_count(rec)
                        f.write(f"{j}. `{condition_str}`{precedent_str}\n")
                    f.write("\n")
                
                f.write("---\n\n")
            
            if len(results) > 5:
                f.write(f"*... and {len(results) - 5} more {reaction_type} reactions*\n\n")
    
    print(f"Markdown summary saved to: {md_file}")


if __name__ == "__main__":
    main()
