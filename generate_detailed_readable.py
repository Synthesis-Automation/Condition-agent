#!/usr/bin/env python3
"""
Create an enhanced readable format with complete condition details.
Shows: Catalyst / Ligand / Base / Solvent / Temperature / Time
"""

import json


def extract_complete_conditions(rec: dict) -> dict:
    """Extract complete condition details from a recommendation."""
    combo = rec.get('combo', {})
    conditions = rec.get('conditions', {})
    summary = rec.get('summary', {})
    
    result = {
        'catalyst': None,
        'ligand': None,
        'base': None,
        'solvent': None,
        'temperature': None,
        'time': None,
        'precedent_count': 0
    }
    
    # Extract from combo (preferred)
    result['catalyst'] = combo.get('catalyst_name') or combo.get('core_name')
    result['ligand'] = combo.get('ligand_name')
    result['base'] = combo.get('base_name')
    result['solvent'] = combo.get('solvent_name')
    
    # Extract temperature and time
    result['temperature'] = conditions.get('temperature')
    result['time'] = conditions.get('time')
    
    # Get precedent count
    support = summary.get('support', {})
    result['precedent_count'] = support.get('count', 0)
    
    # Fallback to summary if combo is empty
    if not result['catalyst'] and summary:
        catalyst = summary.get('catalyst')
        if isinstance(catalyst, dict):
            result['catalyst'] = catalyst.get('name') or catalyst.get('abbreviation')
        elif catalyst:
            result['catalyst'] = str(catalyst)
    
    if not result['ligand'] and summary:
        ligand = summary.get('ligand')
        if isinstance(ligand, dict):
            result['ligand'] = ligand.get('name') or ligand.get('abbreviation')
        elif ligand:
            result['ligand'] = str(ligand)
    
    if not result['base'] and summary:
        base = summary.get('base')
        if isinstance(base, dict):
            result['base'] = base.get('name') or base.get('abbreviation')
        elif base:
            result['base'] = str(base)
    
    if not result['solvent'] and summary:
        solvent = summary.get('solvent')
        if isinstance(solvent, dict):
            result['solvent'] = solvent.get('name') or solvent.get('abbreviation')
        elif solvent:
            result['solvent'] = str(solvent)
    
    return result


def format_condition_line(cond: dict) -> str:
    """Format a condition dict as a readable line."""
    parts = []
    
    if cond['catalyst']:
        parts.append(f"Catalyst: {cond['catalyst']}")
    if cond['ligand']:
        parts.append(f"Ligand: {cond['ligand']}")
    if cond['base']:
        parts.append(f"Base: {cond['base']}")
    if cond['solvent']:
        parts.append(f"Solvent: {cond['solvent']}")
    if cond['temperature']:
        parts.append(f"Temp: {cond['temperature']}")
    if cond['time']:
        parts.append(f"Time: {cond['time']}")
    
    if not parts:
        parts.append("Conditions not fully specified")
    
    line = " | ".join(parts)
    
    if cond['precedent_count'] > 0:
        line += f" [{cond['precedent_count']} precedent(s)]"
    
    return line


def main():
    # Load results
    with open('test_all_sample_reactions_results.json', 'r') as f:
        data = json.load(f)
    
    detailed_results = data['detailed_results']
    
    # Create enhanced output
    output_file = 'test_results_detailed.txt'
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 120 + "\n")
        f.write(" " * 30 + "CROSS-FAMILY RECOMMENDATION TEST RESULTS\n")
        f.write(" " * 35 + "DETAILED CONDITIONS FORMAT\n")
        f.write("=" * 120 + "\n\n")
        
        f.write(f"Total Reactions: {len(detailed_results)}\n")
        f.write(f"Total Recommendations: {sum(r['result']['num_recommendations'] for r in detailed_results)}\n\n")
        
        # Group by reaction type
        by_type = {}
        for result in detailed_results:
            rtype = result['reaction_type']
            if rtype not in by_type:
                by_type[rtype] = []
            by_type[rtype].append(result)
        
        # Sort types
        sorted_types = sorted(by_type.items(), key=lambda x: len(x[1]), reverse=True)
        
        for reaction_type, results in sorted_types:
            f.write(f"\n{'=' * 120}\n")
            f.write(f"  {reaction_type.upper()} ({len(results)} reactions)\n")
            f.write(f"{'=' * 120}\n\n")
            
            for i, result in enumerate(results, 1):
                smiles = result['smiles']
                description = result['description']
                rec_result = result['result']
                
                f.write(f"[{i}] {description}\n")
                f.write(f"{'─' * 120}\n")
                f.write(f"SMILES: {smiles}\n")
                f.write(f"Detected Family: {rec_result.get('reaction_family', 'Unknown')}\n\n")
                
                recommendations = rec_result.get('recommendations', [])
                
                if recommendations:
                    f.write(f"RECOMMENDED CONDITIONS ({len(recommendations)} options):\n\n")
                    
                    for j, rec in enumerate(recommendations, 1):
                        cond = extract_complete_conditions(rec)
                        condition_line = format_condition_line(cond)
                        
                        f.write(f"  Option {j}: {condition_line}\n")
                    
                    f.write("\n")
                else:
                    f.write("No recommendations available\n\n")
                
                f.write("\n")
        
        f.write("=" * 120 + "\n")
        f.write(" " * 50 + "END OF REPORT\n")
        f.write("=" * 120 + "\n")
    
    print(f"✓ Detailed results saved to: {output_file}")
    
    # Also create a concise version showing just top recommendations
    concise_file = 'test_results_concise.txt'
    with open(concise_file, 'w', encoding='utf-8') as f:
        f.write("CROSS-FAMILY RECOMMENDATION TEST - CONCISE RESULTS\n")
        f.write("=" * 100 + "\n\n")
        
        for result in detailed_results:
            description = result['description']
            smiles = result['smiles']
            recommendations = result['result'].get('recommendations', [])
            
            f.write(f"{result['index']}. {description}\n")
            f.write(f"   {smiles}\n")
            
            if recommendations:
                # Show top recommendation
                cond = extract_complete_conditions(recommendations[0])
                condition_line = format_condition_line(cond)
                f.write(f"   → {condition_line}\n")
            else:
                f.write("   → No recommendations\n")
            
            f.write("\n")
    
    print(f"✓ Concise results saved to: {concise_file}")
    print(f"\nTotal reactions processed: {len(detailed_results)}")


if __name__ == "__main__":
    main()
