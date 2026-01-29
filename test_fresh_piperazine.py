#!/usr/bin/env python
"""Fresh test of piperazine C-N coupling detection."""

if __name__ == "__main__":
    from chemtools.featurizers.formatters.reaction import featurize_reaction
    
    rxn = "CN1CCNCC1.Nc1cc(Cl)c(Cl)cc1[N+](=O)[O-]>>CN1CCN(c2cc(N)c([N+](=O)[O-])cc2Cl)CC1"
    result = featurize_reaction(rxn)
    
    print("=" * 70)
    print("REACTION:", rxn)
    print("=" * 70)
    print(f"Detected as: {result['reaction_type']}")
    print(f"Confidence: {result.get('confidence')}")
    print()
    print(f"Primary motif IDs: {result['aggregates']['primary_motif_ids']}")
    print()
    print(f"Reacted: {result['aggregates']['reacted_motifs']}")
    print(f"Formed: {result['aggregates']['formed_motifs']}")
    print(f"Spectators: {result['aggregates']['spectator_motifs']}")
    print()
    print(f"✓ RCH2-NHR in primary? {'RCH2-NHR' in result['aggregates']['primary_motif_ids']}")
    print(f"✓ RCH2-NHR in reacted? {'RCH2-NHR' in result['aggregates']['reacted_motifs']}")
