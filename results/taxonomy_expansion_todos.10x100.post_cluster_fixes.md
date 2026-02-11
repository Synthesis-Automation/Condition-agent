# Taxonomy Expansion TODOs

Generated from unresolved discovery clusters.

## Summary

- Input clusters: 100
- TODO entries: 100
- Expand existing reaction type: 3
- Review constraints/product patterns: 36
- Propose new reaction family: 61

## Priority 1: Cluster none -> none || events:none

- Cluster count: 22
- Priority score: 29
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> none`
- Event signature: `none`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: empty_reaction_key:7, low_reaction_key_quality:22, missing_formed_motifs:22, unclassified_motif:22, unknown_reaction_type:7

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_001_none_reacted_formed",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": []
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 2: Cluster Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Alkyl-N(R)CO2R|Ar-Alkyl|RCH2-CO2R|RCH2-COR|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component

- Cluster count: 9
- Priority score: 11
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Alkyl-N(R)CO2R|Ar-Alkyl|RCH2-CO2R|RCH2-COR|R_acidic-H`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Hydrogenation` (score=0.325)
- Reasons: mapping_warning:9, motif_outside_reaction_taxonomy:9

### Gap Analysis

- Missing reacted-slot motifs: Ar-C=N
- Missing formed-slot motifs: Alkyl-N(R)CO2R, R_acidic-H
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-C=N, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_002_c_c_bond_formation_intermolecular_or_multi_component_alkenyl_alkenyl_alkenyl_or_alkyl_n_r_co2r_ar_alkyl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-Alkenyl",
      "Alkenyl-OR",
      "Alkyl-Si*",
      "Ar-C=N"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-N(R)CO2R",
      "Ar-Alkyl",
      "RCH2-CO2R",
      "RCH2-COR",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 3: Cluster Ar-C=N|Ar-NR2|Bn-OH -> Ar-Ar|Ar-AromN|Bn-Br || events:intermolecular_or_multi_component

- Cluster count: 9
- Priority score: 11
- Action track: `propose_new_reaction_family`
- Motif signature: `Ar-C=N|Ar-NR2|Bn-OH -> Ar-Ar|Ar-AromN|Bn-Br`
- Event signature: `intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3333)
- Reasons: mapping_warning:9, motif_outside_reaction_taxonomy:9

### Gap Analysis

- Missing reacted-slot motifs: Ar-C=N, Ar-NR2
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Ar-C=N, Ar-NR2
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_003_intermolecular_or_multi_component_ar_c_n_ar_nr2_ar_ar_ar_aromn",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Ar-C=N",
      "Ar-NR2",
      "Bn-OH"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-Ar",
      "Ar-AromN",
      "Bn-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 4: Cluster Ar-Alkyl -> Ar-Br|Bn-Br || events:intramolecular_likely

- Cluster count: 8
- Priority score: 11
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkyl -> Ar-Br|Bn-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.65)
- Reasons: unknown_reaction_type:8

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Oxidation_alcohol",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 5: Cluster none -> R2CH-Br || events:intramolecular_likely

- Cluster count: 7
- Priority score: 10
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.35)
- Reasons: unknown_reaction_type:7

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_006_intramolecular_likely_reacted_r2ch_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "R2CH-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 6: Cluster Alkyl-H|Ar-COR|Ar-H -> Ar-C=N || events:c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 7
- Priority score: 9
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-H|Ar-COR|Ar-H -> Ar-C=N`
- Event signature: `c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Ircatalyzed_CH_borylation` (score=0.2167)
- Reasons: mapping_warning:4, motif_outside_reaction_taxonomy:7

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-H
- Missing formed-slot motifs: Ar-C=N
- Existing compound IDs missing in reaction slots: Alkyl-H, Ar-C=N
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_005_c_n_bond_formation_intermolecular_or_multi_component_alkyl_h_ar_cor_ar_c_n",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-H",
      "Ar-COR",
      "Ar-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-C=N"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 7: Cluster Alkyl-N(R)CO2R|CH3-OH -> R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 6
- Priority score: 8
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|CH3-OH -> R2CH-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.675)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:6

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 8: Cluster Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|Ar-F|Ar-OH|R2CH-CO2H|R2CH-CO2R|RCH2-CO2R -> Alkyl-NHCOR|R2CH-CONH2|R2CH-CONHR|R2CH-NH2|RCH2-CO2H|RCH2-CONHR|R_acidic-H || events:amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 6
- Priority score: 8
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|Ar-F|Ar-OH|R2CH-CO2H|R2CH-CO2R|RCH2-CO2R -> Alkyl-NHCOR|R2CH-CONH2|R2CH-CONHR|R2CH-NH2|RCH2-CO2H|RCH2-CONHR|R_acidic-H`
- Event signature: `amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.2812)
- Reasons: mapping_warning:6, motif_outside_reaction_taxonomy:6

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-Ar
- Missing formed-slot motifs: R_acidic-H
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-Ar, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_008_amidation_like_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_alkyl_n_r_co2r_ar_alkyl_alkyl_nhcor_r2ch_conh2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-Alkyl",
      "Ar-Ar",
      "Ar-F",
      "Ar-OH",
      "R2CH-CO2H",
      "R2CH-CO2R",
      "RCH2-CO2R"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-NHCOR",
      "R2CH-CONH2",
      "R2CH-CONHR",
      "R2CH-NH2",
      "RCH2-CO2H",
      "RCH2-CONHR",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 9: Cluster none -> Ar-Br|R2CH-Br || events:intramolecular_likely

- Cluster count: 4
- Priority score: 7
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Ar-Br|R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.175)
- Reasons: unknown_reaction_type:4

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_010_intramolecular_likely_reacted_ar_br_r2ch_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-Br",
      "R2CH-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 10: Cluster Alkyl-H|Ar-CHO|RCH2-OH -> Ar-Alkyl|RCH2-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 4
- Priority score: 6
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-H|Ar-CHO|RCH2-OH -> Ar-Alkyl|RCH2-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3917)
- Reasons: mapping_warning:4, motif_outside_reaction_taxonomy:4

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-H
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_009_c_o_bond_formation_intermolecular_or_multi_component_alkyl_h_ar_cho_ar_alkyl_rch2_or",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-H",
      "Ar-CHO",
      "RCH2-OH"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-Alkyl",
      "RCH2-OR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 11: Cluster Alkenyl-OR|Alkyl-F|Alkyl-Si*|Ar-CHO -> Alkyl-CF3|Ar-Alkyl|R2CH-OH|R3C-CO2R|R3C-OH || events:c_c_bond_formation|intermolecular_or_multi_component

- Cluster count: 4
- Priority score: 6
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-OR|Alkyl-F|Alkyl-Si*|Ar-CHO -> Alkyl-CF3|Ar-Alkyl|R2CH-OH|R3C-CO2R|R3C-OH`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Reduction_carbonyl` (score=0.3725)
- Reasons: mapping_warning:4, motif_outside_reaction_taxonomy:4

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkyl-CF3
- Existing compound IDs missing in reaction slots: Alkyl-CF3
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_011_c_c_bond_formation_intermolecular_or_multi_component_alkenyl_or_alkyl_f_alkyl_cf3_ar_alkyl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-OR",
      "Alkyl-F",
      "Alkyl-Si*",
      "Ar-CHO"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-CF3",
      "Ar-Alkyl",
      "R2CH-OH",
      "R3C-CO2R",
      "R3C-OH"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 12: Cluster Ar-C=N|Ar-NR2|RCH2-OH -> Ar-Ar|Ar-AromN|RCH2-Br || events:intermolecular_or_multi_component

- Cluster count: 4
- Priority score: 6
- Action track: `propose_new_reaction_family`
- Motif signature: `Ar-C=N|Ar-NR2|RCH2-OH -> Ar-Ar|Ar-AromN|RCH2-Br`
- Event signature: `intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3333)
- Reasons: mapping_warning:4, motif_outside_reaction_taxonomy:4

### Gap Analysis

- Missing reacted-slot motifs: Ar-C=N, Ar-NR2
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Ar-C=N, Ar-NR2
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_012_intermolecular_or_multi_component_ar_c_n_ar_nr2_ar_ar_ar_aromn",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Ar-C=N",
      "Ar-NR2",
      "RCH2-OH"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-Ar",
      "Ar-AromN",
      "RCH2-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 13: Cluster Alkyl-Hydrazine|RCH2-CHO|R_acidic-H -> CH3-NR2 || events:c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 3
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-Hydrazine|RCH2-CHO|R_acidic-H -> CH3-NR2`
- Event signature: `c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.35)
- Reasons: mapping_warning:3, motif_outside_reaction_taxonomy:3

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-Hydrazine
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-Hydrazine
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_013_c_n_bond_formation_intermolecular_or_multi_component_alkyl_hydrazine_rch2_cho_ch3_nr2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-Hydrazine",
      "RCH2-CHO",
      "R_acidic-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "CH3-NR2"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 14: Cluster Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-CONHR -> Alkenyl-COR|Alkenyl-NR2|Ar-Hydrazide_NR|R2CH-NR2|R_acidic-H || events:c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

- Cluster count: 3
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-CONHR -> Alkenyl-COR|Alkenyl-NR2|Ar-Hydrazide_NR|R2CH-NR2|R_acidic-H`
- Event signature: `c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation`
- Top taxonomy candidate: `Hydrogenation` (score=0.325)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:3

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-NR2, Ar-Hydrazide_NR, R_acidic-H
- Existing compound IDs missing in reaction slots: Alkenyl-NR2, Ar-Hydrazide_NR, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_014_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_alkenyl_alkenyl_alkenyl_or_alkenyl_cor_alkenyl_nr2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-Alkenyl",
      "Alkenyl-OR",
      "Alkyl-Si*",
      "Ar-CONHR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-COR",
      "Alkenyl-NR2",
      "Ar-Hydrazide_NR",
      "R2CH-NR2",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 15: Cluster Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|R2CH-CO2H|RCH2-Br|RCH2-OH|R_acidic-H -> Alkyl-NRCOR|R2CH-CONR2|RCH2-OR || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 3
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|R2CH-CO2H|RCH2-Br|RCH2-OH|R_acidic-H -> Alkyl-NRCOR|R2CH-CONR2|RCH2-OR`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.419)
- Reasons: mapping_warning:3, motif_outside_reaction_taxonomy:3

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_015_amidation_like_c_n_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkyl_co2r_alkyl_n_r_co2r_alkyl_nrcor_r2ch_conr2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-CO2R",
      "Alkyl-N(R)CO2R",
      "Ar-Alkyl",
      "R2CH-CO2H",
      "RCH2-Br",
      "RCH2-OH",
      "R_acidic-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-NRCOR",
      "R2CH-CONR2",
      "RCH2-OR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 16: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 3
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.435)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:3

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_016_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_co2r",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Inorganic",
      "RCH2-I"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 17: Cluster Alkyl-NHCOR|Alkyl-NRCOR|Ar-Alkyl|Ar-CONHR|Bn-LactamN|R2CH-CONHR|R2CH-CONR2|RCH2-CO2H -> none || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 3
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-NHCOR|Alkyl-NRCOR|Ar-Alkyl|Ar-CONHR|Bn-LactamN|R2CH-CONHR|R2CH-CONR2|RCH2-CO2H -> none`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Amide_formation` (score=0.3187)
- Reasons: mapping_warning:3, missing_formed_motifs:3, motif_outside_reaction_taxonomy:3

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-NRCOR, Bn-LactamN, R2CH-CONHR, R2CH-CONR2
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-NRCOR, Bn-LactamN, R2CH-CONHR, R2CH-CONR2
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_017_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkyl_nhcor_alkyl_nrcor_formed",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-NHCOR",
      "Alkyl-NRCOR",
      "Ar-Alkyl",
      "Ar-CONHR",
      "Bn-LactamN",
      "R2CH-CONHR",
      "R2CH-CONR2",
      "RCH2-CO2H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": []
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 18: Cluster Ar-Alkyl|Inorganic|RCH2-CO2H|Unclassified-Reactant -> RCH2-CO2R || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkyl|Inorganic|RCH2-CO2H|Unclassified-Reactant -> RCH2-CO2R`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Esterification` (score=0.5125)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1, unclassified_motif:1

### Gap Analysis

- Missing reacted-slot motifs: Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Esterification",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 19: Cluster Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Alkyl-N(R)CO2R|RCH2-CO2R|RCH2-COR|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component

- Cluster count: 2
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Alkyl-N(R)CO2R|RCH2-CO2R|RCH2-COR|R_acidic-H`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Hydrogenation` (score=0.325)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: Ar-C=N
- Missing formed-slot motifs: Alkyl-N(R)CO2R, R_acidic-H
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-C=N, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_018_c_c_bond_formation_intermolecular_or_multi_component_alkenyl_alkenyl_alkenyl_or_alkyl_n_r_co2r_rch2_co2r",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-Alkenyl",
      "Alkenyl-OR",
      "Alkyl-Si*",
      "Ar-C=N"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-N(R)CO2R",
      "RCH2-CO2R",
      "RCH2-COR",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 20: Cluster Bn-COR -> Ar-Alkyl|Inorganic|R2CH-Br || events:intramolecular_likely

- Cluster count: 2
- Priority score: 4
- Action track: `expand_existing_reaction_type`
- Motif signature: `Bn-COR -> Ar-Alkyl|Inorganic|R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Reduction_carbonyl` (score=0.7667)
- Reasons: motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Inorganic
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Reduction_carbonyl",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Inorganic"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 21: Cluster Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R|RCH2-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|R2CH-NH2|RCH2-CHO|RCH2-CONHR || events:amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 2
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R|RCH2-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|R2CH-NH2|RCH2-CHO|RCH2-CONHR`
- Event signature: `amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.5062)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-Ar",
      "R2CH-OR"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 22: Cluster Alkyl-N(R)CO2R|Alkyl-NHCOR|Ar-Alkyl|Ar-OR|Ar-SO2NHR|R2CH-CONHR|R3C-OR|RCH2-CO2R -> R2CH-CO2H|R2CH-NH2|R2CH-OH|RCH2-CO2H|RCH2-OH || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 2
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-N(R)CO2R|Alkyl-NHCOR|Ar-Alkyl|Ar-OR|Ar-SO2NHR|R2CH-CONHR|R3C-OR|RCH2-CO2R -> R2CH-CO2H|R2CH-NH2|R2CH-OH|RCH2-CO2H|RCH2-OH`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.2325)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-OR, R2CH-CONHR, R3C-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-OR, R2CH-CONHR, R3C-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_021_c_o_bond_formation_intermolecular_or_multi_component_alkyl_n_r_co2r_alkyl_nhcor_r2ch_co2h_r2ch_nh2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Alkyl-NHCOR",
      "Ar-Alkyl",
      "Ar-OR",
      "Ar-SO2NHR",
      "R2CH-CONHR",
      "R3C-OR",
      "RCH2-CO2R"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "R2CH-CO2H",
      "R2CH-NH2",
      "R2CH-OH",
      "RCH2-CO2H",
      "RCH2-OH"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 23: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Bn-Br|Inorganic -> Alkenyl-CO2R|Ar-Alkyl|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 2
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Bn-Br|Inorganic -> Alkenyl-CO2R|Ar-Alkyl|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Esterification` (score=0.4)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_022_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_alkyl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Bn-Br",
      "Inorganic"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Ar-Alkyl",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 24: Cluster Alkenyl-CO2H|Alkenyl-LactamN|Alkenyl-SR|Alkyl-N(R)CO2R|Alkyl-NRCOR|Ar-CO2H|Ar-NHCOR|Inorganic -> none || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 2
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkenyl-LactamN|Alkenyl-SR|Alkyl-N(R)CO2R|Alkyl-NRCOR|Ar-CO2H|Ar-NHCOR|Inorganic -> none`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Amide_formation` (score=0.3187)
- Reasons: mapping_warning:1, missing_formed_motifs:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Alkyl-NRCOR, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Alkyl-NRCOR
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_023_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkenyl_co2h_alkenyl_lactamn_formed",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkenyl-LactamN",
      "Alkenyl-SR",
      "Alkyl-N(R)CO2R",
      "Alkyl-NRCOR",
      "Ar-CO2H",
      "Ar-NHCOR",
      "Inorganic"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": []
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 25: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Br -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 2
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Br -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Esterification` (score=0.435)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_024_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_co2r",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Inorganic",
      "RCH2-Br"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 26: Cluster Alkyl-N(R)CO2R|Ar-CO2H|Ar-COCl|Ar-OH|Ar-SO2Cl|CH3-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|Ar-CONHR|Ar-CONR2|Ar-OSO2R || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 2
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|Ar-CO2H|Ar-COCl|Ar-OH|Ar-SO2Cl|CH3-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|Ar-CONHR|Ar-CONR2|Ar-OSO2R`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Amide_formation` (score=0.605)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-SO2Cl
- Missing formed-slot motifs: Ar-OSO2R
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-OSO2R, Ar-SO2Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-SO2Cl"
    ]
  },
  "product_slot_additions_todo": {
    "product": [
      "Ar-OSO2R"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 27: Cluster RCH2-CO2R -> Ar-COR|HeteroAr-Br|Indole|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

- Cluster count: 2
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `RCH2-CO2R -> Ar-COR|HeteroAr-Br|Indole|R_acidic-H`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation`
- Top taxonomy candidate: `Claisen_condensation` (score=0.65)
- Reasons: motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Indole, R_acidic-H
- Existing compound IDs missing in reaction slots: Indole, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Claisen_condensation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Indole",
      "R_acidic-H"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 28: Cluster Alkyl-CO2R -> Ar-COR|HeteroAr-Cl|Indole|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

- Cluster count: 2
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-CO2R -> Ar-COR|HeteroAr-Cl|Indole|R_acidic-H`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation`
- Top taxonomy candidate: `Claisen_condensation` (score=0.65)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Indole, R_acidic-H
- Existing compound IDs missing in reaction slots: Indole, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Claisen_condensation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Indole",
      "R_acidic-H"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 29: Cluster Ar-CO2R|AromN-H -> Pyrazole || events:c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 2
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Ar-CO2R|AromN-H -> Pyrazole`
- Event signature: `c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Ircatalyzed_CH_borylation` (score=0.325)
- Reasons: mapping_warning:2, motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Pyrazole
- Existing compound IDs missing in reaction slots: Pyrazole
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_028_c_n_bond_formation_intermolecular_or_multi_component_ar_co2r_aromn_h_pyrazole",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Ar-CO2R",
      "AromN-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Pyrazole"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 30: Cluster R3C-OR -> Alkyl-COR|R_acidic-H || events:none

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `R3C-OR -> Alkyl-COR|R_acidic-H`
- Event signature: `none`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.175)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: R3C-OR
- Missing formed-slot motifs: R_acidic-H
- Existing compound IDs missing in reaction slots: R3C-OR, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_029_none_r3c_or_alkyl_cor_r_acidic_h",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "R3C-OR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-COR",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 31: Cluster none -> Acyl-Cl || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Acyl-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Acyl_Halides_formation` (score=0.35)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_035_intramolecular_likely_reacted_acyl_cl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Acyl-Cl"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 32: Cluster Alkyl-Si* -> RCH2-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-Si* -> RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Hiyama` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Hiyama",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 33: Cluster Ar-Alkenyl -> Ar-Alkyl|Ar-Br|R2CH-Br|R3C-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkenyl -> Ar-Alkyl|Ar-Br|R2CH-Br|R3C-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Hydrogenation` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Hydrogenation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 34: Cluster none -> Ar-Br|R2CH-Br|RCH2-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Ar-Br|R2CH-Br|RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.2333)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_040_intramolecular_likely_reacted_ar_br_r2ch_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-Br",
      "R2CH-Br",
      "RCH2-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 35: Cluster none -> Alkenyl-COR|R2CH-Br|RCH2-Br|R_acidic-H || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-COR|R2CH-Br|RCH2-Br|R_acidic-H`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.175)
- Reasons: mapping_warning:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R_acidic-H
- Existing compound IDs missing in reaction slots: R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_041_intramolecular_likely_reacted_alkenyl_cor_r2ch_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-COR",
      "R2CH-Br",
      "RCH2-Br",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 36: Cluster Ar-Alkenyl -> R2CH-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkenyl -> R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Hydrogenation` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Hydrogenation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 37: Cluster RCH2-OR -> R3C-COR|R_acidic-H || events:none

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `RCH2-OR -> R3C-COR|R_acidic-H`
- Event signature: `none`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.175)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: RCH2-OR
- Missing formed-slot motifs: R_acidic-H
- Existing compound IDs missing in reaction slots: RCH2-OR, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_052_none_rch2_or_r3c_cor_r_acidic_h",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "RCH2-OR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "R3C-COR",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 38: Cluster none -> R2CH-Br || events:intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> R2CH-Br`
- Event signature: `intermolecular_or_multi_component`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.35)
- Reasons: mapping_warning:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_053_intermolecular_or_multi_component_reacted_r2ch_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "R2CH-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 39: Cluster none -> Alkyl-Br|Ar-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkyl-Br|Ar-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.175)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_056_intramolecular_likely_reacted_alkyl_br_ar_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-Br",
      "Ar-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 40: Cluster none -> RCH2-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.35)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_058_intramolecular_likely_reacted_rch2_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "RCH2-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 41: Cluster Ar-Alkyl|R_acidic-H -> Bn-Br|RCH2-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Ar-Alkyl|R_acidic-H -> Bn-Br|RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.35)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_059_intramolecular_likely_ar_alkyl_r_acidic_h_bn_br_rch2_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Ar-Alkyl",
      "R_acidic-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Bn-Br",
      "RCH2-Br"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 42: Cluster Bn-COR -> R2CH-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Bn-COR -> R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Reductive_amination` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Reductive_amination",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 43: Cluster CH3-OH|R2CH-CHO|R_acidic-H -> Ar-Alkyl|Benzofuran|HeteroAr-H|R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `CH3-OH|R2CH-CHO|R_acidic-H -> Ar-Alkyl|Benzofuran|HeteroAr-H|R2CH-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3042)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Benzofuran, HeteroAr-H
- Existing compound IDs missing in reaction slots: Benzofuran, HeteroAr-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_030_c_o_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_ch3_oh_r2ch_cho_ar_alkyl_benzofuran",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "CH3-OH",
      "R2CH-CHO",
      "R_acidic-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-Alkyl",
      "Benzofuran",
      "HeteroAr-H",
      "R2CH-OR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 44: Cluster Alkyl-H|Ar-CHO|R2CH-OH -> Ar-Alkyl|R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-H|Ar-CHO|R2CH-OH -> Ar-Alkyl|R2CH-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3917)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-H
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_031_c_o_bond_formation_intermolecular_or_multi_component_alkyl_h_ar_cho_ar_alkyl_r2ch_or",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-H",
      "Ar-CHO",
      "R2CH-OH"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-Alkyl",
      "R2CH-OR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 45: Cluster Alkyl-H|R2CH-CHO|RCH2-OH|R_acidic-H -> RCH2-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-H|R2CH-CHO|RCH2-OH|R_acidic-H -> RCH2-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.5125)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-H
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-H
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-H"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 46: Cluster Alkyl-H|RCH2-CHO|RCH2-OH|R_acidic-H -> RCH2-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-H|RCH2-CHO|RCH2-OH|R_acidic-H -> RCH2-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.5125)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-H
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-H
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-H"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 47: Cluster Alkenyl-H|R2CH-OH|RCH2-OH|Unclassified-Reactant -> RCH2-CHO|RCH2-OR|R_acidic-H || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-H|R2CH-OH|RCH2-OH|Unclassified-Reactant -> RCH2-CHO|RCH2-OR|R_acidic-H`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.4417)
- Reasons: mapping_warning:1, unclassified_motif:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R_acidic-H
- Existing compound IDs missing in reaction slots: R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_034_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_alkenyl_h_r2ch_oh_rch2_cho_rch2_or",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-H",
      "R2CH-OH",
      "RCH2-OH",
      "Unclassified-Reactant"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "RCH2-CHO",
      "RCH2-OR",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 48: Cluster Ar-CO2H -> Ar-COCl|Inorganic || events:intramolecular_likely

- Cluster count: 1
- Priority score: 3
- Action track: `expand_existing_reaction_type`
- Motif signature: `Ar-CO2H -> Ar-COCl|Inorganic`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Acyl_Halides_formation` (score=0.825)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Inorganic
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Acyl_Halides_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Inorganic"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 49: Cluster Ar-CO2H|Inorganic -> Ar-COCl || events:intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-CO2H|Inorganic -> Ar-COCl`
- Event signature: `intermolecular_or_multi_component`
- Top taxonomy candidate: `Acyl_Halides_formation` (score=0.675)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Acyl_Halides_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 50: Cluster Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Bn-LactamN|R2CH-CO2R|RCH2-CONR2|R_acidic-H || events:c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Bn-LactamN|R2CH-CO2R|RCH2-CONR2|R_acidic-H`
- Event signature: `c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation`
- Top taxonomy candidate: `Hydrogenation` (score=0.325)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Ar-C=N
- Missing formed-slot motifs: Bn-LactamN, R_acidic-H
- Existing compound IDs missing in reaction slots: Ar-C=N, Bn-LactamN, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_043_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_alkenyl_alkenyl_alkenyl_or_bn_lactamn_r2ch_co2r",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-Alkenyl",
      "Alkenyl-OR",
      "Alkyl-Si*",
      "Ar-C=N"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Bn-LactamN",
      "R2CH-CO2R",
      "RCH2-CONR2",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 51: Cluster Acyl-Cl|Alkyl-Si*|RCH2-NR2|RCH2-OR -> Alkyl-N(R)CO2R|R2CH-CO2H|R_acidic-H || events:c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Acyl-Cl|Alkyl-Si*|RCH2-NR2|RCH2-OR -> Alkyl-N(R)CO2R|R2CH-CO2H|R_acidic-H`
- Event signature: `c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Hiyama` (score=0.1625)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: RCH2-NR2, RCH2-OR
- Missing formed-slot motifs: Alkyl-N(R)CO2R, R_acidic-H
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, RCH2-NR2, RCH2-OR, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_044_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_acyl_cl_alkyl_si_alkyl_n_r_co2r_r2ch_co2h",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Acyl-Cl",
      "Alkyl-Si*",
      "RCH2-NR2",
      "RCH2-OR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-N(R)CO2R",
      "R2CH-CO2H",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 52: Cluster Alkyl-Si*|Inorganic|RCH2-NR2|RCH2-OR -> Alkyl-N(R)CO2R|R2CH-CO2H|R_acidic-H || events:c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-Si*|Inorganic|RCH2-NR2|RCH2-OR -> Alkyl-N(R)CO2R|R2CH-CO2H|R_acidic-H`
- Event signature: `c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Hiyama` (score=0.1625)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Inorganic, RCH2-NR2, RCH2-OR
- Missing formed-slot motifs: Alkyl-N(R)CO2R, R_acidic-H
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, RCH2-NR2, RCH2-OR, R_acidic-H
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_045_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_alkyl_si_inorganic_alkyl_n_r_co2r_r2ch_co2h",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-Si*",
      "Inorganic",
      "RCH2-NR2",
      "RCH2-OR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-N(R)CO2R",
      "R2CH-CO2H",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 53: Cluster Alkenyl-OR|Alkyl-F|Alkyl-Si*|Ar-CHO -> Alkyl-CF3|R2CH-OH|R3C-CO2R|R3C-OH || events:c_c_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-OR|Alkyl-F|Alkyl-Si*|Ar-CHO -> Alkyl-CF3|R2CH-OH|R3C-CO2R|R3C-OH`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Reduction_carbonyl` (score=0.3375)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkyl-CF3
- Existing compound IDs missing in reaction slots: Alkyl-CF3
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_046_c_c_bond_formation_intermolecular_or_multi_component_alkenyl_or_alkyl_f_alkyl_cf3_r2ch_oh",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-OR",
      "Alkyl-F",
      "Alkyl-Si*",
      "Ar-CHO"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-CF3",
      "R2CH-OH",
      "R3C-CO2R",
      "R3C-OH"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 54: Cluster Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Alkenyl-COR|Alkenyl-NR2|Ar-Alkyl|Ar-NR2|R2CH-NR2|R_acidic-H || events:c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Alkenyl-COR|Alkenyl-NR2|Ar-Alkyl|Ar-NR2|R2CH-NR2|R_acidic-H`
- Event signature: `c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation`
- Top taxonomy candidate: `Hydrogenation` (score=0.325)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Ar-C=N
- Missing formed-slot motifs: Alkenyl-NR2, R_acidic-H
- Existing compound IDs missing in reaction slots: Alkenyl-NR2, Ar-C=N, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_047_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_alkenyl_alkenyl_alkenyl_or_alkenyl_cor_alkenyl_nr2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-Alkenyl",
      "Alkenyl-OR",
      "Alkyl-Si*",
      "Ar-C=N"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-COR",
      "Alkenyl-NR2",
      "Ar-Alkyl",
      "Ar-NR2",
      "R2CH-NR2",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 55: Cluster R2CH-F -> Alkenyl-H|Alkyl-H|R2CH-I || events:c_c_bond_formation|intramolecular_likely|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `expand_existing_reaction_type`
- Motif signature: `R2CH-F -> Alkenyl-H|Alkyl-H|R2CH-I`
- Event signature: `c_c_bond_formation|intramolecular_likely|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Aliphatic_Halide_Exchange` (score=0.7667)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-H, Alkyl-H
- Existing compound IDs missing in reaction slots: Alkenyl-H, Alkyl-H
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Aliphatic_Halide_Exchange",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Alkenyl-H",
      "Alkyl-H"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 56: Cluster Alkyl-H|HeteroAr-H|RCH2-Br -> Ar-Alkyl|Ar-CHO|Pyrazole|RCH2-Cl || events:c_c_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-H|HeteroAr-H|RCH2-Br -> Ar-Alkyl|Ar-CHO|Pyrazole|RCH2-Cl`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Aliphatic_Halide_Exchange` (score=0.3042)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-H
- Missing formed-slot motifs: Pyrazole
- Existing compound IDs missing in reaction slots: Alkyl-H, Pyrazole
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_049_c_c_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkyl_h_heteroar_h_ar_alkyl_ar_cho",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-H",
      "HeteroAr-H",
      "RCH2-Br"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Ar-Alkyl",
      "Ar-CHO",
      "Pyrazole",
      "RCH2-Cl"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 57: Cluster Alkyl-CF3 -> Alkyl-Cl|Alkyl-F || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-CF3 -> Alkyl-Cl|Alkyl-F`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.35)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-CF3
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-CF3
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_050_intramolecular_likely_leaving_group_displacement_alkyl_cf3_alkyl_cl_alkyl_f",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-CF3"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-Cl",
      "Alkyl-F"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 58: Cluster Alkenyl-Cl|Ar-Alkenyl -> Alkyl-CF3|Ar-Alkyl || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-Cl|Ar-Alkenyl -> Alkyl-CF3|Ar-Alkyl`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Hydrogenation` (score=0.65)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkyl-CF3
- Existing compound IDs missing in reaction slots: Alkyl-CF3
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Hydrogenation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Alkyl-CF3"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 59: Cluster none -> Bn-Cl|HeteroAr-Cl|Quinoline || events:intramolecular_likely|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Bn-Cl|HeteroAr-Cl|Quinoline`
- Event signature: `intramolecular_likely|ring_closure_or_annulation`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.1167)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Quinoline
- Existing compound IDs missing in reaction slots: Quinoline
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_054_intramolecular_likely_ring_closure_or_annulation_reacted_bn_cl_heteroar_cl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Bn-Cl",
      "HeteroAr-Cl",
      "Quinoline"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 60: Cluster Ar-Alkyl -> Bn-Br|HeteroAr-Br|Quinoline || events:intramolecular_likely|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkyl -> Bn-Br|HeteroAr-Br|Quinoline`
- Event signature: `intramolecular_likely|ring_closure_or_annulation`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.65)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Quinoline
- Existing compound IDs missing in reaction slots: Quinoline
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Oxidation_alcohol",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Quinoline"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 61: Cluster none -> Alkyl-Br|Bn-Br|HeteroAr-Br|Pyridine|R_acidic-H || events:c_c_bond_formation|intramolecular_likely|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkyl-Br|Bn-Br|HeteroAr-Br|Pyridine|R_acidic-H`
- Event signature: `c_c_bond_formation|intramolecular_likely|ring_closure_or_annulation`
- Top taxonomy candidate: `Alcohol_to_Alkyl_Halide` (score=0.14)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Pyridine, R_acidic-H
- Existing compound IDs missing in reaction slots: Pyridine, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_057_c_c_bond_formation_intramolecular_likely_ring_closure_or_annulation_reacted_alkyl_br_bn_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-Br",
      "Bn-Br",
      "HeteroAr-Br",
      "Pyridine",
      "R_acidic-H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 62: Cluster Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R3C-OR|RCH2-CO2H -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|R2CH-OH|RCH2-CONHR|RCH2-OH || events:none

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R3C-OR|RCH2-CO2H -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|R2CH-OH|RCH2-CONHR|RCH2-OH`
- Event signature: `none`
- Top taxonomy candidate: `Amide_formation` (score=0.4482)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-Ar, R3C-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-Ar, R3C-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_061_none_alkyl_co2r_alkyl_n_r_co2r_alkyl_nhcor_alkyl_nrcor",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-CO2R",
      "Alkyl-N(R)CO2R",
      "Ar-Alkyl",
      "Ar-Ar",
      "R2CH-CO2H",
      "R3C-OR",
      "RCH2-CO2H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-NHCOR",
      "Alkyl-NRCOR",
      "R2CH-CONH2",
      "R2CH-CONHR",
      "R2CH-CONR2",
      "R2CH-OH",
      "RCH2-CONHR",
      "RCH2-OH"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 63: Cluster Alkyl-N(R)CO2R|Alkyl-Si*|Ar-Alkyl|Ar-Ar|Ar-NHCOR|R2CH-CO2H|R2CH-OR|RCH2-CO2H -> Alkyl-NHCOR|R2CH-CONHR|R2CH-OH|RCH2-CONHR|RCH2-OH || events:none

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|Alkyl-Si*|Ar-Alkyl|Ar-Ar|Ar-NHCOR|R2CH-CO2H|R2CH-OR|RCH2-CO2H -> Alkyl-NHCOR|R2CH-CONHR|R2CH-OH|RCH2-CONHR|RCH2-OH`
- Event signature: `none`
- Top taxonomy candidate: `Amide_formation` (score=0.4537)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-Ar",
      "R2CH-OR"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 64: Cluster Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|RCH2-CHO|RCH2-CONHR || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|RCH2-CHO|RCH2-CONHR`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.4625)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-Ar",
      "R2CH-OR"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 65: Cluster Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R|RCH2-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-NH2|RCH2-CHO|RCH2-CONHR|RCH2-CONR2 || events:amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R|RCH2-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-NH2|RCH2-CHO|RCH2-CONHR|RCH2-CONR2`
- Event signature: `amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.5062)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-Ar",
      "R2CH-OR"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 66: Cluster Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R|RCH2-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|R2CH-NH2|RCH2-CHO|RCH2-CONHR || events:amidation_like|c_c_bond_formation|c_n_bond_formation|c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R|RCH2-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|R2CH-NH2|RCH2-CHO|RCH2-CONHR`
- Event signature: `amidation_like|c_c_bond_formation|c_n_bond_formation|c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.5062)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-Ar, R2CH-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-Ar",
      "R2CH-OR"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 67: Cluster Alkyl-N(R)CO2R|Alkyl-Si*|Ar-Alkyl|Ar-NHCOR|R2CH-CO2R|R3C-OR|RCH2-CO2H|RCH2-CONHR -> R2CH-CO2H|R2CH-NH2|R2CH-OH|RCH2-NH2|RCH2-OH || events:none

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-N(R)CO2R|Alkyl-Si*|Ar-Alkyl|Ar-NHCOR|R2CH-CO2R|R3C-OR|RCH2-CO2H|RCH2-CONHR -> R2CH-CO2H|R2CH-NH2|R2CH-OH|RCH2-NH2|RCH2-OH`
- Event signature: `none`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.2213)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, R3C-OR, RCH2-CONHR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, R3C-OR, RCH2-CONHR
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_066_none_alkyl_n_r_co2r_alkyl_si_r2ch_co2h_r2ch_nh2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Alkyl-Si*",
      "Ar-Alkyl",
      "Ar-NHCOR",
      "R2CH-CO2R",
      "R3C-OR",
      "RCH2-CO2H",
      "RCH2-CONHR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "R2CH-CO2H",
      "R2CH-NH2",
      "R2CH-OH",
      "RCH2-NH2",
      "RCH2-OH"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 68: Cluster Alkyl-AromN|Ar-Alkyl|Ar-CO2H|CH3-NHR|Inorganic|R3C-OR|RCH2-OR -> Alkyl-NRCOR|Ar-CONR2|AromN-H|Imidazole|R2CH-OH|RCH2-OH || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-AromN|Ar-Alkyl|Ar-CO2H|CH3-NHR|Inorganic|R3C-OR|RCH2-OR -> Alkyl-NRCOR|Ar-CONR2|AromN-H|Imidazole|R2CH-OH|RCH2-OH`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.3024)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-AromN, Inorganic, R3C-OR, RCH2-OR
- Missing formed-slot motifs: AromN-H, Imidazole
- Existing compound IDs missing in reaction slots: Alkyl-AromN, AromN-H, Imidazole, R3C-OR, RCH2-OR
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_067_amidation_like_c_n_bond_formation_intermolecular_or_multi_component_alkyl_aromn_ar_alkyl_alkyl_nrcor_ar_conr2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-AromN",
      "Ar-Alkyl",
      "Ar-CO2H",
      "CH3-NHR",
      "Inorganic",
      "R3C-OR",
      "RCH2-OR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-NRCOR",
      "Ar-CONR2",
      "AromN-H",
      "Imidazole",
      "R2CH-OH",
      "RCH2-OH"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 69: Cluster Alkyl-N(R)CO2R|Ar-Alkyl|R2CH-CO2H|R2CH-OR|RCH2-CO2R|RCH2-OR -> Alkyl-NHCOR|R2CH-CONHR|R2CH-OH|RCH2-CO2H|RCH2-OH || events:none

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-N(R)CO2R|Ar-Alkyl|R2CH-CO2H|R2CH-OR|RCH2-CO2R|RCH2-OR -> Alkyl-NHCOR|R2CH-CONHR|R2CH-OH|RCH2-CO2H|RCH2-OH`
- Event signature: `none`
- Top taxonomy candidate: `Amide_formation` (score=0.2483)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, R2CH-OR, RCH2-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, R2CH-OR, RCH2-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_068_none_alkyl_n_r_co2r_ar_alkyl_alkyl_nhcor_r2ch_conhr",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-Alkyl",
      "R2CH-CO2H",
      "R2CH-OR",
      "RCH2-CO2R",
      "RCH2-OR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-NHCOR",
      "R2CH-CONHR",
      "R2CH-OH",
      "RCH2-CO2H",
      "RCH2-OH"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 70: Cluster Alkyl-N(R)CO2R|Bn-NH2|R2CH-CO2H|R3C-OR -> Alkyl-NHCOR|Ar-Alkyl|R2CH-CONHR|R2CH-NH2|R2CH-OH || events:none

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|Bn-NH2|R2CH-CO2H|R3C-OR -> Alkyl-NHCOR|Ar-Alkyl|R2CH-CONHR|R2CH-NH2|R2CH-OH`
- Event signature: `none`
- Top taxonomy candidate: `Amide_formation` (score=0.465)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, R3C-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, R3C-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "R3C-OR"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 71: Cluster Alkyl-AromN|Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R3C-CO2R|RCH2-CO2R -> Alkyl-NHCOR|Alkyl-NRCOR|AromN-H|Imidazole|R2CH-CONHR|R2CH-NH2|R2CH-OH|R3C-CO2H || events:none

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-AromN|Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R3C-CO2R|RCH2-CO2R -> Alkyl-NHCOR|Alkyl-NRCOR|AromN-H|Imidazole|R2CH-CONHR|R2CH-NH2|R2CH-OH|R3C-CO2H`
- Event signature: `none`
- Top taxonomy candidate: `Claisen_condensation` (score=0.2438)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-AromN, Alkyl-N(R)CO2R, Ar-Ar
- Missing formed-slot motifs: AromN-H, Imidazole
- Existing compound IDs missing in reaction slots: Alkyl-AromN, Alkyl-N(R)CO2R, Ar-Ar, AromN-H, Imidazole
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_070_none_alkyl_aromn_alkyl_co2r_alkyl_nhcor_alkyl_nrcor",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-AromN",
      "Alkyl-CO2R",
      "Alkyl-N(R)CO2R",
      "Ar-Alkyl",
      "Ar-Ar",
      "R2CH-CO2H",
      "R3C-CO2R",
      "RCH2-CO2R"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-NHCOR",
      "Alkyl-NRCOR",
      "AromN-H",
      "Imidazole",
      "R2CH-CONHR",
      "R2CH-NH2",
      "R2CH-OH",
      "R3C-CO2H"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 72: Cluster Ar-NH2|HeteroAr-OH|RCH2-Br|RCH2-CO2H -> Ar-NHCOR|Ar-OR || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-NH2|HeteroAr-OH|RCH2-Br|RCH2-CO2H -> Ar-NHCOR|Ar-OR`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.675)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: HeteroAr-OH
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: HeteroAr-OH
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "HeteroAr-OH"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 73: Cluster Alkenyl-CO2H|Ar-OH|Bn-Cl|Bn-NH2|Inorganic -> Alkenyl-CONHR|Alkyl-NHCOR|Ar-OR|Bn-OR || events:amidation_like|benzyl_o_alkylation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-CO2H|Ar-OH|Bn-Cl|Bn-NH2|Inorganic -> Alkenyl-CONHR|Alkyl-NHCOR|Ar-OR|Bn-OR`
- Event signature: `amidation_like|benzyl_o_alkylation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.6525)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 74: Cluster Ar-OH|Bn-Br|Inorganic|R2CH-CO2H|R2CH-NH2 -> Alkyl-NHCOR|Ar-OR|Bn-OR|R2CH-CONHR || events:amidation_like|benzyl_o_alkylation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-OH|Bn-Br|Inorganic|R2CH-CO2H|R2CH-NH2 -> Alkyl-NHCOR|Ar-OR|Bn-OR|R2CH-CONHR`
- Event signature: `amidation_like|benzyl_o_alkylation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.6525)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 75: Cluster Alkyl-F|Inorganic|R2CH-CO2H|R2CH-NH2|R2CH-OH -> Alkyl-N(R)CO2R|Alkyl-NHCOR|R2CH-CONHR || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-F|Inorganic|R2CH-CO2H|R2CH-NH2|R2CH-OH -> Alkyl-N(R)CO2R|Alkyl-NHCOR|R2CH-CONHR`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.5067)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Inorganic
- Missing formed-slot motifs: Alkyl-N(R)CO2R
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": [
      "Alkyl-N(R)CO2R"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 76: Cluster Alkyl-N(R)CO2R|Ar-CO2H|RCH2-I -> Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|Ar-CO2H|RCH2-I -> Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.45)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 77: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Bn-Br|Bn-Cl|Inorganic -> Alkenyl-CO2R|Ar-Alkyl|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Bn-Br|Bn-Cl|Inorganic -> Alkenyl-CO2R|Ar-Alkyl|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3567)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_076_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_leaving_group_displacement_ring_closure_or_annulation_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_alkyl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Bn-Br",
      "Bn-Cl",
      "Inorganic"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Ar-Alkyl",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 78: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.435)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_077_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_co2r",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Inorganic",
      "RCH2-I"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 79: Cluster Alkyl-N(R)CO2R|Alkyl-OTf -> Alkenyl-CO2R|Alkenyl-LactamN|Alkenyl-SR|Alkyl-CF3|Ar-CO2R|Ar-NHCOR|Inorganic|R2CH-CONR2 || events:c_c_bond_formation|c_n_bond_formation|c_o_bond_formation|c_s_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-N(R)CO2R|Alkyl-OTf -> Alkenyl-CO2R|Alkenyl-LactamN|Alkenyl-SR|Alkyl-CF3|Ar-CO2R|Ar-NHCOR|Inorganic|R2CH-CONR2`
- Event signature: `c_c_bond_formation|c_n_bond_formation|c_o_bond_formation|c_s_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3688)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R
- Missing formed-slot motifs: Alkenyl-LactamN, Alkyl-CF3, Inorganic
- Existing compound IDs missing in reaction slots: Alkenyl-LactamN, Alkyl-CF3, Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_078_c_c_bond_formation_c_n_bond_formation_c_o_bond_formation_c_s_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_alkyl_n_r_co2r_alkyl_otf_alkenyl_co2r_alkenyl_lactamn",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Alkyl-OTf"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Alkenyl-LactamN",
      "Alkenyl-SR",
      "Alkyl-CF3",
      "Ar-CO2R",
      "Ar-NHCOR",
      "Inorganic",
      "R2CH-CONR2"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 80: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|R3C-CO2R|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|R3C-CO2R|RCH2-NHR`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.47)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Esterification",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 81: Cluster Ar-CO2H|CH3-I|Inorganic|R2CH-CO2H|R2CH-NH2 -> Alkyl-NHCOR|Ar-CONHR|R2CH-CO2R || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-CO2H|CH3-I|Inorganic|R2CH-CO2H|R2CH-NH2 -> Alkyl-NHCOR|Ar-CONHR|R2CH-CO2R`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Amide_formation` (score=0.6233)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 82: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Bn-I|Inorganic -> Alkenyl-CO2R|Ar-Alkyl|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Bn-I|Inorganic -> Alkenyl-CO2R|Ar-Alkyl|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Esterification` (score=0.4)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_081_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_alkyl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Bn-I",
      "Inorganic"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Ar-Alkyl",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 83: Cluster Alkenyl-CO2H|Alkenyl-LactamN|Alkenyl-SR|Alkyl-N(R)CO2R|Alkyl-NRCOR|Ar-CO2H|Ar-NHCOR|Inorganic -> none || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkenyl-LactamN|Alkenyl-SR|Alkyl-N(R)CO2R|Alkyl-NRCOR|Ar-CO2H|Ar-NHCOR|Inorganic -> none`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Amide_formation` (score=0.3187)
- Reasons: mapping_warning:1, missing_formed_motifs:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Alkyl-NRCOR, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Alkyl-NRCOR
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_082_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkenyl_co2h_alkenyl_lactamn_formed",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkenyl-LactamN",
      "Alkenyl-SR",
      "Alkyl-N(R)CO2R",
      "Alkyl-NRCOR",
      "Ar-CO2H",
      "Ar-NHCOR",
      "Inorganic"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": []
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 84: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|R2CH-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|R2CH-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.435)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_083_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_co2r",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Inorganic",
      "R2CH-I"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 85: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|R2CH-Cl -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|R3C-CO2R|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|R2CH-Cl -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|R3C-CO2R|RCH2-NHR`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.47)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Esterification",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 86: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Cl -> Alkenyl-CO2R|Alkyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Cl -> Alkenyl-CO2R|Alkyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.47)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Esterification",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 87: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Cl -> Alkenyl-CO2R|Ar-CO2R|R2CH-CO2R|R2CH-NHR|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Cl -> Alkenyl-CO2R|Ar-CO2R|R2CH-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.47)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Esterification",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 88: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Cl -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-CO2R|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Cl -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-CO2R|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.47)
- Reasons: motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Esterification",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 89: Cluster Acyl-Br|Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|Ar-CO2H|R2CH-CO2H|RCH2-Br -> Alkyl-Hydrazide_NHR|Alkyl-Hydrazide_NR|Alkyl-NHCOR|Alkyl-NRCOR|Ar-CO2R|R2CH-CONHR || events:amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Acyl-Br|Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|Ar-CO2H|R2CH-CO2H|RCH2-Br -> Alkyl-Hydrazide_NHR|Alkyl-Hydrazide_NR|Alkyl-NHCOR|Alkyl-NRCOR|Ar-CO2R|R2CH-CONHR`
- Event signature: `amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Amide_formation` (score=0.3375)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Acyl-Br, Alkyl-N(R)CO2R, Ar-Ar
- Missing formed-slot motifs: Alkyl-Hydrazide_NHR, Alkyl-Hydrazide_NR
- Existing compound IDs missing in reaction slots: Acyl-Br, Alkyl-Hydrazide_NHR, Alkyl-Hydrazide_NR, Alkyl-N(R)CO2R, Ar-Ar
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_088_amidation_like_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_acyl_br_alkyl_co2r_alkyl_hydrazide_nhr_alkyl_hydrazide_nr",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Acyl-Br",
      "Alkyl-CO2R",
      "Alkyl-N(R)CO2R",
      "Ar-Alkyl",
      "Ar-Ar",
      "Ar-CO2H",
      "R2CH-CO2H",
      "RCH2-Br"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-Hydrazide_NHR",
      "Alkyl-Hydrazide_NR",
      "Alkyl-NHCOR",
      "Alkyl-NRCOR",
      "Ar-CO2R",
      "R2CH-CONHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 90: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Esterification` (score=0.435)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_089_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_co2r",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Inorganic",
      "RCH2-I"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 91: Cluster Alkyl-N(R)CO2R|R2CH-I -> Alkenyl-CO2R|Alkenyl-LactamN|Alkenyl-SR|Alkyl-NRCOR|Ar-CO2R|Ar-NHCOR|Inorganic|R2CH-CONR2 || events:c_c_bond_formation|c_n_bond_formation|c_o_bond_formation|c_s_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-N(R)CO2R|R2CH-I -> Alkenyl-CO2R|Alkenyl-LactamN|Alkenyl-SR|Alkyl-NRCOR|Ar-CO2R|Ar-NHCOR|Inorganic|R2CH-CONR2`
- Event signature: `c_c_bond_formation|c_n_bond_formation|c_o_bond_formation|c_s_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.4125)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R
- Missing formed-slot motifs: Alkenyl-LactamN, Inorganic
- Existing compound IDs missing in reaction slots: Alkenyl-LactamN, Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_090_c_c_bond_formation_c_n_bond_formation_c_o_bond_formation_c_s_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkyl_n_r_co2r_r2ch_i_alkenyl_co2r_alkenyl_lactamn",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "R2CH-I"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Alkenyl-LactamN",
      "Alkenyl-SR",
      "Alkyl-NRCOR",
      "Ar-CO2R",
      "Ar-NHCOR",
      "Inorganic",
      "R2CH-CONR2"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 92: Cluster Alkenyl-CO2H|Alkenyl-LactamN|Alkenyl-SR|Alkyl-N(R)CO2R|Alkyl-NRCOR|Ar-CO2H|Ar-NHCOR|Inorganic -> none || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkenyl-LactamN|Alkenyl-SR|Alkyl-N(R)CO2R|Alkyl-NRCOR|Ar-CO2H|Ar-NHCOR|Inorganic -> none`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Amide_formation` (score=0.3187)
- Reasons: mapping_warning:1, missing_formed_motifs:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Alkyl-NRCOR, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Alkyl-NRCOR
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_091_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkenyl_co2h_alkenyl_lactamn_formed",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkenyl-LactamN",
      "Alkenyl-SR",
      "Alkyl-N(R)CO2R",
      "Alkyl-NRCOR",
      "Ar-CO2H",
      "Ar-NHCOR",
      "Inorganic"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": []
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 93: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Cl -> Alkenyl-CO2R|Ar-CO2R|R2CH-CO2R|R2CH-NHR|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Cl -> Alkenyl-CO2R|Ar-CO2R|R2CH-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Esterification` (score=0.47)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Esterification",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Inorganic"
    ]
  },
  "product_slot_additions_todo": {
    "product": []
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 94: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Alkyl-CF3|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Alkyl-CF3|Ar-CO2R|R2CH-NHR|RCH2-NHR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement`
- Top taxonomy candidate: `Esterification` (score=0.4)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Inorganic
- Missing formed-slot motifs: Alkyl-CF3
- Existing compound IDs missing in reaction slots: Alkyl-CF3, Alkyl-N(R)CO2R
- Unknown compound IDs: Inorganic

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_093_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_alkyl_cf3",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CO2H",
      "Alkyl-N(R)CO2R",
      "Ar-CO2H",
      "Inorganic",
      "RCH2-I"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-CO2R",
      "Alkyl-CF3",
      "Ar-CO2R",
      "R2CH-NHR",
      "RCH2-NHR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```
- `new_compound_entries`
```json
[
  {
    "id": "Inorganic",
    "todo": "define_compound_id_explicitly",
    "notes": "motif has no A-B structure; add explicit id/smarts in organic_compounds.v1.3.json"
  }
]
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 95: Cluster Alkyl-N(R)CO2R|Ar-CO2H|Ar-OH|Ar-SO2Cl|CH3-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|Ar-CONHR|Ar-CONR2|Ar-OSO2R || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|Ar-CO2H|Ar-OH|Ar-SO2Cl|CH3-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|Ar-CONHR|Ar-CONR2|Ar-OSO2R`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.54)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-N(R)CO2R, Ar-SO2Cl
- Missing formed-slot motifs: Ar-OSO2R
- Existing compound IDs missing in reaction slots: Alkyl-N(R)CO2R, Ar-OSO2R, Ar-SO2Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Alkyl-N(R)CO2R",
      "Ar-SO2Cl"
    ]
  },
  "product_slot_additions_todo": {
    "product": [
      "Ar-OSO2R"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 96: Cluster Ar-CO2H|Ar-OH|Ar-SO2Cl|CH3-NH2 -> Alkyl-NHCOR|Ar-CONHR|Ar-OSO2R || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-CO2H|Ar-OH|Ar-SO2Cl|CH3-NH2 -> Alkyl-NHCOR|Ar-CONHR|Ar-OSO2R`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.5583)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Ar-SO2Cl
- Missing formed-slot motifs: Ar-OSO2R
- Existing compound IDs missing in reaction slots: Ar-OSO2R, Ar-SO2Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Ar-SO2Cl"
    ]
  },
  "product_slot_additions_todo": {
    "product": [
      "Ar-OSO2R"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 97: Cluster Ar-CO2H|Ar-OH|Ar-SO2Cl|R2CH-NH2 -> Alkyl-NHCOR|Ar-CONHR|Ar-OSO2R || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-CO2H|Ar-OH|Ar-SO2Cl|R2CH-NH2 -> Alkyl-NHCOR|Ar-CONHR|Ar-OSO2R`
- Event signature: `amidation_like|c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Amide_formation` (score=0.5583)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Ar-SO2Cl
- Missing formed-slot motifs: Ar-OSO2R
- Existing compound IDs missing in reaction slots: Ar-OSO2R, Ar-SO2Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Amide_formation",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": [
      "Ar-SO2Cl"
    ]
  },
  "product_slot_additions_todo": {
    "product": [
      "Ar-OSO2R"
    ]
  },
  "constraints_review_todo": [
    "review include/exclude reacted constraints",
    "review include/exclude formed constraints",
    "confirm min_reactant_slot_matches / min_product_slot_matches"
  ]
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 98: Cluster Alkenyl-Br|Alkenyl-CO2R|Ar-NHR|Bn-NHR|Furan|HeteroAr-H -> Alkyl-NRCOR|Ar-NRCOR|R3C-Br|R3C-CO2H|R3C-Cl|R3C-OR || events:c_c_bond_formation|c_n_bond_formation|ester_hydrolysis_like|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-Br|Alkenyl-CO2R|Ar-NHR|Bn-NHR|Furan|HeteroAr-H -> Alkyl-NRCOR|Ar-NRCOR|R3C-Br|R3C-CO2H|R3C-Cl|R3C-OR`
- Event signature: `c_c_bond_formation|c_n_bond_formation|ester_hydrolysis_like|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3917)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Furan
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Furan
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_098_c_c_bond_formation_c_n_bond_formation_ester_hydrolysis_like_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkenyl_br_alkenyl_co2r_alkyl_nrcor_ar_nrcor",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-Br",
      "Alkenyl-CO2R",
      "Ar-NHR",
      "Bn-NHR",
      "Furan",
      "HeteroAr-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-NRCOR",
      "Ar-NRCOR",
      "R3C-Br",
      "R3C-CO2H",
      "R3C-Cl",
      "R3C-OR"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 99: Cluster AromN-H|Benzimidazole|RCH2-CO2R -> Benzimidazole || events:c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `AromN-H|Benzimidazole|RCH2-CO2R -> Benzimidazole`
- Event signature: `c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Ircatalyzed_CH_borylation` (score=0.2167)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Benzimidazole
- Missing formed-slot motifs: Benzimidazole
- Existing compound IDs missing in reaction slots: Benzimidazole
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_099_c_n_bond_formation_intermolecular_or_multi_component_aromn_h_benzimidazole_benzimidazole",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "AromN-H",
      "Benzimidazole",
      "RCH2-CO2R"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Benzimidazole"
    ]
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.

## Priority 100: Cluster Alkyl-CO2R|AromN-H|Imidazole -> none || events:c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 3
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-CO2R|AromN-H|Imidazole -> none`
- Event signature: `c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Ircatalyzed_CH_borylation` (score=0.2833)
- Reasons: mapping_warning:1, missing_formed_motifs:1, motif_outside_reaction_taxonomy:1

### Gap Analysis

- Missing reacted-slot motifs: Imidazole
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Imidazole
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_100_c_n_bond_formation_intermolecular_or_multi_component_alkyl_co2r_aromn_h_formed",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-CO2R",
      "AromN-H",
      "Imidazole"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": []
  },
  "constraints": {
    "include_reacted": [],
    "exclude_reacted": [],
    "include_formed": [],
    "exclude_formed": [],
    "min_reactant_slot_matches": 0,
    "min_product_slot_matches": 0
  }
}
```

### Validation Checklist

- Confirm representative reactions are true single-step transforms.
- Verify motifs against product/reactant structures with RDKit inspection.
- If updating existing family, ensure constraints still disambiguate close families.
- Run taxonomy validators and reaction reliability report after edits.
