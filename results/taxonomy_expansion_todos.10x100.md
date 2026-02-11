# Taxonomy Expansion TODOs

Generated from unresolved discovery clusters.

## Summary

- Input clusters: 100
- TODO entries: 100
- Expand existing reaction type: 4
- Review constraints/product patterns: 47
- Propose new reaction family: 49

## Priority 1: Cluster Alkyl-Si*|RCH2-NR2|RCH2-OR -> R2CH-CO2H|RCH2-NH2|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component

- Cluster count: 27
- Priority score: 30
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-Si*|RCH2-NR2|RCH2-OR -> R2CH-CO2H|RCH2-NH2|R_acidic-H`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Hiyama` (score=0.2167)
- Reasons: mapping_warning:17, unknown_reaction_type:27

### Gap Analysis

- Missing reacted-slot motifs: RCH2-NR2, RCH2-OR
- Missing formed-slot motifs: R_acidic-H
- Existing compound IDs missing in reaction slots: RCH2-NR2, RCH2-OR, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_001_c_c_bond_formation_intermolecular_or_multi_component_alkyl_si_rch2_nr2_r2ch_co2h_rch2_nh2",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-Si*",
      "RCH2-NR2",
      "RCH2-OR"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "R2CH-CO2H",
      "RCH2-NH2",
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

## Priority 2: Cluster none -> none || events:none

- Cluster count: 22
- Priority score: 29
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> none`
- Event signature: `none`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: empty_reaction_key:7, low_reaction_key_quality:22, missing_formed_motifs:22, unclassified_motif:22, unknown_reaction_type:22

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_003_none_reacted_formed",
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

## Priority 3: Cluster none -> Alkenyl-Br || events:intramolecular_likely

- Cluster count: 22
- Priority score: 25
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:22

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Br
- Existing compound IDs missing in reaction slots: Alkenyl-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_002_intramolecular_likely_reacted_alkenyl_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Br"
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

## Priority 4: Cluster RCH2-OH -> RCH2-Br || events:intramolecular_likely

- Cluster count: 20
- Priority score: 23
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `RCH2-OH -> RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:20

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Br
- Existing compound IDs missing in reaction slots: RCH2-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "RCH2-Br"
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

## Priority 5: Cluster R2CH-OH -> Alkenyl-Alkenyl|RCH2-Br || events:intramolecular_likely

- Cluster count: 19
- Priority score: 22
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R2CH-OH -> Alkenyl-Alkenyl|RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: mapping_warning:1, unknown_reaction_type:19

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Br
- Existing compound IDs missing in reaction slots: RCH2-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "RCH2-Br"
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

## Priority 6: Cluster RCH2-Cl -> RCH2-I || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 14
- Priority score: 17
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `RCH2-Cl -> RCH2-I`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:14

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-I
- Existing compound IDs missing in reaction slots: RCH2-I
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "RCH2-I"
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

## Priority 7: Cluster RCH2-OH -> RCH2-Cl || events:intramolecular_likely

- Cluster count: 12
- Priority score: 15
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `RCH2-OH -> RCH2-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:12

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Cl
- Existing compound IDs missing in reaction slots: RCH2-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "RCH2-Cl"
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

## Priority 8: Cluster R2CH-OH -> RCH2-Br || events:intramolecular_likely

- Cluster count: 11
- Priority score: 14
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R2CH-OH -> RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: mapping_warning:2, unknown_reaction_type:11

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Br
- Existing compound IDs missing in reaction slots: RCH2-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "RCH2-Br"
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

## Priority 9: Cluster Ar-C=N|Ar-NR2|Bn-OH -> Ar-Ar|Ar-AromN|Bn-Br || events:intermolecular_or_multi_component

- Cluster count: 9
- Priority score: 14
- Action track: `propose_new_reaction_family`
- Motif signature: `Ar-C=N|Ar-NR2|Bn-OH -> Ar-Ar|Ar-AromN|Bn-Br`
- Event signature: `intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3333)
- Reasons: mapping_warning:9, motif_outside_reaction_taxonomy:9, unknown_reaction_type:9

### Gap Analysis

- Missing reacted-slot motifs: Ar-C=N, Ar-NR2
- Missing formed-slot motifs: Bn-Br
- Existing compound IDs missing in reaction slots: Ar-C=N, Ar-NR2, Bn-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_014_intermolecular_or_multi_component_ar_c_n_ar_nr2_ar_ar_ar_aromn",
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

## Priority 10: Cluster Ar-Alkyl|R2CH-OH -> Ar-Alkenyl|RCH2-Br || events:intramolecular_likely

- Cluster count: 10
- Priority score: 13
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkyl|R2CH-OH -> Ar-Alkenyl|RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.65)
- Reasons: unknown_reaction_type:10

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Br
- Existing compound IDs missing in reaction slots: RCH2-Br
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
      "RCH2-Br"
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

## Priority 11: Cluster Alkyl-COR|R_acidic-H -> R3C-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 9
- Priority score: 12
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-COR|R_acidic-H -> R3C-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Arylation_acidic_C_H` (score=0.65)
- Reasons: mapping_warning:9, unknown_reaction_type:9

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Arylation_acidic_C_H",
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

## Priority 12: Cluster Bn-OH -> Bn-Cl || events:intramolecular_likely

- Cluster count: 9
- Priority score: 12
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Bn-OH -> Bn-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:9

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Bn-Cl
- Existing compound IDs missing in reaction slots: Bn-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Bn-Cl"
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

## Priority 13: Cluster R2CH-OH -> R2CH-Br || events:intramolecular_likely

- Cluster count: 9
- Priority score: 12
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R2CH-OH -> R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:9

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Br
- Existing compound IDs missing in reaction slots: R2CH-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R2CH-Br"
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

## Priority 14: Cluster Alkyl-H|Ar-COR|Ar-H -> Ar-C=N || events:c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 7
- Priority score: 12
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-H|Ar-COR|Ar-H -> Ar-C=N`
- Event signature: `c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Ircatalyzed_CH_borylation` (score=0.2167)
- Reasons: mapping_warning:4, motif_outside_reaction_taxonomy:7, unknown_reaction_type:7

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-H
- Missing formed-slot motifs: Ar-C=N
- Existing compound IDs missing in reaction slots: Alkyl-H, Ar-C=N
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_019_c_n_bond_formation_intermolecular_or_multi_component_alkyl_h_ar_cor_ar_c_n",
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

## Priority 15: Cluster Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Alkyl-N(R)CO2R|Ar-Alkyl|RCH2-CO2R|RCH2-COR|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component

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
  "id": "todo_family_011_c_c_bond_formation_intermolecular_or_multi_component_alkenyl_alkenyl_alkenyl_or_alkyl_n_r_co2r_ar_alkyl",
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

## Priority 16: Cluster Alkyl-COR|R_acidic-H -> none || events:c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 8
- Priority score: 11
- Action track: `expand_existing_reaction_type`
- Motif signature: `Alkyl-COR|R_acidic-H -> none`
- Event signature: `c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Arylation_acidic_C_H` (score=0.85)
- Reasons: mapping_warning:3, missing_formed_motifs:8, unknown_reaction_type:8

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Arylation_acidic_C_H",
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

## Priority 17: Cluster Bn-Cl -> Bn-I || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 8
- Priority score: 11
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Bn-Cl -> Bn-I`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:8

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Bn-I
- Existing compound IDs missing in reaction slots: Bn-I
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Bn-I"
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

## Priority 18: Cluster Ar-Alkyl -> Ar-Br|Bn-Br || events:intramolecular_likely

- Cluster count: 8
- Priority score: 11
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkyl -> Ar-Br|Bn-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.65)
- Reasons: unknown_reaction_type:8

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Bn-Br
- Existing compound IDs missing in reaction slots: Bn-Br
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
      "Bn-Br"
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

## Priority 19: Cluster Bn-OH -> Bn-Br || events:intramolecular_likely

- Cluster count: 8
- Priority score: 11
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Bn-OH -> Bn-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:8

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Bn-Br
- Existing compound IDs missing in reaction slots: Bn-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Bn-Br"
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

## Priority 20: Cluster Alkyl-N(R)CO2R|CH3-OH -> R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 6
- Priority score: 11
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-N(R)CO2R|CH3-OH -> R2CH-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.675)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:6, unknown_reaction_type:6

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

## Priority 21: Cluster none -> R2CH-Br || events:intramolecular_likely

- Cluster count: 7
- Priority score: 10
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:7

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Br
- Existing compound IDs missing in reaction slots: R2CH-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_020_intramolecular_likely_reacted_r2ch_br",
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

## Priority 22: Cluster Alkyl-Cl -> Alkyl-F || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 6
- Priority score: 9
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-Cl -> Alkyl-F`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:6

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkyl-F
- Existing compound IDs missing in reaction slots: Alkyl-F
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Alkyl-F"
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

## Priority 23: Cluster R2CH-OH -> R2CH-Cl || events:intramolecular_likely

- Cluster count: 6
- Priority score: 9
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R2CH-OH -> R2CH-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:6

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Cl
- Existing compound IDs missing in reaction slots: R2CH-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R2CH-Cl"
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

## Priority 24: Cluster Ar-C=N|Ar-NR2|RCH2-OH -> Ar-Ar|Ar-AromN|RCH2-Br || events:intermolecular_or_multi_component

- Cluster count: 4
- Priority score: 9
- Action track: `propose_new_reaction_family`
- Motif signature: `Ar-C=N|Ar-NR2|RCH2-OH -> Ar-Ar|Ar-AromN|RCH2-Br`
- Event signature: `intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.3333)
- Reasons: mapping_warning:4, motif_outside_reaction_taxonomy:4, unknown_reaction_type:4

### Gap Analysis

- Missing reacted-slot motifs: Ar-C=N, Ar-NR2
- Missing formed-slot motifs: RCH2-Br
- Existing compound IDs missing in reaction slots: Ar-C=N, Ar-NR2, RCH2-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_035_intermolecular_or_multi_component_ar_c_n_ar_nr2_ar_ar_ar_aromn",
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

## Priority 25: Cluster Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|Ar-F|Ar-OH|R2CH-CO2H|R2CH-CO2R|RCH2-CO2R -> Alkyl-NHCOR|R2CH-CONH2|R2CH-CONHR|R2CH-NH2|RCH2-CO2H|RCH2-CONHR|R_acidic-H || events:amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component

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
  "id": "todo_family_023_amidation_like_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_alkyl_n_r_co2r_ar_alkyl_alkyl_nhcor_r2ch_conh2",
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

## Priority 26: Cluster R3C-Br -> R3C-F || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 5
- Priority score: 8
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R3C-Br -> R3C-F`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:5

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R3C-F
- Existing compound IDs missing in reaction slots: R3C-F
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R3C-F"
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

## Priority 27: Cluster Alkenyl-I -> Alkenyl-Br || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 5
- Priority score: 8
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-I -> Alkenyl-Br`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.65)
- Reasons: unknown_reaction_type:5

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Br
- Existing compound IDs missing in reaction slots: Alkenyl-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Suzuki_miyaura",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Alkenyl-Br"
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

## Priority 28: Cluster Alkyl-Hydrazine|RCH2-CHO|R_acidic-H -> CH3-NR2 || events:c_n_bond_formation|intermolecular_or_multi_component

- Cluster count: 3
- Priority score: 8
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-Hydrazine|RCH2-CHO|R_acidic-H -> CH3-NR2`
- Event signature: `c_n_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.35)
- Reasons: mapping_warning:3, motif_outside_reaction_taxonomy:3, unknown_reaction_type:3

### Gap Analysis

- Missing reacted-slot motifs: Alkyl-Hydrazine
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: Alkyl-Hydrazine
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_038_c_n_bond_formation_intermolecular_or_multi_component_alkyl_hydrazine_rch2_cho_ch3_nr2",
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

## Priority 29: Cluster none -> Ar-Br|R2CH-Br || events:intramolecular_likely

- Cluster count: 4
- Priority score: 7
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Ar-Br|R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Sandmeyer` (score=0.175)
- Reasons: unknown_reaction_type:4

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Br
- Existing compound IDs missing in reaction slots: R2CH-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_028_intramolecular_likely_reacted_ar_br_r2ch_br",
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

## Priority 30: Cluster none -> Alkenyl-Br|R2CH-Br|R3C-Br || events:intramolecular_likely

- Cluster count: 4
- Priority score: 7
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Br|R2CH-Br|R3C-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:4

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Br, R2CH-Br, R3C-Br
- Existing compound IDs missing in reaction slots: Alkenyl-Br, R2CH-Br, R3C-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_029_intramolecular_likely_reacted_alkenyl_br_r2ch_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Br",
      "R2CH-Br",
      "R3C-Br"
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

## Priority 31: Cluster Alkenyl-OR|Alkyl-Si* -> RCH2-CO2R || events:c_c_bond_formation|intermolecular_or_multi_component

- Cluster count: 4
- Priority score: 7
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-OR|Alkyl-Si* -> RCH2-CO2R`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Esterification` (score=0.35)
- Reasons: mapping_warning:4, unknown_reaction_type:4

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_030_c_c_bond_formation_intermolecular_or_multi_component_alkenyl_or_alkyl_si_rch2_co2r",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-OR",
      "Alkyl-Si*"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "RCH2-CO2R"
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

## Priority 32: Cluster Bn-Cl -> Bn-Br || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 4
- Priority score: 7
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Bn-Cl -> Bn-Br`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:4

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Bn-Br
- Existing compound IDs missing in reaction slots: Bn-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Bn-Br"
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

## Priority 33: Cluster Alkenyl-I -> Alkenyl-Cl || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 4
- Priority score: 7
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-I -> Alkenyl-Cl`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.65)
- Reasons: unknown_reaction_type:4

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Cl
- Existing compound IDs missing in reaction slots: Alkenyl-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Suzuki_miyaura",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Alkenyl-Cl"
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

## Priority 34: Cluster RCH2-OH -> RCH2-I || events:intramolecular_likely

- Cluster count: 4
- Priority score: 7
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `RCH2-OH -> RCH2-I`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:4

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-I
- Existing compound IDs missing in reaction slots: RCH2-I
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "RCH2-I"
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

## Priority 35: Cluster Alkyl-H|Ar-CHO|RCH2-OH -> Ar-Alkyl|RCH2-OR || events:c_o_bond_formation|intermolecular_or_multi_component

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
  "id": "todo_family_027_c_o_bond_formation_intermolecular_or_multi_component_alkyl_h_ar_cho_ar_alkyl_rch2_or",
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

## Priority 36: Cluster Alkenyl-OR|Alkyl-F|Alkyl-Si*|Ar-CHO -> Alkyl-CF3|Ar-Alkyl|R2CH-OH|R3C-CO2R|R3C-OH || events:c_c_bond_formation|intermolecular_or_multi_component

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
  "id": "todo_family_031_c_c_bond_formation_intermolecular_or_multi_component_alkenyl_or_alkyl_f_alkyl_cf3_ar_alkyl",
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

## Priority 37: Cluster CH3-OH|RCH2-CHO|R_acidic-H -> R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 3
- Priority score: 6
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `CH3-OH|RCH2-CHO|R_acidic-H -> R2CH-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.5667)
- Reasons: unknown_reaction_type:3

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
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

## Priority 38: Cluster none -> Alkenyl-Cl|Alkyl-Cl|R3C-Cl|RCH2-Cl || events:intramolecular_likely

- Cluster count: 3
- Priority score: 6
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Cl|Alkyl-Cl|R3C-Cl|RCH2-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:3

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Cl, Alkyl-Cl, R3C-Cl, RCH2-Cl
- Existing compound IDs missing in reaction slots: Alkenyl-Cl, Alkyl-Cl, R3C-Cl, RCH2-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_037_intramolecular_likely_reacted_alkenyl_cl_alkyl_cl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Cl",
      "Alkyl-Cl",
      "R3C-Cl",
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

## Priority 39: Cluster R3C-Br -> R3C-Cl || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 3
- Priority score: 6
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R3C-Br -> R3C-Cl`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:3

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R3C-Cl
- Existing compound IDs missing in reaction slots: R3C-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R3C-Cl"
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

## Priority 40: Cluster Ar-CF3 -> Alkyl-Br|Ar-Alkyl || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 3
- Priority score: 6
- Action track: `propose_new_reaction_family`
- Motif signature: `Ar-CF3 -> Alkyl-Br|Ar-Alkyl`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.175)
- Reasons: mapping_warning:2, unknown_reaction_type:3

### Gap Analysis

- Missing reacted-slot motifs: Ar-CF3
- Missing formed-slot motifs: Alkyl-Br
- Existing compound IDs missing in reaction slots: Alkyl-Br, Ar-CF3
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_041_intramolecular_likely_leaving_group_displacement_ar_cf3_alkyl_br_ar_alkyl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Ar-CF3"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkyl-Br",
      "Ar-Alkyl"
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

## Priority 41: Cluster Bn-Br -> Bn-I || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 3
- Priority score: 6
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Bn-Br -> Bn-I`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:3

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Bn-I
- Existing compound IDs missing in reaction slots: Bn-I
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Bn-I"
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

## Priority 42: Cluster Propargyl-OH -> Propargyl-Br || events:intramolecular_likely

- Cluster count: 3
- Priority score: 6
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Propargyl-OH -> Propargyl-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:3

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Propargyl-Br
- Existing compound IDs missing in reaction slots: Propargyl-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Propargyl-Br"
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

## Priority 43: Cluster Alkyl-H|R2CH-CHO|RCH2-OH|R_acidic-H -> RCH2-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 6
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-H|R2CH-CHO|RCH2-OH|R_acidic-H -> RCH2-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.5125)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1, unknown_reaction_type:1

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

## Priority 44: Cluster Alkyl-H|RCH2-CHO|RCH2-OH|R_acidic-H -> RCH2-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 6
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-H|RCH2-CHO|RCH2-OH|R_acidic-H -> RCH2-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.5125)
- Reasons: mapping_warning:1, motif_outside_reaction_taxonomy:1, unknown_reaction_type:1

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

## Priority 45: Cluster Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-CONHR -> Alkenyl-COR|Alkenyl-NR2|Ar-Hydrazide_NR|R2CH-NR2|R_acidic-H || events:c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

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
  "id": "todo_family_039_c_c_bond_formation_c_n_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_alkenyl_alkenyl_alkenyl_or_alkenyl_cor_alkenyl_nr2",
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

## Priority 46: Cluster Alkyl-CO2R|Alkyl-N(R)CO2R|Ar-Alkyl|R2CH-CO2H|RCH2-Br|RCH2-OH|R_acidic-H -> Alkyl-NRCOR|R2CH-CONR2|RCH2-OR || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

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
  "id": "todo_family_043_amidation_like_c_n_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkyl_co2r_alkyl_n_r_co2r_alkyl_nrcor_r2ch_conr2",
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

## Priority 47: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-I -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

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
  "id": "todo_family_044_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_co2r",
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

## Priority 48: Cluster Alkyl-NHCOR|Alkyl-NRCOR|Ar-Alkyl|Ar-CONHR|Bn-LactamN|R2CH-CONHR|R2CH-CONR2|RCH2-CO2H -> none || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

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
  "id": "todo_family_045_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkyl_nhcor_alkyl_nrcor_formed",
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

## Priority 49: Cluster Alkyl-COR|R2CH-OR|R_acidic-H -> R3C-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 2
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkyl-COR|R2CH-OR|R_acidic-H -> R3C-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Arylation_acidic_C_H` (score=0.4333)
- Reasons: mapping_warning:2, unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: R2CH-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: R2CH-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_047_c_o_bond_formation_intermolecular_or_multi_component_alkyl_cor_r2ch_or_r3c_or",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkyl-COR",
      "R2CH-OR",
      "R_acidic-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
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

## Priority 50: Cluster Ar-Alkyl|R2CH-OH -> Ar-Alkenyl|RCH2-Cl || events:intramolecular_likely

- Cluster count: 2
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkyl|R2CH-OH -> Ar-Alkenyl|RCH2-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.65)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Cl
- Existing compound IDs missing in reaction slots: RCH2-Cl
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
      "RCH2-Cl"
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

## Priority 51: Cluster none -> Alkenyl-Br|Alkyl-Br|R2CH-Br || events:intramolecular_likely

- Cluster count: 2
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Br|Alkyl-Br|R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Br, Alkyl-Br, R2CH-Br
- Existing compound IDs missing in reaction slots: Alkenyl-Br, Alkyl-Br, R2CH-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_049_intramolecular_likely_reacted_alkenyl_br_alkyl_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Br",
      "Alkyl-Br",
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

## Priority 52: Cluster none -> Alkenyl-Cl|R2CH-Cl|R3C-Cl || events:intramolecular_likely

- Cluster count: 2
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Cl|R2CH-Cl|R3C-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Cl, R2CH-Cl, R3C-Cl
- Existing compound IDs missing in reaction slots: Alkenyl-Cl, R2CH-Cl, R3C-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_050_intramolecular_likely_reacted_alkenyl_cl_r2ch_cl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Cl",
      "R2CH-Cl",
      "R3C-Cl"
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

## Priority 53: Cluster none -> Alkenyl-Cl|Alkyl-Cl|R2CH-Cl || events:intramolecular_likely

- Cluster count: 2
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Cl|Alkyl-Cl|R2CH-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Cl, Alkyl-Cl, R2CH-Cl
- Existing compound IDs missing in reaction slots: Alkenyl-Cl, Alkyl-Cl, R2CH-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_051_intramolecular_likely_reacted_alkenyl_cl_alkyl_cl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Cl",
      "Alkyl-Cl",
      "R2CH-Cl"
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

## Priority 54: Cluster Alkenyl-Br -> Alkenyl-I || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 2
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-Br -> Alkenyl-I`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.65)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-I
- Existing compound IDs missing in reaction slots: Alkenyl-I
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Suzuki_miyaura",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Alkenyl-I"
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

## Priority 55: Cluster R3C-I -> R3C-Br || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 2
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R3C-I -> R3C-Br`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R3C-Br
- Existing compound IDs missing in reaction slots: R3C-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R3C-Br"
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

## Priority 56: Cluster none -> Alkyl-F || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 2
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkyl-F`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkyl-F
- Existing compound IDs missing in reaction slots: Alkyl-F
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_055_intramolecular_likely_leaving_group_displacement_reacted_alkyl_f",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
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

## Priority 57: Cluster R2CH-Br -> R2CH-Cl || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 2
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R2CH-Br -> R2CH-Cl`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Cl
- Existing compound IDs missing in reaction slots: R2CH-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R2CH-Cl"
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

## Priority 58: Cluster R3C-F -> R3C-I || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 2
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R3C-F -> R3C-I`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R3C-I
- Existing compound IDs missing in reaction slots: R3C-I
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R3C-I"
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

## Priority 59: Cluster R2CH-F -> R2CH-Br || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 2
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R2CH-F -> R2CH-Br`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Br
- Existing compound IDs missing in reaction slots: R2CH-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R2CH-Br"
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

## Priority 60: Cluster RCH2-Br -> RCH2-I || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 2
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `RCH2-Br -> RCH2-I`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.65)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-I
- Existing compound IDs missing in reaction slots: RCH2-I
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "RCH2-I"
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

## Priority 61: Cluster none -> Alkenyl-Br|Ar-Br || events:intramolecular_likely

- Cluster count: 2
- Priority score: 5
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Br|Ar-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Sandmeyer` (score=0.175)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Br
- Existing compound IDs missing in reaction slots: Alkenyl-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_060_intramolecular_likely_reacted_alkenyl_br_ar_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Br",
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

## Priority 62: Cluster Propargyl-OH -> Propargyl-Cl || events:intramolecular_likely

- Cluster count: 2
- Priority score: 5
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Propargyl-OH -> Propargyl-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Propargyl-Cl
- Existing compound IDs missing in reaction slots: Propargyl-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "Propargyl-Cl"
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

## Priority 63: Cluster Alkenyl-Alkenyl|Alkenyl-OR|Alkyl-Si*|Ar-C=N -> Alkyl-N(R)CO2R|RCH2-CO2R|RCH2-COR|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component

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
  "id": "todo_family_052_c_c_bond_formation_intermolecular_or_multi_component_alkenyl_alkenyl_alkenyl_or_alkyl_n_r_co2r_rch2_co2r",
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

## Priority 64: Cluster Bn-COR -> Ar-Alkyl|Inorganic|R2CH-Br || events:intramolecular_likely

- Cluster count: 2
- Priority score: 4
- Action track: `expand_existing_reaction_type`
- Motif signature: `Bn-COR -> Ar-Alkyl|Inorganic|R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Reduction_carbonyl` (score=0.7667)
- Reasons: motif_outside_reaction_taxonomy:2

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Inorganic, R2CH-Br
- Existing compound IDs missing in reaction slots: R2CH-Br
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
      "Inorganic",
      "R2CH-Br"
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

## Priority 65: Cluster Alkyl-N(R)CO2R|Ar-Alkyl|Ar-Ar|R2CH-CO2H|R2CH-OR|RCH2-CO2H|RCH2-CO2R|RCH2-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|R2CH-CONH2|R2CH-CONHR|R2CH-CONR2|R2CH-NH2|RCH2-CHO|RCH2-CONHR || events:amidation_like|c_c_bond_formation|c_n_bond_formation|intermolecular_or_multi_component

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

## Priority 66: Cluster Alkyl-N(R)CO2R|Alkyl-NHCOR|Ar-Alkyl|Ar-OR|Ar-SO2NHR|R2CH-CONHR|R3C-OR|RCH2-CO2R -> R2CH-CO2H|R2CH-NH2|R2CH-OH|RCH2-CO2H|RCH2-OH || events:c_o_bond_formation|intermolecular_or_multi_component

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
  "id": "todo_family_063_c_o_bond_formation_intermolecular_or_multi_component_alkyl_n_r_co2r_alkyl_nhcor_r2ch_co2h_r2ch_nh2",
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

## Priority 67: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Bn-Br|Inorganic -> Alkenyl-CO2R|Ar-Alkyl|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

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
  "id": "todo_family_064_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_alkyl",
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

## Priority 68: Cluster Alkenyl-CO2H|Alkenyl-LactamN|Alkenyl-SR|Alkyl-N(R)CO2R|Alkyl-NRCOR|Ar-CO2H|Ar-NHCOR|Inorganic -> none || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

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
  "id": "todo_family_065_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_alkenyl_co2h_alkenyl_lactamn_formed",
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

## Priority 69: Cluster Alkenyl-CO2H|Alkyl-N(R)CO2R|Ar-CO2H|Inorganic|RCH2-Br -> Alkenyl-CO2R|Ar-CO2R|R2CH-NHR|RCH2-NHR || events:c_o_bond_formation|intermolecular_or_multi_component|leaving_group_displacement|ring_closure_or_annulation

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
  "id": "todo_family_066_c_o_bond_formation_intermolecular_or_multi_component_leaving_group_displacement_ring_closure_or_annulation_alkenyl_co2h_alkyl_n_r_co2r_alkenyl_co2r_ar_co2r",
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

## Priority 70: Cluster Alkyl-N(R)CO2R|Ar-CO2H|Ar-COCl|Ar-OH|Ar-SO2Cl|CH3-NH2 -> Alkyl-NHCOR|Alkyl-NRCOR|Ar-CONHR|Ar-CONR2|Ar-OSO2R || events:amidation_like|c_n_bond_formation|intermolecular_or_multi_component|leaving_group_displacement

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

## Priority 71: Cluster RCH2-CO2R -> Ar-COR|HeteroAr-Br|Indole|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

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

## Priority 72: Cluster Alkyl-CO2R -> Ar-COR|HeteroAr-Cl|Indole|R_acidic-H || events:c_c_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

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

## Priority 73: Cluster Ar-CO2R|AromN-H -> Pyrazole || events:c_n_bond_formation|intermolecular_or_multi_component

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
  "id": "todo_family_070_c_n_bond_formation_intermolecular_or_multi_component_ar_co2r_aromn_h_pyrazole",
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

## Priority 74: Cluster R3C-OR -> Alkyl-COR|R_acidic-H || events:none

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
  "id": "todo_family_072_none_r3c_or_alkyl_cor_r_acidic_h",
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

## Priority 75: Cluster Ar-CHO -> R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-CHO -> R2CH-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Minisci_acylation` (score=0.65)
- Reasons: mapping_warning:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Minisci_acylation",
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

## Priority 76: Cluster Ar-COR|R_acidic-H -> R3C-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Ar-COR|R_acidic-H -> R3C-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.35)
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
  "id": "todo_family_074_c_o_bond_formation_intermolecular_or_multi_component_ar_cor_r_acidic_h_r3c_or",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Ar-COR",
      "R_acidic-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
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

## Priority 77: Cluster RCH2-COR|R_acidic-H -> none || events:c_c_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 4
- Action track: `expand_existing_reaction_type`
- Motif signature: `RCH2-COR|R_acidic-H -> none`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Arylation_acidic_C_H` (score=0.85)
- Reasons: mapping_warning:1, missing_formed_motifs:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: (none)
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "Arylation_acidic_C_H",
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

## Priority 78: Cluster Alkenyl-CHO|CH3-OH -> R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-CHO|CH3-OH -> R2CH-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.675)
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
  "reaction_id": "Alkyl_Nucleophilic_Substitution",
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

## Priority 79: Cluster R2CH-OR|RCH2-COR|R_acidic-H -> R3C-OR || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `R2CH-OR|RCH2-COR|R_acidic-H -> R3C-OR`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Arylation_acidic_C_H` (score=0.4333)
- Reasons: mapping_warning:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: R2CH-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: R2CH-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_082_c_o_bond_formation_intermolecular_or_multi_component_r2ch_or_rch2_cor_r3c_or",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "R2CH-OR",
      "RCH2-COR",
      "R_acidic-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
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

## Priority 80: Cluster Alkenyl-CHO|RCH2-OR -> none || events:c_o_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Alkenyl-CHO|RCH2-OR -> none`
- Event signature: `c_o_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Reductive_amination` (score=0.425)
- Reasons: mapping_warning:1, missing_formed_motifs:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: RCH2-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: RCH2-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_083_c_o_bond_formation_intermolecular_or_multi_component_alkenyl_cho_rch2_or_formed",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Alkenyl-CHO",
      "RCH2-OR"
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

## Priority 81: Cluster R2CH-OR|RCH2-COR|R_acidic-H -> R3C-OR || events:c_c_bond_formation|intermolecular_or_multi_component

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `R2CH-OR|RCH2-COR|R_acidic-H -> R3C-OR`
- Event signature: `c_c_bond_formation|intermolecular_or_multi_component`
- Top taxonomy candidate: `Arylation_acidic_C_H` (score=0.4333)
- Reasons: mapping_warning:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: R2CH-OR
- Missing formed-slot motifs: (none)
- Existing compound IDs missing in reaction slots: R2CH-OR
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_084_c_c_bond_formation_intermolecular_or_multi_component_r2ch_or_rch2_cor_r3c_or",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "R2CH-OR",
      "RCH2-COR",
      "R_acidic-H"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
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

## Priority 82: Cluster none -> Acyl-Cl || events:intramolecular_likely

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
  "id": "todo_family_085_intramolecular_likely_reacted_acyl_cl",
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

## Priority 83: Cluster Alkenyl-Br|Ar-Alkenyl -> Alkyl-Cl|Ar-Alkyl || events:intramolecular_likely|leaving_group_displacement

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkenyl-Br|Ar-Alkenyl -> Alkyl-Cl|Ar-Alkyl`
- Event signature: `intramolecular_likely|leaving_group_displacement`
- Top taxonomy candidate: `Hydrogenation` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkyl-Cl
- Existing compound IDs missing in reaction slots: Alkyl-Cl
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
      "Alkyl-Cl"
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

## Priority 84: Cluster Ar-Alkyl|R3C-OH -> Ar-Alkenyl|RCH2-Cl || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkyl|R3C-OH -> Ar-Alkenyl|RCH2-Cl`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Cl
- Existing compound IDs missing in reaction slots: RCH2-Cl
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
      "RCH2-Cl"
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

## Priority 85: Cluster Allyl-OR|RCH2-OH -> R2CH-Br|RCH2-Br|RCH2-Cl|RCH2-OR || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Allyl-OR|RCH2-OH -> R2CH-Br|RCH2-Br|RCH2-Cl|RCH2-OR`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.4125)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: Allyl-OR
- Missing formed-slot motifs: R2CH-Br, RCH2-Br, RCH2-Cl
- Existing compound IDs missing in reaction slots: Allyl-OR, R2CH-Br, RCH2-Br, RCH2-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_090_intramolecular_likely_allyl_or_rch2_oh_r2ch_br_rch2_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Allyl-OR",
      "RCH2-OH"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "R2CH-Br",
      "RCH2-Br",
      "RCH2-Cl",
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

## Priority 86: Cluster R2CH-OH -> Ar-Alkenyl|RCH2-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R2CH-OH -> Ar-Alkenyl|RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: mapping_warning:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Br
- Existing compound IDs missing in reaction slots: RCH2-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "RCH2-Br"
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

## Priority 87: Cluster Allyl-OR|RCH2-OH -> R2CH-Cl|RCH2-Cl|RCH2-OR || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `Allyl-OR|RCH2-OH -> R2CH-Cl|RCH2-Cl|RCH2-OR`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Alkyl_Nucleophilic_Substitution` (score=0.4417)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: Allyl-OR
- Missing formed-slot motifs: R2CH-Cl, RCH2-Cl
- Existing compound IDs missing in reaction slots: Allyl-OR, R2CH-Cl, RCH2-Cl
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_092_intramolecular_likely_allyl_or_rch2_oh_r2ch_cl_rch2_cl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [
      "Allyl-OR",
      "RCH2-OH"
    ],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "R2CH-Cl",
      "RCH2-Cl",
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

## Priority 88: Cluster R3C-OH -> R2CH-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `R3C-OH -> R2CH-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `C_O_Coupling` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Br
- Existing compound IDs missing in reaction slots: R2CH-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `reaction_type_update`
```json
{
  "reaction_id": "C_O_Coupling",
  "reactant_slot_additions_todo": {
    "electrophile_or_substrate": []
  },
  "product_slot_additions_todo": {
    "product": [
      "R2CH-Br"
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

## Priority 89: Cluster Alkyl-Si* -> RCH2-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Alkyl-Si* -> RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Hiyama` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: RCH2-Br
- Existing compound IDs missing in reaction slots: RCH2-Br
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
    "product": [
      "RCH2-Br"
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

## Priority 90: Cluster Ar-Alkenyl -> Ar-Alkyl|Ar-Br|R2CH-Br|R3C-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `review_constraints_or_product_patterns`
- Motif signature: `Ar-Alkenyl -> Ar-Alkyl|Ar-Br|R2CH-Br|R3C-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Hydrogenation` (score=0.65)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Br, R3C-Br
- Existing compound IDs missing in reaction slots: R2CH-Br, R3C-Br
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
      "R2CH-Br",
      "R3C-Br"
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

## Priority 91: Cluster none -> Ar-Br|R2CH-Br|RCH2-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Ar-Br|R2CH-Br|RCH2-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Sandmeyer` (score=0.1167)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Br, RCH2-Br
- Existing compound IDs missing in reaction slots: R2CH-Br, RCH2-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_096_intramolecular_likely_reacted_ar_br_r2ch_br",
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

## Priority 92: Cluster none -> Alkenyl-Br|Alkyl-Br|R3C-Br || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Br|Alkyl-Br|R3C-Br`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Br, Alkyl-Br, R3C-Br
- Existing compound IDs missing in reaction slots: Alkenyl-Br, Alkyl-Br, R3C-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_097_intramolecular_likely_reacted_alkenyl_br_alkyl_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Br",
      "Alkyl-Br",
      "R3C-Br"
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

## Priority 93: Cluster none -> Alkenyl-Cl|R2CH-Cl|R_acidic-H || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Cl|R2CH-Cl|R_acidic-H`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Cl, R2CH-Cl, R_acidic-H
- Existing compound IDs missing in reaction slots: Alkenyl-Cl, R2CH-Cl, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_098_intramolecular_likely_reacted_alkenyl_cl_r2ch_cl",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Cl",
      "R2CH-Cl",
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

## Priority 94: Cluster none -> Alkenyl-Br || events:intramolecular_likely|ring_closure_or_annulation

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-Br`
- Event signature: `intramolecular_likely|ring_closure_or_annulation`
- Top taxonomy candidate: `Suzuki_miyaura` (score=0.0)
- Reasons: mapping_warning:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: Alkenyl-Br
- Existing compound IDs missing in reaction slots: Alkenyl-Br
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_099_intramolecular_likely_ring_closure_or_annulation_reacted_alkenyl_br",
  "aliases": [],
  "description": "TODO: new reaction family candidate from unresolved cluster.",
  "reactants": {
    "electrophile_or_substrate": [],
    "nucleophile_or_partner": []
  },
  "products": {
    "product": [
      "Alkenyl-Br"
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

## Priority 95: Cluster none -> Alkenyl-COR|R2CH-Br|RCH2-Br|R_acidic-H || events:intramolecular_likely

- Cluster count: 1
- Priority score: 4
- Action track: `propose_new_reaction_family`
- Motif signature: `none -> Alkenyl-COR|R2CH-Br|RCH2-Br|R_acidic-H`
- Event signature: `intramolecular_likely`
- Top taxonomy candidate: `Oxidation_alcohol` (score=0.0875)
- Reasons: mapping_warning:1, unknown_reaction_type:1

### Gap Analysis

- Missing reacted-slot motifs: (none)
- Missing formed-slot motifs: R2CH-Br, RCH2-Br, R_acidic-H
- Existing compound IDs missing in reaction slots: R2CH-Br, RCH2-Br, R_acidic-H
- Unknown compound IDs: (none)

### Patch Stubs

- `new_reaction_family`
```json
{
  "id": "todo_family_100_intramolecular_likely_reacted_alkenyl_cor_r2ch_br",
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

## Priority 96: Cluster CH3-OH|R2CH-CHO|R_acidic-H -> Ar-Alkyl|Benzofuran|HeteroAr-H|R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

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
  "id": "todo_family_075_c_o_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_ch3_oh_r2ch_cho_ar_alkyl_benzofuran",
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

## Priority 97: Cluster Alkyl-H|Ar-CHO|R2CH-OH -> Ar-Alkyl|R2CH-OR || events:c_o_bond_formation|intermolecular_or_multi_component

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
  "id": "todo_family_076_c_o_bond_formation_intermolecular_or_multi_component_alkyl_h_ar_cho_ar_alkyl_r2ch_or",
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

## Priority 98: Cluster Alkenyl-H|R2CH-OH|RCH2-OH|Unclassified-Reactant -> RCH2-CHO|RCH2-OR|R_acidic-H || events:c_c_bond_formation|c_o_bond_formation|intermolecular_or_multi_component|ring_closure_or_annulation

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
  "id": "todo_family_080_c_c_bond_formation_c_o_bond_formation_intermolecular_or_multi_component_ring_closure_or_annulation_alkenyl_h_r2ch_oh_rch2_cho_rch2_or",
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

## Priority 99: Cluster Ar-CO2H -> Ar-COCl|Inorganic || events:intramolecular_likely

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

## Priority 100: Cluster Ar-CO2H|Inorganic -> Ar-COCl || events:intermolecular_or_multi_component

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
