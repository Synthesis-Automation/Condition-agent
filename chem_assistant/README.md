# Chem Assistant (Featurization/Analysis)

This package exposes ChemTools featurization and analysis functions as LangChain tools
and provides a lightweight agent + CLI for interactive use.

## Quick start

```powershell
python -m chem_assistant.chemtools_cli
```

```python
from chem_assistant.chemtools_agent import ChemToolsAgent

agent = ChemToolsAgent()
print(agent.run("Featurize molecule: c1ccccc1O"))
```

## Tool groups

### Analysis

- analysis_normalize_smiles
- analysis_normalize_reaction
- analysis_analyze_reaction
- analysis_classify_reactant_smiles
- analysis_classify_reactant_category
- analysis_classify_reactant_group
- analysis_classify_reactant_batch
- analysis_get_reactant_category_matches
- analysis_get_all_reactant_matches
- analysis_normalize_reactant_identifier
- analysis_normalize_reaction_type
- analysis_resolve_reaction_family
- analysis_canonical_family_label
- analysis_classify_reactants_with_context
- analysis_reactant_summary

### Reaction detection

- detection_detect_reaction_types
- detection_detect_reaction_types_from_smiles
- detection_detect_motif_ids_from_smiles

### Featurizers

- unified_featurize_molecule
- unified_featurize_reaction
- motif_featurize_molecule
- motif_featurize_reaction
- reaction_pair_featurize_pair
- reaction_pair_featurize_flat

### Calculable features

- calculable_detect_all_features
- calculable_get_reactant_type_features
- calculable_classify_reactant_smiles

### HTE recommendations

- hte_recommend_conditions
- hte_database_stats

### Knowledge base (RAG)

- rag_search

## Notes

- RDKit is required for SMARTS-driven detection and descriptor computation.
- This package focuses on featurization/analysis; protocol and reagent tools remain out of scope.
- HTE recommendation tools are available when the HTE database is present.
- RAG search reads curated notes from `knowledge_base/`.
