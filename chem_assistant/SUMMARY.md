# Chem Assistant Summary (Featurization/Analysis)

This package exposes ChemTools featurization and analysis tools with a minimal
LangGraph agent and CLI. Recommendation, protocol, HTE, and reagent tooling are
intentionally removed.

Key modules:

- chem_assistant/chemtools_wrapper.py: LangChain tools for featurizers/analysis
- chem_assistant/chemtools_agent.py: LLM agent wrapper
- chem_assistant/chemtools_cli.py: Interactive CLI

Tool categories:

- Analysis and normalization
- Reaction detection
- Unified featurizers
- Motif-based featurizers
- Reaction-pair featurizers
- Calculable features
- MolPipeline descriptors
