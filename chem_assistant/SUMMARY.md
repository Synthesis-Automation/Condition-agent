# Chem Assistant Summary (Featurization/Analysis)

This package exposes ChemTools featurization and analysis tools with a minimal
LangGraph agent and CLI. Protocol and reagent tooling are intentionally removed;
HTE recommendations are now available via the agent tool set.

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
- HTE recommendations
- Knowledge base retrieval (RAG)
