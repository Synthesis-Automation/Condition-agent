# Chem Assistant Summary (Featurization/Analysis)

This package exposes ChemTools featurization, analysis, and reagent registry tools
with a minimal LangGraph agent and CLI. Protocol tooling is intentionally removed;
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
- Reagent registry lookup
- Knowledge base retrieval (RAG)
