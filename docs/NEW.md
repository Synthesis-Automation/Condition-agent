# Strengthen the system prompt with stepwise workflows and “evidence first” output (update chemtools_agent.py)

Add more tools for evidence retrieval (protocols, unified recommender, rule DB, reagent lookup) and expose them in chemtools_wrapper.py so the agent can cross‑validate.
Add a lightweight RAG layer over your protocol/rule/HTE docs so it can cite specifics (e.g., index docs/ + data/protocol_db_v2/ and expose a retriever tool).
Use structured outputs (Pydantic response format) for condition recommendations to reduce hallucination and make responses consistent.

not done yet
Add post‑analysis self‑check: compare detection vs analysis vs HTE for agreement and flag conflicts (extend the preflight in chemtools_agent.py).

# ChemTools Toolkit

If you ever want a clean rebuild, delete results/hte_cache.

# System design

A-S

Ar-CHO, Ar-Cl

Ar-B(OH)2

Ar-COOR vs. ArCO-OR

[Alkyl-NHR]

Ar2N-H

R3N

Ar-SH vs. ArS-H

SCN-H