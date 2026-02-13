# GPT-5.2 Reaction Analysis PoC (Isolated)

This folder is a standalone proof-of-concept for a "thinking-like" reaction analyzer.
It does not modify routing, taxonomy, or recommendation logic in the main system.

## What it does

- Parses single- or multi-component reaction SMILES (`A.B>>C.D`)
- Builds deterministic evidence:
  - core-pair formula deltas and whole-side deltas
  - principal transformed pair selection via pairwise MCS
  - functional signals (aryl halide, N-N motif, aromatic ring nitrogens)
- Ranks mechanistic hypotheses
- Returns visible reasoning steps
- Optionally asks an LLM (`gpt-5.2`) to refine the final call

## Run

```bash
python -m poc_gpt52_reaction.cli "Clc1ncc(-c2ccccc2)cn1>>NN=c1ncc(-c2ccccc2)c[nH]1"
```

JSON output:

```bash
python -m poc_gpt52_reaction.cli "Clc1ncc(-c2ccccc2)cn1>>NN=c1ncc(-c2ccccc2)c[nH]1" --json
```

Optional LLM refinement:

```bash
python -m poc_gpt52_reaction.cli "Clc1ncc(-c2ccccc2)cn1>>NN=c1ncc(-c2ccccc2)c[nH]1" --use-llm --model gpt-5.2
```

Multi-reactant example:

```bash
python -m poc_gpt52_reaction.cli "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1"
```

## Notes

- Deterministic analysis works without API keys.
- LLM refinement requires `OPENAI_API_KEY` (or provider-specific key in `llmtools`).
- This PoC is intentionally simple and can later be migrated into `llmtools` once validated.
