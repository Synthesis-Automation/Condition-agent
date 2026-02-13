# GPT-5.2 General Reaction Analyzer PoC (v2)

This PoC is isolated from the core pipeline and is designed for general-purpose reaction typing.

## Design

1. Deterministic evidence extraction (RDKit + principal pair + stoichiometric deltas)
2. Taxonomy-driven candidate generation using `chemtools.detection.detect_reaction_type()`
3. Optional constrained LLM reranking (`gpt-5.2`) that can only choose from deterministic candidates (or `unknown`)
4. Validator gate to prevent invalid/non-taxonomy outputs

## Why v2

The previous PoC had an SNAr-specific scoring branch. v2 removes family-specific hardcoding and uses taxonomy candidates as the backbone, with LLM as an optional reranker/explainer.

## Run

```bash
python -m poc_gpt52_reaction_v2.cli "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1"
```

JSON output:

```bash
python -m poc_gpt52_reaction_v2.cli "CCO>>CC=O" --json
```

Optional LLM rerank:

```bash
python -m poc_gpt52_reaction_v2.cli "CCO>>CC=O" --use-llm --model gpt-5.2
```

Benchmark v1 vs v2:

```bash
python -m poc_gpt52_reaction_v2.benchmark
```

JSON report:

```bash
python -m poc_gpt52_reaction_v2.benchmark --json
```

## Notes

- Works without API keys in deterministic mode.
- LLM mode requires `OPENAI_API_KEY` (or provider-specific key supported by `llmtools`).
- If deterministic evidence is weak or inconsistent, output falls back to `unknown`.
