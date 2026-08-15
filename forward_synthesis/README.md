# Forward Synthesis

`forward_synthesis` predicts structurally possible products from starting
materials and optionally audits those products against a resolved condition
recipe. It is chemistry-first and type-agnostic: graph operators, reactive
sites, atom correspondence, and normalized edits are authoritative; named
reaction families are optional annotations.

The package provides two workflows:

- `predict_products(...)` performs target-blind product generation.
- `assess_proposed_step(...)` separately performs targeted operator replay and
  target-blind competition analysis for a retrosynthetic proposal.

Operators are projected from the generic operator library but are independently
admitted only after forward source reconstruction and reverse precursor
recovery. Forward execution uses explicit multi-component matching rather than
reversing the single-target RDChiral execution path.

## Python example

```python
from core_retrosynthesis import load_generic_library
from forward_synthesis import build_forward_library, predict_products

generic = load_generic_library("generic_library.json.gz")
library = build_forward_library(generic)
result = predict_products("CCBr.N", library, top_k=10)

for candidate in result.candidates:
    print(candidate.rank, candidate.product_smiles, candidate.score)
```

Scores are deterministic priorities, not calibrated probabilities or yield
predictions. Without supplied conditions, results are possible products rather
than claims about the major product.

## CLI

```powershell
python -m forward_synthesis build-library GENERIC_LIBRARY FORWARD_LIBRARY
python -m forward_synthesis predict FORWARD_LIBRARY "CCBr.N" --top-k 10
python -m forward_synthesis assess-step FORWARD_LIBRARY "CCBr.N" "CCN"
```

Use `--recipe resolved_recipe.json` to apply the condition recommender's
canonical compatibility rules. A recipe compatibility result is not evidence
that the product was experimentally observed.

## Evaluation

`evaluate_product_hidden_replay(...)` hides held-out products and reports exact
product rank. Reference overlap with the operator library is rejected by
default. Enumerated alternatives are not treated as experimental negatives.

Retrosynthesis integration is advisory: route forward assessments currently
have `route_ranking_impact="none_advisory_only"` until broad, leakage-controlled
competition and false-warning benchmarks support stronger use.
