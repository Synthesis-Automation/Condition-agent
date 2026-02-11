# Reaction Reliability Implementation Plan

## Objective
Improve reaction typing reliability while maintaining stable reaction-key quality, using taxonomy-driven updates and deterministic diagnostics.

## Current Baseline (Sampled Benchmark)
- Dataset: `examples/test_reactions_random_sampled.csv`
- Latest report: `results/reaction_key_quality_report.random_sampled.after_taxonomy_patch_v2.json`
- Metrics:
  - `unknown_reaction_type`: `0.4725`
  - `low_reaction_key_quality`: `0.0174`
  - `empty_reaction_key`: `0.0058`

## Phase 1: Patch High-Impact Unknown Clusters
1. Source cluster queue from:
   - `results/unknown_clusters.random_sampled.after_taxonomy_patch_v2.json`
   - `results/taxonomy_expansion_candidates.random_sampled.after_taxonomy_patch_v2.json`
2. Prioritize top unresolved families:
   - oxime/hydrazone formation patterns
   - clemmensen-like reductions
   - curtius-related records
   - decarboxylative arylation edge motifs
   - thiol-ene/thiol-yne residual gaps
3. Apply only taxonomy-driven changes in:
   - `chemtools/taxonomy/data/reaction_types.v4.0.json`
   - `chemtools/taxonomy/data/compound_logic.json`
4. Acceptance criteria:
   - unknown rate decreases on sampled benchmark
   - no increase in `low_reaction_key_quality` or `empty_reaction_key`

## Phase 2: Add Multistep/Underspecified Data Bucket
1. Extend diagnostics to explicitly tag:
   - `none -> none`
   - missing product side
   - malformed reaction records
2. Keep reaction typing conservative for these records.
3. Add separate counts in reports and gates so taxonomy progress is not masked by data quality issues.
4. Acceptance criteria:
   - report contains dedicated multistep/underspecified section
   - gate supports optional exclusion or separate threshold for this bucket

## Phase 3: CI-Style Reliability Gate
1. Standardize gate command with fixed thresholds and deltas.
2. Run on every taxonomy change.
3. Suggested thresholds:
   - max unknown rate: `0.65` (tighten over time)
   - max low-quality rate: `0.03`
   - max empty-key rate: `0.01`
   - max unknown delta vs baseline: `0.20` (tighten over time)
4. Acceptance criteria:
   - gate output persisted in `results/`
   - gate passes before merge

## Phase 4: Curated Benchmark Pack
1. Build two fixed subsets from sampled data:
   - single-step clean subset
   - complex/multistep subset
2. Track metrics separately for both subsets.
3. Acceptance criteria:
   - stable benchmark files committed under `examples/` or `data/` (small only)
   - reports show per-subset metrics

## Phase 5: Iterative Triage Loop
1. Run:
   - `scripts/triage_unknown_clusters.py`
   - taxonomy patch
   - `scripts/reaction_key_quality_report.py`
   - `scripts/reaction_reliability_gate.py`
2. Update candidate ranking each iteration:
   - `scripts/suggest_taxonomy_expansion_candidates.py`
3. Acceptance criteria:
   - measurable unknown-rate drop each iteration
   - no quality metric regressions

## Target Milestone
- Reduce sampled benchmark unknown rate from `0.4725` to `<0.35`
- Keep:
  - `low_reaction_key_quality <= 0.02`
  - `empty_reaction_key <= 0.01`

## Standard Execution Commands
```powershell
# 1) Unknown cluster triage
python scripts/triage_unknown_clusters.py --input examples/test_reactions_random_sampled.csv --output-json results/unknown_clusters.current.json --output-samples-csv results/unknown_cluster_samples.current.csv --top-clusters 80 --sample-per-cluster 5

# 2) Taxonomy candidate suggestion
python scripts/suggest_taxonomy_expansion_candidates.py --input results/unknown_clusters.current.json --output-json results/taxonomy_expansion_candidates.current.json --output-csv results/taxonomy_expansion_candidates.current.csv --top-clusters 80 --top-candidates 5

# 3) Reliability report
python scripts/reaction_key_quality_report.py --input examples/test_reactions_random_sampled.csv --output results/reaction_key_quality_report.current.json

# 4) Reliability gate
python scripts/reaction_reliability_gate.py --report results/reaction_key_quality_report.current.json --baseline results/reaction_key_quality_report.test_reactions.json --max-unknown-rate 0.65 --max-low-quality-rate 0.03 --max-empty-key-rate 0.01 --max-unknown-delta 0.20 --max-low-quality-delta 0.02 --max-empty-key-delta 0.01 --output-json results/reaction_reliability_gate.current.json
```

