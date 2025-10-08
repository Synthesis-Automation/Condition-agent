Input: reaction SMILES, role taxonomy hits, substrate flags (steric/electronic), historical runs

Step 1 (Prune): Embed {substrate + candidate ConditionCore/Tail} with a chemistry model; keep top-K by cosine to literature/rules; optional PCA to reduce dim.
Step 2 (Model): Fit GP/DKL on prior data (features: role one-hots + substrate descriptors) with prior mean from rules.
Step 3 (Acquisition): Weighted-EI or Knowledge-Gradient over the pruned set; select a batch (qKG/qEI) for HTE.
Step 4 (Execute): Orchestrator (function-calling) dispatches runs to liquid handler, reactors, analytics.
Step 5 (Monitor): Camera/HPLC/TLC signals → VLM + heuristics for anomaly tags; if flagged, auto-apply safe “fix templates” (e.g., add 5% MeCN, warm to 60 °C, switch base family).
Step 6 (Update): Log to SQL; retrain/update posterior; iterate.
