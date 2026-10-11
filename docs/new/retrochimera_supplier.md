# Optional RetroChimera single-step supplier

RetroChimera is an experimental additional proposal source. `disconnect_target`
continues to use the existing operator supplier. `predict_retrochimera` calls a
separately installed model; use either operation or both. Their rankings remain
separate. No new multi-step planner or automatic supplier selection is introduced.

## Configure an existing installation

No installation or checkpoint download occurs through this integration. Keep the
ML dependencies in their own environment. The copied checkout is not a runtime
dependency of the main Python environment.

```powershell
python -m chem_coworker.scientific_workspace.adapters.retrochimera_config --checkout retrochimera --python retrochimera/.venv/Scripts/python.exe --model-dir retrochimera/checkpoints --output results/ai_native/retrochimera.local.json
```

This hashes the interpreter, supplier Python sources and all checkpoint assets.
Register the generated file as `retrochimera_config` in the workspace artifact
configuration, alongside any existing processed dataset. Create a **new**
investigation. Baseline capture expands the pinned dependency list, checks its
hashes and subsequently detects changes. The worker independently checks hashes
when loading, records dependency versions and checks file metadata before reuse.
Changing the environment requires a new configuration/investigation. The runtime
records installed versions; it does not provide an environment lockfile.

For the web workspace, add the configuration path to a copy of your normal
artifact JSON, then restart with that file and create a fresh conversation:

```json
{"retrochimera_config": "results/ai_native/retrochimera.local.json"}
```

```powershell
python -m app.web_api --scientific-chat --port 8011 --chat-artifacts PATH_TO_MERGED_ARTIFACT_JSON
```

Retain the existing dataset entries in the merged JSON to compare suppliers and
inspect precedents. Model assets remain local generated dependencies, not files
to commit with the wrapper.

## Use in Tools + guidance or Normal mode

The operation appears in the standard capability catalogue. Tools + guidance is
useful for exploratory supplier comparison; Normal also provides structured
route authoring. An agent may compare the existing `disconnect_target` results
with these predictions, expand a selected precursor, and review critical steps.

```python
event = w.run("predict_retrochimera", {
    "target_smiles": target_smiles,
    "num_results": 10,
    "timeout_seconds": 300,
    "consensus_only": False,
})
print(w.call_summary(event))
```

For successive calls use the persistent scientific console:

```powershell
python -u -m chem_coworker.scientific_workspace.console WORKSPACE
```

Send JSONL requests such as:

```json
{"operation":"predict_retrochimera","arguments":{"target_smiles":"COC(=O)c1ccccc1OC(C)=O","num_results":10,"timeout_seconds":300},"reason":"Compare additional precursor hypotheses with the existing supplier"}
```

The console owns one serial worker and reuses the loaded model. Ordinary Python
`w.run` calls use a temporary worker and close it after the request. Python
clients can explicitly own a `RetroChimeraWorker` context and pass it to
`ScientificOperations(store, retrochimera_worker=worker)`. EOF/quit closes the
console worker. Timeout or failure terminates its process tree; the next call
starts fresh. Deadlines include model loading and hashing; process cleanup can
take additional time. Diagnostics remain under the investigation directory.

## Output and review

Results retain canonical precursor multisets, stereochemistry, supplier ranks,
component ranks (zero-based), model scores, shortlist-normalized probabilities,
checkpoint/configuration identity and runtime versions. Equivalent proposals are
deduplicated within the supplier; duplicate ranks and rejected proposals remain
auditable. Ranking is not renumbered after rejection. Cross-supplier comparison
must retain both sources and must not average their incompatible scores.

`structures_parsed_reaction_unverified` means only that structures parse and the
requested product matches. It does not imply a valid transformation, atom
correspondence, selectivity, conditions, experimental support or purchasability.
No taxonomy operator, observed signature, evidence confidence or precedent ID is
invented. Missing operator coverage does not automatically discard a prediction.
Use `assess_retro_validity` with an explicit proposed reaction (see its help),
then inspect relevant precedents. Supplier candidate IDs are not existing
operator realization IDs. Model probabilities are not success rates or yields.

Prefer 5–10 candidates. The bounded maximum is 50 for deliberate coverage probes;
larger lists require stronger review. Consensus mode intersects the component
models and can remove unique useful hypotheses. Successful execution with no
candidates is distinct from timeout/error and proves no chemical impossibility.

This integration and its development pilot do not establish superiority over the
existing supplier, independent chemist acceptance or untouched evaluation gates.

Upstream references: [RetroChimera](https://github.com/microsoft/retrochimera)
and [Syntheseus](https://microsoft.github.io/syntheseus/stable/).
