import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import type { Capabilities, FragmentGuidedRetroRequest, FragmentGuidedRetroResult, FragmentRetroArm } from '../api/types'
import { FragmentResearch, FragmentResearchOptions, useFragmentResearch } from './FragmentResearch'
import { PrecedentDiscovery, usePrecedentDiscovery } from './PrecedentDiscovery'
import { ReactionEditor } from './ReactionEditor'
import { ReactionImage } from './ReactionImage'
import './fragment-guided-retro.css'

const EXAMPLE = 'Fc(cn1)cc2c1c(c3ccccc3OC)n[nH]2'
const readable = (value: string) => value.replaceAll('_', ' ').toLowerCase()

export function useFragmentGuidedRetro(active: boolean) {
  const [workflow, setWorkflow] = useState<'discovery' | 'manual' | 'automatic'>('discovery')
  const discovery = usePrecedentDiscovery(active && workflow === 'discovery')
  const research = useFragmentResearch(active && workflow === 'manual')
  const [target, setTarget] = useState('')
  const [library, setLibrary] = useState<FragmentGuidedRetroRequest['library_mode']>('full')
  const [queryLimit, setQueryLimit] = useState(3)
  const [focusLimit, setFocusLimit] = useState(3)
  const [topK, setTopK] = useState(3)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const [result, setResult] = useState<FragmentGuidedRetroResult | null>(null)
  const pending = useRef<AbortController | null>(null)
  useEffect(() => () => pending.current?.abort(), [])
  useEffect(() => {
    if (!active) {
      pending.current?.abort()
      setBusy(false)
      setResult(null)
      setError('')
    }
  }, [active])
  const reset = () => { setResult(null); setError(''); research.invalidateTransfer() }
  const changeTarget = (value: string) => {
    pending.current?.abort()
    setBusy(false)
    setTarget(value)
    reset()
  }
  const run = async () => {
    pending.current?.abort()
    const controller = new AbortController()
    pending.current = controller
    reset()
    setBusy(true)
    try {
      const next = await api.fragmentGuidedRetro({ target_smiles: target.trim(),
        library_mode: library, query_limit: queryLimit, max_focus_bonds: focusLimit,
        top_k: topK }, controller.signal)
      if (!controller.signal.aborted) setResult(next)
    } catch (reason) {
      if (!controller.signal.aborted) setError(reason instanceof Error ? reason.message : 'Evaluation failed')
    } finally {
      if (!controller.signal.aborted) setBusy(false)
    }
  }
  const exportResult = () => {
    if (workflow === 'discovery') { discovery.exportResult(); return }
    if (workflow === 'manual') { research.exportHistory(); return }
    if (!result) return
    const url = URL.createObjectURL(new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' }))
    const link = document.createElement('a')
    link.href = url
    link.download = 'fragment_guided_retrosynthesis.json'
    link.click()
    setTimeout(() => URL.revokeObjectURL(url), 1000)
  }
  return { workflow, setWorkflow, discovery, research, target, changeTarget, library, setLibrary, queryLimit, setQueryLimit,
    focusLimit, setFocusLimit, topK, setTopK, busy, error, setError, result, reset, run, exportResult }
}

type State = ReturnType<typeof useFragmentGuidedRetro>

export function FragmentGuidedRetroOptions({ state, capabilities }: { state: State; capabilities: Capabilities | null }) {
  const available = capabilities?.fragment_search && capabilities.retrosynthesis_library_modes?.[state.library]?.library_available
  return <div className="analysis-options fragment-retro-options">
    <div className="option-grid fragment-retro-primary">
    <label><span>Workflow</span><select aria-label="Workflow" value={state.workflow} disabled={state.busy || state.research.busy || state.discovery.busy} onChange={event => { state.setWorkflow(event.target.value as State['workflow']); state.reset() }}><option value="discovery">Find synthesis precedents</option><option value="manual">Refine search / assisted research</option><option value="automatic">Automatic POC comparison</option></select></label>
    <div className="feature-mode-note"><strong>{state.workflow === 'discovery' ? 'Automatic core-first discovery' : state.workflow === 'manual' ? 'Assisted fragment research' : 'Experimental single-step comparison'}</strong><span>{state.workflow === 'discovery' ? 'Enter a target to find related construction precedents. Queries and refinements are selected automatically.' : state.workflow === 'manual' ? 'Search a chosen core, inspect source reactions, and assess construction precedents against your target.' : 'Compare source transfers and fragment-guided proposals with an unrestricted baseline.'}</span></div>
    </div>
    {state.workflow === 'discovery' ? null : state.workflow === 'manual' ? <>
      <FragmentResearchOptions state={state.research}>
        <label><span>Construction bonds</span><select value={state.focusLimit} onChange={event => { state.setFocusLimit(Number(event.target.value)); state.reset() }}>{[1, 2, 3, 4, 5].map(value => <option key={value}>{value}</option>)}</select></label>
        <label><span>Proposals per search arm</span><select value={state.topK} onChange={event => { state.setTopK(Number(event.target.value)); state.reset() }}>{[1, 2, 3, 5, 10].map(value => <option key={value}>{value}</option>)}</select></label>
      </FragmentResearchOptions>
      {capabilities && !capabilities.fragment_search && <p className="alert caution">Fragment index unavailable; suggestions and query editing remain available.</p>}
    </> : <>
    <fieldset disabled={state.busy} className="fragment-fields option-grid fragment-retro-grid">
      <label><span>Fragment queries</span><select value={state.queryLimit} onChange={event => { state.setQueryLimit(Number(event.target.value)); state.reset() }}>{[1, 2, 3, 4, 5].map(value => <option key={value}>{value}</option>)}</select></label>
      <label><span>Construction bonds</span><select value={state.focusLimit} onChange={event => { state.setFocusLimit(Number(event.target.value)); state.reset() }}>{[1, 2, 3, 4, 5].map(value => <option key={value}>{value}</option>)}</select></label>
      <label><span>Proposals per search arm</span><select value={state.topK} onChange={event => { state.setTopK(Number(event.target.value)); state.reset() }}>{[1, 2, 3, 5, 10].map(value => <option key={value}>{value}</option>)}</select></label>
      <button className="button quiet" type="button" onClick={() => state.changeTarget(EXAMPLE)}>Fragment retro example</button>
    </fieldset>
    <details className="advanced-options"><summary>Search details</summary><div><p className="fragment-note fragment-retro-wide">Queries follow deterministic fragment suggestion order. Bonds rank by independent reference support, with atom IDs breaking ties. Chemistry validation is mandatory.</p></div></details>
    {capabilities && !available && <div className="alert warning">Requires a prepared fragment index and the selected operator library. Configure the fragment index with --fragment-index or FRAGMENT_PRECEDENT_INDEX; configure operator libraries with CORE_RETROSYNTHESIS_LIBRARY_ROOT.</div>}
    {state.error && <div className="alert error" role="alert">{state.error}</div>}
    </>}
  </div>
}

function Arm({ title, result, additional }: { title: string; result: FragmentRetroArm; additional: Set<string> }) {
  return <section className="fragment-retro-arm">
    <h3>{title}</h3>
    {!result.candidates.length && <p>{result.status === 'no_admitted_source_operators' ? 'No source operators passed whole-reaction admission.' : 'No verified proposals returned within this search budget.'}</p>}
    {result.candidates.map((candidate, index) => <article className="fragment-retro-candidate" key={`${candidate.proposed_reaction_smiles}:${index}`}>
      <h4>Proposal {index + 1}{additional.has(candidate.precursor_smiles) && <span className="fragment-retro-badge">Additional to baseline</span>}</h4>
      <ReactionImage smiles={candidate.proposed_reaction_smiles} label={`${title} proposal ${index + 1}`} />
      <p>Signature: <strong>{readable(candidate.forward_validation_status)}</strong> · Context: {candidate.abstraction_level}{candidate.bond_focus_check && <> · Construction bond: <strong>{candidate.bond_focus_check.status}</strong></>}</p>
      <p>Operator precedents: {candidate.precedent_reaction_ids.join(', ') || 'Unavailable'}</p>
      <details><summary>Precursor SMILES and validation evidence</summary><code>{candidate.precursor_smiles}</code><pre>{JSON.stringify(candidate, null, 2)}</pre></details>
    </article>)}
    <details><summary>Search work and diagnostics</summary><pre>{JSON.stringify(result.diagnostics, null, 2)}</pre></details>
  </section>
}

function Evidence({ result, showQueries = true }: { result: FragmentGuidedRetroResult; showQueries?: boolean }) {
  const comparison = result.transfers.comparison
  const additional = new Set(comparison.additional_guided_precursor_sets)
  const baselineRequested = comparison.baseline_requested !== false
  return <section className="results-card fragment-results fragment-retro-results" aria-label="Fragment-guided results">
    <div className="section-heading"><h2>{showQueries ? 'Fragment-guided comparison' : 'Selected precedent transfer assessment'}</h2><span>{result.execution.elapsed_seconds.toFixed(1)} s</span></div>
    <div className="fragment-retro-summary">
      <div><strong>{baselineRequested ? comparison.baseline_unique_precursor_count : '—'}</strong><span>Baseline precursor sets</span></div>
      <div><strong>{comparison.guided_unique_precursor_count}</strong><span>Guided precursor sets</span></div>
      <div><strong>{baselineRequested ? additional.size : '—'}</strong><span>Additional to baseline</span></div>
      <div><strong>{result.transfers.compiled_source_template_count}</strong><span>Admitted source templates</span></div>
    </div>
    <div className="fragment-retro-body">
      {baselineRequested ? <p className="alert warning">Guided arms use more total search work. Additional proposals indicate coverage in this bounded experiment, not improved accuracy or experimental feasibility.</p> : <p>Baseline comparison was not requested. Guided results are single-step hypotheses; no additional-coverage claim is made.</p>}
      <table className="fragment-retro-work"><caption>Actual search work across all context levels</caption><thead><tr><th>Search</th><th>Template applications</th><th>Forward validations</th></tr></thead><tbody><tr><th>Baseline</th><td>{baselineRequested ? comparison.baseline_work.template_applications : 'Not requested'}</td><td>{baselineRequested ? comparison.baseline_work.validation_attempts : '—'}</td></tr><tr><th>Guided total</th><td>{comparison.guided_work.template_applications}</td><td>{comparison.guided_work.validation_attempts}</td></tr></tbody></table>
      {showQueries && <details open><summary>Selected fragments and search coverage</summary><div className="fragment-suggestion-grid">
        {result.suggestions.candidates.map((candidate, index) => {
          const search = result.searches.find(item => item.candidate_id === candidate.candidate_id)?.result
          return <article key={candidate.candidate_id}><h3>{index + 1}. {candidate.kind.replaceAll('_', ' ')}</h3>
            {candidate.target_highlight_svg && <img className="fragment-highlight" alt={`Selected fragment ${index + 1} on target`} src={`data:image/svg+xml;charset=utf-8,${encodeURIComponent(candidate.target_highlight_svg)}`} />}
            <code>{candidate.query}</code>
            <p>Search: <strong>{search?.search_status ?? 'unavailable'}</strong> · Returned hits: {search?.returned_count ?? 0}</p>
            {search?.relationship_groups && <p>{Object.entries(search.relationship_groups).map(([name, count]) => `${readable(name)}: ${count.precision === 'at_least' ? '≥' : ''}${count.value}`).join(' · ')}. Distinct observations per group; groups overlap.</p>}
            {search?.source_scope && <p>Source: {readable(search.source_scope)} · Coverage: {search.source_coverage_complete ? 'complete' : 'incomplete'}</p>}
            {search?.stop_reason && <p>Stopped: {search.stop_reason}. Incomplete search excluded from guidance.</p>}
            <p>Selection: {candidate.reasons.join('; ')}</p>
            {candidate.cautions.length > 0 && <p>Cautions: {candidate.cautions.map(readable).join('; ')}</p>}
            <details><summary>Fragment search evidence</summary><pre>{JSON.stringify(search, null, 2)}</pre></details>
          </article>
        })}
      </div></details>}
      {!result.guidance.focus_bonds.length && <p className="alert warning">No eligible construction guidance. Complete searches with resolved internal formed-bond witnesses are required; inspect the exclusions below.</p>}
      {baselineRequested && <Arm title="Unrestricted baseline" result={result.transfers.baseline} additional={new Set()} />}
      {result.transfers.guided.map(branch => {
        const evidence = result.guidance.focus_bonds.find(bond => bond.target_atom_ids.join(',') === branch.target_atom_ids.join(','))
        return <section key={branch.target_atom_ids.join(',')} className="fragment-retro-branch">
          <h2>Construction hypothesis: target atoms {branch.target_atom_ids.join('–')}</h2>
          {evidence?.target_highlight_svg && <img className="fragment-highlight" alt={`Construction bond endpoints ${branch.target_atom_ids.join('–')}`} src={`data:image/svg+xml;charset=utf-8,${encodeURIComponent(evidence.target_highlight_svg)}`} />}
          <p>Canonical target atom IDs start at zero. These sites are projected analogue hypotheses.</p>
          <p>Fragment witness sources: {[...new Set(evidence?.supports.map(item => item.reaction_id))].join(', ')}</p>
          <Arm title="Direct source transfer" result={branch.direct_source_transfer} additional={additional} />
          <Arm title="Witness-directed library" result={branch.witness_directed_library} additional={additional} />
        </section>
      })}
      {result.source_comparisons && <details open><summary>Source products versus complete target</summary>
        <p>Structural differences help assess transfer; they do not establish functional-group tolerance or conditions. Query alignments retained: {result.guidance.target_alignment_count ?? 'not recorded'}{result.guidance.target_alignments_truncated ? ' (truncated; excluded from guidance)' : ''}.</p>
        {result.source_comparisons.map((item, index) => <article key={index}><h3>{item.reaction_id}</h3><p>Comparison: {item.comparison.status}{item.comparison.core_atom_count !== undefined && <> · Common core: {item.comparison.core_atom_count} atoms · Target coverage: {((item.comparison.right_coverage ?? 0) * 100).toFixed(0)}%</>}{item.comparison.alignment_ambiguous && ' · Ambiguous alignment'}</p><pre>{JSON.stringify(item.comparison, null, 2)}</pre></article>)}
      </details>}
      <details open><summary>Source admission and rejection evidence</summary>
        <p>A local construction witness can remain valid even when its complete source reaction fails operator admission.</p>
        {result.transfers.source_admissions.length === 0 && <p>No eligible source observations were selected.</p>}
        {result.transfers.source_admissions.map((item, index) => <p key={index}><strong>{item.reaction_id}</strong>: {item.status} · {readable(item.reason)}</p>)}
        <pre>{JSON.stringify(result.transfers.source_admissions, null, 2)}</pre>
      </details>
      <details><summary>Selected source reactions ({result.guidance.source_records.length})</summary>
        {result.guidance.source_records.map(source => <article key={source.reaction_id}><h3>{source.reaction_id}</h3><p>Reference: {source.reference_id ?? 'Unavailable'}</p><ReactionImage smiles={source.reaction_smiles} label={`Source ${source.reaction_id}`} /><code>{source.reaction_smiles}</code></article>)}
      </details>
      <details><summary>Excluded guidance and bounded selection policy</summary><pre>{JSON.stringify({ exclusions: result.guidance.exclusions, eligible_bonds: result.guidance.eligible_bond_count, eligible_sources: result.guidance.eligible_source_count, policy: result.policy }, null, 2)}</pre></details>
      <h3>Research limits</h3><ul>{result.limitations.map(item => <li key={item}>{item}</li>)}</ul>
    </div>
  </section>
}

export function FragmentGuidedRetro({ state, available, searchAvailable }: { state: State; available: boolean; searchAvailable: boolean }) {
  if (state.workflow === 'discovery') return <PrecedentDiscovery state={state.discovery} available={searchAvailable}
    onRefine={(target, query) => { state.research.changeTarget(target); state.research.changeQuery(query); state.research.setFormat('smarts'); state.research.setTopology('subgraph'); state.setWorkflow('manual') }} />
  if (state.workflow === 'manual') return <FragmentResearch state={state.research} searchAvailable={searchAvailable}
    transferAvailable={available} library={state.library} focusLimit={state.focusLimit} topK={state.topK}
    onRestoreSettings={request => { state.setLibrary(request.library_mode); state.setFocusLimit(request.max_focus_bonds); state.setTopK(request.top_k) }}
    transferView={state.research.transfer && <Evidence result={state.research.transfer} showQueries={false} />} />
  return <div className="fragment-search fragment-retro">
    <div className="editor-action-layout">
      <ReactionEditor value={state.target} onChange={state.changeTarget} onError={state.setError} moleculeOnly moleculePurpose="target" disabled={state.busy} />
      <div className="run-control workbench-action-row">
        <button className="button primary run-button" type="button" disabled={state.busy || !available || !state.target.trim()} onClick={() => void state.run()}>{state.busy ? 'Evaluating…' : 'Evaluate fragment-guided retro'}</button>
        <span role="status" aria-live="polite">{state.busy ? 'Generating fragments, searching precedents and validating proposals. Several searches may take a few minutes.' : state.result ? 'Single-step comparison ready. Export JSON to retain the complete evidence.' : 'Paste or draw one connected target molecule.'}</span>
      </div>
    </div>
    {state.result && <Evidence result={state.result} />}
  </div>
}
