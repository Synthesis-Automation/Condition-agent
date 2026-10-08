import { useEffect, useRef, useState } from 'react'
import type { ReactNode } from 'react'
import { api } from '../api/client'
import type { FragmentGuidedRetroResult, FragmentSearchRequest, FragmentSearchResult, FragmentSuggestionsResult, FragmentTransferRequest, SuggestedSearchFragment } from '../api/types'
import { FragmentSearchResults } from './FragmentSearchResults'
import { ReactionEditor } from './ReactionEditor'
import { ReactionImage } from './ReactionImage'
import type { FragmentQueryAlternatives, FragmentQueryVariant, FragmentInvestigation } from '../api/types'
import { sourceText, sourceCitation } from './FragmentSearchResults'

interface Attempt {
  id: string
  parent_id: string | null
  created_at: string
  kind: 'search' | 'transfer' | 'investigation'
  note: string
  request: FragmentSearchRequest | FragmentTransferRequest | (FragmentSearchRequest & { observation_id: string })
  result?: FragmentSearchResult | FragmentGuidedRetroResult | FragmentInvestigation
  query_revision?: FragmentQueryVariant
  error?: string
}

const HISTORY_LIMIT = 20
const EXAMPLE = 'Fc(cn1)cc2c1c(c3ccccc3OC)n[nH]2'

export function useFragmentResearch(active: boolean) {
  const [target, setTarget] = useState('')
  const [query, setQuery] = useState('')
  const [format, setFormat] = useState<FragmentSearchRequest['query_format']>('smiles')
  const [topology, setTopology] = useState<FragmentSearchRequest['topology']>('preserve_rings')
  const [timeout, setTimeoutBudget] = useState(30)
  const [note, setNote] = useState('')
  const [includeBaseline, setIncludeBaseline] = useState(false)
  const [suggestions, setSuggestions] = useState<FragmentSuggestionsResult | null>(null)
  const [suggestionScope, setSuggestionScope] = useState<'target' | 'query'>('target')
  const [current, setCurrent] = useState<Attempt | null>(null)
  const [parentId, setParentId] = useState<string | null>(null)
  const [history, setHistory] = useState<Attempt[]>([])
  const [selected, setSelected] = useState<string[]>([])
  const [transfer, setTransfer] = useState<FragmentGuidedRetroResult | null>(null)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const [alternatives, setAlternatives] = useState<FragmentQueryAlternatives | null>(null)
  const [aromaticAtoms, setAromaticAtoms] = useState<number[]>([])
  const [revision, setRevision] = useState<FragmentQueryVariant | undefined>()
  const [investigation, setInvestigation] = useState<FragmentInvestigation | null>(null)
  const pending = useRef<AbortController | null>(null)
  const counter = useRef(0)
  useEffect(() => () => pending.current?.abort(), [])
  useEffect(() => { if (!active) { pending.current?.abort(); setBusy(false); setError('') } }, [active])
  const reset = () => { setCurrent(null); setSelected([]); setTransfer(null); setError(''); setInvestigation(null); setAlternatives(null); setAromaticAtoms([]); setRevision(undefined) }
  const changeTarget = (value: string) => {
    pending.current?.abort(); setBusy(false); setTarget(value); setSuggestions(null); setParentId(null); reset()
  }
  const changeQuery = (value: string) => { setQuery(value); reset() }
  const suggest = async (scope: 'target' | 'query') => {
    const controller = new AbortController()
    pending.current?.abort(); pending.current = controller
    setBusy(true); setError(''); setSuggestions(null)
    try {
      const next = await api.suggestFragments(scope === 'target' ? target.trim() : query.trim(), controller.signal)
      if (!controller.signal.aborted) { setSuggestions(next); setSuggestionScope(scope) }
    } catch (reason) {
      if (!controller.signal.aborted) setError(reason instanceof Error ? reason.message : 'Fragment suggestion failed')
    } finally { if (!controller.signal.aborted) setBusy(false) }
  }
  const choose = (candidate: SuggestedSearchFragment) => {
    changeQuery(candidate.query); setFormat('smiles'); setTopology('preserve_rings')
    setNote(suggestionScope === 'query' ? 'Chose a smaller suggested query; removed context is unconstrained.' : 'Chose a target-derived strategic region.')
  }
  const broaden = async (selectedAtoms?: number[]) => {
    const controller = new AbortController()
    pending.current?.abort(); pending.current = controller
    setBusy(true); setError('')
    try {
      const next = await api.proposeFragmentQueries({ target_smiles: target.trim(), query: query.trim(),
        query_format: format, topology, ...(selectedAtoms?.length ? { aromatic_atom_ids: selectedAtoms } : {}) }, controller.signal)
      if (!controller.signal.aborted) setAlternatives(next)
    } catch (reason) {
      if (!controller.signal.aborted) setError(reason instanceof Error ? reason.message : 'Query preview failed')
    } finally { if (!controller.signal.aborted) setBusy(false) }
  }
  const chooseVariant = (variant: FragmentQueryVariant) => {
    changeQuery(variant.query); setFormat(variant.query_format); setTopology(variant.topology)
    setRevision(variant); setNote(variant.reason)
  }
  const execute = async (kind: Attempt['kind'], request: Attempt['request']) => {
    if (history.length >= HISTORY_LIMIT) { setError('History is full. Export and clear it before another attempt.'); return }
    const controller = new AbortController()
    pending.current?.abort(); pending.current = controller
    const attempt: Attempt = { id: `attempt-${++counter.current}`, parent_id: parentId,
      created_at: new Date().toISOString(), kind, note, request, query_revision: revision }
    setBusy(true); setError(''); setTransfer(null); setInvestigation(null)
    if (kind === 'search') { setCurrent(null); setSelected([]) }
    try {
      const result = kind === 'search'
        ? await api.searchFragments(request, controller.signal)
        : kind === 'transfer' ? await api.transferFragmentPrecedents(request as FragmentTransferRequest, controller.signal)
        : await api.investigateFragment(request as FragmentSearchRequest & { observation_id: string }, controller.signal)
      if (!controller.signal.aborted) {
        const completed = { ...attempt, result }
        setHistory(previous => [...previous, completed])
        if (kind === 'search') { setCurrent(completed); setParentId(completed.id) }
        else if (kind === 'transfer') setTransfer(result as FragmentGuidedRetroResult)
        else setInvestigation(result as FragmentInvestigation)
      }
    } catch (reason) {
      if (!controller.signal.aborted) {
        const message = reason instanceof Error ? reason.message : 'Research request failed'
        setError(message); setHistory(previous => [...previous, { ...attempt, error: message }])
      }
    } finally { if (!controller.signal.aborted) setBusy(false) }
  }
  const request = (): FragmentSearchRequest => ({ target_smiles: target.trim(), query: query.trim(),
    query_format: format, topology, limit: 10, timeout_seconds: timeout })
  const search = () => execute('search', request())
  const assess = (library: 'compact' | 'full', maxFocus: number, topK: number) => execute('transfer', {
    ...request(), target_smiles: target.trim(), selected_observation_ids: selected,
    library_mode: library, max_focus_bonds: maxFocus, top_k: topK, include_baseline: includeBaseline,
  })
  const toggleSource = (id: string) => {
    setTransfer(null)
    if (selected.includes(id)) setSelected(selected.filter(value => value !== id))
    else if (selected.length < 6) setSelected([...selected, id])
    else setError('Select at most six source observations.')
  }
  const restore = (attempt: Attempt) => {
    pending.current?.abort(); setBusy(false)
    setTarget(attempt.request.target_smiles ?? ''); setQuery(attempt.request.query)
    setFormat(attempt.request.query_format); setTopology(attempt.request.topology)
    setTimeoutBudget(attempt.request.timeout_seconds); setNote(attempt.note)
    setSuggestions(null); setParentId(attempt.id); reset()
    setRevision(attempt.query_revision)
    if (attempt.kind === 'search' && attempt.result) setCurrent(attempt)
    if (attempt.kind === 'transfer' && attempt.result) {
      const result = attempt.result as FragmentGuidedRetroResult
      const request = attempt.request as FragmentTransferRequest
      setTransfer(result); setSelected(request.selected_observation_ids); setIncludeBaseline(request.include_baseline)
      if (result.query_search) setCurrent({ ...attempt, kind: 'search', result: result.query_search })
    }
    if (attempt.kind === 'investigation' && attempt.result) setInvestigation(attempt.result as FragmentInvestigation)
  }
  const exportHistory = () => {
    const content = { schema_version: 'fragment_research_session.v1', history,
      draft: { ...request(), note, query_revision: revision }, limitations: [
        'Interactive browser history, not a pinned scientific-workspace investigation.',
        'Query revisions record user edits; logical broadening is not automatically established.',
        'Construction witnesses and graph validation do not establish experimental feasibility.',
      ] }
    const url = URL.createObjectURL(new Blob([JSON.stringify(content, null, 2)], { type: 'application/json' }))
    const link = document.createElement('a'); link.href = url; link.download = 'fragment_research_session.json'; link.click()
    setTimeout(() => URL.revokeObjectURL(url), 1000)
  }
  return { target, changeTarget, query, changeQuery, format, setFormat, topology, setTopology,
    timeout, setTimeoutBudget, note, setNote, includeBaseline, setIncludeBaseline, suggestions,
    suggestionScope, suggest, choose, current, history, selected, toggleSource, transfer,
    busy, error, setError, reset, search, assess, restore, exportHistory,
    alternatives, aromaticAtoms, setAromaticAtoms, broaden, chooseVariant, investigation,
    investigate: (observation_id: string) => execute('investigation', { ...request(), observation_id }),
    invalidateTransfer: () => setTransfer(null),
    example: () => { changeTarget(EXAMPLE); setQuery(EXAMPLE); setFormat('smiles'); setTopology('preserve_rings') },
    clearHistory: () => { setHistory([]); setParentId(null); reset() } }
}

type State = ReturnType<typeof useFragmentResearch>

function FragmentInvestigationPanel({ result }: { result: FragmentInvestigation }) {
  const outcomes: Record<string, string> = {
    graph_validated_proposals: 'Source transformation produces graph-validated proposals',
    source_compilation_rejected: 'Source is inspectable; operator compilation rejected it',
    no_resolved_construction: 'No resolved construction witness can seed this transfer',
    incomplete_search: 'Partial search: source is inspectable; transfer was not run',
    product_search_required: 'Product-side evidence is required for transfer',
    no_verified_transfer_within_budget: 'No verified transfer within the tested budget',
  }
  const source = result.source
  const record = source.record
  const show = (value: unknown) => value === null || value === undefined ? 'Not recorded' : sourceText(value) || JSON.stringify(value)
  const conditions = record.conditions ?? record.resolved_recipe
  const hasConditions = conditions && typeof conditions === 'object'
    ? Object.values(conditions).some(value => Array.isArray(value) ? value.length > 0 : value !== null && value !== '')
    : Boolean(conditions)
  return <section className="results-card fragment-results fragment-investigation" aria-label="Precedent investigation">
    <h2>Selected construction precedent</h2>
    <p><strong>{outcomes[result.transfer.status] ?? result.transfer.status}</strong></p>
    <ReactionImage smiles={sourceText(record.reaction_smiles)} label="Investigated source reaction" />
    <p>Reference: <strong>{sourceCitation(record, source.reference_id)}</strong> · Reaction: {show(record.reaction_id)}</p>
    <p>Reported yield (%): {show(record.yield_pct)} · Temperature (°C): {show(record.temperature_c)} · Time (h): {show(record.time_h)}</p>
    <details open><summary>Reported conditions</summary><pre>{hasConditions ? show(conditions) : 'Not recorded'}</pre></details>
    <p>Source evidence: {source.relationships.join(', ')}. Search: {result.search_scope.search_status}.</p>
    {!!result.construction_previews?.length && <div className="fragment-suggestion-grid">{result.construction_previews.map((preview, index) => <figure key={index}>
      <img className="fragment-highlight" alt={`Observed source construction ${index + 1}`} src={`data:image/svg+xml,${encodeURIComponent(preview.svg)}`} />
      <figcaption>Observed formed bond at source product atoms {preview.product_atom_ids.join('–')} (component {preview.component_index}). Up to three distinct witnesses shown.</figcaption>
    </figure>)}</div>}
    <p>Structural comparison: {result.comparison.status}{result.comparison.core_atom_count !== undefined ? ` · ${result.comparison.core_atom_count} common-core atoms` : ''}{result.comparison.alignment_ambiguous ? ' · Multiple alignments' : ''}.</p>
    {result.comparison.left_coverage !== undefined && result.comparison.right_coverage !== undefined && <p>The compared common core covers {Math.round(result.comparison.left_coverage * 100)}% of the source product and {Math.round(result.comparison.right_coverage * 100)}% of the target.</p>}
    {result.comparison.alignments?.[0] && <p>In the first possible alignment, {result.comparison.alignments[0].left_only_atom_ids.length} source atoms and {result.comparison.alignments[0].right_only_atom_ids.length} target atoms lie outside the common core. Inspect the differences below before transferring the chemistry.</p>}
    {result.comparison.warnings?.map(warning => <p key={warning}>{warning.replaceAll('_', ' ')}</p>)}
    <details><summary>Substrate differences and construction witnesses</summary><pre>{JSON.stringify({ comparison: result.comparison, witnesses: source.matches }, null, 2)}</pre></details>
    {source.procedures.map((procedure, index) => <details key={index}><summary>Source procedure {index + 1} · {procedure.link_scope}</summary><pre>{JSON.stringify(procedure.record, (_key, value) => value && typeof value === 'object' && 'chunks' in value ? sourceText(value) + (value.truncated ? '\n[Truncated]' : '') : value, 2)}</pre></details>)}
    {!source.procedures.length && <p>Procedure: {source.procedure_availability.replaceAll('_', ' ')}.</p>}
    {result.transfer.source_admissions?.map((admission, index) => <p key={index}>{admission.reaction_id}: {admission.status} · {admission.reason}</p>)}
    {result.transfer.arms?.map((arm, index) => <section key={index}>
      <h3>Proposed construction at target atoms {arm.target_atom_ids.join('–')}</h3>
      {arm.target_highlight_svg && <img className="fragment-highlight" alt={`Proposed target construction ${index + 1}`} src={`data:image/svg+xml,${encodeURIComponent(arm.target_highlight_svg)}`} />}
      {arm.candidates.map((candidate, i) => <article key={i}><ReactionImage smiles={candidate.proposed_reaction_smiles} label={`Selected source proposal ${index + 1}.${i + 1}`} /><p>{candidate.forward_validation_status} · {candidate.abstraction_level}</p><code>{candidate.precursor_smiles}</code></article>)}
    </section>)}
    <p>Source conditions are observed for the source substrate. Target transfers are hypotheses; graph validation does not establish experimental feasibility.</p>
    <details><summary>Full investigation and transfer diagnostics</summary><pre>{JSON.stringify(result, null, 2)}</pre></details>
  </section>
}

export function FragmentResearchOptions({ state, children }: { state: State; children: ReactNode }) {
  return <>
    <fieldset className="fragment-fields option-grid fragment-retro-grid" disabled={state.busy}>
      <label><span>Query format</span><select aria-label="Query format" value={state.format} onChange={event => { state.setFormat(event.target.value as State['format']); state.reset() }}><option value="smiles">SMILES</option><option value="smarts">SMARTS</option></select></label>
      <label><span>Query topology</span><select aria-label="Query topology" value={state.topology} onChange={event => { state.setTopology(event.target.value as State['topology']); state.reset() }}><option value="preserve_rings">Preserve ring system</option><option value="subgraph">Subgraph: allow extra rings</option></select></label>
      <label><span>Search budget (seconds)</span><input type="number" min={1} max={30} value={state.timeout} onChange={event => { state.setTimeoutBudget(Math.min(30, Math.max(1, Number(event.target.value)))); state.reset() }} /></label>
      {children}
      <button type="button" className="button quiet" onClick={state.example}>Fragment research example</button>
    </fieldset>
    <details className="advanced-options"><summary>Advanced options</summary><div>
      <label className="check-option fragment-retro-wide"><input type="checkbox" disabled={state.busy} checked={state.includeBaseline} onChange={event => { state.setIncludeBaseline(event.target.checked); state.invalidateTransfer() }} /><span>Include bounded unrestricted baseline when assessing transfer</span></label>
      <p className="fragment-note fragment-retro-wide">SMILES searches allow extra substitution unless constrained. Removing a substituent from the query stops requiring it; it does not require its absence. SMARTS changes are explicit. Up to 10 source hits are shown.</p>
    </div></details>
  </>
}

export function FragmentResearch({ state, searchAvailable, transferAvailable, library, focusLimit, topK, transferView, onRestoreSettings }: {
  state: State; searchAvailable: boolean; transferAvailable: boolean; library: 'compact' | 'full';
  focusLimit: number; topK: number; transferView: ReactNode
  onRestoreSettings: (request: FragmentTransferRequest) => void
}) {
  const result = state.current?.result as FragmentSearchResult | undefined
  const suggestions = state.suggestions?.candidates.filter(candidate => state.suggestionScope === 'target' || candidate.query !== state.suggestions?.target_smiles) ?? []
  return <div className="fragment-search fragment-retro fragment-research">
    <ReactionEditor value={state.target} onChange={state.changeTarget} onError={state.setError} moleculeOnly moleculePurpose="target" disabled={state.busy} />
    <div className="research-actions">
      <button className="button secondary" disabled={state.busy || !state.target.trim()} onClick={() => void state.suggest('target')}>Suggest strategic regions</button>
      <button className="button quiet" disabled={state.busy || !state.target.trim()} onClick={() => state.changeQuery(state.target)}>Use full target as query</button>
    </div>
    {state.suggestions && <section className="results-card fragment-results">
      <div className="results-summary"><h2>{state.suggestionScope === 'target' ? 'Choose a search region' : 'Choose a simpler query'}</h2></div>
      <p>Suggestions preserve complete ring systems and required valence/stereo context. No search runs until you choose and submit a query.</p>
      {!suggestions.length && <p>No smaller valid query was suggested. You can edit the query drawing or SMARTS explicitly.</p>}
      <div className="fragment-suggestion-grid">{suggestions.map((candidate, index) => <article key={candidate.candidate_id}>
        <h3>{index + 1}. {candidate.kind.replaceAll('_', ' ')}</h3>
        {candidate.target_highlight_svg && <img className="fragment-highlight" alt={`${state.suggestionScope} region ${index + 1}`} src={`data:image/svg+xml,${encodeURIComponent(candidate.target_highlight_svg)}`} />}
        <code>{candidate.query}</code><p>{candidate.reasons.join(' ')}</p>
        {candidate.cautions.length > 0 && <p className="alert caution">{candidate.cautions.join(' ')}</p>}
        <button className="button secondary" disabled={state.busy} onClick={() => state.choose(candidate)}>Use region {index + 1}</button>
      </article>)}</div>
    </section>}
    <div className="editor-action-layout">
      <ReactionEditor value={state.query} onChange={state.changeQuery} onError={state.setError} moleculeOnly moleculePurpose="fragment" queryFormat={state.format} disabled={state.busy} />
      <div className="run-control workbench-action-row">
        <button className="button primary run-button" disabled={state.busy || !searchAvailable || !state.target.trim() || !state.query.trim()} onClick={() => void state.search()}>Search chosen fragment</button>
        <button className="button quiet" disabled={state.busy || state.format !== 'smiles' || !state.query.trim()} onClick={() => void state.suggest('query')}>Suggest simpler queries</button>
        <button className="button secondary" disabled={state.busy || !state.target.trim() || !state.query.trim()} onClick={() => void state.broaden()}>Preview query alternatives</button>
        <span role="status">{state.busy ? 'Searching or assessing selected sources…' : 'The query is validated against the full target before index search.'}</span>
      </div>
    </div>
    {state.error && <div className="alert error" role="alert">{state.error}</div>}
    {state.alternatives && <section className="results-card fragment-results fragment-query-alternatives" aria-label="Query alternatives">
      <h2>Choose an explicit relaxation</h2>
      <p>Each choice matches the target. Highlighting shows one possible alignment; choosing a query does not run a search.</p>
      <div className="fragment-suggestion-grid">{state.alternatives.variants.map((variant, index) => <article key={variant.variant_id}>
        <h3>{variant.relaxations.map(value => value.replaceAll('_', ' ')).join(' + ')}</h3>
        {variant.target_highlight_svg && <img className="fragment-highlight" alt={`Query alternative ${index + 1}`} src={`data:image/svg+xml,${encodeURIComponent(variant.target_highlight_svg)}`} />}
        <code>{variant.query}</code><p>{variant.reason}</p>
        {variant.alignment_ambiguous && <p>Multiple target alignments are possible{variant.target_alignments_truncated ? '; enumeration was truncated' : ''}.</p>}
        <button className="button secondary" disabled={state.busy} onClick={() => state.chooseVariant(variant)}>Use alternative {index + 1}</button>
      </article>)}</div>
      {!state.alternatives.variants.length && <p>No automatic alternatives are available. Edit custom SMARTS to investigate a different hypothesis.</p>}
      {state.alternatives.query_atoms.some(atom => atom.allows_carbon_nitrogen) && <fieldset disabled={state.busy}>
        <legend>Optional C/N alternatives at selected query atoms</legend>
        {state.alternatives.query_atom_indices_svg && <img className="fragment-highlight query-atom-map" alt="Query with atom indices" src={`data:image/svg+xml,${encodeURIComponent(state.alternatives.query_atom_indices_svg)}`} />}
        <p>Zero-based atom IDs follow the input query SMILES order. Pyrrole [nH], charge, isotope and stereochemical constraints are preserved.</p>
        <div className="research-actions">{state.alternatives.query_atoms.filter(atom => atom.allows_carbon_nitrogen).map(atom => <label className="fragment-atom-choice" key={atom.query_atom_id}>
          <input type="checkbox" checked={state.aromaticAtoms.includes(atom.query_atom_id)} onChange={() => state.setAromaticAtoms(state.aromaticAtoms.includes(atom.query_atom_id) ? state.aromaticAtoms.filter(i => i !== atom.query_atom_id) : [...state.aromaticAtoms, atom.query_atom_id])} />Atom {atom.query_atom_id} ({atom.element})
        </label>)}</div>
        <button className="button quiet" disabled={!state.aromaticAtoms.length} onClick={() => void state.broaden(state.aromaticAtoms)}>Preview selected C/N alternatives</button>
      </fieldset>}
    </section>}
    <label className="research-note"><span>Why this region or query revision?</span><textarea disabled={state.busy} maxLength={5000} value={state.note} onChange={event => state.setNote(event.target.value)} placeholder="For example: retain the fused core, omit the peripheral methoxy group." /></label>
    {result && <>
      <section className="results-card research-selection">
        <h2>Select source observations for transfer</h2>
        <p>Select up to six observations. Construction, retention and unresolved evidence remain distinct. Only resolved internal construction witnesses can seed proposals.</p>
        {result.search_status !== 'complete' && <p className="alert caution">This search is {result.search_status}; inspect its results, then refine or rerun it. It cannot seed transfer.</p>}
        {result.hits.map(hit => <label className="check-option" key={hit.hit_id}>
          <input type="checkbox" checked={state.selected.includes(hit.observation_id)} disabled={state.busy || result.search_status !== 'complete'} onChange={() => state.toggleSource(hit.observation_id)} />
          <span>{String(hit.record.reaction_id ?? hit.observation_id)} · {hit.relationships.map(value => value.replaceAll('_', ' ')).join(', ')}</span>
        </label>)}
        <button className="button primary" disabled={state.busy || !transferAvailable || !state.selected.length || result.search_status !== 'complete'} onClick={() => void state.assess(library, focusLimit, topK)}>Assess selected precedents on target</button>
        {!transferAvailable && <p>Transfer requires the selected operator library; source discovery remains available independently.</p>}
      </section>
      <FragmentSearchResults result={result} onInvestigate={id => void state.investigate(id)} busy={state.busy} />
    </>}
    {state.investigation && <FragmentInvestigationPanel result={state.investigation} />}
    {transferView}
    <section className="results-card research-history">
      <div className="results-summary"><h2>Research history ({state.history.length}/{HISTORY_LIMIT})</h2><div className="research-actions">
        <button className="button quiet" disabled={!state.history.length} onClick={state.exportHistory}>Export research JSON</button>
        <button className="button quiet" disabled={state.busy || !state.history.length} onClick={state.clearHistory}>Clear history</button>
      </div></div>
      <p>Completed searches and errors remain here while you revise queries. History stays in this browser session; export before closing the page. Query edits are not automatically certified as broader searches.</p>
      {state.history.map(attempt => <details key={attempt.id} open={Boolean(attempt.error)}><summary>{attempt.id} · {attempt.kind} · {attempt.error ? 'error' : attempt.kind === 'search' ? (attempt.result as FragmentSearchResult).search_status : 'assessed'} · {attempt.request.query}</summary>
        <p>Target: <code>{attempt.request.target_smiles}</code></p><p>Note: {attempt.note || 'Not recorded'} · Parent: {attempt.parent_id ?? 'none'}</p>
        {attempt.error && <p>{attempt.error}</p>}
        <button className="button quiet" disabled={state.busy} onClick={() => {
          state.restore(attempt)
          if (attempt.kind === 'transfer') onRestoreSettings(attempt.request as FragmentTransferRequest)
        }}>Restore {attempt.id}</button>
        <pre>{JSON.stringify(attempt, null, 2)}</pre>
      </details>)}
    </section>
  </div>
}
