import { useEffect, useRef, useState } from 'react'
import type { ReactNode } from 'react'
import { api } from '../api/client'
import type { FragmentGuidedRetroResult, FragmentSearchRequest, FragmentSearchResult, FragmentSuggestionsResult, FragmentTransferRequest, SuggestedSearchFragment } from '../api/types'
import { FragmentSearchResults } from './FragmentSearchResults'
import { ReactionEditor } from './ReactionEditor'

interface Attempt {
  id: string
  parent_id: string | null
  created_at: string
  kind: 'search' | 'transfer'
  note: string
  request: FragmentSearchRequest | FragmentTransferRequest
  result?: FragmentSearchResult | FragmentGuidedRetroResult
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
  const pending = useRef<AbortController | null>(null)
  const counter = useRef(0)
  useEffect(() => () => pending.current?.abort(), [])
  useEffect(() => { if (!active) { pending.current?.abort(); setBusy(false); setError('') } }, [active])
  const reset = () => { setCurrent(null); setSelected([]); setTransfer(null); setError('') }
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
  const execute = async (kind: Attempt['kind'], request: Attempt['request']) => {
    if (history.length >= HISTORY_LIMIT) { setError('History is full. Export and clear it before another attempt.'); return }
    const controller = new AbortController()
    pending.current?.abort(); pending.current = controller
    const attempt: Attempt = { id: `attempt-${++counter.current}`, parent_id: parentId,
      created_at: new Date().toISOString(), kind, note, request }
    setBusy(true); setError(''); setTransfer(null)
    if (kind === 'search') { setCurrent(null); setSelected([]) }
    try {
      const result = kind === 'search'
        ? await api.searchFragments(request, controller.signal)
        : await api.transferFragmentPrecedents(request as FragmentTransferRequest, controller.signal)
      if (!controller.signal.aborted) {
        const completed = { ...attempt, result }
        setHistory(previous => [...previous, completed])
        if (kind === 'search') { setCurrent(completed); setParentId(completed.id) }
        else setTransfer(result as FragmentGuidedRetroResult)
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
    if (attempt.kind === 'search' && attempt.result) setCurrent(attempt)
    if (attempt.kind === 'transfer' && attempt.result) {
      const result = attempt.result as FragmentGuidedRetroResult
      const request = attempt.request as FragmentTransferRequest
      setTransfer(result); setSelected(request.selected_observation_ids); setIncludeBaseline(request.include_baseline)
      if (result.query_search) setCurrent({ ...attempt, kind: 'search', result: result.query_search })
    }
  }
  const exportHistory = () => {
    const content = { schema_version: 'fragment_research_session.v1', history,
      draft: { ...request(), note }, limitations: [
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
    invalidateTransfer: () => setTransfer(null),
    example: () => { changeTarget(EXAMPLE); setQuery(EXAMPLE); setFormat('smiles'); setTopology('preserve_rings') },
    clearHistory: () => { setHistory([]); setParentId(null); reset() } }
}

type State = ReturnType<typeof useFragmentResearch>

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
    {state.error && <div className="alert error" role="alert">{state.error}</div>}
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
        <span role="status">{state.busy ? 'Searching or assessing selected sources…' : 'The query is validated against the full target before index search.'}</span>
      </div>
    </div>
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
      <FragmentSearchResults result={result} />
    </>}
    {transferView}
    <section className="results-card research-history">
      <div className="results-summary"><h2>Research history ({state.history.length}/{HISTORY_LIMIT})</h2><div className="research-actions">
        <button className="button quiet" disabled={!state.history.length} onClick={state.exportHistory}>Export research JSON</button>
        <button className="button quiet" disabled={state.busy || !state.history.length} onClick={state.clearHistory}>Clear history</button>
      </div></div>
      <p>Completed searches and errors remain here while you revise queries. History stays in this browser session; export before closing the page. Query edits are not automatically certified as broader searches.</p>
      {state.history.map(attempt => <details key={attempt.id}><summary>{attempt.id} · {attempt.kind} · {attempt.error ? 'error' : attempt.kind === 'search' ? (attempt.result as FragmentSearchResult).search_status : 'assessed'} · {attempt.request.query}</summary>
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
