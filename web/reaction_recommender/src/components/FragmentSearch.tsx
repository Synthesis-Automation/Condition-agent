import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import type { FragmentCount, FragmentSearchRequest, FragmentSearchResult, FragmentSuggestionsResult, SuggestedSearchFragment } from '../api/types'
import { ReactionImage } from './ReactionImage'
import { ReactionEditor } from './ReactionEditor'
import './fragment-search.css'

const CORE = 'c1ccc2c(c1)COc1ccccc1-2'
const label = (text: string) => text.replaceAll('_', ' ')
const count = (value: FragmentCount) => `${value.precision === 'at_least' ? '≥ ' : ''}${value.value.toLocaleString()}`

// Preserve truncation notices when displaying the domain's chunked source text.
function sourceText(value: unknown): string {
  if (typeof value === 'string') return value
  if (value && typeof value === 'object' && 'chunks' in value && Array.isArray(value.chunks)) {
    return value.chunks.map(chunk => String(chunk.text ?? '')).join('')
  }
  return ''
}

function readableEvidence(_key: string, value: unknown): unknown {
  if (value && typeof value === 'object' && 'chunks' in value) {
    return sourceText(value) + ('truncated' in value && value.truncated ? '\n[Source text truncated]' : '')
  }
  return value
}

export function useFragmentSearch(active: boolean) {
  const [query, setQuery] = useState('')
  const [format, setFormat] = useState<FragmentSearchRequest['query_format']>('smiles')
  const [topology, setTopology] = useState<FragmentSearchRequest['topology']>('preserve_rings')
  const [searchSide, setSearchSide] = useState<NonNullable<FragmentSearchRequest['search_side']>>('product')
  const [limit, setLimit] = useState(5)
  const [timeout, setBudget] = useState(30)
  const [busy, setBusy] = useState(false)
  const [suggesting, setSuggesting] = useState(false)
  const [suggestions, setSuggestions] = useState<FragmentSuggestionsResult | null>(null)
  const [error, setError] = useState('')
  const [result, setResult] = useState<FragmentSearchResult | null>(null)
  const pending = useRef<AbortController | null>(null)
  useEffect(() => () => pending.current?.abort(), [])
  useEffect(() => {
    if (!active) {
      pending.current?.abort()
      setBusy(false)
      setSuggesting(false)
      setResult(null)
      setError('')
    }
  }, [active])

  const reset = () => { setResult(null); setError('') }
  const changeQuery = (value: string) => { setQuery(value); setSuggestions(null) }
  const suggest = async () => {
    const controller = new AbortController()
    pending.current = controller
    setBusy(true)
    setSuggesting(true)
    setSuggestions(null)
    reset()
    try {
      const next = await api.suggestFragments(query.trim(), controller.signal)
      if (!controller.signal.aborted) setSuggestions(next)
    } catch (reason) {
      if (!controller.signal.aborted) setError(reason instanceof Error ? reason.message : 'Suggestion failed')
    } finally {
      if (!controller.signal.aborted) { setBusy(false); setSuggesting(false) }
    }
  }
  const choose = (candidate: SuggestedSearchFragment) => {
    setQuery(candidate.query)
    setFormat(candidate.query_format)
    setTopology(candidate.topology)
    reset()
  }
  const search = async () => {
    const controller = new AbortController()
    pending.current = controller
    setBusy(true)
    reset()
    try {
      const next = await api.searchFragments({ query: query.trim(), query_format: format,
        topology, search_side: searchSide, limit, timeout_seconds: timeout }, controller.signal)
      if (!controller.signal.aborted) setResult(next)
    } catch (reason) {
      if (!controller.signal.aborted) setError(reason instanceof Error ? reason.message : 'Search failed')
    } finally {
      if (!controller.signal.aborted) setBusy(false)
    }
  }

  const exportResult = () => {
    const url = URL.createObjectURL(new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' }))
    const link = document.createElement('a')
    link.href = url
    link.download = 'fragment_precedents.json'
    link.click()
    setTimeout(() => URL.revokeObjectURL(url), 1000)
  }

  return { query, setQuery: changeQuery, format, setFormat, topology, setTopology, searchSide, setSearchSide, limit, setLimit,
    timeout, setBudget, busy, error, setError, result, reset, search, exportResult,
    suggestions, suggesting, suggest, choose }
}

type FragmentSearchState = ReturnType<typeof useFragmentSearch>

export function FragmentSearchOptions({ state, available }: { state: FragmentSearchState; available?: boolean }) {
  const { format, setFormat, topology, setTopology, limit, setLimit, timeout, setBudget,
    busy, error, reset, setQuery, query, suggest, suggesting, searchSide, setSearchSide } = state
  return <div className="analysis-options">
    <fieldset disabled={busy} className="fragment-fields">
      <div className="option-grid">
        <label><span>Query format</span><select aria-label="Query format" value={format} onChange={event => { setFormat(event.target.value as typeof format); reset() }}>
          <option value="smiles">SMILES</option><option value="smarts">SMARTS</option>
        </select></label>
        <label><span>Reaction side</span><select aria-label="Reaction side" value={searchSide} onChange={event => { setSearchSide(event.target.value as typeof searchSide); reset() }}>
          <option value="product">Products</option><option value="reactant">Starting materials</option><option value="either">Either side</option>
        </select></label>
        <label><span>Topology</span><select aria-label="Topology" value={topology} onChange={event => { setTopology(event.target.value as typeof topology); reset() }}>
          <option value="preserve_rings">Preserve ring system</option><option value="subgraph">Subgraph (allow extra rings)</option>
        </select></label>
        <label><span>Top results</span><input type="number" min={1} max={10} value={limit} onChange={event => { setLimit(Math.min(10, Math.max(1, Number(event.target.value)))); reset() }} /></label>
      </div>
      <div className="feature-mode-note"><strong>Fragment precedent search</strong><span>Find reactions containing your core and inspect construction, modification or retention evidence.</span></div>
      <div className="inline-checks"><button className="button quiet" type="button" onClick={() => { setQuery(CORE); setFormat('smiles'); setTopology('preserve_rings'); reset() }}>Cyclic ether example</button>
        <button className="button secondary" type="button" disabled={format !== 'smiles' || !query.trim()} onClick={() => void suggest()}>{suggesting ? 'Suggesting…' : 'Suggest fragments'}</button>
      </div>
      <small>Enter or draw a whole target, then optionally suggest search fragments. Choose a candidate before searching; no index is needed for suggestions.</small>
      <details className="advanced-options"><summary>Advanced options</summary><div>
        <label><span>Search budget (seconds)</span><input type="number" min={1} max={30} value={timeout} onChange={event => { setBudget(Math.min(30, Math.max(1, Number(event.target.value)))); reset() }} /></label>
        <div className="feature-mode-note"><span>Preserve ring system excludes additional fused, bridged or spiro rings. The search does not broaden automatically.</span></div>
      </div></details>
    </fieldset>
    {available === false && <div className="alert caution" role="alert">Fragment index unavailable. Configure a prepared index with --fragment-index and restart the server.</div>}
    {error && <div className="alert error" role="alert">{error}</div>}
  </div>
}

export function FragmentSearch({ state, available }: { state: FragmentSearchState; available?: boolean }) {
  const { query, setQuery, format, busy, result, reset, search, setError, suggestions, suggesting, choose } = state
  return <section className="fragment-search" aria-label="Fragment precedent search">
    <form className="editor-action-layout" onSubmit={event => { event.preventDefault(); void search() }}>
      <ReactionEditor value={query} onChange={value => { setQuery(value); reset() }} onError={setError}
        moleculeOnly moleculePurpose="fragment" queryFormat={format} disabled={busy} />
      <div className="run-control workbench-action-row" aria-label="Analysis action">
        <button className="button primary run-button" type="submit" disabled={busy || available === false || !query.trim()}>{busy && !suggesting ? 'Searching…' : 'Search fragments'}</button>
        <span role="status" aria-live="polite">{suggesting ? 'Extracting target-derived search regions…' : busy ? 'Searching the local index and checking reaction evidence…' : result ? 'Search finished' : 'Ready'}</span>
      </div>
    </form>
    {suggestions && <section className="results-card fragment-results" aria-label="Suggested fragments">
      <div className="results-summary"><div><span className="eyebrow">OPTIONAL SEARCH REGIONS</span><h2>Suggested search fragments</h2></div></div>
      <p>Highlighted atoms belong to the original target. These overlapping queries are suggestions, not precursors or a ranked synthesis plan. Choose one, then search.</p>
      <details><summary>Original target and generation scope</summary><code>{suggestions.target_smiles}</code><p>{suggestions.generated_count} distinct candidates · {suggestions.rejected_count} rejected extractions</p>
        {(suggestions.generation_truncated || suggestions.output_truncated) && <p>Only a bounded subset is shown; this is not an exhaustive fragmentation.</p>}
      </details>
      {!suggestions.candidates.length && <p>No valid candidate within the current limits. You can enter your own connected core.</p>}
      <div className="fragment-suggestion-grid">{suggestions.candidates.map((candidate, index) => <article key={candidate.candidate_id}>
        <h3>{index + 1}. {label(candidate.kind)}</h3>
        {candidate.target_highlight_svg && <img className="fragment-highlight" src={`data:image/svg+xml,${encodeURIComponent(candidate.target_highlight_svg)}`} alt={`Target atoms for candidate ${index + 1}`} />}
        <ReactionImage smiles={candidate.query} kind="molecule" label={`Suggested fragment ${index + 1}`} />
        <p>{candidate.reasons.join(' ')}</p>
        {candidate.cautions.includes('LOW_STRUCTURAL_SPECIFICITY_MAY_BE_BROAD') && <p>Simple core; precedent search may be broad.</p>}
        <details><summary>Query and boundaries</summary><code>{candidate.query}</code><p>Target atom IDs: {candidate.target_atom_ids.join(', ')}</p><pre>{JSON.stringify(candidate.boundaries, null, 2)}</pre><p>{candidate.cautions.map(label).join('; ')}</p></details>
        <button className="button secondary" type="button" disabled={busy} onClick={() => choose(candidate)}>Use candidate {index + 1}</button>
      </article>)}</div>
    </section>}
    {result && <FragmentSearchResults result={result} />}
    {!result && <section className="empty-state"><span>3</span><div><h2>Inspect fragment precedents</h2><p>Matching reactions, source references and graph evidence will appear here.</p></div></section>}
  </section>
}

export function FragmentSearchResults({ result }: { result: FragmentSearchResult }) {
  return <section className="results-card fragment-results" aria-label="Fragment search results">
      <div className="results-summary"><div><span className="eyebrow">FRAGMENT SEARCH RESULT</span><h2>{result.search_status === 'too_broad' ? 'Query too broad' : result.search_status === 'partial' ? 'Partial search results' : 'Fragment precedents'}</h2></div>
        <div className="metric-strip">
          <div><strong>{count(result.counts.molecules ?? result.counts.products)}</strong><span>matching compounds</span></div>
          <div><strong>{count(result.counts.observations)}</strong><span>observations</span></div>
          <div><strong>{count(result.counts.known_references)}</strong><span>references</span></div>
        </div>
      </div>
      <p>{result.returned_count} hits shown · {result.execution.elapsed_seconds.toFixed(1)} s</p>
      <p className="fragment-note">{result.source_scope === 'prefix_pilot' ? 'Prefix pilot corpus' : 'Indexed corpus'} · {result.source_coverage_complete ? 'Source coverage complete' : 'Source coverage incomplete'}</p>
      {result.stop_reason && <p>Search stopped: {label(result.stop_reason)}. Counts marked ≥ are lower bounds; evidence may not yet have been examined.</p>}
      {result.refinement_hints.map(hint => <p key={hint}>{hint}</p>)}
      {result.search_status === 'complete' && (result.counts.molecules ?? result.counts.products).value === 0 && <p>No matching compounds on the selected reaction side.</p>}
      {result.output_truncated && <p>Some hits were omitted to keep the response within its size limit.</p>}
      <p className="fragment-note">{Object.entries(result.relationship_groups).map(([name, value]) => `${label(name)}: ${count(value)}`).join(' · ')}. Groups can overlap.</p>
      {result.hits.map((hit, index) => {
        const reaction = sourceText(hit.record.reaction_smiles)
        const identity = hit.record.reference_identity
        const citation = identity && typeof identity === 'object'
          ? ['raw_reference', 'patent_number', 'doi', 'normalized_citation']
            .map(key => sourceText((identity as Record<string, unknown>)[key])).find(Boolean)
          : ''
        return <article className="fragment-hit" key={hit.hit_id}>
          <h3>{index + 1}. {hit.relationships.map(label).join(' · ')}</h3>
          <ReactionImage smiles={reaction || hit.matched_molecule_smiles || hit.product_smiles || ''} kind={reaction ? 'reaction' : 'molecule'} compact label={`Precedent ${index + 1}`} />
          {hit.matches.map((match, i) => <p key={i}><strong>{match.matched_side === 'reactant' ? 'Reported as starting material' : 'Reported product'}</strong> · {match.match_extent === 'whole_molecule' ? 'Exact compound match' : 'Substructure match'}</p>)}
          <p><strong>Reference:</strong> {citation || hit.reference_id || 'Not recorded'} · {label(hit.citation_availability)}</p>
          <p><strong>Procedure:</strong> {label(hit.procedure_availability)}{hit.procedure_match_scope ? ` · ${label(hit.procedure_match_scope)}` : ''}</p>
          <p className="fragment-note">Observation: {hit.observation_id} · Admission: {hit.admission_tier ?? 'not recorded'}</p>
          {hit.warnings.length > 0 && <p>{hit.warnings.map(label).join('; ')}</p>}
          <details><summary>Source record and conditions</summary><pre>{JSON.stringify(hit.record, readableEvidence, 2)}</pre></details>
          <details><summary>Graph evidence and atom correspondence</summary><pre>{JSON.stringify(hit.matches, null, 2)}</pre></details>
          {hit.procedures.map((procedure, i) => <details key={i}><summary>Procedure {i + 1} · {label(procedure.link_scope)}</summary><pre>{JSON.stringify(procedure.record, readableEvidence, 2)}</pre></details>)}
          {hit.procedure_records_truncated && <p>Additional procedure records were omitted.</p>}
        </article>
      })}
      <details><summary>Search limitations</summary>{result.limitations.map(text => <p key={text}>{text}</p>)}</details>
  </section>
}
