import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import type { SynthesisPrecedentResult } from '../api/types'
import { ReactionEditor } from './ReactionEditor'
import { ReactionImage } from './ReactionImage'
import { sourceCitation, sourceText } from './FragmentSearchResults'

export function usePrecedentDiscovery(active: boolean) {
  const [target, setTarget] = useState('')
  const [result, setResult] = useState<SynthesisPrecedentResult | null>(null)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const pending = useRef<AbortController | null>(null)
  useEffect(() => () => pending.current?.abort(), [])
  useEffect(() => { if (!active) { pending.current?.abort(); setBusy(false) } }, [active])
  const changeTarget = (value: string) => {
    pending.current?.abort(); setTarget(value); setBusy(false); setResult(null); setError('')
  }
  const run = async () => {
    const controller = new AbortController()
    pending.current?.abort(); pending.current = controller
    setBusy(true); setError(''); setResult(null)
    try {
      const next = await api.findSynthesisPrecedents(target.trim(), controller.signal)
      if (!controller.signal.aborted) setResult(next)
    } catch (reason) {
      if (!controller.signal.aborted) setError(reason instanceof Error ? reason.message : 'Precedent discovery failed')
    } finally { if (!controller.signal.aborted) setBusy(false) }
  }
  const exportResult = () => {
    if (!result) return
    const url = URL.createObjectURL(new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' }))
    const link = document.createElement('a'); link.href = url; link.download = 'synthesis_precedents.json'; link.click()
    setTimeout(() => URL.revokeObjectURL(url), 1000)
  }
  return { target, changeTarget, result, busy, error, setError, run, exportResult }
}

export function PrecedentDiscovery({ state, available, onRefine }: {
  state: ReturnType<typeof usePrecedentDiscovery>; available: boolean
  onRefine: (target: string, query: string) => void
}) {
  const result = state.result
  return <div className="fragment-search fragment-research">
    <ReactionEditor value={state.target} onChange={state.changeTarget} onError={state.setError} moleculeOnly disabled={state.busy} />
    <div className="research-actions">
      <button className="button primary" disabled={!available || state.busy || !state.target.trim()} onClick={() => void state.run()}>{state.busy ? 'Finding precedents…' : 'Find synthesis precedents'}</button>
      {result && <button className="button quiet" onClick={state.exportResult}>Export search evidence</button>}
    </div>
    <p role="status">{state.busy ? 'Searching distinctive cores and adding context where useful. This may take up to two minutes.' : 'Start with a target. Core selection and search refinement are automatic.'}</p>
    {!available && <p className="alert caution">A prepared fragment index is required.</p>}
    {state.error && <p role="alert" className="alert error">{state.error}</p>}
    {result && <section className="results-card fragment-results" aria-label="Synthesis precedents">
      <h2>Synthesis precedents</h2>
      <p>{result.returned_count} source observations · {result.execution.elapsed_seconds.toFixed(1)} s · {result.search_status === 'partial' ? 'Bounded or partial search' : 'Search policy completed'}</p>
      <p>Exact target lookup: {result.exact_target.status.replaceAll('_', ' ')}. Results are analogues to investigate, not validated routes.</p>
      {!result.hits.length && <p>No precedents returned within this search scope and budget. This does not establish absence from the literature.</p>}
      {result.hits.map((hit, index) => <article className="fragment-hit" key={hit.observation_id}>
        <h3>{index + 1}. {hit.discovery.core_relationship.replaceAll('_', ' ')}</h3>
        <p>{hit.discovery.explanation}</p>
        <ReactionImage smiles={sourceText(hit.record.reaction_smiles) || hit.matched_molecule_smiles || ''} kind={sourceText(hit.record.reaction_smiles) ? 'reaction' : 'molecule'} compact label={`Discovery precedent ${index + 1}`} />
        <p><strong>Reference:</strong> {sourceCitation(hit.record, hit.reference_id)}</p>
        {hit.discovery.nitrogen_hydrogen_differences.map(d => <p key={d.target_atom_id}>Nitrogen substitution differs: target has {d.target_hydrogens} N-bound H; source match has {d.source_hydrogen_counts.join(' or ')}.</p>)}
        {hit.discovery.alignment_ambiguous && <p>More than one core alignment is possible.</p>}
        <details><summary>Inspect source and conditions</summary><pre>{JSON.stringify(hit.record, null, 2)}</pre>
          {hit.procedures.map((p, i) => <div key={i}><h4>Procedure {i + 1} · {p.link_scope}</h4><pre>{JSON.stringify(p.record, null, 2)}</pre></div>)}
        </details>
        <details><summary>Refine search or investigate transfer</summary><p>Open the retained query in the detailed research workflow. Transfer still requires independent graph checks.</p><button className="button secondary" onClick={() => onRefine(result.target_smiles, hit.discovery.query)}>Use this core in detailed research</button></details>
      </article>)}
      <details><summary>Search history and relaxations</summary><pre>{JSON.stringify(result.attempts, null, 2)}</pre></details>
      <details><summary>Search limitations</summary>{result.limitations.map(text => <p key={text}>{text}</p>)}</details>
    </section>}
  </div>
}
