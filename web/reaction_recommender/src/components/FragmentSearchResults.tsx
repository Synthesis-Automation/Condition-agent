import { useState } from 'react'
import type { FragmentCount, FragmentSearchResult } from '../api/types'
import { ReactionImage } from './ReactionImage'
import './fragment-search.css'

const label = (text: string) => text.replaceAll('_', ' ')
const count = (value: FragmentCount) => `${value.precision === 'at_least' ? '≥ ' : ''}${value.value.toLocaleString()}`

// Preserve truncation notices when displaying the domain's chunked source text.
export function sourceText(value: unknown): string {
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

export function sourceCitation(record: Record<string, unknown>, fallback?: string | null): string {
  const identity = record.reference_identity
  return (identity && typeof identity === 'object'
    ? ['raw_reference', 'patent_number', 'doi', 'normalized_citation']
      .map(key => sourceText((identity as Record<string, unknown>)[key])).find(Boolean) : '') || fallback || 'Not recorded'
}

export function FragmentSearchResults({ result, onInvestigate, busy = false }: {
  result: FragmentSearchResult; onInvestigate?: (observationId: string) => void; busy?: boolean
}) {
  const [relationship, setRelationship] = useState('all')
  const [representativesOnly, setRepresentativesOnly] = useState(false)
  const seen = new Set<string>()
  const hits = result.hits.filter(hit => {
    if (relationship !== 'all' && !hit.relationships.includes(relationship)) return false
    const key = hit.transformation_key ?? hit.observation_id
    if (representativesOnly && seen.has(key)) return false
    seen.add(key)
    return true
  })
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
      <div className="research-actions">
        <label>Evidence to inspect <select value={relationship} onChange={event => setRelationship(event.target.value)}><option value="all">All evidence</option>{Object.keys(result.relationship_groups).map(value => <option key={value} value={value}>{label(value)}</option>)}</select></label>
        <label><input type="checkbox" checked={representativesOnly} onChange={event => setRepresentativesOnly(event.target.checked)} />One representative per observed transformation</label>
      </div>
      <p>{hits.length} of {result.hits.length} returned records displayed. Grouping applies to the returned sample.</p>
      {hits.map((hit, index) => {
        const reaction = sourceText(hit.record.reaction_smiles)
        const citation = sourceCitation(hit.record, hit.reference_id)
        return <article className="fragment-hit" key={hit.hit_id}>
          <h3>{index + 1}. {hit.relationships.map(label).join(' · ')}</h3>
          <ReactionImage smiles={reaction || hit.matched_molecule_smiles || hit.product_smiles || ''} kind={reaction ? 'reaction' : 'molecule'} compact label={`Precedent ${index + 1}`} />
          {hit.matches.map((match, i) => <p key={i}><strong>{match.matched_side === 'reactant' ? 'Reported as starting material' : 'Reported product'}</strong> · {match.match_extent === 'whole_molecule' ? 'Exact compound match' : 'Substructure match'}</p>)}
          <p><strong>Reference:</strong> {citation || hit.reference_id || 'Not recorded'} · {label(hit.citation_availability)}</p>
          <p><strong>Procedure:</strong> {label(hit.procedure_availability)}{hit.procedure_match_scope ? ` · ${label(hit.procedure_match_scope)}` : ''}</p>
          <p className="fragment-note">Observation: {hit.observation_id} · Admission: {hit.admission_tier ?? 'not recorded'}</p>
          {onInvestigate && <button className="button secondary" disabled={busy} onClick={() => onInvestigate(hit.observation_id)}>Investigate this precedent</button>}
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
