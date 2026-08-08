import { useEffect, useState } from 'react'
import type {
  DiscoveryHit,
  DiscoveryResult,
  RecipeComponent,
  Recommendation,
  RecommendationResult,
  ResolvedRecipe,
  SynthesisProtocolDraft,
} from '../api/types'
import { ReactionImage } from './ReactionImage'

const ROLE_FIELDS: Array<[keyof ResolvedRecipe, string]> = [
  ['catalysts', 'Catalyst'],
  ['ligands', 'Ligand'],
  ['bases', 'Base'],
  ['condensation_agents', 'Condensation agent'],
  ['oxidants', 'Oxidant'],
  ['reductants', 'Reductant'],
  ['acids', 'Acid'],
  ['additives', 'Additive'],
  ['solvents', 'Solvent'],
  ['other_components', 'Other'],
]

const FACTOR_LABELS: Record<string, string> = {
  similarity: 'Structural similarity',
  partner_category: 'Reactant category',
  functional_group_tolerance: 'Functional-group tolerance',
  yield: 'Expected yield',
  independent_support: 'Independent support',
  reaction_breadth: 'Reaction breadth',
  dataset_diversity: 'Dataset diversity',
  compatibility: 'Condition compatibility',
  condition_certainty: 'Procedure completeness',
  edit_similarity: 'Bond-edit similarity',
  reaction_center: 'Reaction center',
  local_environment: 'Local environment',
}

export function displayName(value: unknown): string {
  return String(value ?? '')
    .replaceAll('_', ' ')
    .replace(/\b\w/g, (letter) => letter.toUpperCase())
}

function componentName(component: RecipeComponent): string {
  return String(
    component.canonical_name ??
      component.display_name ??
      component.name ??
      component.raw_identifier ??
      component.substance_id ??
      'Unresolved',
  )
}

export function recipeSummary(recipe: ResolvedRecipe): string {
  const parts: string[] = []
  for (const [field, label] of ROLE_FIELDS) {
    const values = recipe[field]
    if (Array.isArray(values) && values.length) {
      parts.push(`${label}: ${values.map(componentName).join(', ')}`)
    }
  }
  if (recipe.temperature_c != null) parts.push(`${recipe.temperature_c} °C`)
  if (recipe.time_h != null) parts.push(`${recipe.time_h} h`)
  return parts.join(' · ') || 'Resolved recipe details unavailable'
}

function Conditions({ recipe }: { recipe: ResolvedRecipe }) {
  const rows: Array<[string, string]> = []
  for (const [field, label] of ROLE_FIELDS) {
    const values = recipe[field]
    if (Array.isArray(values) && values.length) {
      rows.push([label, values.map(componentName).join(', ')])
    }
  }
  for (const [field, label, unit] of [
    ['temperature_c', 'Temperature', '°C'],
    ['time_h', 'Time', 'h'],
    ['concentration_m', 'Concentration', 'M'],
    ['pressure_bar', 'Pressure', 'bar'],
  ] as const) {
    const value = recipe[field]
    if (value != null) rows.push([label, `${value} ${unit}`])
  }
  if (recipe.atmosphere) rows.push(['Atmosphere', recipe.atmosphere])
  return (
    <dl className="detail-list">
      {rows.map(([label, value]) => (
        <div key={label}><dt>{label}</dt><dd>{value}</dd></div>
      ))}
    </dl>
  )
}

function MessageList({ title, values, tone = '' }: { title: string; values: string[]; tone?: string }) {
  if (!values.length) return null
  return (
    <div className={`message-group ${tone}`}>
      <h4>{title}</h4>
      <ul>{values.map((value, index) => <li key={`${value}-${index}`}>{value}</li>)}</ul>
    </div>
  )
}

function protocolFilename(rank: number, recipeId: string): string {
  const recipeToken = recipeId.replace(/[^a-z0-9]+/gi, '_').replace(/^_|_$/g, '').slice(0, 48)
  return `condition_protocol_rank_${rank}_${recipeToken || 'recipe'}.json`
}

function downloadProtocol(protocol: SynthesisProtocolDraft, rank: number) {
  const blob = new Blob([JSON.stringify(protocol, null, 2)], { type: 'application/json' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = protocolFilename(rank, protocol.recipe_id)
  document.body.appendChild(link)
  link.click()
  link.remove()
  URL.revokeObjectURL(url)
}

function ProtocolPanel({ protocol, rank }: { protocol: SynthesisProtocolDraft; rank: number }) {
  return (
    <details className="trace-panel protocol-panel" open>
      <summary>
        <span>Condition protocol JSON</span>
        <span className="readiness-badge">{displayName(protocol.execution_readiness)}</span>
      </summary>
      <div className="protocol-toolbar">
        <span>{protocol.materials.length} materials · {protocol.missing_required_fields.length} missing fields</span>
        <button className="button secondary" type="button" onClick={() => downloadProtocol(protocol, rank)}>Download condition JSON</button>
      </div>
      <pre>{JSON.stringify(protocol, null, 2)}</pre>
    </details>
  )
}

function RecommendationDetails({ item }: { item: Recommendation }) {
  const precedent = item.precedent_reaction_smiles[0] ?? ''
  const protocol = item.synthesis_protocol
  return (
    <article className="result-detail">
      <div className="detail-title-row">
        <div>
          <span className="eyebrow">SELECTED RECIPE</span>
          <h3>Rank {item.rank} · {item.recipe_id}</h3>
          <p>{recipeSummary(item.resolved_recipe)}</p>
        </div>
        <div className="score-orbit"><strong>{item.score.toFixed(3)}</strong><span>score</span></div>
      </div>
      <ReactionImage smiles={precedent} label="Selected precedent reaction" compact />
      <div className="detail-columns">
        <section><h4>Conditions</h4><Conditions recipe={item.resolved_recipe} /></section>
        <section>
          <h4>Evidence</h4>
          <dl className="detail-list">
            <div><dt>Similarity</dt><dd>{item.similarity_score.toFixed(3)}</dd></div>
            <div><dt>Compatibility</dt><dd>{item.compatibility_score.toFixed(3)}</dd></div>
            <div><dt>Reaction support</dt><dd>{item.support}</dd></div>
            <div><dt>Reference support</dt><dd>{item.reference_support}</dd></div>
            <div><dt>Expected yield</dt><dd>{item.expected_yield_pct == null ? 'Unreported' : `${item.expected_yield_pct.toFixed(1)}%`}</dd></div>
          </dl>
        </section>
      </div>
      {protocol && <ProtocolPanel protocol={protocol} rank={item.rank} />}
      <MessageList title="Why this recipe" values={item.explanation} />
      <MessageList title="Compatibility evidence" values={item.compatibility_evidence} />
      <MessageList title="Cautions" values={item.cautions} tone="caution" />
      <details className="trace-panel">
        <summary>Ranking factor trace</summary>
        <div className="factor-grid">
          {Object.entries(item.score_trace.applied_ranking_weights).map(([name, weight]) => {
            const value = item.score_trace.ranking_components[name]
            const contribution = item.score_trace.ranking_contributions[name] ?? 0
            return (
              <div key={name}>
                <strong>{FACTOR_LABELS[name] ?? displayName(name)}</strong>
                <span>Value {value == null ? '—' : value.toFixed(3)}</span>
                <span>Weight {weight.toFixed(3)}</span>
                <span>Contribution {contribution.toFixed(3)}</span>
              </div>
            )
          })}
        </div>
      </details>
      {(item.precedent_reference_ids.length > 0 || item.precedent_reaction_ids.length > 0) && (
        <details className="trace-panel"><summary>Precedent provenance</summary>
          <p className="mono-wrap">References: {item.precedent_reference_ids.join(', ') || 'Unavailable'}</p>
          <p className="mono-wrap">Reactions: {item.precedent_reaction_ids.join(', ') || 'Unavailable'}</p>
        </details>
      )}
    </article>
  )
}

export function RecommendationResults({ result }: { result: RecommendationResult }) {
  const [selected, setSelected] = useState(0)
  useEffect(() => setSelected(0), [result])
  const active = result.recommendations[selected]
  const activeProtocol = active?.synthesis_protocol
  return (
    <section className="results-card">
      <div className="results-summary">
        <div><span className="eyebrow">RECOMMENDATION RESULT</span><h2>{result.recommendations.length} ranked recipe{result.recommendations.length === 1 ? '' : 's'}</h2></div>
        <div className="results-summary-actions">
          {activeProtocol && active && (
            <button className="button secondary" type="button" onClick={() => downloadProtocol(activeProtocol, active.rank)}>Download condition JSON</button>
          )}
          <div className="metric-strip">
            <div><strong>{result.candidate_count}</strong><span>candidates</span></div>
            <div><strong>{result.compatible_candidate_count}</strong><span>compatible</span></div>
            <div><strong>{displayName(result.retrieval_level ?? 'none')}</strong><span>fallback level</span></div>
          </div>
        </div>
      </div>
      <ReactionImage smiles={result.effective_query_reaction_smiles ?? result.query_reaction_smiles} label="Analyzed query reaction" compact />
      {!result.valid && <div className="alert error">{displayName(result.error ?? 'No recommendation')}</div>}
      <MessageList title="Query warnings" values={result.warnings} tone="caution" />
      {result.recommendations.length > 0 && (
        <div className="results-layout">
          <div className="table-scroll">
            <table>
              <thead><tr><th>Rank</th><th>Score</th><th>Similarity</th><th>Compatibility</th><th>Yield</th><th>Support</th><th>Conditions</th></tr></thead>
              <tbody>
                {result.recommendations.map((item, index) => (
                  <tr key={item.recipe_id} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}>
                    <td><strong>{item.rank}</strong>{item.rank_change !== 0 && <span className="rank-change">{item.rank_change > 0 ? '+' : ''}{item.rank_change}</span>}</td>
                    <td>{item.score.toFixed(3)}</td><td>{item.similarity_score.toFixed(3)}</td><td>{item.compatibility_score.toFixed(3)}</td>
                    <td>{item.expected_yield_pct == null ? '—' : `${item.expected_yield_pct.toFixed(1)}%`}</td><td>{item.support}</td><td>{recipeSummary(item.resolved_recipe)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          {active && <RecommendationDetails item={active} />}
        </div>
      )}
    </section>
  )
}

function DiscoveryDetails({ hit }: { hit: DiscoveryHit }) {
  return (
    <article className="result-detail">
      <div className="detail-title-row"><div><span className="eyebrow">SELECTED PRECEDENT</span><h3>Rank {hit.rank} · {displayName(hit.relation_class)}</h3><p>{hit.reference_id || hit.source_dataset}</p></div><div className="score-orbit"><strong>{hit.discovery_score.toFixed(3)}</strong><span>score</span></div></div>
      <ReactionImage smiles={hit.reaction_smiles} label="Selected precedent reaction" compact />
      <div className="detail-columns"><section><h4>Observed conditions</h4><Conditions recipe={hit.resolved_recipe} /></section><section><h4>Provenance</h4><dl className="detail-list"><div><dt>Evidence tier</dt><dd>{displayName(hit.evidence_tier)}</dd></div><div><dt>Chemistry</dt><dd>{displayName(hit.chemistry_status)}</dd></div><div><dt>Outcome</dt><dd>{displayName(hit.outcome_status)}</dd></div><div><dt>Observed yield</dt><dd>{hit.yield_pct == null ? 'Unreported' : `${hit.yield_pct.toFixed(1)}%`}</dd></div></dl></section></div>
      <MessageList title="Why it is related" values={hit.score_trace.matches} />
      <MessageList title="Structural differences" values={hit.score_trace.mismatches} tone="caution" />
      <MessageList title="Insights" values={hit.insights} />
      <MessageList title="Cautions" values={hit.cautions} tone="caution" />
      <details className="trace-panel"><summary>Discovery score trace</summary><div className="factor-grid">
        {Object.entries(hit.score_trace.configured_weights).map(([name, configured]) => {
          const value = hit.score_trace.components[name]
          const effective = hit.score_trace.effective_weights[name]
          const contribution = hit.score_trace.contributions[name]
          return <div key={name}><strong>{FACTOR_LABELS[name] ?? displayName(name)}</strong><span>Value {value == null ? '—' : value.toFixed(3)}</span><span>Configured {configured.toFixed(3)}</span><span>Effective {effective == null ? '—' : effective.toFixed(3)}</span><span>Contribution {contribution == null ? '—' : contribution.toFixed(3)}</span></div>
        })}
      </div></details>
    </article>
  )
}

export function DiscoveryResults({ result }: { result: DiscoveryResult }) {
  const [selected, setSelected] = useState(0)
  useEffect(() => setSelected(0), [result])
  const active = result.hits[selected]
  return (
    <section className="results-card">
      <div className="results-summary"><div><span className="eyebrow">DISCOVERY RESULT</span><h2>{result.hits.length} related precedent{result.hits.length === 1 ? '' : 's'}</h2></div><div className="metric-strip"><div><strong>{result.candidate_count}</strong><span>candidates</span></div><div><strong>{displayName(result.discovery_view)}</strong><span>view</span></div><div><strong>{Object.keys(result.relation_counts).length}</strong><span>relation classes</span></div></div></div>
      <ReactionImage smiles={result.query_reaction_smiles} label="Discovery query reaction" compact />
      {!result.valid && <div className="alert error">{displayName(result.error ?? 'No discovery result')}</div>}
      <MessageList title="Discovery warnings" values={result.warnings} tone="caution" />
      {result.hits.length > 0 && <div className="results-layout"><div className="table-scroll"><table><thead><tr><th>Rank</th><th>Score</th><th>Relationship</th><th>Yield</th><th>Evidence</th><th>Observed conditions</th></tr></thead><tbody>
        {result.hits.map((hit, index) => <tr key={hit.observation_id} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}><td><strong>{hit.rank}</strong></td><td>{hit.discovery_score.toFixed(3)}</td><td>{displayName(hit.relation_class)}</td><td>{hit.yield_pct == null ? '—' : `${hit.yield_pct.toFixed(1)}%`}</td><td>{displayName(hit.evidence_tier)}</td><td>{recipeSummary(hit.resolved_recipe)}</td></tr>)}
      </tbody></table></div>{active && <DiscoveryDetails hit={active} />}</div>}
    </section>
  )
}
