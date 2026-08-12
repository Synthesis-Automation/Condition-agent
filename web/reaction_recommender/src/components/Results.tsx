import { useEffect, useState } from 'react'
import type {
  DiscoveryHit,
  DiscoveryResult,
  ExperimentalDetail,
  MultistepRetrosynthesisResult,
  MultistepRetrosynthesisRoute,
  RecipeComponent,
  Recommendation,
  RecommendationResult,
  ReferenceRecord,
  ResolvedRecipe,
  RetrosynthesisCandidate,
  RetrosynthesisResult,
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

export function compactRecipeSummary(recipe: ResolvedRecipe): string {
  const parts: string[] = []
  for (const [field] of ROLE_FIELDS) {
    const values = recipe[field]
    if (Array.isArray(values)) parts.push(...values.map(componentName))
  }
  if (recipe.temperature_c != null) parts.push(`${recipe.temperature_c} °C`)
  if (recipe.time_h != null) parts.push(`${recipe.time_h} h`)
  return parts.join(' · ') || 'Conditions unavailable'
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

function CitationText({ record }: { record: ReferenceRecord }) {
  const rawCitation = record.raw_reference || record.normalized_citation || ''
  const sections = rawCitation.split('|').map((value) => value.trim()).filter(Boolean)
  const title = sections.length >= 3 ? sections[0] : ''
  const authors = sections.length >= 3 ? sections[1] : ''
  const publication = sections.length >= 3 ? sections.slice(2).join(' | ') : rawCitation
  const publicationMatch = publication.match(/^(.*)\s+\(((?:18|19|20|21)\d{2})\)(.*)$/)
  const terminalPeriod = (value: string) => /[.!?]$/.test(value) ? value : `${value}.`

  if (!rawCitation) {
    return (
      <span className="reference-citation">
        {record.doi || record.patent_number || 'Reference details unavailable'}
        {record.publication_year && <> <strong>{record.publication_year}</strong>.</>}
      </span>
    )
  }

  return (
    <span className="reference-citation">
      {authors && <>{terminalPeriod(authors)} </>}
      {title && <>{terminalPeriod(title)} </>}
      {publicationMatch ? (
        <>
          <em>{publicationMatch[1].trim()}</em>{' '}
          <strong>{publicationMatch[2]}</strong>{publicationMatch[3]}.
        </>
      ) : (
        <>
          {publication}
          {record.publication_year && <> <strong>{record.publication_year}</strong></>}
          {!/[.!?]$/.test(publication) && '.'}
        </>
      )}
    </span>
  )
}

function ReferenceRecords({ records }: { records: ReferenceRecord[] }) {
  if (!records.length) return null
  return (
    <div className="reference-records">
      <h5>Dataset references</h5>
      <ul>
        {records.map((record) => {
          const metadata = [
            record.doi ? `DOI ${record.doi}` : '',
            record.patent_number ? `Patent ${record.patent_number}` : '',
          ].filter(Boolean)
          return (
            <li key={record.reference_id}>
              <CitationText record={record} />
              {metadata.length > 0 && <small>{metadata.join(' · ')}</small>}
            </li>
          )
        })}
      </ul>
    </div>
  )
}

function ExperimentalDetails({ details }: { details: ExperimentalDetail[] }) {
  if (!details.length) return null
  return (
    <details className="trace-panel experimental-panel">
      <summary>
        Experimental details{details.length > 1 ? ` (${details.length} precedents)` : ''}
      </summary>
      <div className="experimental-detail-list">
        {details.map((detail, index) => (
          <article key={detail.observation_id || detail.reaction_id || index}>
            {details.length > 1 && <h5>Precedent {index + 1}</h5>}
            <p>{detail.procedure_text}</p>
            {detail.notes && <small>Notes: {detail.notes}</small>}
          </article>
        ))}
      </div>
    </details>
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
    <details className="trace-panel protocol-panel">
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
          <h3>Rank {item.rank}</h3>
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
            <div><dt>Dataset support</dt><dd>{item.dataset_support}</dd></div>
            <div><dt>Expected yield</dt><dd>{item.expected_yield_pct == null ? 'Unreported' : `${item.expected_yield_pct.toFixed(1)}%`}</dd></div>
          </dl>
          <ReferenceRecords records={item.precedent_references ?? []} />
        </section>
      </div>
      <ExperimentalDetails details={item.precedent_experimental_details ?? []} />
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
      {!result.valid && <div className="alert error">{displayName(result.error ?? 'No recommendation')}</div>}
      <MessageList title="Query warnings" values={result.warnings} tone="caution" />
      {result.recommendations.length > 0 && (
        <div className="results-layout">
          <div className="table-scroll">
            <table>
              <thead><tr><th>Rank</th><th>Score</th><th>Similarity</th><th>Yield</th><th>Conditions</th></tr></thead>
              <tbody>
                {result.recommendations.map((item, index) => (
                  <tr key={item.recipe_id} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}>
                    <td><strong>{item.rank}</strong>{item.rank_change !== 0 && <span className="rank-change">{item.rank_change > 0 ? '+' : ''}{item.rank_change}</span>}</td>
                    <td>{item.score.toFixed(3)}</td><td>{item.similarity_score.toFixed(3)}</td>
                    <td>{item.expected_yield_pct == null ? '—' : `${item.expected_yield_pct.toFixed(1)}%`}</td><td>{compactRecipeSummary(item.resolved_recipe)}</td>
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
      <div className="detail-title-row"><div><span className="eyebrow">SELECTED PRECEDENT</span><h3>Rank {hit.rank} · {displayName(hit.relation_class)}</h3><p>Source: {displayName(hit.source_dataset)}</p></div><div className="score-orbit"><strong>{hit.discovery_score.toFixed(3)}</strong><span>score</span></div></div>
      <ReactionImage smiles={hit.reaction_smiles} label="Selected precedent reaction" compact />
      <div className="detail-columns"><section><h4>Observed conditions</h4><Conditions recipe={hit.resolved_recipe} /></section><section><h4>Provenance</h4><dl className="detail-list"><div><dt>Dataset</dt><dd>{displayName(hit.source_dataset)}</dd></div><div><dt>Evidence tier</dt><dd>{displayName(hit.evidence_tier)}</dd></div><div><dt>Chemistry</dt><dd>{displayName(hit.chemistry_status)}</dd></div><div><dt>Outcome</dt><dd>{displayName(hit.outcome_status)}</dd></div><div><dt>Observed yield</dt><dd>{hit.yield_pct == null ? 'Unreported' : `${hit.yield_pct.toFixed(1)}%`}</dd></div></dl><ReferenceRecords records={hit.reference_record ? [hit.reference_record] : []} /></section></div>
      <ExperimentalDetails details={hit.experimental_detail ? [hit.experimental_detail] : []} />
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
      {result.hits.length > 0 && <div className="results-layout"><div className="table-scroll"><table><thead><tr><th>Rank</th><th>Score</th><th>Relationship</th><th>Yield</th><th>Conditions</th></tr></thead><tbody>
        {result.hits.map((hit, index) => <tr key={hit.observation_id} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}><td><strong>{hit.rank}</strong></td><td>{hit.discovery_score.toFixed(3)}</td><td>{displayName(hit.relation_class)}</td><td>{hit.yield_pct == null ? '—' : `${hit.yield_pct.toFixed(1)}%`}</td><td>{compactRecipeSummary(hit.resolved_recipe)}</td></tr>)}
      </tbody></table></div>{active && <DiscoveryDetails hit={active} />}</div>}
    </section>
  )
}

function RetrosynthesisDetails({
  candidate,
  scopeWarnings,
}: {
  candidate: RetrosynthesisCandidate
  scopeWarnings: string[]
}) {
  const conditionRecommendations = candidate.condition_evidence?.recommendations ?? []
  const combinedCautions = Array.from(new Set([
    ...scopeWarnings,
    ...(candidate.condition_evidence?.warnings ?? []),
  ]))
  return (
    <article className="result-detail retrosynthesis-detail">
      <div className="detail-title-row">
        <div>
          <span className="eyebrow">SELECTED DISCONNECTION</span>
          <h3>Rank {candidate.rank} · {displayName(candidate.transformation_kind ?? 'graph operator')}</h3>
          <p>{candidate.precursor_smiles}</p>
        </div>
        <div className="score-orbit"><strong>{candidate.score.toFixed(3)}</strong><span>score</span></div>
      </div>
      <ReactionImage smiles={candidate.proposed_reaction_smiles} label="Proposed single-step retrosynthesis" compact />
      <MessageList
        title="Functional-group competition"
        values={(candidate.selectivity_warnings ?? []).map((warning) => warning.message)}
        tone="caution"
      />
      {(candidate.selectivity_warnings ?? []).map((warning) => (
        <details className="trace-panel" key={warning.code}>
          <summary>Review competing structural outcomes ({warning.competing_outcomes.length})</summary>
          <div className="supporting-precedent-list">
            {warning.competing_outcomes.map((outcome) => (
              <article key={outcome.candidate_id}>
                <strong>{outcome.element} endpoint · atom {outcome.atom_index}</strong>
                <ReactionImage smiles={outcome.product_smiles} label={`Possible competing ${outcome.element} outcome`} compact />
                <small>{displayName(outcome.site_type)} · {displayName(outcome.availability)}</small>
              </article>
            ))}
          </div>
        </details>
      ))}
      <div className="detail-columns">
        <section>
          <h4>Structural evidence</h4>
          <dl className="detail-list">
            <div><dt>Abstraction level</dt><dd>{candidate.abstraction_level}</dd></div>
            <div><dt>Product similarity</dt><dd>{candidate.product_similarity.toFixed(3)}</dd></div>
            <div><dt>Precursor similarity</dt><dd>{candidate.precursor_similarity.toFixed(3)}</dd></div>
            <div><dt>Context similarity</dt><dd>{candidate.context_similarity.toFixed(3)}</dd></div>
            <div><dt>Template specificity</dt><dd>{candidate.template_specificity.toFixed(3)}</dd></div>
            <div><dt>Independent references</dt><dd>{candidate.independent_reference_support}</dd></div>
            <div><dt>Forward validation</dt><dd>{displayName(candidate.forward_validation_status)}</dd></div>
          </dl>
        </section>
      </div>
      <details className="trace-panel retrosynthesis-conditions" open>
        <summary>
          Recommended conditions
          {conditionRecommendations.length > 0 ? ` (${conditionRecommendations.length})` : ''}
        </summary>
        {conditionRecommendations.length > 0 ? (
          <div className="retrosynthesis-condition-list">
            {conditionRecommendations.map((recommendation) => (
              <article key={recommendation.recipe_id}>
                <div className="retrosynthesis-condition-heading">
                  <div><strong>Recipe rank {recommendation.rank}</strong><span>{compactRecipeSummary(recommendation.resolved_recipe)}</span></div>
                  <small>Score {recommendation.score.toFixed(3)}{recommendation.expected_yield_pct == null ? '' : ` · ${recommendation.expected_yield_pct.toFixed(1)}% expected yield`}</small>
                </div>
                <Conditions recipe={recommendation.resolved_recipe} />
                <ReferenceRecords records={recommendation.precedent_references ?? []} />
                <ExperimentalDetails details={recommendation.precedent_experimental_details ?? []} />
              </article>
            ))}
          </div>
        ) : (
          <p className="retrosynthesis-empty-evidence">
            {candidate.condition_evidence?.status === 'pending'
              ? 'Loading condition recommendations…'
              : 'No compatible condition recipe was found for this proposed reaction.'}
          </p>
        )}
      </details>
      {candidate.supporting_precedents.length > 0 && (
        <details className="trace-panel supporting-precedents">
          <summary>Supporting precedent reactions ({candidate.supporting_precedents.length})</summary>
          <div className="supporting-precedent-list">
            {candidate.supporting_precedents.map((precedent, index) => (
              <article key={`${precedent.reaction_smiles}:${index}`}>
                <ReactionImage smiles={precedent.reaction_smiles} label={`Supporting precedent reaction ${index + 1}`} compact />
                <ReferenceRecords records={precedent.reference_record ? [precedent.reference_record] : []} />
                {!precedent.reference_record && <p className="retrosynthesis-empty-evidence">Publication details unavailable for this precedent.</p>}
              </article>
            ))}
          </div>
        </details>
      )}
      <details className="trace-panel">
        <summary>Ranking trace</summary>
        <div className="factor-grid">
          <div><strong>Policy</strong><span>{candidate.ranking_policy_definition_id || 'Default structural ranking'}</span></div>
          <div><strong>Structural rank</strong><span>{candidate.pre_diversity_rank || candidate.rank}</span></div>
          <div><strong>Diversity rank</strong><span>{candidate.diversity_rank || candidate.rank}</span></div>
          <div><strong>Score band</strong><span>{candidate.structural_score_band}</span></div>
        </div>
      </details>
      <details className="trace-panel machine-code-panel">
        <summary>Operator identity</summary>
        <dl className="detail-list">
          <div><dt>Operator</dt><dd className="mono-value">{candidate.operator_id || 'Unavailable'}</dd></div>
          <div><dt>Realization</dt><dd className="mono-value">{candidate.realization_id || 'Unavailable'}</dd></div>
          <div><dt>Site</dt><dd className="mono-value">{candidate.disconnection_site_key || 'Unavailable'}</dd></div>
          <div><dt>Synthon</dt><dd className="mono-value">{candidate.synthon_signature || 'Unavailable'}</dd></div>
        </dl>
      </details>
      <MessageList title="Scope and condition cautions" values={combinedCautions} tone="caution" />
    </article>
  )
}

export function RetrosynthesisResults({ result }: { result: RetrosynthesisResult }) {
  const [selected, setSelected] = useState(0)
  useEffect(() => setSelected(0), [result])
  const active = result.candidates[selected]
  return (
    <section className="results-card">
      <div className="results-summary">
        <div><span className="eyebrow">RETROSYNTHESIS RESULT</span><h2>{result.candidates.length} validated disconnection{result.candidates.length === 1 ? '' : 's'}</h2></div>
        <div className="metric-strip">
          <div><strong>{result.candidate_count}</strong><span>candidates</span></div>
          <div><strong>{result.library_operator_count}</strong><span>operators</span></div>
          <div><strong>{result.library_template_count}</strong><span>templates</span></div>
        </div>
      </div>
      {!result.valid && <div className="alert error">{displayName(result.error ?? 'No retrosynthesis candidates')}</div>}
      {result.candidates.length === 0 && <MessageList title="Scope and cautions" values={result.warnings} tone="caution" />}
      {result.candidates.length > 0 && (
        <div className="results-layout">
          <div className="table-scroll">
            <table>
              <thead><tr><th>Rank</th><th>Score</th><th>Level</th><th>Transformation</th><th>Conditions</th></tr></thead>
              <tbody>
                {result.candidates.map((candidate, index) => (
                  <tr key={`${candidate.template_id}:${candidate.precursor_smiles}`} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}>
                    <td><strong>{candidate.rank}</strong></td>
                    <td>{candidate.score.toFixed(3)}</td>
                    <td>{candidate.abstraction_level}</td>
                    <td>{displayName(candidate.transformation_kind ?? 'graph operator')}</td>
                    <td>{candidate.condition_evidence?.recommendations?.[0] ? compactRecipeSummary(candidate.condition_evidence.recommendations[0].resolved_recipe) : 'No compatible conditions'}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          {active && <RetrosynthesisDetails candidate={active} scopeWarnings={result.warnings} />}
        </div>
      )}
    </section>
  )
}

function MultistepRouteDetails({
  route,
}: {
  route: MultistepRetrosynthesisRoute
}) {
  return (
    <article className="result-detail multistep-route-detail">
      <div className="detail-title-row">
        <div>
          <span className="eyebrow">{route.solved ? 'SOLVED ROUTE' : 'PARTIAL ROUTE'}</span>
          <h3>{route.reaction_count} reaction{route.reaction_count === 1 ? '' : 's'} · depth {route.maximum_depth}</h3>
          <p>{route.solved ? 'Every leaf satisfies the configured terminal rule.' : 'At least one branch remains unresolved within the search limits.'}</p>
        </div>
        <div className="score-orbit"><strong>{route.route_cost.toFixed(3)}</strong><span>cost</span></div>
      </div>

      <section className="route-step-list" aria-label="Retrosynthetic route steps">
        {route.steps.map((step, index) => (
          <article className="route-step-card" key={step.step_id}>
            <div className="route-step-heading">
              <div><span>Step {index + 1} · depth {step.depth}</span><strong>{displayName(step.candidate.transformation_kind ?? 'graph operator')}</strong></div>
              <small>{step.candidate.abstraction_level} · score {step.candidate.score.toFixed(3)}</small>
            </div>
            <ReactionImage smiles={step.candidate.proposed_reaction_smiles} label={`Route reaction ${index + 1}`} compact />
            <dl className="detail-list">
              <div><dt>Precursors</dt><dd>{step.precursor_smiles.join(' · ')}</dd></div>
              <div><dt>Forward validation</dt><dd>{displayName(step.candidate.forward_validation_status)}</dd></div>
              <div><dt>Independent support</dt><dd>{step.candidate.independent_reference_support}</dd></div>
              <div><dt>Operator</dt><dd className="mono-value">{step.candidate.operator_id}</dd></div>
            </dl>
          </article>
        ))}
      </section>

      <section className="route-leaf-section">
        <h4>Starting-material leaves</h4>
        <div className="route-leaf-grid">
          {route.leaves.map((leaf, index) => (
            <article className={leaf.terminal ? 'route-leaf terminal' : 'route-leaf unresolved'} key={`${leaf.canonical_smiles}:${leaf.depth}:${index}`}>
              <div className="route-leaf-heading">
                <strong>{leaf.canonical_smiles}</strong>
                <span>{leaf.terminal ? 'Terminal' : 'Unresolved'}</span>
              </div>
              <small>MW {leaf.molecular_weight == null ? '—' : leaf.molecular_weight.toFixed(2)} · depth {leaf.depth}</small>
              {leaf.terminal_reasons.length > 0 && <p>{leaf.terminal_reasons.map(displayName).join(' · ')}</p>}
              {leaf.unresolved_reason && <p>{displayName(leaf.unresolved_reason)}</p>}
              {(leaf.literature_match?.source_records ?? []).map((record, recordIndex) => (
                <div className="route-leaf-source" key={`${record.reaction_id ?? 'source'}:${recordIndex}`}>
                  {record.cas_no && <strong>CAS {record.cas_no}</strong>}
                  {record.citation && <span>{record.citation}</span>}
                  {record.reaction_id && <small>{record.reaction_id}</small>}
                </div>
              ))}
            </article>
          ))}
        </div>
      </section>
      <MessageList title="Route cautions" values={route.warnings} tone="caution" />
    </article>
  )
}

export function MultistepRetrosynthesisResults({ result }: { result: MultistepRetrosynthesisResult }) {
  const [selected, setSelected] = useState(0)
  useEffect(() => setSelected(0), [result])
  const availableRoutes = [...result.routes, ...result.partial_routes]
  const active = availableRoutes[selected]
  return (
    <section className="results-card">
      <div className="results-summary">
        <div><span className="eyebrow">MULTI-STEP RETROSYNTHESIS</span><h2>{result.route_count} solved route{result.route_count === 1 ? '' : 's'}</h2></div>
        <div className="metric-strip">
          <div><strong>{result.route_count}</strong><span>solved</span></div>
          <div><strong>{result.partial_route_count}</strong><span>partial</span></div>
          <div><strong>{result.diagnostics.expanded_states}</strong><span>expanded</span></div>
          <div><strong>{result.max_depth}</strong><span>max depth</span></div>
        </div>
      </div>
      {!result.valid && <div className="alert error">{displayName(result.error ?? 'No solved multistep routes')}</div>}
      {availableRoutes.length === 0 && <MessageList title="Scope and cautions" values={result.warnings} tone="caution" />}
      {availableRoutes.length > 0 && (
        <div className="results-layout multistep-results-layout">
          <div className="table-scroll">
            <table>
              <thead><tr><th>Route</th><th>Status</th><th>Steps</th><th>Depth</th><th>Cost</th><th>Terminal leaves</th></tr></thead>
              <tbody>
                {availableRoutes.map((route, index) => (
                  <tr key={route.route_id} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}>
                    <td><strong>{index + 1}</strong></td>
                    <td>{route.solved ? 'Solved' : 'Partial'}</td>
                    <td>{route.reaction_count}</td>
                    <td>{route.maximum_depth}</td>
                    <td>{route.route_cost.toFixed(3)}</td>
                    <td>{route.leaves.filter((leaf) => leaf.terminal).length}/{route.leaves.length}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          {active && <MultistepRouteDetails route={active} />}
        </div>
      )}
      <MessageList title="Search cautions" values={result.warnings} tone="caution" />
    </section>
  )
}
