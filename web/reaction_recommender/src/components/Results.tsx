import { useEffect, useState } from 'react'
import type {
  ConditionPrecedent,
  DiscoveryHit,
  DiscoveryResult,
  ExperimentalDetail,
  ForwardProductCandidate,
  ForwardSynthesisResult,
  MultistepRetrosynthesisResult,
  MultistepRetrosynthesisRoute,
  RecipeComponent,
  Recommendation,
  RecommendationResult,
  ReferenceRecord,
  ResolvedRecipe,
  RetrosynthesisCandidate,
  RetrosynthesisForwardAssessment,
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
      component.cas ??
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

function ForwardValidityAudit({ audit }: { audit?: RetrosynthesisForwardAssessment | null }) {
  if (!audit?.evaluated) return null
  const caution = audit.validity !== 'structurally_supported'
  return (
    <details className="trace-panel" open>
      <summary>Independent forward validity audit · {displayName(audit.validity)}</summary>
      <div className={`alert ${caution ? 'caution' : ''}`}>
        This audit is advisory. It independently replays the proposed operator and runs a target-blind competing-product search.
      </div>
      <dl className="detail-list">
        <div><dt>Targeted operator replay</dt><dd>{displayName(audit.targeted_replay_status)}</dd></div>
        <div><dt>Blind target recovery</dt><dd>{audit.intended_product_rank == null ? 'Not recovered' : `Rank ${audit.intended_product_rank} · ${displayName(audit.intended_match)}`}</dd></div>
        <div><dt>Validated pathways</dt><dd>{audit.blind_prediction_summary.valid_pathway_count}</dd></div>
        <div><dt>Generated products</dt><dd>{audit.blind_prediction_summary.candidate_count}</dd></div>
        <div><dt>Condition basis</dt><dd>{audit.blind_prediction_summary.conditions_supplied ? 'Canonical recipe applied' : audit.blind_prediction_summary.condition_profile_supplied ? 'Coarse condition profile applied' : 'Unconditioned structural possibilities'}</dd></div>
        <div><dt>Best competitor</dt><dd className="mono-value">{audit.best_competitor_product ?? 'None generated'}</dd></div>
        <div><dt>Target margin</dt><dd>{audit.score_margin == null ? 'Not available' : audit.score_margin.toFixed(3)}</dd></div>
      </dl>
      <div className="table-scroll">
        <table>
          <thead><tr><th>Check</th><th>Status</th><th>Evidence</th></tr></thead>
          <tbody>{audit.checks.map((check) => (
            <tr key={check.check_id}>
              <td>{displayName(check.check_id)}</td>
              <td>{displayName(check.status)}</td>
              <td>{check.detail}</td>
            </tr>
          ))}</tbody>
        </table>
      </div>
      {audit.blind_prediction_summary.top_products.length > 0 && (
        <details className="trace-panel">
          <summary>Blind forward products ({audit.blind_prediction_summary.top_products.length} shown)</summary>
          <div className="table-scroll"><table>
            <thead><tr><th>Rank</th><th>Score</th><th>Product</th><th>Role</th></tr></thead>
            <tbody>{audit.blind_prediction_summary.top_products.map((product) => (
              <tr key={`${product.rank}:${product.product_smiles}`}>
                <td>{product.rank}</td><td>{product.score.toFixed(3)}</td><td className="mono-value">{product.product_smiles}</td><td>{product.is_intended ? 'Intended' : 'Competitor'}</td>
              </tr>
            ))}</tbody>
          </table></div>
        </details>
      )}
      <MessageList title="Forward-audit cautions" values={[...audit.warnings, ...audit.blind_prediction_summary.warnings]} tone="caution" />
    </details>
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

function ConditionPrecedents({ precedents }: { precedents: ConditionPrecedent[] }) {
  if (!precedents.length) return null
  return (
    <section className="condition-precedents" aria-label="Condition precedent reactions">
      <h5>Condition precedent reactions</h5>
      <div className="condition-precedent-list">
        {precedents.map((precedent, index) => (
          <article key={precedent.observation_id || precedent.reaction_id || `${precedent.reaction_smiles}:${index}`}>
            <div className="condition-precedent-heading">
              <strong>Precedent {index + 1}</strong>
              {precedent.reaction_id && <small>{precedent.reaction_id}</small>}
            </div>
            <ReactionImage
              smiles={precedent.reaction_smiles}
              label={`Condition precedent reaction ${index + 1}`}
              compact
            />
            <ReferenceRecords records={precedent.reference_record ? [precedent.reference_record] : []} />
            {!precedent.reference_record && (
              <p className="retrosynthesis-empty-evidence">
                Publication details unavailable for this condition precedent.
              </p>
            )}
            <ExperimentalDetails details={precedent.experimental_detail ? [precedent.experimental_detail] : []} />
          </article>
        ))}
      </div>
    </section>
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
      <ForwardValidityAudit audit={candidate.forward_assessment} />
      {candidate.strategic_complexity && (
        <details className="trace-panel" open>
          <summary>Strategic complexity reduction · {(100 * candidate.strategic_complexity_score).toFixed(1)}/100 · {displayName(candidate.strategic_class)}</summary>
          <dl className="detail-list">
            <div><dt>Strategic scaffold step</dt><dd>{candidate.strategic_candidate ? 'Yes' : 'No — tactical or insufficient simplification'}</dd></div>
            <div><dt>Retained target core</dt><dd>{(100 * candidate.strategic_complexity.largest_retained_core_fraction).toFixed(1)}%</dd></div>
            <div><dt>Product-derived fragments</dt><dd>{candidate.strategic_complexity.product_derived_component_heavy_atom_counts.join(' + ') || 'Mapping unavailable'} heavy atoms</dd></div>
            <div><dt>Core fragmentation</dt><dd>{candidate.strategic_complexity.core_fragmentation_score.toFixed(3)}</dd></div>
            <div><dt>Ring-topology reduction</dt><dd>{candidate.strategic_complexity.ring_topology_reduction_score.toFixed(3)}</dd></div>
            <div><dt>Graph-complexity reduction</dt><dd>{candidate.strategic_complexity.graph_complexity_reduction_fraction.toFixed(3)}</dd></div>
            <div><dt>Convergency</dt><dd>{candidate.strategic_complexity.convergency_score.toFixed(3)}</dd></div>
            <div><dt>Tactical penalty</dt><dd>−{candidate.strategic_complexity.tactical_penalty.toFixed(3)}</dd></div>
            <div><dt>Target graph</dt><dd>{candidate.strategic_complexity.target.heavy_atom_count} heavy atoms · cycle rank {candidate.strategic_complexity.target.cycle_rank} · raw complexity {candidate.strategic_complexity.target.raw_complexity.toFixed(2)}</dd></div>
          </dl>
          <MessageList title="Strategic audit cautions" values={candidate.strategic_complexity.warnings} tone="caution" />
        </details>
      )}
      <MessageList
        title="Strong intramolecular compatibility warning"
        values={(candidate.precursor_compatibility_assessments ?? []).map((assessment) => assessment.message)}
        tone="caution"
      />
      {(candidate.precursor_compatibility_assessments ?? []).length > 0 && (
        <details className="trace-panel" open>
          <summary>Review incompatible sites in one molecular component</summary>
          {(candidate.precursor_compatibility_assessments ?? []).map((assessment) => (
            <dl className="detail-list" key={assessment.assessment_id}>
              <div><dt>Rule</dt><dd className="mono-value">{assessment.rule_id}</dd></div>
              <div><dt>Component</dt><dd className="mono-value">{assessment.component_smiles}</dd></div>
              <div><dt>Reactive pair</dt><dd>{assessment.left_site.chemist_label} + {assessment.right_site.chemist_label}</dd></div>
              <div><dt>Severity</dt><dd>{displayName(assessment.intrinsic_severity)} · {displayName(assessment.warning_strength)} warning</dd></div>
              {assessment.graph_distance != null && <div><dt>Graph distance</dt><dd>{assessment.graph_distance} bonds</dd></div>}
              {assessment.potential_closure_ring_size != null && <div><dt>Possible closure</dt><dd>{assessment.potential_closure_ring_size}-membered ring</dd></div>}
              <div><dt>Scope</dt><dd>Same molecular component</dd></div>
            </dl>
          ))}
        </details>
      )}
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
            <div><dt>Signature sanity check</dt><dd>{displayName(candidate.forward_validation_status)}</dd></div>
            <div><dt>Strategic reduction</dt><dd>{(100 * candidate.strategic_complexity_score).toFixed(1)}/100 · {displayName(candidate.strategic_class)}</dd></div>
            {candidate.precursor_realism_score != null && <div><dt>Precursor realism</dt><dd>{candidate.precursor_realism_score.toFixed(3)}</dd></div>}
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
                <ConditionPrecedents precedents={recommendation.condition_precedents ?? []} />
                {recommendation.condition_precedents == null && (
                  <>
                    <ReferenceRecords records={recommendation.precedent_references ?? []} />
                    <ExperimentalDetails details={recommendation.precedent_experimental_details ?? []} />
                  </>
                )}
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
          <summary>Disconnection-template precedents ({candidate.supporting_precedents.length})</summary>
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
          {candidate.precursor_compatibility_policy_definition_id && <div><strong>Compatibility policy</strong><span>{candidate.precursor_compatibility_policy_definition_id}</span></div>}
          {candidate.strategic_complexity && <div><strong>Strategic-complexity policy</strong><span>{candidate.strategic_complexity.definition_id}</span></div>}
          {candidate.strategic_reserve_selected && <div><strong>Strategic reserve</strong><span>Retained as a bounded scaffold-level alternative</span></div>}
          {(candidate.precursor_compatibility_band_penalty ?? 0) > 0 && <div><strong>Compatibility action</strong><span>{displayName(candidate.precursor_compatibility_disposition ?? 'demote')} · band +{candidate.precursor_compatibility_band_penalty}</span></div>}
          {candidate.hierarchical_ranking_definition_id && <div><strong>Hierarchical policy</strong><span>{candidate.hierarchical_ranking_definition_id}</span></div>}
          {candidate.hierarchical_rank > 0 && <div><strong>Hierarchical rerank</strong><span>{candidate.pre_hierarchical_rank} → {candidate.hierarchical_rank}</span></div>}
          {candidate.hierarchical_site_rank > 0 && <div><strong>SITE1 rank / score</strong><span>{candidate.hierarchical_site_rank} / {candidate.hierarchical_site_score.toFixed(3)}</span></div>}
          {candidate.hierarchical_synthon_rank > 0 && <div><strong>SYN1 rank / score</strong><span>{candidate.hierarchical_synthon_rank} / {candidate.hierarchical_synthon_score.toFixed(3)}</span></div>}
          {candidate.hierarchical_realization_rank > 0 && <div><strong>REAL1 rank / score</strong><span>{candidate.hierarchical_realization_rank} / {candidate.hierarchical_realization_score.toFixed(3)}</span></div>}
          {candidate.completion_group_id && <div><strong>Completion group</strong><span>{candidate.completion_group_id}</span></div>}
          {candidate.completion_group_id && <div><strong>Completion prior</strong><span>{candidate.completion_prior_probability == null ? 'Unavailable' : candidate.completion_prior_probability.toFixed(3)} ({displayName(candidate.completion_prior_backoff_level)}, {candidate.completion_prior_independent_support}/{candidate.completion_prior_total_support} support)</span></div>}
          {candidate.precursor_realism_score != null && <div><strong>Realism rerank</strong><span>{candidate.pre_realism_rank} → {candidate.precursor_realism_rank}</span></div>}
          {candidate.precursor_realism_aggregation && <div><strong>Weakest precursor</strong><span>{candidate.precursor_realism_aggregation.weakest_component_score.toFixed(3)}</span></div>}
          {candidate.precursor_realism_aggregation && <div><strong>Known substantial component bonus</strong><span>+{candidate.precursor_realism_aggregation.known_substantial_component_bonus.toFixed(3)}</span></div>}
          {candidate.precursor_realism_score != null && <div><strong>Realism band penalty</strong><span>+{candidate.precursor_realism_band_penalty}</span></div>}
          {(candidate.precursor_realism_score != null || (candidate.precursor_compatibility_band_penalty ?? 0) > 0) && <div><strong>Effective band</strong><span>{candidate.effective_structural_score_band}</span></div>}
        </div>
        {(candidate.precursor_realism_assessments ?? []).map((assessment) => (
          <dl className="detail-list" key={assessment.canonical_smiles}>
            <div><dt>Precursor</dt><dd className="mono-value">{assessment.canonical_smiles}</dd></div>
            <div><dt>Realism tier</dt><dd>{displayName(assessment.evidence_tier)}</dd></div>
            <div><dt>Molecular weight</dt><dd>{assessment.molecular_weight.toFixed(2)} Da</dd></div>
            <div><dt>MW penalty</dt><dd>−{assessment.molecular_weight_penalty.toFixed(3)}</dd></div>
            <div><dt>Evidence</dt><dd>{[
              assessment.evidence.buyable && 'buyable',
              assessment.evidence.in_compound_registry && 'registry',
              assessment.evidence.in_literature && 'literature',
            ].filter(Boolean).join(', ') || 'none'}</dd></div>
          </dl>
        ))}
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
          <div><strong>{result.strategic_candidate_count}</strong><span>strategic</span></div>
          <div><strong>{(result.forward_validity_counts.structurally_supported ?? 0) + (result.forward_validity_counts.structurally_supported_with_competition ?? 0)}</strong><span>forward supported</span></div>
        </div>
      </div>
      {!result.valid && <div className="alert error">{displayName(result.error ?? 'No retrosynthesis candidates')}</div>}
      {result.candidates.length === 0 && <MessageList title="Scope and cautions" values={result.warnings} tone="caution" />}
      {result.candidates.length > 0 && (
        <div className="results-layout">
          <div className="table-scroll">
            <table>
              <thead><tr><th>Rank</th><th>Score</th><th>Strategic</th>{result.precursor_realism_enabled && <th>Realism</th>}<th>Forward audit</th><th>Level</th><th>Transformation</th><th>Conditions</th></tr></thead>
              <tbody>
                {result.candidates.map((candidate, index) => (
                  <tr key={`${candidate.template_id}:${candidate.precursor_smiles}`} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}>
                    <td><strong>{candidate.rank}</strong></td>
                    <td>{candidate.score.toFixed(3)}</td>
                    <td>{(100 * candidate.strategic_complexity_score).toFixed(1)} · {displayName(candidate.strategic_class)}</td>
                    {result.precursor_realism_enabled && <td>{candidate.precursor_realism_score?.toFixed(3) ?? '—'}</td>}
                    <td>{candidate.forward_assessment ? displayName(candidate.forward_assessment.validity) : 'Not evaluated'}</td>
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

function ForwardProductDetails({
  candidate,
}: {
  candidate: ForwardProductCandidate
}) {
  const alternativeCount = candidate.alternative_pathway_ids.length
  return (
    <article className="result-detail">
      <div className="detail-title-row">
        <div>
          <span className="eyebrow">PREDICTED PRODUCT</span>
          <h3>Rank {candidate.rank} · {candidate.abstraction_level}</h3>
          <p className="mono-value">{candidate.product_smiles}</p>
        </div>
        <div className="score-orbit"><strong>{candidate.score.toFixed(3)}</strong><span>score</span></div>
      </div>
      <ReactionImage smiles={candidate.reaction_smiles} label={`Forward product rank ${candidate.rank}`} />
      <dl className="detail-list">
        <div><dt>Reverse round trip</dt><dd>{candidate.reverse_round_trip ? 'Passed' : 'Failed'}</dd></div>
        <div><dt>Observed edit agreement</dt><dd>{candidate.operator_edit_agreement ? 'Passed' : 'Failed'}</dd></div>
        <div><dt>Independent references</dt><dd>{candidate.independent_reference_support}</dd></div>
        <div><dt>Observation support</dt><dd>{candidate.observation_support}</dd></div>
        <div><dt>Alternative pathways</dt><dd>{alternativeCount}</dd></div>
        <div><dt>Structural score band</dt><dd>{candidate.structural_score_band}</dd></div>
        <div><dt>Condition compatibility</dt><dd>{candidate.recipe_evidence.evaluated ? `${candidate.recipe_evidence.compatible ? 'Compatible' : 'Conflict'}${candidate.recipe_evidence.score == null ? '' : ` · ${candidate.recipe_evidence.score.toFixed(3)}`}` : 'Not evaluated'}</dd></div>
        <div><dt>Condition-profile adjustment</dt><dd>{candidate.condition_profile_evidence.evaluated ? `${candidate.condition_profile_evidence.score_adjustment >= 0 ? '+' : ''}${candidate.condition_profile_evidence.score_adjustment.toFixed(3)}` : 'Not applied'}</dd></div>
        <div><dt>Participating components</dt><dd>{candidate.participating_component_indices.map((value) => value + 1).join(', ')}</dd></div>
        <div><dt>Stoichiometric input</dt><dd>{candidate.reactant_stoichiometry.map(([index, count]) => `component ${index + 1} × ${count}`).join(' · ')}</dd></div>
        <div><dt>Self-reaction assumption</dt><dd>{candidate.uses_virtual_copies ? 'Yes — multiple equivalents assumed' : 'No'}</dd></div>
      </dl>
      <details className="trace-panel" open>
        <summary>Ranking evidence</summary>
        <dl className="detail-list">
          {Object.entries(candidate.score_components).map(([name, value]) => (
            <div key={name}><dt>{displayName(name)}</dt><dd>{value.toFixed(3)}</dd></div>
          ))}
        </dl>
      </details>
      {candidate.condition_profile_evidence.evaluated && (
        <details className="trace-panel" open>
          <summary>Condition-profile evidence</summary>
          <dl className="detail-list">
            <div><dt>Strategy</dt><dd>{displayName(candidate.condition_profile_evidence.profile.strategy)}</dd></div>
            <div><dt>Catalyst family</dt><dd>{displayName(candidate.condition_profile_evidence.profile.catalyst_family)}</dd></div>
            <div><dt>Redox environment</dt><dd>{displayName(candidate.condition_profile_evidence.profile.redox_mode)}</dd></div>
            <div><dt>Medium</dt><dd>{displayName(candidate.condition_profile_evidence.profile.medium)}</dd></div>
            <div><dt>Matched rules</dt><dd>{candidate.condition_profile_evidence.matched_rules.map(displayName).join(' · ') || 'No structure-only ranking rule'}</dd></div>
            <div><dt>Score adjustment</dt><dd>{candidate.condition_profile_evidence.score_adjustment >= 0 ? '+' : ''}{candidate.condition_profile_evidence.score_adjustment.toFixed(3)}</dd></div>
          </dl>
        </details>
      )}
      <details className="trace-panel">
        <summary>Graph and operator trace</summary>
        <dl className="detail-list">
          <div><dt>Operator</dt><dd className="mono-value">{candidate.operator_id}</dd></div>
          <div><dt>Forward operator</dt><dd className="mono-value">{candidate.forward_operator_id}</dd></div>
          <div><dt>Realization</dt><dd className="mono-value">{candidate.realization_id || 'Unavailable'}</dd></div>
          <div><dt>Template</dt><dd className="mono-value">{candidate.template_id}</dd></div>
          <div><dt>Reaction signature</dt><dd className="mono-value">{candidate.reaction_signature_id || 'Unavailable'}</dd></div>
          <div><dt>Observed edits</dt><dd className="mono-value">{candidate.observed_edit_tokens.join(' · ') || 'Unavailable'}</dd></div>
          <div><dt>Atom correspondence</dt><dd>{candidate.atom_correspondence.length} product atoms traced to precursors</dd></div>
        </dl>
      </details>
      {alternativeCount > 0 && (
        <details className="trace-panel">
          <summary>Alternative pathways ({alternativeCount})</summary>
          <dl className="detail-list">
            <div><dt>Operators</dt><dd className="mono-value">{candidate.alternative_operator_ids.join(' · ') || 'Same operator, different assignment'}</dd></div>
            <div><dt>Templates</dt><dd className="mono-value">{candidate.alternative_template_ids.join(' · ') || 'Same template'}</dd></div>
          </dl>
        </details>
      )}
      <MessageList title="Condition cautions" values={[...candidate.recipe_evidence.hard_conflicts, ...candidate.recipe_evidence.cautions]} tone="caution" />
      <MessageList title="Condition-profile cautions" values={candidate.condition_profile_evidence.cautions} tone="caution" />
      <MessageList title="Prediction warnings" values={candidate.warnings} tone="caution" />
    </article>
  )
}

export function ForwardSynthesisResults({ result }: { result: ForwardSynthesisResult }) {
  const [selected, setSelected] = useState(0)
  useEffect(() => setSelected(0), [result])
  const active = result.prediction.candidates[selected]
  const diagnostics = result.prediction.diagnostics
  const assessment = result.assessment
  return (
    <section className="results-card">
      <div className="results-summary">
        <div>
          <span className="eyebrow">{assessment ? 'FORWARD ROUTE AUDIT' : 'FORWARD SYNTHESIS'}</span>
          <h2>{result.prediction.candidates.length} possible product{result.prediction.candidates.length === 1 ? '' : 's'}</h2>
        </div>
        <div className="metric-strip">
          <div><strong>{diagnostics.unique_product_count}</strong><span>products</span></div>
          <div><strong>{diagnostics.valid_pathway_count}</strong><span>pathways</span></div>
          <div><strong>{diagnostics.self_reaction_pathway_count}</strong><span>self pathways</span></div>
          <div><strong>{diagnostics.condition_profile_conflict_count}</strong><span>conditions excluded</span></div>
          <div><strong>{diagnostics.applied_operator_count}</strong><span>operators applied</span></div>
        </div>
      </div>
      {assessment && (
        <div className={`alert ${assessment.disposition === 'clear' ? '' : 'caution'}`}>
          <strong>{displayName(assessment.disposition)} route-step assessment.</strong>{' '}
          Intended product: {displayName(assessment.intended_match)}
          {assessment.intended_product_rank != null && ` at blind rank ${assessment.intended_product_rank}`}.
          {' '}Targeted replay: {displayName(assessment.targeted_replay_status)}.
          {assessment.best_competitor_product && <> Best competitor: <code>{assessment.best_competitor_product}</code>.</>}
          {assessment.score_margin != null && <> Score margin: {assessment.score_margin.toFixed(3)}.</>}
        </div>
      )}
      {!result.valid && <div className="alert error">{displayName(result.prediction.error ?? result.prediction.status)}</div>}
      <MessageList title="Assessment warnings" values={assessment?.warnings ?? []} tone="caution" />
      {result.prediction.candidates.length === 0 && <MessageList title="Scope and cautions" values={result.prediction.warnings} tone="caution" />}
      {result.prediction.candidates.length > 0 && (
        <div className="results-layout">
          <div className="table-scroll">
            <table>
              <thead><tr><th>Rank</th><th>Score</th><th>Product</th><th>Level</th><th>Support</th><th>Pathways</th><th>Conditions</th></tr></thead>
              <tbody>
                {result.prediction.candidates.map((candidate, index) => (
                  <tr key={candidate.pathway_id} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}>
                    <td><strong>{candidate.rank}</strong></td>
                    <td>{candidate.score.toFixed(3)}</td>
                    <td className="mono-value">{candidate.product_smiles}</td>
                    <td>{candidate.abstraction_level}</td>
                    <td>{candidate.independent_reference_support}</td>
                    <td>{1 + candidate.alternative_pathway_ids.length}</td>
                    <td>{candidate.recipe_evidence.evaluated ? (candidate.recipe_evidence.compatible ? 'Compatible' : 'Conflict') : 'Not supplied'}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          {active && <ForwardProductDetails candidate={active} />}
        </div>
      )}
      {result.prediction.competition_groups.length > 0 && (
        <details className="trace-panel">
          <summary>Competing pathway groups ({result.prediction.competition_groups.length})</summary>
          <div className="table-scroll">
            <table>
              <thead><tr><th>Level</th><th>Candidate ranks</th><th>Products</th><th>Operators</th></tr></thead>
              <tbody>{result.prediction.competition_groups.map((group) => (
                <tr key={`${group.competition_level}:${group.group_key}`}>
                  <td>{displayName(group.competition_level)}</td>
                  <td>{group.candidate_ranks.join(', ')}</td>
                  <td>{group.product_smiles.length}</td>
                  <td>{group.operator_ids.length}</td>
                </tr>
              ))}</tbody>
            </table>
          </div>
        </details>
      )}
      <MessageList title="Prediction scope" values={result.prediction.warnings} tone="caution" />
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

      <dl className="detail-list">
        <div><dt>Strategic steps</dt><dd>{route.evidence_summary.strategic_step_count}/{route.evidence_summary.strategic_complexity_assessed_step_count}</dd></div>
        <div><dt>Tactical steps</dt><dd>{route.evidence_summary.tactical_step_count}</dd></div>
        <div><dt>Mean strategic reduction</dt><dd>{route.evidence_summary.mean_strategic_complexity_reduction_score == null ? '—' : `${(100 * route.evidence_summary.mean_strategic_complexity_reduction_score).toFixed(1)}/100`}</dd></div>
        <div><dt>Target → frontier reduction</dt><dd>{route.evidence_summary.target_to_frontier_complexity_reduction_fraction == null ? '—' : `${(100 * route.evidence_summary.target_to_frontier_complexity_reduction_fraction).toFixed(1)}%`}</dd></div>
      </dl>

      <section className="route-step-list" aria-label="Retrosynthetic route steps">
        {route.steps.map((step, index) => (
          <article className="route-step-card" key={step.step_id}>
            <div className="route-step-heading">
              <div><span>Step {index + 1} · depth {step.depth}</span><strong>{displayName(step.candidate.transformation_kind ?? 'graph operator')}</strong></div>
              <small>{step.candidate.abstraction_level} · score {step.candidate.score.toFixed(3)}</small>
            </div>
            <ReactionImage smiles={step.candidate.proposed_reaction_smiles} label={`Route reaction ${index + 1}`} compact />
            <ForwardValidityAudit audit={step.forward_assessment} />
            <MessageList
              title="Strong intramolecular compatibility warning"
              values={(step.candidate.precursor_compatibility_assessments ?? []).map((assessment) => assessment.message)}
              tone="caution"
            />
            <dl className="detail-list">
              <div><dt>Precursors</dt><dd>{step.precursor_smiles.join(' · ')}</dd></div>
              <div><dt>Step cost</dt><dd>{step.step_cost.toFixed(3)} · {Object.entries(step.step_cost_components).filter(([, value]) => value > 0).map(([name, value]) => `${displayName(name)} ${value.toFixed(3)}`).join(' · ')}</dd></div>
              <div><dt>Signature sanity check</dt><dd>{displayName(step.candidate.forward_validation_status)}</dd></div>
              <div><dt>Independent support</dt><dd>{step.candidate.independent_reference_support}</dd></div>
              {step.candidate.precursor_realism_score != null && <div><dt>Precursor realism</dt><dd>{step.candidate.precursor_realism_score.toFixed(3)} · band +{step.candidate.precursor_realism_band_penalty}</dd></div>}
              <div><dt>Strategic reduction</dt><dd>{(100 * step.candidate.strategic_complexity_score).toFixed(1)}/100 · {displayName(step.candidate.strategic_class)}</dd></div>
              {(step.candidate.precursor_compatibility_band_penalty ?? 0) > 0 && <div><dt>Compatibility audit</dt><dd>{displayName(step.candidate.precursor_compatibility_disposition ?? 'warn')} · band +{step.candidate.precursor_compatibility_band_penalty}</dd></div>}
              {step.condition_evidence && <div><dt>Condition availability</dt><dd>{displayName(step.condition_evidence.status)} · {step.condition_evidence.compatible_candidate_count} compatible precedents</dd></div>}
              <div><dt>Operator</dt><dd className="mono-value">{step.candidate.operator_id}</dd></div>
            </dl>
            {(step.candidate.precursor_compatibility_assessments ?? []).length > 0 && (
              <details className="trace-panel">
                <summary>Compatibility audit trace</summary>
                {(step.candidate.precursor_compatibility_assessments ?? []).map((assessment) => (
                  <dl className="detail-list" key={assessment.assessment_id}>
                    <div><dt>Rule</dt><dd className="mono-value">{assessment.rule_id}</dd></div>
                    <div><dt>Component</dt><dd className="mono-value">{assessment.component_smiles}</dd></div>
                    <div><dt>Reactive pair</dt><dd>{assessment.left_site.chemist_label} + {assessment.right_site.chemist_label}</dd></div>
                    <div><dt>Warning</dt><dd>{displayName(assessment.warning_strength)} · same molecular component</dd></div>
                  </dl>
                ))}
              </details>
            )}
            {(step.candidate.precursor_realism_assessments ?? []).length > 0 && (
              <details className="trace-panel">
                <summary>Precursor realism trace</summary>
                {(step.candidate.precursor_realism_assessments ?? []).map((assessment) => (
                  <dl className="detail-list" key={assessment.canonical_smiles}>
                    <div><dt>Precursor</dt><dd className="mono-value">{assessment.canonical_smiles}</dd></div>
                    <div><dt>Tier / score</dt><dd>{displayName(assessment.evidence_tier)} · {assessment.score.toFixed(3)}</dd></div>
                    <div><dt>Evidence</dt><dd>{[assessment.evidence.buyable && 'buyable', assessment.evidence.in_compound_registry && 'registry', assessment.evidence.in_literature && 'literature'].filter(Boolean).join(' · ') || 'none'}</dd></div>
                  </dl>
                ))}
              </details>
            )}
            {step.condition_evidence && (
              <details className="trace-panel retrosynthesis-conditions" open>
                <summary>Recommended conditions ({step.condition_evidence.recommendations.length})</summary>
                {step.condition_evidence.recommendations.length > 0 ? (
                  <div className="retrosynthesis-condition-list">
                    {step.condition_evidence.recommendations.map((recommendation) => (
                      <article key={recommendation.recipe_id}>
                        <div className="retrosynthesis-condition-heading">
                          <div><strong>Recipe rank {recommendation.rank}</strong><span>{compactRecipeSummary(recommendation.resolved_recipe)}</span></div>
                          <small>Score {recommendation.score.toFixed(3)}</small>
                        </div>
                        <Conditions recipe={recommendation.resolved_recipe} />
                        <ConditionPrecedents precedents={recommendation.condition_precedents ?? []} />
                        {recommendation.condition_precedents == null && (
                          <>
                            <ReferenceRecords records={recommendation.precedent_references ?? []} />
                            <ExperimentalDetails details={recommendation.precedent_experimental_details ?? []} />
                          </>
                        )}
                      </article>
                    ))}
                  </div>
                ) : <p className="retrosynthesis-empty-evidence">No compatible condition recipe was found for this reaction.</p>}
                <MessageList title="Condition cautions" values={step.condition_evidence.warnings} tone="caution" />
              </details>
            )}
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
              <p>{displayName(leaf.terminal_evidence)} · {displayName(leaf.catalog_role_status)}</p>
              {leaf.terminal_reasons.length > 0 && <p>{leaf.terminal_reasons.map(displayName).join(' · ')}</p>}
              {leaf.unresolved_reason && <p>{displayName(leaf.unresolved_reason)}</p>}
              {(leaf.literature_match?.source_records ?? []).map((record, recordIndex) => (
                <div className="route-leaf-source" key={`${record.supplier_record_id ?? record.reaction_id ?? 'source'}:${recordIndex}`}>
                  {record.supplier && <strong>{record.supplier}</strong>}
                  {record.source_collection && <span>{record.source_collection}</span>}
                  {record.snapshot_date && <small>Snapshot {record.snapshot_date}</small>}
                  {record.supplier_record_id && <small>{record.supplier_record_id}</small>}
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
          <div><strong>{result.diagnostics.validation_attempts}</strong><span>validated</span></div>
          <div><strong>{result.search_elapsed_seconds.toFixed(1)}s</strong><span>search time</span></div>
        </div>
      </div>
      {!result.valid && <div className="alert error">{displayName(result.error ?? 'No solved multistep routes')}</div>}
      {availableRoutes.length === 0 && <MessageList title="Scope and cautions" values={result.warnings} tone="caution" />}
      {availableRoutes.length > 0 && (
        <div className="results-layout multistep-results-layout">
          <div className="table-scroll">
            <table>
              <thead><tr><th>Route</th><th>Status</th><th>Steps</th><th>Depth</th><th>Cost</th><th>Strategic steps</th><th>Frontier reduction</th><th>Realism</th><th>Conditions</th><th>Strong incompat.</th><th>Terminal leaves</th></tr></thead>
              <tbody>
                {availableRoutes.map((route, index) => (
                  <tr key={route.route_id} className={selected === index ? 'selected' : ''} onClick={() => setSelected(index)}>
                    <td><strong>{index + 1}</strong></td>
                    <td>{route.solved ? 'Solved' : 'Partial'}</td>
                    <td>{route.reaction_count}</td>
                    <td>{route.maximum_depth}</td>
                    <td>{route.route_cost.toFixed(3)}</td>
                    <td>{route.evidence_summary.strategic_step_count}/{route.evidence_summary.strategic_complexity_assessed_step_count}</td>
                    <td>{route.evidence_summary.target_to_frontier_complexity_reduction_fraction == null ? '—' : `${(100 * route.evidence_summary.target_to_frontier_complexity_reduction_fraction).toFixed(1)}%`}</td>
                    <td>{route.evidence_summary.weakest_precursor_realism_score == null ? '—' : route.evidence_summary.weakest_precursor_realism_score.toFixed(3)}</td>
                    <td>{route.evidence_summary.condition_assessed_step_count === 0 ? '—' : `${route.evidence_summary.condition_supported_step_count}/${route.evidence_summary.condition_assessed_step_count}`}</td>
                    <td>{route.evidence_summary.strong_compatibility_warning_step_count}</td>
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
