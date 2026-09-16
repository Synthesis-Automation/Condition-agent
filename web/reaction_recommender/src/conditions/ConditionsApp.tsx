import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import type { Capabilities, CompletionChoice, CompletionProposal, RecipeComponent, ResolvedRecipe } from '../api/types'
import { ReactionEditor } from '../components/ReactionEditor'
import { ReactionImage } from '../components/ReactionImage'
import { CompletionDialog } from '../components/CompletionDialog'
import type { ConditionOption, ConditionSearchScope, ConditionsResult } from './types'
import examples from './examples.json'
import './conditions.css'

const ROLES = [
  ['catalysts', 'Catalyst'], ['ligands', 'Ligand'], ['bases', 'Base'],
  ['condensation_agents', 'Coupling reagent'], ['oxidants', 'Oxidant'],
  ['reductants', 'Reductant'], ['acids', 'Acid'], ['additives', 'Additive'],
  ['solvents', 'Solvent'], ['other_components', 'Other material'],
] as const

const SOURCES = ['generic', 'weak_label'] as const
type ConditionSource = typeof SOURCES[number]
const SOURCE_LABELS = { generic: 'Literature-based', weak_label: 'Screening suggestions' }

function materialName(item: RecipeComponent): string {
  return String(item.canonical_name || item.display_name || item.name || item.raw_identifier || item.cas || item.substance_id || 'Identity not resolved')
}

function groups(recipe: ResolvedRecipe): Array<{ label: string; values: RecipeComponent[] }> {
  return ROLES.flatMap(([key, label]) => {
    const values = recipe[key]
    return Array.isArray(values) && values.length ? [{ label, values }] : []
  })
}

function setup(recipe: ResolvedRecipe): string[] {
  return [
    recipe.temperature_c != null ? `${recipe.temperature_c} °C` : 'Temperature not reported',
    recipe.time_h != null ? `${recipe.time_h} h` : 'Time not reported',
    recipe.concentration_m != null ? `${recipe.concentration_m} M` : '',
    recipe.pressure_bar != null ? `${recipe.pressure_bar} bar` : '',
    recipe.atmosphere || '',
  ].filter(Boolean)
}

function download(filename: string, value: unknown): void {
  const url = URL.createObjectURL(new Blob([JSON.stringify(value, null, 2) + '\n'], { type: 'application/json' }))
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  link.click()
  window.setTimeout(() => URL.revokeObjectURL(url), 1000)
}

function readable(value: string): string {
  return value.replaceAll('_', ' ').replaceAll('.', ' › ').toLowerCase()
}

function ChemistrySummary({ result }: { result: ConditionsResult }) {
  const structural = result.sources.find(source => source.source === 'generic')?.result
  const screening = result.sources.find(source => source.source === 'weak_label')?.result
  const label = structural?.reaction_label?.text
  const participants = screening?.query_participants || []
  const nearby = [...new Set((structural?.reaction_partners || []).flatMap(partner =>
    (partner.nearby_groups || []).map(group => group.label)))]
  if (!label && !participants.length && !nearby.length) return null
  return <section className="condition-chemistry" aria-label="Reaction chemistry">
    <h3>{label || 'Reactive groups in your drawing'}</h3>
    {participants.length > 0 && <div className="condition-chemistry-tags" aria-label="Reactive groups">{participants.map((partner, index) => <span key={`${partner.site_id}-${index}`}>{partner.chemist_label}</span>)}</div>}
    {nearby.length > 0 && <details className="condition-nearby"><summary>Nearby groups</summary><p>{nearby.join(' · ')}</p></details>}
  </section>
}

function Procedure({ recipe }: { recipe: ResolvedRecipe }) {
  const stages = Array.isArray(recipe.stages) ? recipe.stages as Array<Record<string, unknown>> : []
  return <>
    <dl className="condition-materials">
      {groups(recipe).map(group => <div key={group.label}>
        <dt>{group.label}</dt><dd>{group.values.map((item, index) => <span key={index}>
          {materialName(item)}
          {item.amount != null && <small> · {String(item.amount)} {String(item.amount_unit || '(unit not reported)')}</small>}
          {item.identity_status && item.identity_status !== 'resolved' ? <small> · identity unresolved</small> : null}
        </span>)}</dd>
      </div>)}
    </dl>
    {stages.length > 1 ? <ol className="condition-stages">{stages.map((stage, i) => <li key={i}>
      Stage {String(stage.stage_index ?? i + 1)}: {setup(stage as ResolvedRecipe).join(' · ')}
    </li>)}</ol> : <div className="condition-setpoints">{setup(stages.length ? stages[0] as ResolvedRecipe : recipe).map(text => <span key={text}>{text}</span>)}</div>}
  </>
}

function Evidence({ option }: { option: ConditionOption }) {
  return <details className="condition-evidence" open>
    <summary>Details & references</summary>
    {option.cautions.length > 0 && <section><h4>Review notes</h4><ul>{option.cautions.map((caution, i) => <li key={i}>{caution.includes(' ') ? caution : readable(caution)}</li>)}</ul></section>}
    {option.evidence.map((evidence, index) => {
      const item = evidence.recommendation
      const refs = item.precedent_references || []
      const precedents = item.condition_precedents || []
      return <section key={index}>
        <h4>{evidence.source === 'generic' ? 'Reaction precedent evidence' : 'Label-based screening evidence'}</h4>
        {item.match_level && <p>{item.match_namespace === 'shared_reaction_core.v2' ? item.match_label : `Level ${item.match_level} · ${item.match_label}`} — evidence distance, not a success probability.</p>}
        {item.evidence_relation === 'analogue_evidence' && <p>Related precedent: inspect its reported source and substrate differences.</p>}
        {item.match_details?.map((detail) => <p key={detail}>{detail}</p>)}
        <p>{evidence.source === 'generic'
          ? `${item.reference_support ?? 0} independent reference(s) · ${item.support ?? 0} observation(s)`
          : `${item.support ?? 0} source observation(s). Precedent reaction structures are not verified.`}</p>
        {(item.explanation || []).map((line, i) => <p key={i}>{line}</p>)}
        {item.historical_yield_pct != null && <p className="condition-note">Historical yield summary: {item.historical_yield_pct}%. This aggregates matched source observations; it is not a yield prediction for this reaction.</p>}
        {precedents.slice(0, 3).map((precedent, i) => <div className="condition-precedent" key={i}>
          {precedent.reaction_smiles && <ReactionImage smiles={precedent.reaction_smiles} label={`Precedent reaction ${i + 1}`} compact />}
          <p>{precedent.reference_record?.raw_reference || precedent.reference_record?.normalized_citation || precedent.reference_id || 'Citation not available'}</p>
          {precedent.experimental_detail?.procedure_text && <details><summary>Reported procedure</summary><p>{precedent.experimental_detail.procedure_text}</p></details>}
        </div>)}
        {!precedents.length && refs.map((ref, i) => <p key={i}>{ref.raw_reference || ref.normalized_citation || ref.reference_id}</p>)}
        {evidence.source === 'weak_label' && <p>Source rows: {(item.source_row_numbers || []).join(', ') || 'See full result JSON'}</p>}
        {evidence.source === 'weak_label' && item.source_matches?.slice(0, 3).map((match, i) => <p className="condition-label-match" key={i}><strong>{match.participant_display_labels.join(' + ')}</strong><br />{match.source_reaction_type} · source row {match.source_row_number}</p>)}
        {(item.compatibility_evidence || []).length > 0 && <details><summary>Compatibility checks</summary><ul>{item.compatibility_evidence?.map((line, i) => <li key={i}>{line}</li>)}</ul></details>}
        {evidence.warnings.length > 0 && <details><summary>Evidence qualifications ({evidence.warnings.length})</summary><ul>{evidence.warnings.map((line, i) => <li key={i}>{readable(line)}</li>)}</ul></details>}
      </section>
    })}
    <section>
      <h4>Experiment preparation</h4>
      <p>Confirm unreported quantities and procedure details before use.</p>
      <ul>{option.synthesis_protocol.missing_required_fields.map(field => <li key={field}>{readable(field)}</li>)}</ul>
    </section>
  </details>
}

function ConditionDetails({ option, source, number, onDownload }: {
  option: ConditionOption; source: ConditionSource; number: number
  onDownload: () => void
}) {
  const [copyStatus, setCopyStatus] = useState('Copy conditions')
  const recipe = option.resolved_recipe
  const primary = option.evidence.find(item => item.source === source)!.recommendation
  const screening = source === 'weak_label'
  const copy = async () => {
    const lines = groups(recipe).map(group => `${group.label}: ${group.values.map(item => `${materialName(item)}${item.amount != null ? ` (${item.amount} ${item.amount_unit || 'unit not reported'})` : ' (amount not reported)'}`).join(', ')}`)
    const stages = Array.isArray(recipe.stages) ? recipe.stages as ResolvedRecipe[] : []
    lines.push(...(stages.length ? stages.map((stage, i) => `Stage ${i + 1}: ${setup(stage).join(' · ')}`) : [setup(recipe).join(' · ')]))
    lines.push(screening ? 'Weak-label screening suggestion; source reactions are not structure-verified.' : option.evidence_label, ...option.cautions, 'Protocol requires preparation and review before execution.')
    try { await navigator.clipboard.writeText(lines.join('\n')); setCopyStatus('Copied') }
    catch { setCopyStatus('Copy unavailable — download JSON') }
  }
  return <article className={`condition-option ${screening ? 'weak_label' : option.evidence_kind}`}>
    <div className="condition-option-heading">
      <h3>{screening ? 'Screen' : 'Precedent'} {number}</h3>
    </div>
    <div className="condition-badges"><span className={`condition-badge ${screening ? 'weak_label' : option.evidence_kind}`}>{screening ? 'Weak-label screening' : primary.match_label || option.evidence_label}</span>
      <span>{!screening ? `${primary.reference_support ?? 0} reference(s)` : `${primary.support ?? 0} source observation(s)`}</span>
      {primary.historical_yield_pct != null && <span>Reported mean yield {primary.historical_yield_pct}%</span>}
      {new Set(option.evidence.map(item => item.source)).size > 1 && <span>Supported by both sources</span>}
    </div>
    <Procedure recipe={recipe} />
    <div className="condition-card-actions">
      <button className="button secondary" onClick={() => void copy()}>{copyStatus}</button>
      <button className="button quiet" onClick={onDownload}>Export recipe</button>
    </div>
    <Evidence option={option} />
  </article>
}

function ConditionResults({ result, selected, onSelect, onDownload }: {
  result: ConditionsResult; selected: string[]
  onSelect: (optionId: string) => void; onDownload: (option: ConditionOption) => void
}) {
  const optionsFor = (source: ConditionSource) => result.recommendations
    .filter(option => option.evidence.some(item => item.source === source))
    .sort((left, right) => (left.evidence.find(item => item.source === source)?.recommendation.rank ?? left.rank)
      - (right.evidence.find(item => item.source === source)?.recommendation.rank ?? right.rank))
  const [activeSource, setActiveSource] = useState<ConditionSource>(() => optionsFor('generic').length ? 'generic' : optionsFor('weak_label').length ? 'weak_label' : 'generic')
  const [viewed, setViewed] = useState<Partial<Record<ConditionSource, string>>>({})
  const detailRef = useRef<HTMLElement>(null)

  const viewOption = (source: ConditionSource, optionId: string) => {
    setViewed(current => ({ ...current, [source]: optionId }))
    if (window.matchMedia('(max-width: 900px)').matches) {
      window.requestAnimationFrame(() => detailRef.current?.scrollIntoView({ block: 'start', behavior: 'smooth' }))
    }
  }

  return <>
    <div className="condition-tabs" role="tablist" aria-label="Condition sources">
      {SOURCES.map(source => <button key={source} type="button" role="tab"
        id={`condition-tab-${source}`} aria-controls={`condition-panel-${source}`}
        aria-selected={activeSource === source} tabIndex={activeSource === source ? 0 : -1}
        onClick={() => setActiveSource(source)} onKeyDown={event => {
          if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return
          event.preventDefault()
          const next = event.key === 'Home' ? SOURCES[0] : event.key === 'End' ? SOURCES[1] : SOURCES.find(item => item !== source)!
          setActiveSource(next)
          document.getElementById(`condition-tab-${next}`)?.focus()
        }}>{SOURCE_LABELS[source]} <span>{optionsFor(source).length}</span></button>)}
    </div>
    {result.sources.filter(source => source.status !== 'ok').map(source => <p className="condition-note" key={source.source}>{SOURCE_LABELS[source.source]}: {source.message}</p>)}
    {SOURCES.map(source => {
      const options = optionsFor(source)
      const current = options.find(option => option.option_id === viewed[source]) || options[0]
      const screening = source === 'weak_label'
      const noun = screening ? 'screen' : 'precedent'
      return <section key={source} role="tabpanel" id={`condition-panel-${source}`}
        aria-labelledby={`condition-tab-${source}`} hidden={activeSource !== source} tabIndex={0}>
        {activeSource === source && <>
          <p className="conditions-tab-description">{screening ? 'Weak-label suggestions. Source reaction structures are not verified.' : 'Literature precedents, closest structural matches first.'}</p>
          {!current ? <div className="conditions-empty">{screening ? 'No supported screening suggestions for this reaction.' : 'No qualifying literature precedents for this search.'}</div> : <div className="conditions-browser">
            <div className="condition-table-panel">
              <p className="condition-table-hint">Select a row for details. Check recipes to build a screening set.</p>
              <div className="condition-table-scroll" role="region" aria-label={`${SOURCE_LABELS[source]} table`} tabIndex={0}>
                <table className="condition-table" aria-label={`${SOURCE_LABELS[source]} suggestions`}>
                  <colgroup><col className="condition-check-col" /><col className="condition-number-col" /><col /><col className="condition-solvent-col" /><col className="condition-yield-col" /></colgroup>
                  <thead><tr><th scope="col"><span className="condition-sr-only">Add to screening set</span></th><th scope="col">No.</th><th scope="col">Catalyst / reagents</th><th scope="col">Solvent</th><th scope="col">Reported yield</th></tr></thead>
                  <tbody>{options.map((option, index) => {
                    const materials = groups(option.resolved_recipe)
                    const reagents = materials.filter(group => group.label !== 'Solvent').map(group => group.values.map(materialName).join(', ')).join(' · ') || 'Not reported'
                    const solvents = materials.find(group => group.label === 'Solvent')?.values.map(materialName).join(', ') || 'Not reported'
                    const primary = option.evidence.find(item => item.source === source)!.recommendation
                    const active = option.option_id === current.option_id
                    return <tr key={option.option_id} data-active={active} onClick={() => viewOption(source, option.option_id)}>
                      <td onClick={event => event.stopPropagation()}><input type="checkbox" aria-label={`Select ${noun} ${index + 1}`} checked={selected.includes(option.option_id)} onChange={() => onSelect(option.option_id)} /></td>
                      <td><button type="button" className="condition-row-link" aria-label={`View ${noun} ${index + 1}`} aria-current={active ? 'true' : undefined} aria-controls={`condition-details-${source}`}>{index + 1}</button></td>
                      <td><span className="condition-table-material" title={reagents}>{reagents}</span></td>
                      <td><span className="condition-table-material" title={solvents}>{solvents}</span></td>
                      <td>{primary.historical_yield_pct != null ? `${primary.historical_yield_pct}%` : <span aria-label="Not reported">—</span>}</td>
                    </tr>
                  })}</tbody>
                </table>
              </div>
            </div>
            <section key={`${source}-${current.option_id}`} ref={detailRef} className="condition-detail-panel" id={`condition-details-${source}`} aria-label="Selected suggestion details" tabIndex={0}>
              <ConditionDetails option={current} source={source} number={options.indexOf(current) + 1} onDownload={() => onDownload(current)} />
            </section>
          </div>}
        </>}
      </section>
    })}
  </>
}

export default function ConditionsApp() {
  const [reaction, setReaction] = useState('')
  const [exampleLabel, setExampleLabel] = useState('')
  const [capabilities, setCapabilities] = useState<Capabilities | null>(null)
  const [result, setResult] = useState<ConditionsResult | null>(null)
  const [error, setError] = useState('')
  const [busy, setBusy] = useState(false)
  const [status, setStatus] = useState('')
  const [proposal, setProposal] = useState<CompletionProposal | null>(null)
  const [topK, setTopK] = useState(10)
  const [searchScope, setSearchScope] = useState<ConditionSearchScope>('automatic')
  const [selected, setSelected] = useState<string[]>([])
  const requestId = useRef(0)
  const resultsRef = useRef<HTMLElement>(null)

  useEffect(() => {
    api.capabilities().then(setCapabilities).catch(() => setError('The recommendation service is unavailable. Check that the Python API is running, then reload.'))
  }, [])

  const changeReaction = (value: string) => {
    requestId.current += 1
    setReaction(value); setExampleLabel(''); setResult(null); setSelected([]); setProposal(null); setStatus(''); setBusy(false); setError('')
  }

  const loadExample = () => {
    const candidates = examples.filter(example => example.reaction_smiles !== reaction.trim())
    const example = candidates[Math.floor(Math.random() * candidates.length)]
    changeReaction(example.reaction_smiles)
    setExampleLabel(example.label)
  }

  const search = async (choices?: CompletionChoice[]) => {
    const currentId = ++requestId.current
    const query = reaction.trim()
    setBusy(true); setError(''); setResult(null); setSelected([]); setProposal(null)
    setStatus('Checking the reaction and searching both condition libraries…')
    try {
      if (!choices) {
        const prepared = await api.prepareReaction(query)
        if (currentId !== requestId.current) return
        if (!prepared.valid) throw new Error('The reaction could not be validated. Check the starting materials and product.')
        if (prepared.completion_proposal.requirements.length) {
          setProposal(prepared.completion_proposal); setStatus('Confirm the source of the missing fragment.'); return
        }
      }
      const response = await fetch('/api/v1/conditions/recommend', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ reaction_smiles: query, completion_choices: choices || [], top_k: topK, search_scope: searchScope }),
      })
      const body = await response.json()
      if (!response.ok) throw new Error(typeof body.detail?.message === 'string' ? body.detail.message : 'The search could not complete. Check the reaction and try again.')
      if (currentId !== requestId.current) return
      const next = body.data as ConditionsResult
      setResult(next); setStatus('')
      window.setTimeout(() => resultsRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' }), 50)
    } catch (caught) {
      if (currentId === requestId.current) { setError(caught instanceof Error ? caught.message : 'Search failed.'); setStatus('') }
    } finally {
      if (currentId === requestId.current) setBusy(false)
    }
  }

  const exportScreen = () => {
    if (!result || !selected.length) return
    download('condition-screening-selection.json', {
      artifact_type: 'condition_screening_selection', schema_version: '1.0', execution_ready: false,
      query_reaction_smiles: result.query_reaction_smiles,
      selection_method: 'chemist_selected_intact_recipes',
      handoffs: result.recommendations.filter(option => selected.includes(option.option_id)).map(option => result.automation_exports[option.option_id]),
    })
  }

  const available = capabilities?.recommendation || capabilities?.weak_label_recommendation
  return <div className="conditions-app">
    <header className="conditions-header"><a href="/" className="conditions-brand"><span className="conditions-logo" aria-hidden="true">ZBS</span><h1>ZBS chemistry recommender</h1></a><span className="conditions-service"><i className={available ? 'available' : ''} />{capabilities ? result?.sources.some(source => source.status === 'unavailable') ? 'Some libraries unavailable' : available ? 'Connected' : 'Libraries unavailable' : 'Connecting…'}</span></header>
    <main className="conditions-main">
      <section className="conditions-query" aria-label="Reaction search">
        <ReactionEditor value={reaction} onChange={changeReaction} onError={setError} />
        <div className="conditions-search-bar">
          <label>Max recipes per tab<select aria-label="Max recipes per tab" value={topK} disabled={busy} onChange={event => setTopK(Number(event.target.value))}><option value={5}>5</option><option value={10}>10</option><option value={20}>20</option><option value={50}>50</option></select></label>
          <label>Reaction matches<select aria-label="Reaction matches" value={searchScope} disabled={busy} onChange={event => setSearchScope(event.target.value as ConditionSearchScope)}><option value="automatic">Automatic broadening</option><option value="same_handle">Same reactive handle</option><option value="broad">Broader analogues</option></select></label>
          <div className="conditions-example"><button className="button quiet" onClick={loadExample}>Try an example</button>{exampleLabel && <span role="status">{exampleLabel}</span>}</div>
          <button className="button primary" disabled={busy || !reaction.trim() || !available} onClick={() => void search()}>{busy ? 'Finding conditions…' : 'Find conditions'}</button>
        </div>
      </section>
      {error && <div className="condition-error" role="alert">{error}</div>}
      {(busy || status) && <div className="condition-progress" role="status" aria-live="polite">{busy && <span className="condition-spinner" />}{status}</div>}
      {result && <section ref={resultsRef} className="conditions-results" aria-label="Recommended conditions">
        <div className="conditions-result-heading"><h2>{result.valid ? 'Recommended conditions' : 'No supported conditions found'}</h2><button className="button quiet" onClick={() => download('condition-recommendations.json', result)}>Export all results</button></div>
        <ChemistrySummary result={result} />
        <ConditionResults result={result} selected={selected}
          onSelect={optionId => setSelected(current => current.includes(optionId) ? current.filter(id => id !== optionId) : [...current, optionId])}
          onDownload={option => download(`condition-option-${option.rank}-automation.json`, result.automation_exports[option.option_id])} />
        {selected.length > 0 && <div className="condition-screen-bar has-selection"><strong>{selected.length} condition{selected.length === 1 ? '' : 's'} selected</strong><button className="button primary" onClick={exportScreen}>Export screening set</button></div>}
      </section>}
    </main>
    {proposal && <CompletionDialog proposal={proposal} onCancel={() => { setProposal(null); setStatus('') }} onConfirm={choices => void search(choices)} />}
  </div>
}
