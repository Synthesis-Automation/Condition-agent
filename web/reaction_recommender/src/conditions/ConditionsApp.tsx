import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import type { Capabilities, CompletionChoice, CompletionProposal, RecipeComponent, ResolvedRecipe } from '../api/types'
import { ReactionEditor } from '../components/ReactionEditor'
import { ReactionImage } from '../components/ReactionImage'
import { CompletionDialog } from '../components/CompletionDialog'
import type { ConditionOption, ConditionsResult } from './types'
import './conditions.css'

const EXAMPLE = 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1'
const ROLES = [
  ['catalysts', 'Catalyst'], ['ligands', 'Ligand'], ['bases', 'Base'],
  ['condensation_agents', 'Coupling reagent'], ['oxidants', 'Oxidant'],
  ['reductants', 'Reductant'], ['acids', 'Acid'], ['additives', 'Additive'],
  ['solvents', 'Solvent'], ['other_components', 'Other material'],
] as const

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

function Procedure({ recipe }: { recipe: ResolvedRecipe }) {
  const stages = Array.isArray(recipe.stages) ? recipe.stages as Array<Record<string, unknown>> : []
  return <>
    <dl className="condition-materials">
      {groups(recipe).map(group => <div key={group.label}>
        <dt>{group.label}</dt><dd>{group.values.map((item, index) => <span key={index}>
          {materialName(item)}
          <small>{item.amount != null ? ` · ${item.amount} ${item.amount_unit || '(unit not reported)'}` : ' · amount not reported'}</small>
          {item.identity_status && item.identity_status !== 'resolved' ? <small> · identity unresolved</small> : null}
        </span>)}</dd>
      </div>)}
    </dl>
    {stages.length > 0 ? <ol className="condition-stages">{stages.map((stage, i) => <li key={i}>
      Stage {String(stage.stage_index ?? i + 1)}: {setup(stage as ResolvedRecipe).join(' · ')}
    </li>)}</ol> : <div className="condition-setpoints">{setup(recipe).map(text => <span key={text}>{text}</span>)}</div>}
  </>
}

function Evidence({ option }: { option: ConditionOption }) {
  return <details className="condition-evidence">
    <summary>Precedents & details</summary>
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
        {(item.compatibility_evidence || []).length > 0 && <details><summary>Compatibility checks</summary><ul>{item.compatibility_evidence?.map((line, i) => <li key={i}>{line}</li>)}</ul></details>}
        {evidence.warnings.length > 0 && <details><summary>Evidence qualifications ({evidence.warnings.length})</summary><ul>{evidence.warnings.map((line, i) => <li key={i}>{readable(line)}</li>)}</ul></details>}
      </section>
    })}
  </details>
}

function ConditionCard({ option, onDownload, selected, onSelect }: {
  option: ConditionOption; onDownload: () => void; selected: boolean; onSelect: () => void
}) {
  const [copyStatus, setCopyStatus] = useState('Copy conditions')
  const recipe = option.resolved_recipe
  const primary = option.evidence[0].recommendation
  const specificCautions = option.cautions.filter(caution =>
    !['Source reactions are not structure-verified', 'Reaction-type labels are weak evidence'].includes(caution))
  const why = option.evidence_kind === 'weak_label'
    ? 'Compatible reaction handles in label-based records; use as a screening starting point.'
    : option.evidence_kind === 'structure_review'
      ? 'Retrieved through broader structural evidence; inspect the reaction match before use.'
      : 'Matched reaction changes and compatible conditions from structural precedents.'
  const copy = async () => {
    const lines = groups(recipe).map(group => `${group.label}: ${group.values.map(item => `${materialName(item)}${item.amount != null ? ` (${item.amount} ${item.amount_unit || 'unit not reported'})` : ' (amount not reported)'}`).join(', ')}`)
    const stages = Array.isArray(recipe.stages) ? recipe.stages as ResolvedRecipe[] : []
    lines.push(...(stages.length ? stages.map((stage, i) => `Stage ${i + 1}: ${setup(stage).join(' · ')}`) : [setup(recipe).join(' · ')]))
    lines.push(option.evidence_label, ...option.cautions, 'Protocol requires preparation and review before execution.')
    try { await navigator.clipboard.writeText(lines.join('\n')); setCopyStatus('Copied') }
    catch { setCopyStatus('Copy unavailable — download JSON') }
  }
  return <article className={`condition-option ${option.evidence_kind}`}>
    <div className="condition-option-heading">
      <div><span className="condition-rank">{String(option.rank).padStart(2, '0')}</span><h3>Condition option {option.rank}</h3></div>
      <label className="condition-select"><input type="checkbox" checked={selected} onChange={onSelect} />Select for screen</label>
    </div>
    <div className="condition-badges"><span className={`condition-badge ${option.evidence_kind}`}>{option.evidence_label}</span>
      <span>{option.evidence[0].source === 'generic' ? `${primary.reference_support ?? 0} reference(s)` : `${primary.support ?? 0} label observation(s)`}</span>
      {new Set(option.evidence.map(item => item.source)).size > 1 && <span>Supported by both sources</span>}
    </div>
    <Procedure recipe={recipe} />
    <p className="condition-why">{why}</p>
    {specificCautions.length > 0 && <div className="condition-caution"><strong>Check before use</strong><ul>{specificCautions.map((caution, i) => <li key={i}>{caution}</li>)}</ul></div>}
    <Evidence option={option} />
    <details className="condition-automation">
      <summary>Automation preparation · {option.synthesis_protocol.missing_required_fields.length} missing fields</summary>
      <p>JSON includes material identities, reported amounts and units, ordered source stages, and provenance. A robot adapter must resolve the gaps and validate the setup before execution.</p>
      <ul>{option.synthesis_protocol.missing_required_fields.map(field => <li key={field}>{readable(field)}</li>)}</ul>
      <p>No dispensing order, scale, or workup is invented.</p>
    </details>
    <div className="condition-card-actions">
      <button className="button secondary" onClick={onDownload}>Download automation JSON</button>
      <button className="button quiet" onClick={() => void copy()}>{copyStatus}</button>
      <span>Protocol draft · review required</span>
    </div>
  </article>
}

export default function ConditionsApp() {
  const [reaction, setReaction] = useState('')
  const [capabilities, setCapabilities] = useState<Capabilities | null>(null)
  const [result, setResult] = useState<ConditionsResult | null>(null)
  const [error, setError] = useState('')
  const [busy, setBusy] = useState(false)
  const [status, setStatus] = useState('')
  const [proposal, setProposal] = useState<CompletionProposal | null>(null)
  const [visible, setVisible] = useState(3)
  const [selected, setSelected] = useState<string[]>([])
  const requestId = useRef(0)
  const resultsRef = useRef<HTMLElement>(null)

  useEffect(() => {
    api.capabilities().then(setCapabilities).catch(() => setError('The recommendation service is unavailable. Check that the Python API is running, then reload.'))
  }, [])

  const changeReaction = (value: string) => {
    requestId.current += 1
    setReaction(value); setResult(null); setSelected([]); setProposal(null); setStatus(''); setBusy(false); setError('')
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
        body: JSON.stringify({ reaction_smiles: query, completion_choices: choices || [] }),
      })
      const body = await response.json()
      if (!response.ok) throw new Error(typeof body.detail?.message === 'string' ? body.detail.message : 'The search could not complete. Check the reaction and try again.')
      if (currentId !== requestId.current) return
      const next = body.data as ConditionsResult
      setResult(next); setVisible(next.shortlist_size); setStatus('')
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
    <header className="conditions-header"><a href="/" className="conditions-brand"><span className="conditions-logo">C</span>Condition Desk</a><span className="conditions-service"><i className={available ? 'available' : ''} />{capabilities ? result?.sources.some(source => source.status === 'unavailable') ? 'Some libraries unavailable' : available ? 'Service connected' : 'Libraries unavailable' : 'Connecting…'}</span></header>
    <main className="conditions-main">
      <section className="conditions-intro"><span className="conditions-kicker">FROM REACTION TO EXPERIMENT</span><h1>Find your starting conditions.</h1><p>Draw your reaction. Compare precedent conditions and screening suggestions in one place.</p></section>
      <section className="conditions-query" aria-label="Reaction search">
        <ReactionEditor value={reaction} onChange={changeReaction} onError={setError} />
        <div className="conditions-search-bar"><button className="button quiet" onClick={() => changeReaction(EXAMPLE)}>Try a Suzuki example</button><span>Include the intended product.</span><button className="button primary" disabled={busy || !reaction.trim() || !available} onClick={() => void search()}>{busy ? 'Finding conditions…' : 'Find conditions →'}</button></div>
      </section>
      {error && <div className="condition-error" role="alert">{error}</div>}
      <div className="condition-progress" role="status" aria-live="polite">{busy && <span className="condition-spinner" />}{status}</div>
      {!result && !busy && <div className="conditions-introduction"><div><b>01</b><h3>Start with precedent</h3><p>Conditions matched to reaction changes and compatible functional groups.</p></div><div><b>02</b><h3>Explore a screen</h3><p>Additional label-based suggestions, with their evidence limits visible.</p></div><div><b>03</b><h3>Prepare the experiment</h3><p>Copy conditions or export a structured protocol for automation planning.</p></div></div>}
      {result && <section ref={resultsRef} className="conditions-results" aria-label="Recommended conditions">
        <div className="conditions-result-heading"><div><span className="conditions-kicker">YOUR STARTING POINTS</span><h2>{result.valid ? 'Conditions to consider' : 'No supported conditions found'}</h2><p>{result.valid ? 'Reaction precedents first; additional screening suggestions follow. Expand a card to inspect its evidence.' : 'Check the reaction drawing or try a supported transformation. No conditions have been invented.'}</p></div><button className="button secondary" onClick={() => download('condition-recommendations.json', result)}>Download full results JSON</button></div>
        <div className="condition-source-status">{result.sources.map(source => <span key={source.source} className={source.status === 'ok' ? 'ok' : 'limited'}>{source.source === 'generic' ? 'Reaction library' : 'Screening library'} · {source.status === 'ok' ? 'matches found' : source.status === 'abstained' ? 'no supported matches' : source.status}</span>)}</div>
        {result.sources.filter(source => source.status !== 'ok').map(source => <p className="condition-note" key={source.source}>{source.source === 'generic' ? 'Reaction library' : 'Screening library'}: {source.message}</p>)}
        {result.recommendations.some(option => option.evidence_kind === 'weak_label') && <p className="condition-note">Screening suggestions use source labels; their precedent reaction structures are not verified.</p>}
        {result.recommendations.slice(0, visible).map(option => <ConditionCard key={option.option_id} option={option} selected={selected.includes(option.option_id)} onSelect={() => setSelected(current => current.includes(option.option_id) ? current.filter(id => id !== option.option_id) : [...current, option.option_id])} onDownload={() => download(`condition-option-${option.rank}-automation.json`, result.automation_exports[option.option_id])} />)}
        {visible < result.recommendations.length && <button className="button secondary conditions-more" onClick={() => setVisible(count => count + 3)}>Show 3 more options ({result.recommendations.length - visible} remaining)</button>}
        {result.recommendations.length > 0 && <div className="condition-screen-bar"><div><strong>{selected.length ? `${selected.length} condition${selected.length === 1 ? '' : 's'} selected` : 'Build your screening set'}</strong><p>Select intact recipes above. Each export retains its own evidence and preparation gaps.</p></div><button className="button primary" disabled={!selected.length} onClick={exportScreen}>Download screening JSON</button></div>}
      </section>}
      <footer className="conditions-footer">Condition Desk <span>Structure-based recommendations. Explicit uncertainty. Traceable recipes.</span></footer>
    </main>
    {proposal && <CompletionDialog proposal={proposal} onCancel={() => { setProposal(null); setStatus('') }} onConfirm={choices => void search(choices)} />}
  </div>
}
