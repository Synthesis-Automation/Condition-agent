import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import type { Capabilities, RetrosynthesisCandidate } from '../api/types'
import type { PlanningChoice, PlanningNode, PlanningRequest, PlanningResponse, PlanningSearch, PlanningSession, PlanningSettings } from '../api/planning'
import { ReactionEditor } from './ReactionEditor'
import { ReactionImage } from './ReactionImage'
import { compactRecipeSummary, displayName } from './Results'
import './interactive-planner.css'

const STORAGE_KEY = 'reaction-workbench.interactive-planning.v1'
const DEFAULT_SETTINGS: PlanningSettings = {
  library_mode: 'full', top_k: 5, include_l0: true, use_context: true,
  diversify: true, use_precursor_realism: true, use_forward_validation: true,
}

function nodes(root: PlanningNode): PlanningNode[] {
  return [root, ...root.children.flatMap(nodes)]
}

function choiceKey(choice: PlanningChoice): string {
  return `${choice.search_id}:${choice.strategy_index}:${choice.realization_index}`
}

function selectedCandidate(session: PlanningSession, choice: PlanningChoice | null): RetrosynthesisCandidate | undefined {
  if (!choice) return undefined
  const strategy = session.searches.find(search => search.search_id === choice.search_id)?.result.strategies[choice.strategy_index]
  return strategy && [strategy.representative, ...strategy.alternate_realizations][choice.realization_index]
}

function download(value: unknown, filename: string) {
  const url = URL.createObjectURL(new Blob([JSON.stringify(value, null, 2)], { type: 'application/json' }))
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  link.click()
  setTimeout(() => URL.revokeObjectURL(url), 1000)
}

function savedSettings(value: unknown): PlanningSettings {
  if (!value || typeof value !== 'object') return DEFAULT_SETTINGS
  const result = { ...DEFAULT_SETTINGS }
  for (const key of Object.keys(result) as Array<keyof PlanningSettings>) {
    const item = (value as Record<string, unknown>)[key]
    if (key === 'library_mode' && (item === 'full' || item === 'compact')) result[key] = item
    else if (key === 'top_k' && typeof item === 'number' && Number.isInteger(item) && item >= 1 && item <= 50) result[key] = item
    else if (key !== 'library_mode' && key !== 'top_k' && typeof item === 'boolean') result[key] = item
  }
  return result
}

export function useInteractivePlanner(active: boolean) {
  const [data, setData] = useState<PlanningResponse | null>(null)
  const [settings, setSettings] = useState<PlanningSettings>(DEFAULT_SETTINGS)
  const [target, setTarget] = useState('')
  const [editingTarget, setEditingTarget] = useState(true)
  const [selectedId, setSelectedId] = useState('root')
  const [busy, setBusy] = useState('')
  const [error, setError] = useState('')
  const [notice, setNotice] = useState('')
  const [saveError, setSaveError] = useState('')
  const pending = useRef<AbortController | null>(null)
  const started = useRef(false)
  const session = data?.session
  const selected = session ? nodes(session.root).find(node => node.node_id === selectedId) ?? session.root : null

  const cancel = () => { pending.current?.abort(); pending.current = null; setBusy('') }
  const run = async (request: PlanningRequest, message: string): Promise<PlanningResponse | null> => {
    cancel()
    const controller = new AbortController()
    pending.current = controller
    setBusy(message)
    setError('')
    try {
      const next = await api.planningAction(request, controller.signal)
      if (controller.signal.aborted) return null
      setData(next)
      setEditingTarget(false)
      setSelectedId(previous => nodes(next.session.root).some(node => node.node_id === previous) ? previous : 'root')
      return next
    } catch (reason) {
      if (!controller.signal.aborted) setError(reason instanceof Error ? reason.message : 'Planning request failed')
      return null
    } finally {
      if (pending.current === controller) { pending.current = null; setBusy('') }
    }
  }

  const restore = async (raw: unknown) => {
    if (!raw || typeof raw !== 'object') { setError('Choose a planning session JSON file.'); return }
    const saved = raw as { schema_version?: string; session?: PlanningSession; settings?: unknown; selected_id?: string }
    if (saved.schema_version !== 'interactive_planning_browser.v1' || !saved.session) {
      setError('Unsupported session file. Import an exported planning session, not a route-only export.')
      return
    }
    const next = await run({ action: 'restore', session: saved.session }, 'Checking saved route…')
    if (next) {
      setSettings(savedSettings(saved.settings))
      setSelectedId(nodes(next.session.root).some(node => node.node_id === saved.selected_id) ? saved.selected_id! : 'root')
      setNotice('Session restored. Structures and route connections were checked; saved search, condition, and stock evidence has not been rerun.')
    }
  }

  useEffect(() => {
    if (!active) { cancel(); if (!data) started.current = false; return }
    if (started.current) return
    started.current = true
    try {
      const saved = localStorage.getItem(STORAGE_KEY)
      if (saved) void restore(JSON.parse(saved))
    } catch { setError('The saved session could not be read. You can import a session or start a new target.') }
  }, [active]) // The planner owns its session across mode changes.
  useEffect(() => () => pending.current?.abort(), [])
  useEffect(() => {
    if (!data) return
    try {
      localStorage.setItem(STORAGE_KEY, JSON.stringify({
        schema_version: 'interactive_planning_browser.v1', session: data.session, settings, selected_id: selectedId,
      }))
      setSaveError('')
    } catch { setSaveError('Browser autosave is unavailable or full. Export the session to keep your work.') }
  }, [data, settings, selectedId])

  const act = (action: PlanningRequest['action'], extra: Partial<PlanningRequest> = {}) => session
    ? run({ action, session, node_id: selected?.node_id, settings, ...extra }, action === 'search' ? 'Finding disconnections…' : action === 'conditions' ? 'Finding conditions…' : 'Updating route…')
    : Promise.resolve(null)
  const start = async () => {
    const next = await run({ action: 'start', target_smiles: target }, 'Starting plan…')
    if (next) { setSelectedId('root'); setNotice(''); setTarget(next.session.root.smiles) }
  }
  const select = (nodeId: string) => { cancel(); setSelectedId(nodeId); setError('') }
  const updateSettings = (patch: Partial<PlanningSettings>) => { cancel(); setSettings(current => ({ ...current, ...patch })) }
  const importFile = async (file: File) => {
    if (file.size > 20_000_000) { setError('Session files must be smaller than 20 MB.'); return }
    try { await restore(JSON.parse(await file.text())) } catch { setError('The file is not valid JSON.') }
  }
  const exportSession = () => data && download({ schema_version: 'interactive_planning_browser.v1',
    session: data.session, settings, selected_id: selectedId }, 'interactive_planning_session.json')
  return { data, settings, target, setTarget, editingTarget, setEditingTarget, selected, busy,
    error, setError, notice, saveError, cancel, act, start, select, updateSettings, importFile, exportSession }
}

type State = ReturnType<typeof useInteractivePlanner>

export function InteractivePlannerOptions({ state }: { state: State }) {
  return <div className="analysis-options planner-options">
    <div className="option-grid feature-options">
      <label><span>Strategies per expansion</span><input aria-label="Strategies per expansion" type="number" min={1} max={50} value={state.settings.top_k} disabled={Boolean(state.busy)} onChange={event => state.updateSettings({ top_k: Math.min(50, Math.max(1, Number(event.target.value))) })} /></label>
      <div className="feature-mode-note"><strong>Build a route one choice at a time</strong><span>Select a molecule, find disconnections, and choose the precursor set to continue. Every precursor becomes its own branch.</span></div>
    </div>
    <details className="advanced-options"><summary>Advanced options</summary><div>
      <label><span>Operator library</span><select value={state.settings.library_mode} disabled={Boolean(state.busy)} onChange={event => state.updateSettings({ library_mode: event.target.value as 'full' | 'compact' })}><option value="full">Full</option><option value="compact">Compact</option></select></label>
      {([['use_context', 'Rank with local reaction context'], ['diversify', 'Prioritize distinct disconnections'], ['use_precursor_realism', 'Consider precursor realism'], ['use_forward_validation', 'Audit forward products and competing pathways'], ['include_l0', 'Include broad L0 fallback operators']] as const).map(([key, label]) => <label className="check-option" key={key}><input type="checkbox" checked={state.settings[key]} disabled={Boolean(state.busy)} onChange={event => state.updateSettings({ [key]: event.target.checked })} /><span>{label}</span></label>)}
    </div></details>
  </div>
}

function RouteNode({ node, path, selectedId, onSelect, stock }: { node: PlanningNode; path: string; selectedId: string; onSelect: (id: string) => void; stock: PlanningSession['stock'] }) {
  const [collapsed, setCollapsed] = useState(false)
  const title = path === '0' ? 'Target' : `Molecule ${path}`
  return <li className="planner-branch">
    <div className={`planner-molecule ${node.node_id === selectedId ? 'selected' : ''}`}>
      <button className="planner-molecule-select" type="button" aria-label={`Select ${title.toLowerCase()}`} aria-pressed={node.node_id === selectedId} onClick={() => onSelect(node.node_id)}>
        <strong>{title}</strong>
        <ReactionImage smiles={node.smiles} kind="molecule" label={`${title} structure`} focusable={false} />
        <span className={`planner-status ${node.stopped ? 'stopped' : ''}`}>{node.stopped ? 'User-designated starting material' : node.choice ? 'Step selected' : 'Needs expansion'}</span>
        {stock[node.smiles]?.status === 'verified_stock_match' && <span className="planner-status stopped">Stock match · last check</span>}
        <code>{node.smiles}</code>
      </button>
      {node.children.length > 0 && <button className="button quiet" type="button" aria-expanded={!collapsed} onClick={() => setCollapsed(value => !value)}>{collapsed ? 'Expand branch' : 'Collapse branch'}</button>}
    </div>
    {node.children.length > 0 && !collapsed && <><div className="planner-connector"><span>Chosen reaction<br />All {node.children.length} precursors required</span></div><ul className="planner-children">
      {node.children.map((child, index) => <RouteNode key={child.node_id} node={child} path={path === '0' ? `${index + 1}` : `${path}.${index + 1}`} selectedId={selectedId} onSelect={onSelect} stock={stock} />)}
    </ul></>}
  </li>
}

function CandidateEvidence({ candidate }: { candidate: RetrosynthesisCandidate }) {
  return <>
    <div className="planner-metrics"><span>Score <strong>{Number(candidate.score ?? 0).toFixed(3)}</strong></span><span>References <strong>{candidate.independent_reference_support ?? 0}</strong></span><span>Context <strong>{candidate.abstraction_level}</strong></span></div>
    <p className="planner-note">Signature: {displayName(candidate.forward_validation_status)} · Forward audit: {displayName(candidate.forward_assessment?.validity ?? 'not evaluated')}</p>
    {(candidate.selectivity_warnings?.length || candidate.precursor_compatibility_assessments?.length) ? <p className="alert caution">This proposal has selectivity or precursor compatibility cautions. Review the evidence before proceeding.</p> : null}
    <details className="planner-evidence"><summary>Validation, precedents, and cautions</summary><pre>{JSON.stringify(candidate, null, 2)}</pre></details>
  </>
}

function StrategyCard({ search, index, state }: { search: PlanningSearch; index: number; state: State }) {
  const strategy = search.result.strategies[index]
  const [variant, setVariant] = useState(() => state.selected?.choice?.search_id === search.search_id && state.selected.choice.strategy_index === index ? state.selected.choice.realization_index : 0)
  const variants = [strategy.representative, ...strategy.alternate_realizations]
  const candidate = variants[variant] ?? variants[0]
  const choice: PlanningChoice = { search_id: search.search_id, strategy_index: index, realization_index: variant }
  const chosen = state.selected?.choice && choiceKey(state.selected.choice) === choiceKey(choice)
  return <article className={`planner-candidate ${chosen ? 'chosen' : ''}`}>
    <div className="planner-candidate-heading"><h3>Strategy {index + 1}</h3><span>{displayName(candidate.transformation_kind ?? 'graph transformation')}</span></div>
    {variants.length > 1 && <label className="planner-variant">Precursor choice<select aria-label={`Strategy ${index + 1} precursor choice`} value={variant} onChange={event => setVariant(Number(event.target.value))}>{variants.map((item, number) => <option key={number} value={number}>{number + 1}. {item.precursor_smiles}</option>)}</select></label>}
    <ReactionImage smiles={candidate.proposed_reaction_smiles} label={`Strategy ${index + 1} reaction`} />
    <CandidateEvidence candidate={candidate} />
    <button className="button primary" type="button" disabled={Boolean(state.busy || chosen)} onClick={() => void state.act('select', choice)}>{chosen ? 'Step selected' : 'Use this step'}</button>
  </article>
}

export function InteractivePlanner({ state, capabilities }: { state: State; capabilities: Capabilities | null }) {
  const fileInput = useRef<HTMLInputElement>(null)
  const { data, selected, settings, busy } = state
  const available = Boolean(capabilities?.retrosynthesis_library_modes?.[settings.library_mode]?.library_available)
  const session = data?.session
  const matchingSearches = session?.searches.filter(search => search.target_smiles === selected?.smiles) ?? []
  const search = [...matchingSearches].reverse().find(item => Object.entries(settings).every(([key, value]) => item.settings[key as keyof PlanningSettings] === value))
  const candidate = session && selected ? selectedCandidate(session, selected.choice) : undefined
  const conditions = selected?.choice ? session?.conditions[choiceKey(selected.choice)] : undefined
  const stock = selected ? session?.stock[selected.smiles] : undefined
  return <section className="interactive-planner" aria-label="Interactive route planning">
    <div className="planner-toolbar">
      <div className="button-row"><button className="button quiet" type="button" disabled={Boolean(busy)} onClick={() => { state.setTarget(''); state.setEditingTarget(true) }}>New target</button>
        <button className="button quiet" type="button" disabled={Boolean(busy || !session?.past.length)} onClick={() => void state.act('undo')}>Undo</button>
        <button className="button quiet" type="button" disabled={Boolean(busy || !session?.future.length)} onClick={() => void state.act('redo')}>Redo</button></div>
      <div className="button-row"><button className="button quiet" type="button" disabled={!data || Boolean(busy)} onClick={state.exportSession}>Export session</button>
        <button className="button quiet" type="button" disabled={!data || Boolean(busy)} onClick={() => data && download(data.route_tree, 'selected_route.json')}>Export route</button>
        <button className="button quiet" type="button" disabled={Boolean(busy)} onClick={() => fileInput.current?.click()}>Import session</button>
        <input ref={fileInput} className="sr-only" aria-label="Import planning session" type="file" accept=".json,application/json" onChange={event => { const file = event.target.files?.[0]; if (file) void state.importFile(file); event.target.value = '' }} /></div>
    </div>
    {state.error && <p className="alert error" role="alert">{state.error}</p>}
    {state.saveError && <p className="alert caution" role="alert">{state.saveError}</p>}
    {state.notice && <p className="planner-note">{state.notice}</p>}
    {busy && <div className="planner-progress"><span role="status" aria-live="polite">{busy} {busy.startsWith('Finding') && 'Bounded searches may take a few minutes.'}</span><button className="button quiet" type="button" onClick={state.cancel}>Cancel</button></div>}
    {state.editingTarget && <div className="planner-target">
      <ReactionEditor value={state.target} onChange={state.setTarget} onError={state.setError} moleculeOnly moleculePurpose="target" disabled={Boolean(busy)} />
      {data && <p className="planner-note">Starting a new plan replaces the saved session. Export the current session to keep both.</p>}
      <div className="button-row"><button className="button primary" type="button" disabled={!state.target.trim() || Boolean(busy)} onClick={() => void state.start()}>Start planning</button>{data && <button className="button quiet" type="button" onClick={() => state.setEditingTarget(false)}>Keep current plan</button>}</div>
    </div>}
    {data && selected && <>
      <div className="planner-summary"><strong>{data.summary.reaction_count} selected steps</strong><span>{data.summary.unresolved_count} unresolved molecules</span><span>{data.summary.starting_material_count} user-designated starting materials</span><span>Depth {data.summary.maximum_depth}</span><small>{state.saveError ? 'Autosave unavailable' : 'Saved in this browser'}</small></div>
      {!data.summary.unresolved_count && <p className="planner-note">All branches end at user-designated starting materials. Stock availability and route feasibility are not established by these choices.</p>}
      <div className="planner-layout">
        <section className="planner-tree-panel" aria-label="Selected route"><h2>Selected route</h2><p className="planner-note" id="planner-tree-help">Target → precursor branches. Select any molecule to continue or revise its step. Scroll to explore the tree.</p><div className="planner-tree" role="region" aria-label="Route tree" aria-describedby="planner-tree-help" tabIndex={0}><ul><RouteNode node={data.session.root} path="0" selectedId={selected.node_id} onSelect={state.select} stock={data.session.stock} /></ul></div></section>
        <section className="planner-inspector" aria-label="Selected molecule planning">
          <div className="planner-inspector-heading"><h2>{selected.node_id === 'root' ? 'Plan the target' : 'Plan this precursor'}</h2><code>{selected.smiles}</code></div>
          <div className="planner-actions"><button className="button primary" type="button" disabled={Boolean(busy || !available || selected.stopped)} onClick={() => void state.act('search')}>{search ? 'Search again' : 'Find disconnections'}</button>
            {selected.choice ? <button className="button quiet" type="button" disabled={Boolean(busy)} onClick={() => void state.act('remove')}>Remove step and its branches</button> : <button className="button quiet" type="button" disabled={Boolean(busy)} onClick={() => void state.act(selected.stopped ? 'reopen' : 'stop')}>{selected.stopped ? 'Reopen for planning' : 'Use as starting material'}</button>}
          </div>
          {!available && capabilities && <p className="alert caution">The selected operator library is unavailable. You can inspect, save, and revise existing choices.</p>}
          {(capabilities?.stock_portfolio_available || stock) && <div className="planner-stock"><button className="button quiet" type="button" disabled={Boolean(busy)} onClick={() => void state.act('stock')}>{stock ? 'Refresh stock check' : 'Check supplier stock'}</button>{stock && <details className="planner-evidence"><summary>{displayName(stock.status)} · last check</summary><p>{stock.status === 'verified_stock_match' ? 'Exact match in the local supplier snapshot. See dates and availability evidence below; a saved check is not live inventory.' : 'No verified availability established by this lookup. See the retained source evidence below.'}</p><pre>{JSON.stringify(stock, null, 2)}</pre></details>}</div>}
          {selected.stopped && <p className="planner-note">You designated this molecule as a starting material. This is not a verified stock match. Reopen it to search further.</p>}
          {candidate && <details className="planner-chosen" open><summary>Chosen step and conditions</summary>
            <ReactionImage smiles={candidate.proposed_reaction_smiles} label="Chosen route step" />
            <CandidateEvidence candidate={candidate} />
            <button className="button secondary" type="button" disabled={Boolean(busy)} onClick={() => void state.act('conditions')}>{conditions ? 'Refresh step conditions' : 'Find step conditions'}</button>
            {conditions && <div className="planner-conditions"><p>{displayName(conditions.status)}</p>{conditions.recommendations?.map((recipe, index) => <p key={index}><strong>{index + 1}.</strong> {compactRecipeSummary(recipe.resolved_recipe)}</p>)}{conditions.warnings?.map((warning, index) => <p key={index}>{displayName(warning)}</p>)}<details className="planner-evidence"><summary>Condition evidence and precedents</summary><pre>{JSON.stringify(conditions, null, 2)}</pre></details></div>}
          </details>}
          {!selected.stopped && <>
            {search ? <><div className="planner-search-heading"><h3>{search.result.strategies.length} candidate strategies</h3><p className="planner-note">{displayName(search.settings.library_mode)} library · up to {search.settings.top_k} strategies requested. Saved results are retained until you search again.</p></div>
              {selected.children.length > 0 && <p className="planner-note">Choosing another step replaces this molecule’s downstream branches. Undo restores them.</p>}
              {!search.result.strategies.length && <p className="alert caution">No candidates found within this search scope. This molecule remains unresolved; try a different scope or designate it as a starting material.</p>}
              {search.result.warnings?.length > 0 && <details className="planner-evidence"><summary>Search coverage and cautions</summary><ul>{search.result.warnings.map((warning, index) => <li key={index}>{warning}</li>)}</ul><pre>{JSON.stringify(search.result.search_diagnostics, null, 2)}</pre></details>}
              {search.result.strategies.map((_, index) => <StrategyCard key={`${selected.node_id}:${search.search_id}:${index}`} search={search} index={index} state={state} />)}
            </> : <p className="planner-note">{matchingSearches.length ? 'Search settings changed. Find disconnections to generate candidates with these settings. Previously chosen steps retain their original evidence.' : 'Find disconnections to see precursor choices for this molecule.'}</p>}
          </>}
        </section>
      </div>
    </>}
  </section>
}
