import { useCallback, useEffect, useRef, useState } from 'react'
import type { PlanningChoice, PlanningNode, PlanningSession } from '../api/planning'
import { ReactionImage } from './ReactionImage'

export function planningChoiceKey(choice: PlanningChoice): string {
  return `${choice.search_id}:${choice.strategy_index}:${choice.realization_index}`
}

interface TreeProps {
  session: PlanningSession
  selectedId: string
  busy: boolean
  onSelect: (id: string) => void
  onChoose: (id: string, choice: PlanningChoice) => void
}

function MoleculeNode({ node, path, active, ...props }: TreeProps & {
  node: PlanningNode; path: string; active: boolean
}) {
  const [collapsed, setCollapsed] = useState(false)
  const title = path ? `Molecule ${path}` : 'Target'
  return <li className={`planner-branch molecule-branch ${active ? 'on-route' : ''}`}>
    <div className={`planner-molecule ${node.node_id === props.selectedId ? 'selected' : ''} ${node.stopped ? 'stopped' : ''}`} data-node-id={node.node_id}>
      <button className="planner-molecule-select" type="button" aria-label={`Select ${title.toLowerCase()}`} aria-pressed={node.node_id === props.selectedId} onClick={() => props.onSelect(node.node_id)}>
        <strong>{title}</strong>
        <ReactionImage smiles={node.smiles} kind="molecule" label={`${title} structure`} focusable={false} />
        <span className={`planner-status ${node.stopped ? 'stopped' : ''}`}>{node.stopped ? 'User-designated starting material' : node.choice ? 'Step selected' : node.alternatives.length ? 'Choose a reaction' : 'Needs expansion'}</span>
        {props.session.stock[node.smiles]?.status === 'verified_stock_match' && <span className="planner-status stopped">Stock match · last check</span>}
        <code>{node.smiles}</code>
      </button>
      {node.alternatives.length > 0 && <button className="button quiet" type="button" aria-expanded={!collapsed} onClick={() => setCollapsed(value => !value)}>{collapsed ? `Expand ${node.alternatives.length} alternatives` : 'Collapse alternatives'}</button>}
    </div>
    {node.alternatives.length > 0 && !collapsed && <ul className="planner-alternatives" aria-label={`Alternative reactions for ${title.toLowerCase()}`}>
      {node.alternatives.map((option, index) => {
        const key = planningChoiceKey(option.choice)
        const chosen = node.choice !== null && planningChoiceKey(node.choice) === key
        const onRoute = active && chosen
        const reactionPath = path ? `${path}/${index + 1}` : `${index + 1}`
        return <li key={key} className={`planner-branch reaction-branch ${onRoute ? 'on-route' : ''}`} data-choice-id={key}>
          <button className={`planner-reaction ${onRoute ? 'on-route' : ''}`} type="button" aria-label={`Use reaction ${reactionPath}`} aria-pressed={onRoute} disabled={props.busy} onClick={() => props.onChoose(node.node_id, option.choice)} title={`Strategy ${option.choice.strategy_index + 1}, precursor choice ${option.choice.realization_index + 1}. Select this reaction and its path to the target.`}>
            <strong>#{index + 1}</strong><span>{onRoute ? 'Selected route' : chosen ? 'Selected here' : 'Alternative'}</span>
            <small>All {option.children.length} precursors</small>
          </button>
          <ul className="planner-precursors" aria-label={`Required precursors for reaction ${reactionPath}`}>
            {option.children.map((child, childIndex) => <MoleculeNode key={child.node_id} {...props} node={child} path={`${reactionPath}.${childIndex + 1}`} active={onRoute} />)}
          </ul>
        </li>
      })}
    </ul>}
  </li>
}

export function PlanningTree(props: TreeProps) {
  const viewport = useRef<HTMLDivElement>(null)
  const content = useRef<HTMLUListElement>(null)
  const focus = useCallback(() => {
    const selected = content.current?.querySelector<HTMLElement>('.planner-molecule.selected')
    const panel = viewport.current
    if (!selected || !panel) return
    const box = selected.getBoundingClientRect()
    const bounds = panel.getBoundingClientRect()
    panel.scrollBy({ left: box.left - bounds.left - (panel.clientWidth - box.width) / 2,
      top: box.top - bounds.top - (panel.clientHeight - box.height) / 2 })
  }, [])
  useEffect(() => {
    if (!viewport.current || !content.current) return
    // Keep the current molecule in view as branches grow, without shrinking it.
    const observer = new ResizeObserver(focus)
    observer.observe(viewport.current)
    observer.observe(content.current)
    focus()
    return () => observer.disconnect()
  }, [focus, props.selectedId])
  return <section className="planner-tree-panel" aria-label="Route alternatives">
    <div className="planner-tree-heading"><h2>Route alternatives</h2><span className="planner-tree-legend">Green path: selected route · Blue outline: selected molecule</span></div>
    <p className="planner-note" id="planner-tree-help">Select a molecule to find candidates, then use Add to tree for the reaction you want. Each numbered reaction is an added alternative; its precursors are all required. Other added branches stay saved.</p>
    <div className="planner-tree-tools" aria-label="Tree view controls">
      <span className="planner-note">Fixed molecule size · Scroll to explore branches</span>
      <button className="button quiet" type="button" onClick={focus}>Focus molecule</button>
    </div>
    <div ref={viewport} className="planner-tree" role="region" aria-label="Route tree" aria-describedby="planner-tree-help" tabIndex={0}>
      <div className="planner-tree-canvas"><ul ref={content} className="planner-tree-content">
        <MoleculeNode {...props} node={props.session.root} path="" active />
      </ul></div>
    </div>
  </section>
}
