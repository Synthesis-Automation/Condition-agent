import type { CanonicalRouteTree, RetrosynthesisConditionEvidence, RetrosynthesisRequest, RetrosynthesisResult } from './types'

export type PlanningSettings = Omit<RetrosynthesisRequest, 'target_smiles'>
export interface PlanningChoice { search_id: string; strategy_index: number; realization_index: number }
export interface PlanningAlternative {
  choice: PlanningChoice
  children: PlanningNode[]
}
export interface PlanningNode {
  node_id: string
  smiles: string
  stopped: boolean
  choice: PlanningChoice | null
  alternatives: PlanningAlternative[]
  expansion_warnings: string[]
}
export interface PlanningSearch {
  search_id: string
  target_smiles: string
  settings: PlanningSettings
  result: RetrosynthesisResult
}
export interface PlanningSession {
  schema_version: 'interactive_planning.v2'
  root: PlanningNode
  past: PlanningNode[]
  future: PlanningNode[]
  searches: PlanningSearch[]
  conditions: Record<string, RetrosynthesisConditionEvidence>
  stock: Record<string, { status: string; checked_at: string; source_records: Array<Record<string, string>> }>
}
export interface PlanningResponse {
  session: PlanningSession
  route_tree: CanonicalRouteTree
  summary: { reaction_count: number; unresolved_count: number; starting_material_count: number; maximum_depth: number }
}
export interface PlanningRequest {
  action: 'start' | 'restore' | 'search' | 'select' | 'clear' | 'remove' | 'stop' | 'reopen' | 'undo' | 'redo' | 'conditions' | 'stock'
  session?: PlanningSession
  target_smiles?: string
  node_id?: string
  settings?: PlanningSettings
  search_id?: string
  strategy_index?: number
  realization_index?: number
}
