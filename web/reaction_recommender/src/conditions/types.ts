import type { JsonObject, ResolvedRecipe, SynthesisProtocolDraft, Recommendation, WeakLabelRecommendation } from '../api/types'

export interface ConditionEvidence {
  source: 'generic' | 'weak_label'
  recommendation_mode: string
  warnings: string[]
  recommendation: Partial<Recommendation & WeakLabelRecommendation>
}

export interface ConditionOption {
  option_id: string
  rank: number
  evidence_kind: 'verified_signature' | 'structure_review' | 'weak_label'
  evidence_label: string
  resolved_recipe: ResolvedRecipe
  synthesis_protocol: SynthesisProtocolDraft
  evidence: ConditionEvidence[]
  cautions: string[]
}

export interface ConditionsResult {
  query_reaction_smiles: string
  valid: boolean
  recommendations: ConditionOption[]
  sources: Array<{
    source: 'generic' | 'weak_label'
    status: 'ok' | 'abstained' | 'unavailable' | 'skipped'
    message: string
    result: JsonObject
  }>
  warnings: string[]
  shortlist_size: number
  schema_version: string
  automation_exports: Record<string, JsonObject>
}
