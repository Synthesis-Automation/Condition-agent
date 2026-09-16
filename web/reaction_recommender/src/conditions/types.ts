import type { JsonObject, ResolvedRecipe, SynthesisProtocolDraft, Recommendation, WeakLabelRecommendation, RecommendationResult, WeakLabelRecommendationResult } from '../api/types'

export type ConditionSearchScope = 'same_handle' | 'automatic' | 'broad'

export interface ConditionQueryChemistry {
  reaction_label?: RecommendationResult['reaction_label']
  query_participants?: WeakLabelRecommendationResult['query_participants']
  reaction_partners?: Array<{
    component_index: number
    site_type: string
    nearby_groups?: Array<{ label: string }>
  }>
}

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
    result: ConditionQueryChemistry
  }>
  warnings: string[]
  shortlist_size: number
  schema_version: string
  automation_exports: Record<string, JsonObject>
}
