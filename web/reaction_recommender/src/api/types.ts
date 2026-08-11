export type JsonObject = Record<string, unknown>

export interface ApiEnvelope<T> {
  api_schema_version: string
  data: T
}

export interface Capabilities {
  service: string
  index_name: string
  index_available: boolean
  loaded_runtime_variants: number
  rxnmapper_available: boolean
  recommendation: boolean
  discovery: boolean
  featurization: boolean
  reaction_rendering: boolean
  retrosynthesis?: boolean
  local_only: boolean
  default_library_mode?: 'full' | 'compact'
  library_modes?: Record<string, {
    label: string
    index_name: string
    index_available: boolean
  }>
  default_retrosynthesis_library_mode?: 'full' | 'compact'
  retrosynthesis_library_modes?: Record<string, {
    label: string
    library_name: string
    library_available: boolean
  }>
}

export interface RankingProfile {
  profile_id: string
  label: string
  description: string
  weights: Record<string, number>
}

export interface CompletionOption {
  option_id: string
  option_kind: string
  display_name: string
  capability_id?: string | null
  substance_id?: string | null
}

export interface CompletionRequirement {
  requirement_id: string
  fragment_key: string
  canonical_fragment_smiles: string
  attachment_element: string
  options: CompletionOption[]
}

export interface CompletionProposal {
  query_reaction_smiles: string
  proposal_id: string
  status: string
  requirements: CompletionRequirement[]
  warnings: string[]
}

export interface PrepareReactionResult {
  valid: boolean
  input_reaction_smiles: string
  completion_proposal: CompletionProposal
  warnings: string[]
}

export interface CompletionChoice {
  requirement_id: string
  option_id?: string
  custom_identifier?: string
}

export interface RecipeComponent {
  cas?: string | null
  canonical_name?: string
  display_name?: string
  name?: string
  raw_identifier?: string
  substance_id?: string
  [key: string]: unknown
}

export interface ProtocolMaterial {
  material_id: string
  category: 'reaction_input' | 'condition' | 'reaction_output'
  identity_status: string
  substance_id?: string | null
  canonical_name?: string | null
  cas?: string | null
  smiles?: string | null
  role?: string | null
  role_status: string
  amount?: number | null
  amount_unit?: string | null
  source_field?: string | null
  warnings: string[]
}

export interface ProtocolOperatingConditions {
  temperature_c?: number | null
  time_h?: number | null
  concentration_m?: number | null
  atmosphere?: string | null
}

export interface ProtocolOperation {
  operation_id: string
  operation_type: 'maintain_conditions'
  stage_index: number
  temperature_c?: number | null
  time_h?: number | null
  atmosphere?: string | null
}

export interface SynthesisProtocolDraft {
  protocol_id: string
  recipe_id: string
  recipe_core_id: string
  reaction_smiles?: string | null
  materials: ProtocolMaterial[]
  operating_conditions: ProtocolOperatingConditions
  operations: ProtocolOperation[]
  execution_readiness: 'review_required'
  missing_required_fields: string[]
  warnings: string[]
  schema_version: string
}

export interface ResolvedRecipe {
  recipe_id?: string
  recipe_core_id?: string
  catalysts?: RecipeComponent[]
  ligands?: RecipeComponent[]
  bases?: RecipeComponent[]
  condensation_agents?: RecipeComponent[]
  oxidants?: RecipeComponent[]
  reductants?: RecipeComponent[]
  acids?: RecipeComponent[]
  additives?: RecipeComponent[]
  solvents?: RecipeComponent[]
  other_components?: RecipeComponent[]
  temperature_c?: number | null
  time_h?: number | null
  concentration_m?: number | null
  pressure_bar?: number | null
  atmosphere?: string | null
  [key: string]: unknown
}

export interface ReferenceRecord {
  reference_id: string
  raw_reference?: string | null
  normalized_citation?: string | null
  doi?: string | null
  patent_number?: string | null
  publication_year?: number | null
  resolution_status?: string | null
}

export interface ExperimentalDetail {
  observation_id?: string
  reaction_id?: string
  source_dataset?: string
  reference_id?: string
  procedure_text: string
  notes?: string | null
  stages?: unknown
  steps?: unknown
}

export interface ScoreTrace {
  ranking_profile: string
  ranking_components: Record<string, number | null>
  ranking_contributions: Record<string, number>
  applied_ranking_weights: Record<string, number>
  independent_evidence_count: number
  observed_outcome_count: number
  definition_versions: Record<string, string>
}

export interface Recommendation {
  rank: number
  default_rank?: number | null
  rank_change: number
  recipe_id: string
  recipe_core_id: string
  score: number
  default_score?: number | null
  similarity_score: number
  compatibility_score: number
  expected_yield_pct?: number | null
  support: number
  reference_support: number
  dataset_support: number
  retrieval_level: string
  resolved_recipe: ResolvedRecipe
  synthesis_protocol?: SynthesisProtocolDraft
  explanation: string[]
  compatibility_evidence: string[]
  cautions: string[]
  precedent_reaction_smiles: string[]
  precedent_reaction_ids: string[]
  precedent_reference_ids: string[]
  precedent_references?: ReferenceRecord[]
  precedent_experimental_details?: ExperimentalDetail[]
  score_trace: ScoreTrace
  factor_evidence: JsonObject
}

export interface RetrievalTrace {
  level: string
  candidate_count: number
  independent_candidate_count: number
  compatible_candidate_count: number
  independent_compatible_candidate_count: number
  excluded_candidate_count: number
  status: string
}

export interface RecommendationResult {
  query_reaction_smiles: string
  effective_query_reaction_smiles?: string | null
  valid: boolean
  error?: string | null
  schema_version: string
  query_signature_id?: string | null
  query_reaction_core_id?: string | null
  named_family?: string | null
  transformation_class?: string | null
  reaction_label?: { text?: string; concise?: string; detailed?: string } | null
  recommendation_mode: string
  retrieval_level?: string | null
  candidate_count: number
  independent_candidate_count: number
  compatible_candidate_count: number
  warnings: string[]
  retrieval_trace: RetrievalTrace[]
  recommendations: Recommendation[]
  ranking_preferences: JsonObject
}

export interface DiscoveryScoreTrace {
  components: Record<string, number | null>
  contributions: Record<string, number>
  configured_weights: Record<string, number>
  effective_weights: Record<string, number>
  matches: string[]
  mismatches: string[]
}

export interface DiscoveryHit {
  rank: number
  reaction_id: string
  observation_id: string
  reaction_smiles: string
  relation_class: string
  relation_tiers: string[]
  discovery_score: number
  yield_pct?: number | null
  outcome_status: string
  evidence_tier: string
  chemistry_status: string
  source_dataset: string
  reference_id: string
  reference_record?: ReferenceRecord | null
  experimental_detail?: ExperimentalDetail | null
  resolved_recipe: ResolvedRecipe
  recipe_id: string
  recipe_core_id: string
  hypothesis_id?: string | null
  score_trace: DiscoveryScoreTrace
  insights: string[]
  cautions: string[]
}

export interface DiscoveryResult {
  query_reaction_smiles: string
  valid: boolean
  error?: string | null
  schema_version: string
  query_signature_id?: string | null
  query_reaction_core_id?: string | null
  named_family?: string | null
  transformation_class?: string | null
  discovery_view: string
  candidate_count: number
  relation_counts: Record<string, number>
  warnings: string[]
  hits: DiscoveryHit[]
}

export interface RecommendationRequest {
  reaction_smiles: string
  library_mode: 'full' | 'compact'
  top_k: number
  minimum_pool_size: number | null
  unrestricted_fallback: boolean
  use_rxnmapper: boolean
  ranking_preferences: {
    profile_id: string
    weights: Record<string, number>
  }
  completion_choices: CompletionChoice[]
}

export interface DiscoveryRequest {
  reaction_smiles: string
  library_mode: 'full' | 'compact'
  top_k: number
  view: string
  include_low_yield: boolean
  include_unreported_outcomes: boolean
  use_rxnmapper: boolean
  include_review: boolean
}

export interface FeatureMotif {
  motif_id: string
  label: string
  component_index: number
  side?: string
  atom_indices: number[]
  tags: string[]
  confidence: number
}

export interface FeatureReactiveSite {
  hypothesis_id: string
  label: string
  site_type: string
  availability: string
  component_index: number
  side?: string
  atom_indices: number[]
  context_kind?: string | null
}

export interface ReactionCoreEvent {
  event_id: string
  edit_tokens: string[]
  reactant_component_indices: number[]
  product_component_indices: number[]
}

export interface ReactionCoreSummary {
  core_id: string
  evidence: string
  evidence_status: string
  confidence: number
  active_atom_count: number
  event_count: number
  events: ReactionCoreEvent[]
  quality: {
    status: string
    review_reasons: string[]
    blocking_reasons: string[]
    checked_edit_fraction: number
    active_atom_mapping_coverage: number
  }
  warnings: string[]
}

export interface FeaturePartner {
  role?: string | null
  component_index: number
  label: string
  role_confidence: number
  anchor_contexts: string[]
}

export interface FeatureAnalysisResult {
  input_kind: 'molecule' | 'reaction'
  input_smiles: string
  valid: boolean
  schema_version: string
  overview: Record<string, unknown>
  motifs: FeatureMotif[]
  reactive_sites: FeatureReactiveSite[]
  reaction_core?: ReactionCoreSummary | null
  partners: FeaturePartner[]
  mapping?: Record<string, unknown> | null
  core_graphic_svg?: string | null
  core_graphic_error?: string | null
  core_projection?: Record<string, unknown> | null
  warnings: string[]
  error?: string | null
  analysis: JsonObject
}

export interface FeatureAnalysisRequest {
  input_smiles: string
  use_rxnmapper: boolean
  force_resolved_mapping: boolean
}

export interface RetrosynthesisRequest {
  target_smiles: string
  library_mode: 'full' | 'compact'
  top_k: number
  include_l0: boolean
  use_context: boolean
  diversify: boolean
}

export interface RetrosynthesisCandidate {
  rank: number
  target_smiles: string
  precursor_smiles: string
  proposed_reaction_smiles: string
  transformation_kind?: string | null
  abstraction_level: string
  compiler_engine: string
  template_id: string
  score: number
  context_similarity: number
  product_similarity: number
  precursor_similarity: number
  template_specificity: number
  independent_reference_support: number
  forward_validation_status: string
  center_transition_key: string
  disconnection_site_key: string
  precedent_reaction_ids: string[]
  operator_id: string
  realization_id: string
  operator_signature: string
  synthon_signature: string
  pre_diversity_rank: number
  diversity_rank: number
  diversity_group_key: string[]
  structural_score_band: number
  ranking_policy_definition_id: string
}

export interface RetrosynthesisResult {
  target_smiles: string
  library_mode: 'full' | 'compact'
  valid: boolean
  error?: string | null
  schema_version: string
  candidate_count: number
  library_operator_count: number
  library_template_count: number
  warnings: string[]
  candidates: RetrosynthesisCandidate[]
}
