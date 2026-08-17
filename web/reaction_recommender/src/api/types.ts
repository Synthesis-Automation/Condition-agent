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
  weak_label_recommendation?: boolean
  weak_label_dataset_name?: string
  featurization: boolean
  reaction_rendering: boolean
  retrosynthesis?: boolean
  forward_synthesis?: boolean
  multistep_retrosynthesis?: boolean
  literature_molecule_index_available?: boolean
  literature_molecule_index_name?: string
  stock_portfolio_available?: boolean
  stock_portfolio_name?: string
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
  forward_library_modes?: Record<string, {
    label: string
    library_name: string
    library_available: boolean
    prepared: boolean
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
  condition_precedents?: ConditionPrecedent[]
  score_trace: ScoreTrace
  factor_evidence: JsonObject
}

export interface ConditionPrecedent {
  reaction_id: string
  observation_id: string
  reaction_smiles: string
  reference_id: string
  reference_record?: ReferenceRecord | null
  experimental_detail?: ExperimentalDetail | null
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

export interface WeakLabelParticipant {
  component_index: number
  site_id: string
  site_type: string
  canonical_signature: string
  chemist_label: string
  active_atom_indices: number[]
  role?: string | null
  role_confidence: number
  center_substitution_class?: string | null
  attachment_carbon_classes: string[]
  alpha_branched?: boolean | null
  evidence: string[]
  schema_version: string
}

export interface WeakLabelRecommendation {
  rank: number
  recipe_id: string
  resolved_recipe: ResolvedRecipe
  score: number
  label_similarity: number
  signature_similarity: number
  qualifier_similarity: number
  compatibility_score: number
  expected_yield_pct?: number | null
  mean_z_score?: number | null
  support: number
  source_reaction_types: string[]
  source_row_numbers: number[]
  source_matches: WeakLabelSourceMatch[]
  explanation: string[]
  compatibility_evidence: string[]
  cautions: string[]
}

export interface WeakLabelSourceMatch {
  source_row_number: number
  source_reaction_type: string
  participant_roles: [string, string]
  participant_display_labels: [string, string]
  participant_signatures: [string, string]
  label_similarity: number
  signature_similarity: number
  qualifier_similarity: number
}

export interface WeakLabelRecommendationResult {
  query_reaction_smiles: string
  valid: boolean
  recommendation_mode: 'weak_label_fallback' | 'weak_label_screening'
  reaction_type_hint_id?: string | null
  reaction_type_id?: string | null
  source_reaction_type_candidates: string[]
  query_participants: WeakLabelParticipant[]
  candidate_count: number
  compatible_candidate_count: number
  excluded_candidate_count: number
  recipe_count: number
  recommendations: WeakLabelRecommendation[]
  warnings: string[]
  error?: string | null
  schema_version: string
}

export type RecommendationApiResult = RecommendationResult | WeakLabelRecommendationResult

export interface RecommendationRequest {
  reaction_smiles: string
  recommendation_mode: 'generic' | 'weak_label_fallback' | 'weak_label_screening'
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
  use_precursor_realism: boolean
  use_forward_validation: boolean
}

export interface ForwardSynthesisRequest {
  starting_materials: string
  intended_product?: string | null
  operator_hint?: string | null
  recipe?: JsonObject | null
  condition_profile?: ForwardConditionProfile | null
  library_mode: 'full' | 'compact'
  top_k: number
  include_l0: boolean
  include_self_reactions: boolean
}

export interface ForwardConditionProfileOption {
  id: string
  label: string
  description?: string
}

export interface ForwardConditionProfileCatalog {
  definition_id: string
  definition_version: string
  schema_version: string
  strategies: ForwardConditionProfileOption[]
  redox_modes: ForwardConditionProfileOption[]
  media: ForwardConditionProfileOption[]
  catalyst_families: ForwardConditionProfileOption[]
}

export interface ForwardConditionProfile {
  strategy: string
  redox_mode: string
  medium: string
  catalyst_family: string
  source?: string
  schema_version?: string
}

export interface ForwardConditionProfileEvidence {
  evaluated: boolean
  profile: ForwardConditionProfile
  score_adjustment: number
  compatible?: boolean | null
  hard_conflicts: string[]
  matched_rules: string[]
  cautions: string[]
  definition_id: string
  definition_version: string
}

export interface ForwardRecipeEvidence {
  evaluated: boolean
  compatible?: boolean | null
  score?: number | null
  hard_conflicts: string[]
  cautions: string[]
  definition_id: string
  definition_version: string
}

export interface OperatorAtomCorrespondence {
  product_component_index: number
  product_atom_index: number
  precursor_component_index: number
  precursor_atom_index: number
  atom_map_number: number
  operator_map_number?: number | null
  precursor_instance_index: number
}

export interface ForwardProductCandidate {
  rank: number
  product_smiles: string
  reaction_smiles: string
  mapped_reaction_smiles: string
  pathway_id: string
  operator_id: string
  forward_operator_id: string
  realization_id: string
  template_id: string
  abstraction_level: string
  participating_component_indices: number[]
  participating_precursor_smiles: string
  assignment: number[]
  reactant_stoichiometry: Array<[number, number]>
  uses_virtual_copies: boolean
  atom_correspondence: OperatorAtomCorrespondence[]
  score: number
  score_components: Record<string, number>
  structural_score_band: number
  reverse_round_trip: boolean
  reaction_signature_id?: string | null
  reaction_signature_schema_version?: string | null
  operator_edit_agreement: boolean
  observed_edit_tokens: string[]
  independent_reference_support: number
  observation_support: number
  precedent_reaction_ids: string[]
  precedent_reference_ids: string[]
  recipe_evidence: ForwardRecipeEvidence
  condition_profile_evidence: ForwardConditionProfileEvidence
  alternative_pathway_ids: string[]
  alternative_operator_ids: string[]
  alternative_template_ids: string[]
  named_annotations: string[]
  warnings: string[]
  schema_version: string
}

export interface ForwardSearchDiagnostics {
  library_operator_count: number
  indexed_operator_count: number
  applied_operator_count: number
  generated_outcome_count: number
  reverse_round_trip_failure_count: number
  invalid_reaction_count: number
  missing_signature_count: number
  operator_edit_mismatch_count: number
  recipe_conflict_count: number
  condition_profile_conflict_count: number
  valid_pathway_count: number
  self_reaction_pathway_count: number
  unique_product_count: number
}

export interface ForwardCompetitionGroup {
  competition_level: 'operator' | 'site' | 'product'
  group_key: string
  candidate_ranks: number[]
  product_smiles: string[]
  operator_ids: string[]
}

export interface ForwardPrediction {
  query_starting_materials: string
  canonical_starting_materials: string
  conditions_supplied: boolean
  condition_profile_supplied: boolean
  condition_profile: ForwardConditionProfile
  self_reactions_considered: boolean
  valid: boolean
  status: string
  candidates: ForwardProductCandidate[]
  competition_groups: ForwardCompetitionGroup[]
  diagnostics: ForwardSearchDiagnostics
  warnings: string[]
  error?: string | null
  ranking_definition_id: string
  schema_version: string
}

export interface ForwardStepAssessment {
  starting_materials: string
  intended_product: string
  intended_match: 'exact' | 'stereo_relaxed' | 'connectivity_only' | 'absent' | 'invalid'
  targeted_replay_status: string
  intended_product_rank?: number | null
  best_competitor_product?: string | null
  score_margin?: number | null
  disposition: 'clear' | 'competitive' | 'unsupported' | 'structurally_inconsistent' | 'condition_incompatible' | 'out_of_scope'
  validity: 'structurally_supported' | 'structurally_supported_with_competition' | 'inconclusive' | 'contradicted' | 'out_of_scope'
  checks: ForwardValidityCheck[]
  advisory_only: boolean
  operator_hint?: string | null
  warnings: string[]
  schema_version: string
}

export interface ForwardValidityCheck {
  check_id: string
  status: 'pass' | 'warning' | 'fail' | 'not_evaluated'
  detail: string
}

export interface ForwardBlindPredictionSummary {
  valid: boolean
  status: string
  conditions_supplied: boolean
  condition_profile_supplied: boolean
  candidate_count: number
  valid_pathway_count: number
  top_products: Array<{
    rank: number
    product_smiles: string
    score: number
    is_intended: boolean
  }>
  warnings: string[]
}

export interface RetrosynthesisForwardAssessment extends ForwardStepAssessment {
  evaluated: boolean
  blind_prediction_summary: ForwardBlindPredictionSummary
}

export interface ForwardSynthesisResult {
  analysis_mode: 'blind_prediction' | 'step_assessment'
  library_mode: 'full' | 'compact'
  valid: boolean
  schema_version: string
  forward_library_operator_count: number
  prediction: ForwardPrediction
  assessment?: ForwardStepAssessment | null
}

export interface MultistepRetrosynthesisRequest {
  target_smiles: string
  library_mode: 'full' | 'compact'
  top_k_routes: number
  max_depth: 2 | 3
  molecular_weight_threshold: number
  include_l0: boolean
  use_context: boolean
  diversify: boolean
  use_precursor_realism: boolean
  use_condition_availability: boolean
  use_forward_validation: boolean
}

export interface RetrosynthesisConditionsRequest {
  reaction_smiles: string
  library_mode: 'full' | 'compact'
  top_k: number
  starting_materials?: string | null
  intended_product?: string | null
  operator_hint?: string | null
  use_forward_validation: boolean
  include_l0: boolean
}

export interface RetrosynthesisConditionEvidence {
  status: 'pending' | 'recommended_direct' | 'recommended_fallback' | 'insufficient_evidence'
  query_reaction_smiles: string
  recommender_valid: boolean
  recommendation_mode: string
  retrieval_level?: string | null
  uses_type_agnostic_fallback: boolean
  candidate_count: number
  independent_candidate_count: number
  compatible_candidate_count: number
  independent_compatible_candidate_count: number
  excluded_candidate_count: number
  best_recipe_score?: number | null
  best_recipe_compatibility_score?: number | null
  best_recipe_reference_support: number
  recommendations: Recommendation[]
  warnings: string[]
  error?: string | null
  forward_assessment?: RetrosynthesisForwardAssessment | null
}

export interface RetrosynthesisSupportingPrecedent {
  reaction_smiles: string
  reference_record?: ReferenceRecord | null
}

export interface FunctionalGroupCompetitionOutcome {
  candidate_id: string
  component_index: number
  atom_index: number
  element: string
  site_type: string
  site_signature: string
  availability: string
  product_smiles: string
}

export interface FunctionalGroupCompetitionWarning {
  code: string
  message: string
  selected_outcome: FunctionalGroupCompetitionOutcome
  competing_outcomes: FunctionalGroupCompetitionOutcome[]
  assessment_mode: string
  ranking_impact: 'none'
  conditions_evaluated: false
  audit_id: string
  schema_version: string
}

export interface ReactivePairSiteReference {
  hypothesis_id: string
  component_index: number
  atom_index: number
  site_type: string
  canonical_signature: string
  chemist_label: string
  availability: string
}

export interface ReactivePairInteractionAssessment {
  assessment_id: string
  rule_id: string
  interaction_class: string
  scope: 'same_component'
  component_index: number
  component_smiles: string
  left_site: ReactivePairSiteReference
  right_site: ReactivePairSiteReference
  graph_distance?: number | null
  potential_closure_ring_size?: number | null
  intrinsic_severity: string
  warning_strength: 'advisory' | 'strong'
  message: string
  definition_id: string
  schema_version: string
}

export interface MoleculeComplexityTrace {
  canonical_smiles: string
  heavy_atom_count: number
  heavy_bond_count: number
  cycle_rank: number
  ring_system_count: number
  largest_ring_system_atom_count: number
  fused_ring_bond_count: number
  branching_excess: number
  assigned_stereocenter_count: number
  atom_environment_count: number
  raw_complexity: number
  definition_id: string
  schema_version: string
}

export interface RetrosyntheticComplexityReduction {
  score: number
  strategic_class: string
  is_strategic: boolean
  evidence: string
  target: MoleculeComplexityTrace
  precursors: MoleculeComplexityTrace[]
  product_derived_component_heavy_atom_counts: number[]
  largest_retained_core_fraction: number
  core_fragmentation_score: number
  ring_topology_reduction_score: number
  graph_complexity_reduction_score: number
  graph_complexity_reduction_fraction: number
  stereochemical_simplification_score: number
  convergency_score: number
  tactical_penalty: number
  formed_product_bond_count: number
  formed_product_ring_bond_count: number
  product_cycle_rank_reduction: number
  extra_precursor_heavy_atom_count: number
  warnings: string[]
  definition_id: string
  schema_version: string
}

export interface RetrosynthesisCandidate {
  rank: number
  target_smiles: string
  precursor_smiles: string
  proposed_reaction_smiles: string
  condition_query_reaction_smiles?: string
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
  forward_assessment?: RetrosynthesisForwardAssessment | null
  center_transition_key: string
  disconnection_site_key: string
  supporting_precedents: RetrosynthesisSupportingPrecedent[]
  condition_evidence: RetrosynthesisConditionEvidence
  operator_id: string
  realization_id: string
  operator_signature: string
  synthon_signature: string
  pre_diversity_rank: number
  diversity_rank: number
  diversity_group_key: string[]
  structural_score_band: number
  ranking_policy_definition_id: string
  selectivity_warnings?: FunctionalGroupCompetitionWarning[]
  precursor_compatibility_assessments?: ReactivePairInteractionAssessment[]
  precursor_compatibility_disposition?: 'pass' | 'warn' | 'demote' | 'reject'
  precursor_compatibility_warning_strength?: 'advisory' | 'strong' | null
  precursor_compatibility_band_penalty?: number
  precursor_compatibility_policy_definition_id?: string
  strategic_complexity?: RetrosyntheticComplexityReduction | null
  strategic_complexity_score: number
  strategic_class: string
  strategic_candidate: boolean
  strategic_reserve_selected: boolean
  precursor_realism_score?: number | null
  precursor_realism_assessments?: PrecursorRealismAssessment[]
  precursor_realism_aggregation?: PrecursorRealismAggregation | null
  pre_realism_rank: number
  precursor_realism_rank: number
  precursor_realism_band_penalty: number
  effective_structural_score_band: number
  completion_group_id: string
  completion_prior_probability?: number | null
  completion_prior_backoff_level: string
  completion_prior_independent_support: number
  completion_prior_total_support: number
  completion_prior_alternative_count: number
  hierarchical_site_score: number
  hierarchical_synthon_score: number
  hierarchical_realization_score: number
  hierarchical_partition_key: [string, number] | []
  hierarchical_site_rank: number
  hierarchical_synthon_rank: number
  hierarchical_realization_rank: number
  pre_hierarchical_rank: number
  hierarchical_rank: number
  hierarchical_ranking_definition_id: string
}

export interface PrecursorRealismAssessment {
  canonical_smiles: string
  inchi_key?: string | null
  molecular_weight: number
  evidence: {
    buyable: boolean
    in_compound_registry: boolean
    in_literature: boolean
  }
  evidence_tier: string
  base_score: number
  molecular_weight_smallness: number
  molecular_weight_penalty: number
  score: number
  definition_id: string
  schema_version: string
}

export interface PrecursorRealismAggregation {
  weakest_component_score: number
  known_substantial_component_bonus: number
  score: number
  supporting_component_smiles: string | null
  supporting_evidence_tier: string | null
  substantive_component_molecular_weight_threshold_da: number
  definition_id: string
  schema_version: string
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
  strategic_complexity_definition_id: string
  strategic_candidate_count: number
  precursor_realism_enabled: boolean
  forward_validation_enabled: boolean
  forward_validity_counts: Record<string, number>
  precursor_realism_sources: {
    buyable: boolean
    compound_registry: boolean
    literature: boolean
  }
  warnings: string[]
  candidates: RetrosynthesisCandidate[]
}

export interface MoleculeIndexMatch {
  canonical_smiles: string
  inchi_key?: string | null
  occurrence_count: number
  source_records: Array<Record<string, string>>
}

export interface StartingMaterialAssessment {
  route_node_id: string
  smiles: string
  canonical_smiles: string
  depth: number
  molecular_weight?: number | null
  terminal: boolean
  terminal_reasons: string[]
  terminal_evidence: string
  catalog_role_status: string
  unresolved_reason?: string | null
  literature_match?: MoleculeIndexMatch | null
}

export interface MultistepRouteCandidate {
  proposed_reaction_smiles: string
  precursor_smiles: string
  transformation_kind?: string | null
  abstraction_level: string
  score: number
  independent_reference_support: number
  forward_validation_status: string
  operator_id: string
  realization_id: string
  disconnection_site_key: string
  synthon_signature: string
  condition_query_reaction_smiles?: string
  precursor_compatibility_assessments?: ReactivePairInteractionAssessment[]
  precursor_compatibility_disposition?: 'pass' | 'warn' | 'demote' | 'reject'
  precursor_compatibility_warning_strength?: 'advisory' | 'strong' | null
  precursor_compatibility_band_penalty?: number
  precursor_compatibility_policy_definition_id?: string
  strategic_complexity?: RetrosyntheticComplexityReduction | null
  strategic_complexity_score: number
  strategic_class: string
  strategic_candidate: boolean
  strategic_reserve_selected: boolean
  precursor_realism_score?: number | null
  precursor_realism_assessments?: PrecursorRealismAssessment[]
  precursor_realism_aggregation?: PrecursorRealismAggregation | null
  precursor_realism_band_penalty: number
  effective_structural_score_band: number
}

export interface MultistepRouteStep {
  step_id: string
  depth: number
  product_smiles: string
  precursor_smiles: string[]
  step_cost: number
  step_cost_components: Record<string, number>
  candidate: MultistepRouteCandidate
  product_node_id: string
  precursor_node_ids: string[]
  condition_evidence?: RetrosynthesisConditionEvidence | null
  forward_assessment?: RetrosynthesisForwardAssessment | null
}

export interface MultistepRouteEvidenceSummary {
  compatibility_warning_step_count: number
  strong_compatibility_warning_step_count: number
  realism_assessed_step_count: number
  weakest_precursor_realism_score?: number | null
  strategic_complexity_assessed_step_count: number
  strategic_step_count: number
  tactical_step_count: number
  complexity_increasing_step_count: number
  mean_strategic_complexity_reduction_score?: number | null
  target_to_frontier_complexity_reduction_fraction?: number | null
  condition_assessed_step_count: number
  condition_supported_step_count: number
  condition_direct_step_count: number
  condition_fallback_step_count: number
  condition_insufficient_step_count: number
  condition_coverage_fraction?: number | null
  condition_cost_penalty: number
}

export interface RouteTreeMoleculeNode {
  molecule_node_id: string
  smiles: string
  depth: number
  terminal: boolean
  terminal_evidence: string
  unresolved_reason?: string | null
  reaction?: RouteTreeReactionNode | null
}

export interface RouteTreeReactionNode {
  reaction_node_id: string
  step_id: string
  depth: number
  proposed_reaction_smiles: string
  operator_id: string
  disconnection_site_key: string
  children: RouteTreeMoleculeNode[]
}

export interface CanonicalRouteTree {
  tree_id: string
  target_smiles: string
  root: RouteTreeMoleculeNode
  reaction_count: number
  maximum_depth: number
  fingerprint_tokens: string[]
  schema_version: string
}

export interface MultistepRetrosynthesisRoute {
  route_id: string
  target_smiles: string
  solved: boolean
  route_cost: number
  reaction_count: number
  maximum_depth: number
  steps: MultistepRouteStep[]
  leaves: StartingMaterialAssessment[]
  route_tree: CanonicalRouteTree
  evidence_summary: MultistepRouteEvidenceSummary
  forward_validity_counts?: Record<string, number>
  warnings: string[]
}

export interface MultistepSearchDiagnostics {
  expanded_states: number
  one_step_calls: number
  one_step_cache_hits: number
  generated_candidates: number
  rejected_cycles: number
  rejected_invalid_candidates: number
  duplicate_states: number
  retained_alternative_paths: number
  frontier_states: number
  solved_routes_found: number
  partial_routes_found: number
  stopped_by_expansion_limit: boolean
  proposed_actions: number
  validation_attempts: number
  valid_actions: number
  first_solution_expansion?: number | null
  beam_pruned_states: number
  dead_end_states: number
  expansion_level_calls: Array<[string, number]>
}

export interface MultistepRetrosynthesisResult {
  target_smiles: string
  library_mode: 'full' | 'compact'
  valid: boolean
  error?: string | null
  schema_version: string
  max_depth: 2 | 3
  molecular_weight_threshold: number
  ranking_policy_definition_id: string
  precursor_realism_enabled: boolean
  condition_availability_enabled: boolean
  forward_validation_requested: boolean
  forward_validity_counts: Record<string, number>
  precursor_realism_requested: boolean
  condition_availability_requested: boolean
  strategic_complexity_definition_id: string
  precursor_realism_sources: {
    buyable: boolean
    compound_registry: boolean
    literature: boolean
  }
  route_count: number
  partial_route_count: number
  library_operator_count: number
  library_template_count: number
  search_elapsed_seconds: number
  terminal_stock_source: {
    type: 'supplier_stock_portfolio' | 'literature_molecule_index'
    name: string
  }
  search_budget: Record<string, number>
  route_postprocessing: {
    definition_id: string
    tree_ids: string[]
    distance_matrix: number[][]
  }
  routes: MultistepRetrosynthesisRoute[]
  partial_routes: MultistepRetrosynthesisRoute[]
  diagnostics: MultistepSearchDiagnostics
  warnings: string[]
}
