import type {
  ApiEnvelope,
  Capabilities,
  CoupledStrategyRetrosynthesisRequest,
  CoupledStrategyRetrosynthesisResult,
  FeatureAnalysisRequest,
  FeatureAnalysisResult,
  FragmentSearchRequest,
  FragmentQueryAlternatives,
  FragmentInvestigation,
  FragmentTransferRequest,
  FragmentSearchResult,
  SynthesisPrecedentResult,
  FragmentSuggestionsResult,
  FragmentGuidedRetroRequest,
  FragmentGuidedRetroResult,
  ForwardSynthesisRequest,
  ForwardSynthesisResult,
  ForwardConditionProfileCatalog,
  MultistepRetrosynthesisRequest,
  MultistepRetrosynthesisResult,
  PrepareReactionResult,
  RankingProfile,
  RecommendationRequest,
  ReactionContextResult,
  RecommendationApiResult,
  RetrosynthesisRequest,
  RetrosynthesisConditionsRequest,
  RetrosynthesisConditionEvidence,
  RetrosynthesisResult,
} from './types'

import type { PlanningRequest, PlanningResponse } from './planning'

const API_ROOT = '/api/v1'

export class ApiError extends Error {
  readonly code: string
  readonly status: number

  constructor(message: string, code: string, status: number) {
    super(message)
    this.name = 'ApiError'
    this.code = code
    this.status = status
  }
}

async function jsonRequest<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(`${API_ROOT}${path}`, {
    ...init,
    headers: {
      'Content-Type': 'application/json',
      ...(init?.headers ?? {}),
    },
  })
  const payload = await response.json()
  if (!response.ok) {
    const detail = payload.detail ?? payload
    const message = Array.isArray(detail)
      ? detail.map(issue => `${(issue.loc ?? []).filter((part: unknown) => part !== 'body').join('.')}: ${issue.msg ?? 'Invalid value'}`).join('; ')
      : typeof detail === 'string' ? detail : detail.message
    throw new ApiError(
      message || `Request failed with status ${response.status}`,
      detail.code ?? 'REQUEST_FAILED',
      response.status,
    )
  }
  return (payload as ApiEnvelope<T>).data
}

export const api = {
  findSynthesisPrecedents: (target: string, signal?: AbortSignal) =>
    jsonRequest<SynthesisPrecedentResult>('/fragments/discover', {
      method: 'POST', body: JSON.stringify({ target_smiles: target }), signal,
    }),
  planningAction: (request: PlanningRequest, signal?: AbortSignal) =>
    jsonRequest<PlanningResponse>('/retrosynthesis/planner', {
      method: 'POST', body: JSON.stringify(request), signal,
    }),
  transferFragmentPrecedents: (request: FragmentTransferRequest, signal?: AbortSignal) =>
    jsonRequest<FragmentGuidedRetroResult>('/retrosynthesis/fragment-transfer', {
      method: 'POST', body: JSON.stringify(request), signal,
    }),
  fragmentGuidedRetro: (request: FragmentGuidedRetroRequest, signal?: AbortSignal) =>
    jsonRequest<FragmentGuidedRetroResult>('/retrosynthesis/fragment-guided', {
      method: 'POST', body: JSON.stringify(request), signal,
    }),
  proposeFragmentQueries: (request: Pick<FragmentSearchRequest, 'query' | 'query_format' | 'topology' | 'target_smiles'> & { aromatic_atom_ids?: number[] }, signal?: AbortSignal) =>
    jsonRequest<FragmentQueryAlternatives>('/fragments/query-alternatives', {
      method: 'POST', body: JSON.stringify(request), signal,
    }),
  investigateFragment: (request: FragmentSearchRequest & { observation_id: string }, signal?: AbortSignal) =>
    jsonRequest<FragmentInvestigation>('/fragments/investigate', {
      method: 'POST', body: JSON.stringify(request), signal,
    }),
  suggestFragments: (target: string, signal?: AbortSignal) =>
    jsonRequest<FragmentSuggestionsResult>('/fragments/suggest', {
      method: 'POST', body: JSON.stringify({ target_smiles: target }), signal,
    }),
  searchFragments: (request: FragmentSearchRequest, signal?: AbortSignal) =>
    jsonRequest<FragmentSearchResult>('/fragments/search', {
      method: 'POST', body: JSON.stringify(request), signal,
    }),
  reactionContext: (reactionSmiles: string, libraryMode: 'full' | 'compact') =>
    jsonRequest<ReactionContextResult>('/recommendations/context', {
      method: 'POST',
      body: JSON.stringify({ reaction_smiles: reactionSmiles, library_mode: libraryMode }),
    }),
  capabilities: () => jsonRequest<Capabilities>('/capabilities'),

  rankingProfiles: async () => {
    const result = await jsonRequest<{ profiles: RankingProfile[] }>(
      '/ranking-profiles',
    )
    return result.profiles
  },

  prepareReaction: (reactionSmiles: string) =>
    jsonRequest<PrepareReactionResult>('/reactions/prepare', {
      method: 'POST',
      body: JSON.stringify({ reaction_smiles: reactionSmiles }),
    }),

  recommend: (request: RecommendationRequest) =>
    jsonRequest<RecommendationApiResult>('/recommendations', {
      method: 'POST',
      body: JSON.stringify(request),
    }),

  analyzeFeatures: (request: FeatureAnalysisRequest) =>
    jsonRequest<FeatureAnalysisResult>('/features/analyze', {
      method: 'POST',
      body: JSON.stringify(request),
    }),

  forwardSynthesize: (request: ForwardSynthesisRequest) =>
    jsonRequest<ForwardSynthesisResult>('/forward-synthesis', {
      method: 'POST',
      body: JSON.stringify(request),
    }),

  forwardConditionProfiles: () =>
    jsonRequest<ForwardConditionProfileCatalog>('/forward-synthesis/condition-profiles'),

  retrosynthesize: (request: RetrosynthesisRequest) =>
    jsonRequest<RetrosynthesisResult>('/retrosynthesis', {
      method: 'POST',
      body: JSON.stringify(request),
    }),

  multistepRetrosynthesize: (request: MultistepRetrosynthesisRequest) =>
    jsonRequest<MultistepRetrosynthesisResult>('/retrosynthesis/routes', {
      method: 'POST',
      body: JSON.stringify(request),
    }),
  coupledStrategyRetrosynthesize: (request: CoupledStrategyRetrosynthesisRequest) =>
    jsonRequest<CoupledStrategyRetrosynthesisResult>('/retrosynthesis/coupled-strategies', {
      method: 'POST',
      body: JSON.stringify(request),
    }),

  retrosynthesisConditions: (request: RetrosynthesisConditionsRequest) =>
    jsonRequest<RetrosynthesisConditionEvidence>('/retrosynthesis/conditions', {
      method: 'POST',
      body: JSON.stringify(request),
    }),

  renderReaction: async (
    reactionSmiles: string,
    width = 760,
    height = 220,
  ) => {
    const response = await fetch(`${API_ROOT}/render/reaction`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        reaction_smiles: reactionSmiles,
        width,
        height,
      }),
    })
    if (!response.ok) {
      throw new ApiError('Reaction rendering failed.', 'RENDER_FAILED', response.status)
    }
    return response.blob()
  },

  renderMolecule: async (
    moleculeSmiles: string,
    width = 760,
    height = 220,
  ) => {
    const response = await fetch(`${API_ROOT}/render/molecule`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        molecule_smiles: moleculeSmiles,
        width,
        height,
      }),
    })
    if (!response.ok) {
      throw new ApiError('Molecule rendering failed.', 'RENDER_FAILED', response.status)
    }
    return response.blob()
  },
}
