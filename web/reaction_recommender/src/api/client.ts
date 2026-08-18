import type {
  ApiEnvelope,
  Capabilities,
  CoupledStrategyRetrosynthesisRequest,
  CoupledStrategyRetrosynthesisResult,
  FeatureAnalysisRequest,
  FeatureAnalysisResult,
  ForwardSynthesisRequest,
  ForwardSynthesisResult,
  ForwardConditionProfileCatalog,
  MultistepRetrosynthesisRequest,
  MultistepRetrosynthesisResult,
  PrepareReactionResult,
  RankingProfile,
  RecommendationRequest,
  RecommendationApiResult,
  RetrosynthesisRequest,
  RetrosynthesisConditionsRequest,
  RetrosynthesisConditionEvidence,
  RetrosynthesisResult,
} from './types'

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
    throw new ApiError(
      detail.message ?? `Request failed with status ${response.status}`,
      detail.code ?? 'REQUEST_FAILED',
      response.status,
    )
  }
  return (payload as ApiEnvelope<T>).data
}

export const api = {
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
