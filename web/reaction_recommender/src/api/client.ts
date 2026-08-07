import type {
  ApiEnvelope,
  Capabilities,
  DiscoveryRequest,
  DiscoveryResult,
  PrepareReactionResult,
  RankingProfile,
  RecommendationRequest,
  RecommendationResult,
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
    jsonRequest<RecommendationResult>('/recommendations', {
      method: 'POST',
      body: JSON.stringify(request),
    }),

  discover: (request: DiscoveryRequest) =>
    jsonRequest<DiscoveryResult>('/discovery', {
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
}

