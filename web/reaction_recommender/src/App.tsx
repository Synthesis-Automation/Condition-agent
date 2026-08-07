import { useEffect, useMemo, useState } from 'react'
import { ApiError, api } from './api/client'
import type {
  Capabilities,
  CompletionChoice,
  CompletionProposal,
  DiscoveryResult,
  RankingProfile,
  RecommendationResult,
} from './api/types'
import { CompletionDialog } from './components/CompletionDialog'
import { RankingDialog } from './components/RankingDialog'
import { ReactionEditor } from './components/ReactionEditor'
import { DiscoveryResults, RecommendationResults } from './components/Results'

type Mode = 'recommendation' | 'discovery'

const ERROR_MESSAGES: Record<string, string> = {
  INVALID_REACTION: 'The reaction could not be parsed. Check both sides and the reaction arrow.',
  RXNMAPPER_UNAVAILABLE: 'RXNMapper is not installed. Turn off mapping or install the mapping requirements.',
  REACTION_COMPLETION_CHOICES_INCOMPLETE: 'Confirm a source for every missing product fragment.',
  NO_COMPATIBLE_PRECEDENTS: 'No chemically compatible precedents were found.',
}

function friendlyError(error: unknown): string {
  if (error instanceof ApiError) return ERROR_MESSAGES[error.code] ?? error.message
  if (error instanceof Error) return error.message
  return 'An unexpected error occurred.'
}

function App() {
  const [reactionSmiles, setReactionSmiles] = useState('')
  const [mode, setMode] = useState<Mode>('recommendation')
  const [capabilities, setCapabilities] = useState<Capabilities | null>(null)
  const [profiles, setProfiles] = useState<RankingProfile[]>([])
  const [profileId, setProfileId] = useState('default')
  const [customWeights, setCustomWeights] = useState<Record<string, number> | null>(null)
  const [rankingOpen, setRankingOpen] = useState(false)
  const [completionProposal, setCompletionProposal] = useState<CompletionProposal | null>(null)
  const [busy, setBusy] = useState(false)
  const [status, setStatus] = useState('Ready')
  const [error, setError] = useState('')
  const [result, setResult] = useState<RecommendationResult | DiscoveryResult | null>(null)
  const [topK, setTopK] = useState(5)
  const [minimumPoolSize, setMinimumPoolSize] = useState<number | null>(null)
  const [unrestrictedFallback, setUnrestrictedFallback] = useState(false)
  const [useRxnmapper, setUseRxnmapper] = useState(true)
  const [discoveryView, setDiscoveryView] = useState('closest_chemistry')
  const [includeLowYield, setIncludeLowYield] = useState(true)
  const [includeUnreported, setIncludeUnreported] = useState(true)

  useEffect(() => {
    Promise.all([api.capabilities(), api.rankingProfiles()])
      .then(([nextCapabilities, nextProfiles]) => {
        setCapabilities(nextCapabilities)
        setProfiles(nextProfiles)
        setUseRxnmapper(nextCapabilities.rxnmapper_available)
      })
      .catch((nextError) => setError(`Local API unavailable: ${friendlyError(nextError)}`))
  }, [])

  useEffect(() => {
    setResult(null)
    setCompletionProposal(null)
  }, [reactionSmiles])

  const selectedProfile = useMemo(
    () => profiles.find((profile) => profile.profile_id === profileId),
    [profiles, profileId],
  )

  const runRecommendation = async (completionChoices: CompletionChoice[] = []) => {
    setBusy(true)
    setError('')
    setStatus('Analyzing reaction and ranking compatible recipes…')
    try {
      const next = await api.recommend({
        reaction_smiles: reactionSmiles.trim(),
        top_k: topK,
        minimum_pool_size: minimumPoolSize,
        unrestricted_fallback: unrestrictedFallback,
        use_rxnmapper: useRxnmapper,
        ranking_preferences: {
          profile_id: profileId,
          weights: customWeights ?? {},
        },
        completion_choices: completionChoices,
      })
      setResult(next)
      setStatus(next.valid ? `Done — ${next.recommendations.length} recipe(s)` : 'No recommendation')
      setCompletionProposal(null)
    } catch (nextError) {
      setError(friendlyError(nextError))
      setStatus('Recommendation failed')
    } finally {
      setBusy(false)
    }
  }

  const startRecommendation = async () => {
    setBusy(true)
    setError('')
    setStatus('Validating reaction and checking fragment sources…')
    try {
      const prepared = await api.prepareReaction(reactionSmiles.trim())
      if (prepared.completion_proposal.requirements.length > 0) {
        setCompletionProposal(prepared.completion_proposal)
        setStatus('Chemist confirmation required')
        setBusy(false)
        return
      }
      await runRecommendation([])
    } catch (nextError) {
      setError(friendlyError(nextError))
      setStatus('Reaction validation failed')
      setBusy(false)
    }
  }

  const runDiscovery = async () => {
    setBusy(true)
    setError('')
    setStatus('Comparing bond edits and local environments…')
    try {
      const next = await api.discover({
        reaction_smiles: reactionSmiles.trim(),
        top_k: topK,
        view: discoveryView,
        include_low_yield: includeLowYield,
        include_unreported_outcomes: includeUnreported,
        use_rxnmapper: useRxnmapper,
        include_review: unrestrictedFallback,
      })
      setResult(next)
      setStatus(next.valid ? `Done — ${next.hits.length} analogue(s)` : 'No structural analogue')
    } catch (nextError) {
      setError(friendlyError(nextError))
      setStatus('Discovery failed')
    } finally {
      setBusy(false)
    }
  }

  const run = () => {
    if (!reactionSmiles.trim()) {
      setError('Draw or paste a reaction before running the analysis.')
      return
    }
    if (mode === 'recommendation') void startRecommendation()
    else void runDiscovery()
  }

  const exportResult = () => {
    if (!result) return
    const blob = new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' })
    const url = URL.createObjectURL(blob)
    const link = document.createElement('a')
    link.href = url
    link.download = mode === 'recommendation' ? 'generic_recommendation.json' : 'reaction_discovery.json'
    link.click()
    URL.revokeObjectURL(url)
  }

  const recommendationResult = result && 'recommendations' in result ? result : null
  const discoveryResult = result && 'hits' in result ? result : null

  return (
    <main className="app-shell">
      <header className="topbar">
        <div><span className="eyebrow">LOCAL · STRUCTURE-FIRST</span><h1>Reaction Condition Recommender</h1><p>Draw a transformation, retrieve compatible precedent chemistry, and inspect every ranking decision.</p></div>
        <div className="service-status">
          <span className={`status-dot ${capabilities?.index_available ? '' : 'offline'}`} />
          <div><strong>{capabilities?.index_available ? 'Index ready' : capabilities ? 'Index unavailable' : 'Connecting…'}</strong><span>{capabilities?.index_name ?? 'Local recommendation service'}</span></div>
        </div>
      </header>

      <ReactionEditor value={reactionSmiles} onChange={setReactionSmiles} onError={setError} />

      <section className="control-card" aria-labelledby="analysis-title">
        <div className="section-heading">
          <div><span className="step-number">2</span><div><h2 id="analysis-title">Choose the analysis</h2><p>Hard chemistry filters remain active in every mode.</p></div></div>
          {result && <button className="button quiet" type="button" onClick={exportResult}>Export JSON</button>}
        </div>
        <div className="mode-switch" role="tablist" aria-label="Analysis mode">
          <button type="button" className={mode === 'recommendation' ? 'active' : ''} onClick={() => { setMode('recommendation'); setResult(null) }}><strong>Condition recommendation</strong><span>Rank compatible recipes</span></button>
          <button type="button" className={mode === 'discovery' ? 'active' : ''} onClick={() => { setMode('discovery'); setResult(null) }}><strong>Reaction discovery</strong><span>Explore structural precedents</span></button>
        </div>

        <div className="option-grid">
          <label><span>Top results</span><input type="number" min="1" max="50" value={topK} onChange={(event) => setTopK(Math.min(50, Math.max(1, Number(event.target.value))))} /></label>
          {mode === 'recommendation' ? (
            <label className="wide-option"><span>Ranking profile</span><div className="joined-control"><select value={profileId} onChange={(event) => { setProfileId(event.target.value); setCustomWeights(null) }}>{profiles.map((profile) => <option key={profile.profile_id} value={profile.profile_id}>{profile.label}</option>)}</select><button type="button" className="button quiet" onClick={() => setRankingOpen(true)} disabled={!selectedProfile}>Customize</button></div><small>{selectedProfile?.description}</small></label>
          ) : (
            <label className="wide-option"><span>Discovery view</span><select value={discoveryView} onChange={(event) => setDiscoveryView(event.target.value)}><option value="closest_chemistry">Closest chemistry</option><option value="diverse_strategies">Diverse strategies</option><option value="successful_precedents">Successful precedents</option><option value="failure_informed">Failure-informed</option></select></label>
          )}
          <div className="run-control"><button className="button primary run-button" type="button" onClick={run} disabled={busy || !capabilities?.index_available}>{busy ? 'Working…' : mode === 'recommendation' ? 'Recommend conditions' : 'Discover precedents'}</button><span>{status}</span></div>
        </div>

        {mode === 'discovery' && <div className="inline-checks"><label><input type="checkbox" checked={includeLowYield} onChange={(event) => setIncludeLowYield(event.target.checked)} /> Include low-yield precedents</label><label><input type="checkbox" checked={includeUnreported} onChange={(event) => setIncludeUnreported(event.target.checked)} /> Include unreported outcomes</label></div>}

        <details className="advanced-options"><summary>Advanced options</summary><div>
          {mode === 'recommendation' && <label><span>Minimum precedent pool</span><input type="number" min="1" max="100" placeholder="Definition default" value={minimumPoolSize ?? ''} onChange={(event) => setMinimumPoolSize(event.target.value ? Number(event.target.value) : null)} /></label>}
          <label className="check-option"><input type="checkbox" checked={useRxnmapper} disabled={!capabilities?.rxnmapper_available} onChange={(event) => setUseRxnmapper(event.target.checked)} /><span>Use RXNMapper for unresolved or ambiguous queries</span></label>
          <label className="check-option"><input type="checkbox" checked={unrestrictedFallback} onChange={(event) => setUnrestrictedFallback(event.target.checked)} /><span>Review-core and unrestricted fallback (expert review required)</span></label>
        </div></details>
        {error && <div className="alert error" role="alert">{error}</div>}
      </section>

      {recommendationResult && <RecommendationResults result={recommendationResult} />}
      {discoveryResult && <DiscoveryResults result={discoveryResult} />}

      {!result && <section className="empty-state"><span>3</span><div><h2>Inspect ranked evidence</h2><p>Recommendations, discovery hits, reaction drawings, conditions, score traces, cautions, and precedent provenance will appear here.</p></div></section>}

      <footer>All chemistry and data remain on this machine. Molecular structure is the source of truth.</footer>

      {completionProposal && <CompletionDialog proposal={completionProposal} onCancel={() => setCompletionProposal(null)} onConfirm={(choices) => void runRecommendation(choices)} />}
      {rankingOpen && selectedProfile && <RankingDialog weights={customWeights ?? selectedProfile.weights} onCancel={() => setRankingOpen(false)} onSave={(weights) => { setCustomWeights(weights); setRankingOpen(false) }} />}
    </main>
  )
}

export default App

