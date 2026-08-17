import { useEffect, useMemo, useRef, useState } from 'react'
import { ApiError, api } from './api/client'
import type {
  Capabilities,
  CompletionChoice,
  CompletionProposal,
  FeatureAnalysisResult,
  ForwardSynthesisResult,
  ForwardConditionProfileCatalog,
  MultistepRetrosynthesisResult,
  RankingProfile,
  RecommendationApiResult,
  RecommendationResult,
  RetrosynthesisResult,
  WeakLabelRecommendationResult,
} from './api/types'
import { CompletionDialog } from './components/CompletionDialog'
import { FeatureResults } from './components/FeatureResults'
import { RankingDialog } from './components/RankingDialog'
import { ReactionEditor } from './components/ReactionEditor'
import { ForwardSynthesisResults, MultistepRetrosynthesisResults, RecommendationResults, RetrosynthesisResults, WeakLabelRecommendationResults } from './components/Results'

type Mode = 'recommendation' | 'weak_label' | 'forward_synthesis' | 'retrosynthesis' | 'multistep_retrosynthesis' | 'features'
type LibraryMode = 'full' | 'compact'
type WeakLabelOutput = 'weak_label_fallback' | 'weak_label_screening'

const ERROR_MESSAGES: Record<string, string> = {
  INVALID_REACTION: 'The reaction could not be parsed. Check both sides and the reaction arrow.',
  RXNMAPPER_UNAVAILABLE: 'RXNMapper is not installed. Turn off mapping or install the mapping requirements.',
  REACTION_COMPLETION_CHOICES_INCOMPLETE: 'Confirm a source for every missing product fragment.',
  NO_COMPATIBLE_PRECEDENTS: 'No chemically compatible precedents were found.',
  QUERY_NOT_SUPPORTED_BY_WEAK_LABEL_DATASET: 'The graph-derived reaction type or reactive-site pair is not supported by the weak-label dataset.',
  NO_COMPATIBLE_WEAK_LABEL_RECIPE: 'No weak-label recipes passed participant and chemistry compatibility checks.',
  NO_RETROSYNTHESIS_CANDIDATES: 'No structurally validated single-step disconnections were found.',
  NO_MULTISTEP_ROUTES: 'No fully terminated route was found within the selected depth and search limits.',
  INVALID_STARTING_MATERIALS: 'The starting materials could not be parsed. Enter dot-separated molecule SMILES.',
}

function friendlyError(error: unknown): string {
  if (error instanceof ApiError) return ERROR_MESSAGES[error.code] ?? error.message
  if (error instanceof Error) return error.message
  return 'An unexpected error occurred.'
}

function App() {
  const [reactionSmiles, setReactionSmiles] = useState('')
  const [mode, setMode] = useState<Mode>('recommendation')
  const [libraryMode, setLibraryMode] = useState<LibraryMode>('full')
  const [weakLabelOutput, setWeakLabelOutput] = useState<WeakLabelOutput>('weak_label_fallback')
  const [capabilities, setCapabilities] = useState<Capabilities | null>(null)
  const [profiles, setProfiles] = useState<RankingProfile[]>([])
  const [profileId, setProfileId] = useState('default')
  const [customWeights, setCustomWeights] = useState<Record<string, number> | null>(null)
  const [rankingOpen, setRankingOpen] = useState(false)
  const [completionProposal, setCompletionProposal] = useState<CompletionProposal | null>(null)
  const [busy, setBusy] = useState(false)
  const [status, setStatus] = useState('Ready')
  const [error, setError] = useState('')
  const [result, setResult] = useState<RecommendationApiResult | ForwardSynthesisResult | RetrosynthesisResult | MultistepRetrosynthesisResult | FeatureAnalysisResult | null>(null)
  const [topK, setTopK] = useState(5)
  const [minimumPoolSize, setMinimumPoolSize] = useState<number | null>(null)
  const [unrestrictedFallback, setUnrestrictedFallback] = useState(false)
  const [useRxnmapper, setUseRxnmapper] = useState(true)
  const [forceResolvedMapping, setForceResolvedMapping] = useState(false)
  const [includeL0, setIncludeL0] = useState(true)
  const [useRetrosynthesisContext, setUseRetrosynthesisContext] = useState(true)
  const [diversifyRetrosynthesis, setDiversifyRetrosynthesis] = useState(true)
  const [usePrecursorRealism, setUsePrecursorRealism] = useState(true)
  const [useForwardValidation, setUseForwardValidation] = useState(true)
  const [useConditionAvailability, setUseConditionAvailability] = useState(true)
  const [multistepDepth, setMultistepDepth] = useState<2 | 3>(3)
  const [molecularWeightThreshold, setMolecularWeightThreshold] = useState(150)
  const [forwardRecipeText, setForwardRecipeText] = useState('')
  const [includeSelfReactions, setIncludeSelfReactions] = useState(true)
  const [forwardConditionCatalog, setForwardConditionCatalog] = useState<ForwardConditionProfileCatalog | null>(null)
  const [forwardStrategy, setForwardStrategy] = useState('unspecified')
  const [forwardRedoxMode, setForwardRedoxMode] = useState('neutral')
  const [forwardMedium, setForwardMedium] = useState('neutral')
  const [forwardCatalystFamily, setForwardCatalystFamily] = useState('unspecified')
  const retrosynthesisRun = useRef(0)

  useEffect(() => {
    Promise.all([api.capabilities(), api.rankingProfiles(), api.forwardConditionProfiles()])
      .then(([nextCapabilities, nextProfiles, conditionCatalog]) => {
        setCapabilities(nextCapabilities)
        setProfiles(nextProfiles)
        setForwardConditionCatalog(conditionCatalog)
        setUseRxnmapper(nextCapabilities.rxnmapper_available)
      })
      .catch((nextError) => setError(`Local API unavailable: ${friendlyError(nextError)}`))
  }, [])

  useEffect(() => {
    retrosynthesisRun.current += 1
    setBusy(false)
    setResult(null)
    setCompletionProposal(null)
  }, [reactionSmiles])

  const selectedProfile = useMemo(
    () => profiles.find((profile) => profile.profile_id === profileId),
    [profiles, profileId],
  )

  const changeMode = (nextMode: Mode) => {
    retrosynthesisRun.current += 1
    setBusy(false)
    setMode(nextMode)
    setResult(null)
    setError('')
    if ((nextMode === 'forward_synthesis' || nextMode === 'retrosynthesis' || nextMode === 'multistep_retrosynthesis') && capabilities) {
      setLibraryMode(capabilities.default_retrosynthesis_library_mode ?? 'compact')
    }
    if (nextMode === 'multistep_retrosynthesis') {
      setTopK((current) => Math.min(10, current))
    }
  }

  const loadRetrosynthesisConditions = async (
    baseResult: RetrosynthesisResult,
    runId: number,
    activeLibraryMode: LibraryMode,
  ) => {
    for (let index = 0; index < baseResult.candidates.length; index += 1) {
      if (retrosynthesisRun.current !== runId) return
      const candidate = baseResult.candidates[index]
      setStatus(`Loading conditions for hit ${index + 1} of ${baseResult.candidates.length}…`)
      try {
        const conditionEvidence = await api.retrosynthesisConditions({
          reaction_smiles: candidate.condition_query_reaction_smiles
            || candidate.proposed_reaction_smiles,
          library_mode: activeLibraryMode,
          top_k: 3,
          preferred_reaction_ids: candidate.condition_precedent_reaction_ids,
          starting_materials: candidate.precursor_smiles,
          intended_product: candidate.target_smiles,
          operator_hint: candidate.operator_id || null,
          use_forward_validation: useForwardValidation,
          include_l0: includeL0,
        })
        if (retrosynthesisRun.current !== runId) return
        setResult((current) => {
          if (!current || !('candidates' in current)) return current
          return {
            ...current,
            candidates: current.candidates.map((value) =>
              value.template_id === candidate.template_id
              && value.proposed_reaction_smiles === candidate.proposed_reaction_smiles
                ? {
                    ...value,
                    condition_evidence: conditionEvidence,
                    forward_assessment: conditionEvidence.forward_assessment
                      ?? value.forward_assessment,
                  }
                : value,
            ),
          }
        })
      } catch (conditionError) {
        if (retrosynthesisRun.current !== runId) return
        setResult((current) => {
          if (!current || !('candidates' in current)) return current
          return {
            ...current,
            candidates: current.candidates.map((value) =>
              value.template_id === candidate.template_id
              && value.proposed_reaction_smiles === candidate.proposed_reaction_smiles
                ? {
                    ...value,
                    condition_evidence: {
                      ...value.condition_evidence,
                      status: 'insufficient_evidence' as const,
                      recommendation_mode: 'unavailable',
                      warnings: ['CONDITION_RECOMMENDATION_UNAVAILABLE'],
                      error: friendlyError(conditionError),
                    },
                  }
                : value,
            ),
          }
        })
      }
    }
    if (retrosynthesisRun.current === runId) {
      setStatus(`Done — ${baseResult.candidates.length} disconnection(s), conditions loaded`)
    }
  }

  const runRecommendation = async (completionChoices: CompletionChoice[] = []) => {
    setBusy(true)
    setError('')
    setStatus(mode === 'weak_label'
      ? weakLabelOutput === 'weak_label_screening'
        ? 'Matching graph hints and diversifying screening conditions…'
        : 'Matching graph hints against weak-label recipes…'
      : 'Analyzing reaction and ranking compatible recipes…')
    try {
      const next = await api.recommend({
        reaction_smiles: reactionSmiles.trim(),
        recommendation_mode: mode === 'weak_label' ? weakLabelOutput : 'generic',
        library_mode: libraryMode,
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
    if (mode === 'weak_label') {
      await runRecommendation([])
      return
    }
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

  const runFeatureAnalysis = async () => {
    setBusy(true)
    setError('')
    setStatus('Reading the molecular graph and identifying reactive features…')
    try {
      const next = await api.analyzeFeatures({
        input_smiles: reactionSmiles.trim(),
        use_rxnmapper: useRxnmapper,
        force_resolved_mapping: forceResolvedMapping,
      })
      setResult(next)
      setStatus(`Done — ${next.input_kind} features analyzed`)
    } catch (nextError) {
      setError(friendlyError(nextError))
      setStatus('Reaction analysis failed')
    } finally {
      setBusy(false)
    }
  }

  const runRetrosynthesis = async () => {
    const runId = retrosynthesisRun.current + 1
    retrosynthesisRun.current = runId
    const activeLibraryMode = libraryMode
    const startedAt = Date.now()
    setBusy(true)
    setError('')
    setStatus('Applying graph operators and validating generated precursors… 0s')
    const timer = window.setInterval(() => {
      if (retrosynthesisRun.current === runId) {
        const elapsed = Math.floor((Date.now() - startedAt) / 1000)
        setStatus(`Applying graph operators and validating generated precursors… ${elapsed}s`)
      }
    }, 1000)
    try {
      const next = await api.retrosynthesize({
        target_smiles: reactionSmiles.trim(),
        library_mode: activeLibraryMode,
        top_k: topK,
        include_l0: includeL0,
        use_context: useRetrosynthesisContext,
        diversify: diversifyRetrosynthesis,
        use_precursor_realism: usePrecursorRealism,
        use_forward_validation: useForwardValidation,
      })
      if (retrosynthesisRun.current !== runId) return
      setResult(next)
      if (next.valid) {
        setStatus(`Found ${next.candidates.length} disconnection(s); loading conditions…`)
        void loadRetrosynthesisConditions(next, runId, activeLibraryMode)
      } else {
        setStatus('No validated disconnection')
      }
    } catch (nextError) {
      if (retrosynthesisRun.current !== runId) return
      setError(friendlyError(nextError))
      setStatus('Retrosynthesis failed')
    } finally {
      window.clearInterval(timer)
      if (retrosynthesisRun.current === runId) setBusy(false)
    }
  }

  const runForwardSynthesis = async () => {
    const runId = retrosynthesisRun.current + 1
    retrosynthesisRun.current = runId
    let recipe: Record<string, unknown> | null = null
    if (forwardRecipeText.trim()) {
      try {
        const parsed: unknown = JSON.parse(forwardRecipeText)
        if (!parsed || typeof parsed !== 'object' || Array.isArray(parsed)) {
          throw new Error('Recipe must be a JSON object.')
        }
        recipe = parsed as Record<string, unknown>
      } catch (recipeError) {
        setError(recipeError instanceof Error ? `Condition recipe: ${recipeError.message}` : 'Condition recipe JSON is invalid.')
        return
      }
    }
    const startedAt = Date.now()
    setBusy(true)
    setError('')
    setStatus('Applying forward graph operators and validating possible products… 0s')
    const timer = window.setInterval(() => {
      if (retrosynthesisRun.current === runId) {
        const elapsed = Math.floor((Date.now() - startedAt) / 1000)
        setStatus(`Applying forward graph operators and validating possible products… ${elapsed}s`)
      }
    }, 1000)
    try {
      const next = await api.forwardSynthesize({
        starting_materials: reactionSmiles.trim(),
        recipe,
        library_mode: libraryMode,
        top_k: topK,
        include_l0: includeL0,
        include_self_reactions: includeSelfReactions,
        condition_profile: {
          strategy: forwardStrategy,
          redox_mode: forwardRedoxMode,
          medium: forwardMedium,
          catalyst_family: forwardCatalystFamily,
        },
      })
      if (retrosynthesisRun.current !== runId) return
      setResult(next)
      setStatus(next.prediction.candidates.length
        ? `Done — ${next.prediction.candidates.length} possible product(s), ${next.prediction.diagnostics.valid_pathway_count} pathway(s)`
        : 'No supported forward product generated')
    } catch (nextError) {
      if (retrosynthesisRun.current !== runId) return
      setError(friendlyError(nextError))
      setStatus('Forward synthesis failed')
    } finally {
      window.clearInterval(timer)
      if (retrosynthesisRun.current === runId) setBusy(false)
    }
  }

  const runMultistepRetrosynthesis = async () => {
    const runId = retrosynthesisRun.current + 1
    retrosynthesisRun.current = runId
    const startedAt = Date.now()
    setBusy(true)
    setError('')
    setStatus('Searching routes and auditing realism, compatibility, and conditions… 0s')
    const timer = window.setInterval(() => {
      if (retrosynthesisRun.current === runId) {
        const elapsed = Math.floor((Date.now() - startedAt) / 1000)
        setStatus(`Searching routes and auditing realism, compatibility, and conditions… ${elapsed}s`)
      }
    }, 1000)
    try {
      const next = await api.multistepRetrosynthesize({
        target_smiles: reactionSmiles.trim(),
        library_mode: libraryMode,
        top_k_routes: Math.min(10, topK),
        max_depth: multistepDepth,
        molecular_weight_threshold: molecularWeightThreshold,
        include_l0: includeL0,
        use_context: useRetrosynthesisContext,
        diversify: diversifyRetrosynthesis,
        use_precursor_realism: usePrecursorRealism,
        use_condition_availability: useConditionAvailability,
        use_forward_validation: useForwardValidation,
      })
      if (retrosynthesisRun.current !== runId) return
      setResult(next)
      setStatus(next.valid
        ? `Done — ${next.route_count} solved route(s), ${next.partial_route_count} partial`
        : `No solved route — ${next.partial_route_count} partial route(s) retained`)
    } catch (nextError) {
      if (retrosynthesisRun.current !== runId) return
      setError(friendlyError(nextError))
      setStatus('Multi-step retrosynthesis failed')
    } finally {
      window.clearInterval(timer)
      if (retrosynthesisRun.current === runId) setBusy(false)
    }
  }

  const run = () => {
    if (!reactionSmiles.trim()) {
      setError(mode === 'features' ? 'Enter or draw a molecule or reaction before analyzing features.' : mode === 'forward_synthesis' ? 'Draw or paste the starting materials before predicting products.' : mode === 'retrosynthesis' || mode === 'multistep_retrosynthesis' ? 'Draw or paste a target molecule before running retrosynthesis.' : 'Draw or paste a reaction before running the analysis.')
      return
    }
    if ((mode === 'recommendation' || mode === 'weak_label') && !reactionSmiles.includes('>')) {
      setError('Condition recommendation requires reaction SMILES.')
      return
    }
    if ((mode === 'retrosynthesis' || mode === 'multistep_retrosynthesis') && reactionSmiles.includes('>')) {
      setError('Retrosynthesis requires one target molecule, not a reaction.')
      return
    }
    if (mode === 'forward_synthesis' && reactionSmiles.includes('>')) {
      setError('Forward synthesis requires starting materials only, not a reaction arrow.')
      return
    }
    if (mode === 'recommendation' || mode === 'weak_label') void startRecommendation()
    else if (mode === 'forward_synthesis') void runForwardSynthesis()
    else if (mode === 'retrosynthesis') void runRetrosynthesis()
    else if (mode === 'multistep_retrosynthesis') void runMultistepRetrosynthesis()
    else void runFeatureAnalysis()
  }

  const exportResult = () => {
    if (!result) return
    const blob = new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' })
    const url = URL.createObjectURL(blob)
    const link = document.createElement('a')
    link.href = url
    link.download = mode === 'recommendation' || mode === 'weak_label'
      ? mode === 'weak_label' ? 'weak_label_conditions.json' : 'generic_recommendation.json'
      : mode === 'forward_synthesis'
          ? 'forward_synthesis_products.json'
        : mode === 'retrosynthesis'
          ? 'retrosynthesis_candidates.json'
          : mode === 'multistep_retrosynthesis'
            ? 'multistep_retrosynthesis_routes.json'
          : 'structure_features.json'
    link.click()
    URL.revokeObjectURL(url)
  }

  const weakLabelResult = result && 'recommendations' in result && result.recommendation_mode.startsWith('weak_label')
    ? result as WeakLabelRecommendationResult
    : null
  const recommendationResult = result && 'recommendations' in result && !result.recommendation_mode.startsWith('weak_label')
    ? result as RecommendationResult
    : null
  const forwardSynthesisResult = result && 'prediction' in result ? result : null
  const retrosynthesisResult = result && 'candidates' in result ? result : null
  const multistepRetrosynthesisResult = result && 'routes' in result && 'diagnostics' in result ? result : null
  const featureResult = result && 'input_kind' in result ? result : null
  const isRetrosynthesisMode = mode === 'retrosynthesis' || mode === 'multistep_retrosynthesis'
  const isForwardMode = mode === 'forward_synthesis'
  const isOperatorMode = isRetrosynthesisMode || isForwardMode
  const selectedLibraryAvailable = mode === 'weak_label'
    ? capabilities?.weak_label_recommendation ?? false
    : isForwardMode
    ? capabilities?.forward_library_modes?.[libraryMode]?.library_available ?? false
    : isRetrosynthesisMode
    ? (capabilities?.retrosynthesis_library_modes?.[libraryMode]?.library_available ?? false)
      && (mode !== 'multistep_retrosynthesis'
        || capabilities?.stock_portfolio_available === true
        || capabilities?.literature_molecule_index_available === true)
    : capabilities?.library_modes?.[libraryMode]?.index_available
      ?? capabilities?.index_available
      ?? false
  const libraryKind = isOperatorMode ? 'operator library' : 'index'
  const serviceStatus = mode === 'weak_label'
    ? `Weak-label dataset ${selectedLibraryAvailable ? 'ready' : 'unavailable'}`
    : `${libraryMode === 'full' ? 'Full' : 'Compact'} ${libraryKind} ${selectedLibraryAvailable ? 'ready' : 'unavailable'}`

  return (
    <main className="app-shell">
      <header className="topbar">
        <div><h1>Reaction Chemistry Workbench</h1><p>Predict products, audit competing pathways, plan disconnections, and find compatible conditions.</p></div>
        <div className="service-status">
          <span className={`status-dot ${selectedLibraryAvailable ? '' : 'offline'}`} />
          <strong>{capabilities ? serviceStatus : 'Connecting…'}</strong>
        </div>
      </header>

      <div className="analysis-workbench">
        <section className="control-card" aria-labelledby="analysis-title">
        <div className="section-heading">
          <div><span className="step-number">1</span><h2 id="analysis-title">Analysis mode</h2></div>
          {result && <button className="button quiet" type="button" onClick={exportResult}>Export JSON</button>}
        </div>
        <div className="analysis-control-layout">
          <fieldset className="mode-switch" aria-labelledby="analysis-title">
            <label className={mode === 'features' ? 'active' : ''}><input type="radio" name="analysis-mode" value="features" checked={mode === 'features'} onChange={() => changeMode('features')} /><strong>Analyze reactions</strong></label>
            <label className={mode === 'recommendation' ? 'active' : ''}><input type="radio" name="analysis-mode" value="recommendation" checked={mode === 'recommendation'} onChange={() => changeMode('recommendation')} /><strong>Condition recommendation</strong></label>
            <label className={mode === 'weak_label' ? 'active weak-label-mode' : 'weak-label-mode'}><input type="radio" name="analysis-mode" value="weak_label" checked={mode === 'weak_label'} onChange={() => changeMode('weak_label')} /><strong>Weak-label conditions</strong></label>
            <label className={mode === 'forward_synthesis' ? 'active' : ''}><input type="radio" name="analysis-mode" value="forward_synthesis" checked={mode === 'forward_synthesis'} onChange={() => changeMode('forward_synthesis')} /><strong>Forward synthesis</strong></label>
            <label className={mode === 'retrosynthesis' ? 'active' : ''}><input type="radio" name="analysis-mode" value="retrosynthesis" checked={mode === 'retrosynthesis'} onChange={() => changeMode('retrosynthesis')} /><strong>Single-step retrosynthesis</strong></label>
            <label className={mode === 'multistep_retrosynthesis' ? 'active' : ''}><input type="radio" name="analysis-mode" value="multistep_retrosynthesis" checked={mode === 'multistep_retrosynthesis'} onChange={() => changeMode('multistep_retrosynthesis')} /><strong>Multi-step retrosynthesis</strong></label>
          </fieldset>

          <div className="analysis-options">
            <div className={`option-grid ${mode === 'features' ? 'feature-options' : ''}`}>
              {mode !== 'features' && mode !== 'weak_label' && <label className="library-option"><span>{isOperatorMode ? 'Operator library' : 'Precedent library'}</span><select aria-label={isOperatorMode ? 'Operator library' : 'Precedent library'} value={libraryMode} onChange={(event) => { retrosynthesisRun.current += 1; setBusy(false); setLibraryMode(event.target.value as LibraryMode); setResult(null) }}><option value="full">Full — complete</option><option value="compact">Compact — faster</option></select></label>}
              {mode !== 'features' && <label><span>{mode === 'multistep_retrosynthesis' ? 'Top routes' : 'Top results'}</span><input type="number" min="1" max={mode === 'multistep_retrosynthesis' ? 10 : 50} value={topK} onChange={(event) => setTopK(Math.min(mode === 'multistep_retrosynthesis' ? 10 : 50, Math.max(1, Number(event.target.value))))} /></label>}
              {mode === 'recommendation' ? (
                <label className="wide-option"><span>Ranking profile</span><div className="joined-control"><select value={profileId} onChange={(event) => { setProfileId(event.target.value); setCustomWeights(null) }}>{profiles.map((profile) => <option key={profile.profile_id} value={profile.profile_id}>{profile.label}</option>)}</select><button type="button" className="button quiet" onClick={() => setRankingOpen(true)} disabled={!selectedProfile}>Customize</button></div></label>
              ) : mode === 'weak_label' ? (
                <label className="wide-option"><span>Weak-label output</span><select value={weakLabelOutput} onChange={(event) => { setWeakLabelOutput(event.target.value as WeakLabelOutput); setResult(null) }}><option value="weak_label_fallback">Ranked fallback recipes</option><option value="weak_label_screening">Diverse screening array</option></select><small>Graph-derived reaction hints gate an unverified, label-only recipe dataset.</small></label>
              ) : mode === 'forward_synthesis' ? (
                <div className="feature-mode-note"><strong>Product and pathway prediction</strong><span>Forward operators generate products without seeing a target, then validate graph edits and reverse recovery.</span></div>
              ) : mode === 'retrosynthesis' ? (
                <div className="feature-mode-note"><strong>Single-step operator search</strong><span>Structure-derived operators propose and forward-validate precursor sets.</span></div>
              ) : mode === 'multistep_retrosynthesis' ? (
                <div className="multistep-primary-options"><label><span>Maximum depth</span><select value={multistepDepth} onChange={(event) => setMultistepDepth(Number(event.target.value) as 2 | 3)}><option value={2}>2 steps</option><option value={3}>3 steps</option></select></label><label><span>Terminal MW threshold</span><input type="number" min="1" max="500" value={molecularWeightThreshold} onChange={(event) => setMolecularWeightThreshold(Math.min(500, Math.max(1, Number(event.target.value))))} /></label></div>
              ) : (
                <div className="feature-mode-note"><strong>Automatic input detection</strong><span>Molecules show motifs and sites; reactions also include edits and mapping.</span></div>
              )}
            </div>

            {mode === 'forward_synthesis' && forwardConditionCatalog && <div className="forward-condition-panel"><div className="forward-condition-heading"><strong>Condition profile</strong><span>Coarse categories guide ranking without pretending to be a complete recipe.</span></div><div className="forward-condition-grid"><label><span>Reaction strategy</span><select value={forwardStrategy} onChange={(event) => { const value = event.target.value; setForwardStrategy(value); if (value !== 'transition_metal_catalysis') setForwardCatalystFamily('unspecified') }}>{forwardConditionCatalog.strategies.map((option) => <option key={option.id} value={option.id}>{option.label}</option>)}</select></label><label><span>Catalyst family</span><select value={forwardCatalystFamily} disabled={forwardStrategy !== 'transition_metal_catalysis'} onChange={(event) => setForwardCatalystFamily(event.target.value)}>{forwardConditionCatalog.catalyst_families.map((option) => <option key={option.id} value={option.id}>{option.label}</option>)}</select></label><label><span>Redox environment</span><select value={forwardRedoxMode} onChange={(event) => setForwardRedoxMode(event.target.value)}>{forwardConditionCatalog.redox_modes.map((option) => <option key={option.id} value={option.id}>{option.label}</option>)}</select></label><label><span>Medium</span><select value={forwardMedium} onChange={(event) => setForwardMedium(event.target.value)}>{forwardConditionCatalog.media.map((option) => <option key={option.id} value={option.id}>{option.label}</option>)}</select></label></div><p>{forwardConditionCatalog.strategies.find((option) => option.id === forwardStrategy)?.description}</p></div>}

            {mode !== 'weak_label' && <details className="advanced-options"><summary>Advanced options</summary><div>
              {mode === 'recommendation' && <label><span>Minimum precedent pool</span><input type="number" min="1" max="100" placeholder="Definition default" value={minimumPoolSize ?? ''} onChange={(event) => setMinimumPoolSize(event.target.value ? Number(event.target.value) : null)} /></label>}
              {!isOperatorMode && <label className="check-option"><input type="checkbox" checked={useRxnmapper} disabled={!capabilities?.rxnmapper_available || (mode === 'features' && !reactionSmiles.includes('>'))} onChange={(event) => { setUseRxnmapper(event.target.checked); if (!event.target.checked) setForceResolvedMapping(false) }} /><span>Use RXNMapper for unresolved or ambiguous reactions</span></label>}
              {isForwardMode ? (
                <><label className="check-option"><input type="checkbox" checked={includeSelfReactions} onChange={(event) => setIncludeSelfReactions(event.target.checked)} /><span>Include intermolecular self-reactions by allowing multiple equivalents of one input</span></label><label className="check-option"><input type="checkbox" checked={includeL0} onChange={(event) => setIncludeL0(event.target.checked)} /><span>Use broad L0 operators as the final fallback tier</span></label><label className="wide-option"><span>Expert override: canonical recipe JSON <small>(optional)</small></span><textarea value={forwardRecipeText} onChange={(event) => setForwardRecipeText(event.target.value)} placeholder='{"bases": [{"substance_id": "…"}], "temperature_c": 80}' spellCheck={false} /></label></>
              ) : isRetrosynthesisMode ? (
                <><label className="check-option"><input type="checkbox" checked={useRetrosynthesisContext} onChange={(event) => setUseRetrosynthesisContext(event.target.checked)} /><span>Rank with local reaction-context similarity</span></label><label className="check-option"><input type="checkbox" checked={diversifyRetrosynthesis} onChange={(event) => setDiversifyRetrosynthesis(event.target.checked)} /><span>Rank SITE1 → SYN1/REAL1 with completion priors and diversify within score bands</span></label><label className="check-option"><input type="checkbox" checked={usePrecursorRealism} onChange={(event) => setUsePrecursorRealism(event.target.checked)} /><span>De-rank unlikely precursors using stock, registry, literature, and molecular weight</span></label><label className="check-option"><input type="checkbox" checked={useForwardValidation} onChange={(event) => setUseForwardValidation(event.target.checked)} /><span>Independently replay each proposed reaction forward and audit competing products</span></label>{mode === 'multistep_retrosynthesis' && <label className="check-option"><input type="checkbox" checked={useConditionAvailability} onChange={(event) => setUseConditionAvailability(event.target.checked)} /><span>Audit condition availability for each retained reaction and rerank routes</span></label>}<label className="check-option"><input type="checkbox" checked={includeL0} onChange={(event) => setIncludeL0(event.target.checked)} /><span>Use broad L0 operators as the final fallback tier</span></label></>
              ) : mode === 'features' ? (
                <label className="check-option"><input type="checkbox" checked={forceResolvedMapping} disabled={!useRxnmapper || !reactionSmiles.includes('>')} onChange={(event) => setForceResolvedMapping(event.target.checked)} /><span>Map resolved reactions too, for additional atom-mapping evidence</span></label>
              ) : (
                <label className="check-option"><input type="checkbox" checked={unrestrictedFallback} onChange={(event) => setUnrestrictedFallback(event.target.checked)} /><span>Review-core and unrestricted fallback (expert review required)</span></label>
              )}
            </div></details>}
            {error && <div className="alert error" role="alert">{error}</div>}
          </div>
        </div>
        </section>

        <div className="editor-action-layout">
          <ReactionEditor
            value={reactionSmiles}
            onChange={setReactionSmiles}
            onError={setError}
            allowMolecule={mode === 'features' || isRetrosynthesisMode}
            moleculeOnly={isOperatorMode}
            moleculePurpose={isForwardMode ? 'starting_materials' : 'target'}
          />

          <div className="run-control workbench-action-row" aria-label="Analysis action">
            <button className="button primary run-button" type="button" onClick={run} disabled={busy || (mode !== 'features' && !selectedLibraryAvailable) || (mode === 'features' && !capabilities)}>{busy ? 'Working…' : mode === 'recommendation' ? 'Recommend conditions' : mode === 'weak_label' ? weakLabelOutput === 'weak_label_screening' ? 'Build screening array' : 'Find weak-label recipes' : mode === 'forward_synthesis' ? 'Predict products' : mode === 'retrosynthesis' ? 'Plan one step' : mode === 'multistep_retrosynthesis' ? 'Plan multi-step routes' : 'Analyze reactions'}</button>
            <span role="status" aria-live="polite">{status}</span>
          </div>
        </div>
      </div>

      {recommendationResult && <RecommendationResults result={recommendationResult} />}
      {weakLabelResult && <WeakLabelRecommendationResults result={weakLabelResult} />}
      {forwardSynthesisResult && <ForwardSynthesisResults result={forwardSynthesisResult} />}
      {retrosynthesisResult && <RetrosynthesisResults result={retrosynthesisResult} />}
      {multistepRetrosynthesisResult && <MultistepRetrosynthesisResults result={multistepRetrosynthesisResult} />}
      {featureResult && <FeatureResults result={featureResult} />}

      {!result && <section className="empty-state"><span>3</span><div><h2>{mode === 'features' ? 'Inspect graph-derived features' : mode === 'weak_label' ? 'Inspect weak-label condition hypotheses' : mode === 'forward_synthesis' ? 'Inspect possible products and competing pathways' : mode === 'retrosynthesis' ? 'Inspect proposed disconnections' : mode === 'multistep_retrosynthesis' ? 'Inspect solved and partial routes' : 'Inspect ranked evidence'}</h2><p>{mode === 'features' ? 'Structure summaries, motifs, reactive sites, reaction-core events, mapping evidence, and the canonical analysis will appear here.' : mode === 'weak_label' ? 'The graph-derived reaction type, matched reactive sites, label-only support, canonical recipes, and unverified-evidence cautions will appear here.' : mode === 'forward_synthesis' ? 'Validated products, blind ranks, pathway alternatives, operator identities, graph correspondence, and condition compatibility will appear here.' : mode === 'retrosynthesis' ? 'Validated precursor proposals, operator identities, structural scores, support, and ranking traces will appear here.' : mode === 'multistep_retrosynthesis' ? 'Each route shows validated reaction steps, terminal starting materials, supplier-stock provenance, and unresolved stopping reasons.' : 'Recommendations, reaction drawings, conditions, score traces, cautions, and precedent provenance will appear here.'}</p></div></section>}

      <footer>All chemistry and data remain on this machine. Molecular structure is the source of truth.</footer>

      {completionProposal && <CompletionDialog proposal={completionProposal} onCancel={() => setCompletionProposal(null)} onConfirm={(choices) => void runRecommendation(choices)} />}
      {rankingOpen && selectedProfile && <RankingDialog weights={customWeights ?? selectedProfile.weights} onCancel={() => setRankingOpen(false)} onSave={(weights) => { setCustomWeights(weights); setRankingOpen(false) }} />}
    </main>
  )
}

export default App
