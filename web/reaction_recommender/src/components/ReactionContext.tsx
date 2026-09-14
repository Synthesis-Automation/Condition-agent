import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import type { ReactionContextResult, ResolvedRecipe } from '../api/types'
import { ReactionImage } from './ReactionImage'

const RELATIONS = {
  supplied_reaction: 'Supplied reaction recovered',
  precursor_alternative: 'Alternative precursors for the shared transformation',
  route_alternative: 'Alternative synthesis route',
  different_product: 'Different product',
  unresolved: 'Relationship to your reaction is unresolved',
}

export function ReactionContext({ reaction, libraryMode, summarizeRecipe }: {
  reaction: string
  libraryMode: 'full' | 'compact'
  summarizeRecipe: (recipe: ResolvedRecipe) => string
}) {
  const [result, setResult] = useState<ReactionContextResult | null>(null)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const generation = useRef(0)
  useEffect(() => {
    generation.current += 1
    setResult(null)
    setError('')
    setBusy(false)
    return () => { generation.current += 1 }
  }, [reaction, libraryMode])

  async function explore() {
    const run = ++generation.current
    setBusy(true)
    setError('')
    try {
      const next = await api.reactionContext(reaction, libraryMode)
      if (run === generation.current) setResult(next)
    } catch (err) {
      if (run === generation.current) setError(err instanceof Error ? err.message : String(err))
    } finally {
      if (run === generation.current) setBusy(false)
    }
  }

  const forward = result?.reactant_analysis
  const retro = result?.product_analysis
  return <details className="reaction-context">
    <summary>Reaction context: possible products and alternative precursors</summary>
    <p>Explore your starting materials and product with a bounded graph search. Generated structures are hypotheses; alternative conditions apply to the displayed alternative reaction.</p>
    {!result && <button type="button" className="button secondary" disabled={busy} onClick={() => void explore()}>{busy ? 'Exploring reaction context…' : 'Explore reaction context'}</button>}
    {error && <p role="alert">{error}</p>}
    {forward && <section aria-label="Reactant-side analysis">
      <h3>Possible products from your reactants</h3>
      {forward.status !== 'complete' ? <p>Analysis {forward.status}: {forward.error}</p> : <>
        <p>Target match: {forward.intended_match?.replaceAll('_', ' ')}. Assessment: {forward.validity?.replaceAll('_', ' ')}.</p>
        <p>These are structural possibilities. Their order does not predict experimental selectivity under a recipe.</p>
        {forward.blind_prediction?.candidates.map(candidate => <details key={candidate.rank}>
          <summary>Proposed product {candidate.rank}: {candidate.product_smiles}</summary>
          <ReactionImage smiles={candidate.reaction_smiles} label="Proposed forward reaction" />
          <p>Source precedent IDs: {candidate.precedent_reaction_ids.join(', ') || 'Unavailable'}</p>
        </details>)}
        {!forward.blind_prediction?.candidates.length && <p>No products found within this search. This does not establish that the requested reaction is impossible.</p>}
        {Boolean(forward.blind_prediction?.competition_groups.length) && <details>
          <summary>Potential competing sites and pathways</summary>
          <ul>{forward.blind_prediction?.competition_groups.map(group => <li key={group.group_key}>
            {group.competition_level === 'site' ? 'Site alternatives' : group.competition_level === 'operator' ? 'Operator alternatives' : 'Product alternatives'}: proposed products {group.candidate_ranks.join(', ')}.
          </li>)}</ul>
          <p>These groups describe generated graph alternatives; their relative experimental rates are unknown.</p>
        </details>}
        <details><summary>Forward search evidence</summary><pre>{JSON.stringify({ checks: forward.checks, warnings: forward.warnings, diagnostics: forward.blind_prediction?.diagnostics, search_warnings: forward.blind_prediction?.warnings }, null, 2)}</pre></details>
      </>}
    </section>}
    {retro && <section aria-label="Product-side analysis">
      <h3>Precursor choices for your product</h3>
      {retro.status !== 'complete' && <p>Analysis {retro.status}: {retro.error}</p>}
      {retro.alternatives.map((item, index) => <details key={`${item.reaction_smiles}-${index}`}>
        <summary>{RELATIONS[item.relation.kind]}{item.relation.inputs_changed ? ' — inputs change' : ''}</summary>
        <ReactionImage smiles={item.reaction_smiles} label="Proposed precursor alternative" />
        <p>{item.relation.differences.join('. ')}</p>
        <p>Source precedent IDs: {item.precedent_reaction_ids.join(', ') || 'Unavailable'}</p>
        {item.relation.kind === 'supplied_reaction' && <p>Use the recommendations above for your supplied reaction. This recovery adds no independent experimental support.</p>}
        {item.conditions && <>
          <h4>Condition evidence for this alternative</h4>
          {item.conditions.recommendations.map(recipe => <div key={recipe.recipe_id}>
            <p>{summarizeRecipe(recipe.resolved_recipe)}</p>
            <p>{recipe.match_label}. Precedents: {recipe.precedent_reaction_ids.join(', ')}</p>
            {recipe.cautions.length > 0 && <details><summary>Cautions for this recipe</summary><ul>{recipe.cautions.map((caution, i) => <li key={i}>{caution}</li>)}</ul></details>}
          </div>)}
          {!item.conditions.recommendations.length && <p>No qualified condition recipe was found.</p>}
          <details><summary>Condition evidence and cautions</summary><pre>{JSON.stringify(item.conditions, null, 2)}</pre></details>
        </>}
        {item.condition_error && <p>Condition evidence unavailable: {item.condition_error}</p>}
        {item.relation.reasons.length > 0 && <p>{item.relation.reasons.join('; ')}</p>}
      </details>)}
      {retro.status === 'complete' && !retro.alternatives.length && <p>No precursor proposals found within this search.</p>}
      {(retro.search_diagnostics?.budget_limited || retro.display_limit_reached) && <p>Search or display limits were reached; further alternatives may exist.</p>}
    </section>}
    {result && <details><summary>Analysis limitations</summary><ul>{result.warnings.map(w => <li key={w}>{w}</li>)}</ul></details>}
  </details>
}
