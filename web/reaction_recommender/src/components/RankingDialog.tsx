import { useState } from 'react'

const LABELS: Record<string, string> = {
  similarity: 'Structural similarity',
  partner_category: 'Reactant category',
  functional_group_tolerance: 'Functional-group tolerance',
  yield: 'Expected yield',
  independent_support: 'Independent support',
  reaction_breadth: 'Reaction breadth',
  dataset_diversity: 'Dataset diversity',
  compatibility: 'Condition compatibility',
  condition_certainty: 'Procedure completeness',
}

interface RankingDialogProps {
  weights: Record<string, number>
  onCancel: () => void
  onSave: (weights: Record<string, number>) => void
}

export function RankingDialog({ weights, onCancel, onSave }: RankingDialogProps) {
  const [values, setValues] = useState<Record<string, number>>(weights)
  const total = Object.values(values).reduce((sum, value) => sum + value, 0)

  return (
    <div className="modal-backdrop" role="presentation">
      <section className="modal-card ranking-modal" role="dialog" aria-modal="true" aria-labelledby="ranking-title">
        <div className="modal-heading">
          <div>
            <span className="eyebrow">TRANSPARENT RERANKING</span>
            <h2 id="ranking-title">Customize priorities</h2>
          </div>
          <button className="icon-button" type="button" onClick={onCancel} aria-label="Close">×</button>
        </div>
        <p className="modal-intro">
          Values are normalized automatically. Chemistry admission and hard compatibility gates remain locked.
        </p>
        <div className="weight-grid">
          {Object.entries(values).map(([name, value]) => (
            <label key={name}>
              <span>{LABELS[name] ?? name.replaceAll('_', ' ')}</span>
              <input
                type="number"
                min="0"
                max="100"
                step="1"
                value={Number((value * 100).toFixed(1))}
                onChange={(event) =>
                  setValues((current) => ({
                    ...current,
                    [name]: Math.max(0, Number(event.target.value) / 100),
                  }))
                }
              />
            </label>
          ))}
        </div>
        {total <= 0 && <div className="alert error">At least one priority must be positive.</div>}
        <div className="modal-actions">
          <button className="button quiet" type="button" onClick={onCancel}>Cancel</button>
          <button className="button primary" type="button" disabled={total <= 0} onClick={() => onSave(values)}>
            Apply priorities
          </button>
        </div>
      </section>
    </div>
  )
}

