import { useMemo, useState } from 'react'
import type { CompletionChoice, CompletionProposal } from '../api/types'

interface CompletionDialogProps {
  proposal: CompletionProposal
  onCancel: () => void
  onConfirm: (choices: CompletionChoice[]) => void
}

export function CompletionDialog({
  proposal,
  onCancel,
  onConfirm,
}: CompletionDialogProps) {
  const defaults = useMemo(
    () =>
      Object.fromEntries(
        proposal.requirements.map((requirement) => [
          requirement.requirement_id,
          requirement.options.find((option) => option.option_kind !== 'unresolved')
            ?.option_id ??
            requirement.options[0]?.option_id ??
            '__custom__',
        ]),
      ),
    [proposal],
  )
  const [selections, setSelections] = useState<Record<string, string>>(defaults)
  const [custom, setCustom] = useState<Record<string, string>>({})
  const [error, setError] = useState('')

  const confirm = () => {
    const choices: CompletionChoice[] = []
    for (const requirement of proposal.requirements) {
      const selected = selections[requirement.requirement_id]
      if (selected === '__custom__') {
        const identifier = (custom[requirement.requirement_id] ?? '').trim()
        if (!identifier) {
          setError('Enter a custom identifier for every custom source choice.')
          return
        }
        choices.push({
          requirement_id: requirement.requirement_id,
          custom_identifier: identifier,
        })
      } else {
        choices.push({
          requirement_id: requirement.requirement_id,
          option_id: selected,
        })
      }
    }
    onConfirm(choices)
  }

  return (
    <div className="modal-backdrop" role="presentation">
      <section className="modal-card" role="dialog" aria-modal="true" aria-labelledby="completion-title">
        <div className="modal-heading">
          <div>
            <span className="eyebrow">CHEMIST CONFIRMATION</span>
            <h2 id="completion-title">Confirm missing reaction source</h2>
          </div>
          <button className="icon-button" type="button" onClick={onCancel} aria-label="Close">×</button>
        </div>
        <p className="modal-intro">
          The product contains a fragment whose source is not explicit in the drawing.
          Choose a curated source or enter the actual material used.
        </p>
        <div className="completion-list">
          {proposal.requirements.map((requirement) => {
            const selected = selections[requirement.requirement_id]
            return (
              <div className="completion-item" key={requirement.requirement_id}>
                <div>
                  <strong>Fragment {requirement.canonical_fragment_smiles}</strong>
                  <span>Attachment element: {requirement.attachment_element}</span>
                </div>
                <select
                  value={selected}
                  onChange={(event) =>
                    setSelections((current) => ({
                      ...current,
                      [requirement.requirement_id]: event.target.value,
                    }))
                  }
                >
                  {requirement.options.map((option) => (
                    <option value={option.option_id} key={option.option_id}>
                      {option.display_name}
                    </option>
                  ))}
                  <option value="__custom__">Enter custom identifier…</option>
                </select>
                {selected === '__custom__' && (
                  <input
                    value={custom[requirement.requirement_id] ?? ''}
                    onChange={(event) =>
                      setCustom((current) => ({
                        ...current,
                        [requirement.requirement_id]: event.target.value,
                      }))
                    }
                    placeholder="CAS, name, registry ID, or SMILES"
                  />
                )}
              </div>
            )
          })}
        </div>
        {error && <div className="alert error">{error}</div>}
        <div className="modal-actions">
          <button className="button quiet" type="button" onClick={onCancel}>Cancel</button>
          <button className="button primary" type="button" onClick={confirm}>Confirm and recommend</button>
        </div>
      </section>
    </div>
  )
}

