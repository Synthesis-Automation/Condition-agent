import { useMemo, useState } from 'react'
import type { Ketcher } from 'ketcher-core'
import { Editor } from 'ketcher-react'
import { StandaloneStructServiceProvider } from 'ketcher-standalone'

const EXAMPLE_REACTION =
  'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1'

interface ReactionEditorProps {
  value: string
  onChange: (value: string) => void
  onError: (message: string) => void
}

function formatError(smiles: string): string | null {
  const sections = smiles.split('>')
  if (sections.length !== 3) return 'Draw a reaction arrow between both sides.'
  if (!sections[0].trim()) return 'The reactant side is empty.'
  if (!sections[2].trim()) return 'The product side is empty.'
  return null
}

export function ReactionEditor({ value, onChange, onError }: ReactionEditorProps) {
  const provider = useMemo(() => new StandaloneStructServiceProvider(), [])
  const [ketcher, setKetcher] = useState<Ketcher | null>(null)
  const [editorStatus, setEditorStatus] = useState('Loading editor…')

  const generate = async () => {
    if (!ketcher) return
    try {
      const smiles = (await ketcher.getSmiles()).trim()
      const error = formatError(smiles)
      onChange(smiles)
      setEditorStatus(error ?? 'Reaction SMILES generated from the drawing.')
      if (error) onError(error)
    } catch (error) {
      onError(error instanceof Error ? error.message : 'Could not export drawing.')
    }
  }

  const load = async (smiles: string) => {
    if (!ketcher) return
    try {
      await ketcher.setMolecule(smiles)
      setEditorStatus('Reaction loaded into the drawing canvas.')
    } catch {
      onError('Ketcher could not load this reaction SMILES.')
    }
  }

  const clear = async () => {
    if (!ketcher) return
    await ketcher.setMolecule('')
    onChange('')
    setEditorStatus('Canvas cleared.')
  }

  return (
    <section className="editor-card" aria-labelledby="editor-title">
      <div className="section-heading">
        <div>
          <span className="step-number">1</span>
          <div>
            <h2 id="editor-title">Define the reaction</h2>
            <p>Draw it or paste an existing reaction SMILES.</p>
          </div>
        </div>
        <div className="button-row">
          <button
            className="button quiet"
            type="button"
            onClick={() => {
              onChange(EXAMPLE_REACTION)
              void load(EXAMPLE_REACTION)
            }}
            disabled={!ketcher}
          >
            Load example
          </button>
          <button className="button quiet" type="button" onClick={clear} disabled={!ketcher}>
            Clear
          </button>
        </div>
      </div>

      <div className="editor-frame">
        {!ketcher && <div className="editor-loading">Loading editor…</div>}
        <Editor
          staticResourcesUrl="/"
          structServiceProvider={provider}
          onInit={(instance) => {
            setKetcher(instance)
            setEditorStatus('Editor ready.')
          }}
          errorHandler={(error) => onError(String(error))}
          disableMacromoleculesEditor
        />
      </div>

      <div className="smiles-workbench">
        <div className="smiles-field">
          <label htmlFor="reaction-smiles">Reaction SMILES</label>
          <textarea
            id="reaction-smiles"
            value={value}
            onChange={(event) => onChange(event.target.value)}
            placeholder="reactants>>products"
            spellCheck={false}
          />
        </div>
        <div className="smiles-actions">
          <span className="inline-status">{editorStatus}</span>
          <div className="button-row">
            <button
              className="button quiet"
              type="button"
              onClick={() => void load(value)}
              disabled={!ketcher || !value.trim()}
            >
              Load SMILES into editor
            </button>
            <button
              className="button secondary"
              type="button"
              onClick={generate}
              disabled={!ketcher}
            >
              Generate from drawing
            </button>
          </div>
        </div>
      </div>
    </section>
  )
}
