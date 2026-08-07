import { useMemo, useState } from 'react'
import type { Ketcher } from 'ketcher-core'
import { Editor } from 'ketcher-react'
import { StandaloneStructServiceProvider } from 'ketcher-standalone'

const EXAMPLE_REACTION = 'CCBr.N>>CCN'

type Message = {
  kind: 'idle' | 'success' | 'error'
  text: string
}

function reactionError(smiles: string): string | null {
  const sections = smiles.split('>')
  if (sections.length !== 3) {
    return 'Draw a reaction arrow so the structure has reactant and product sides.'
  }
  if (!sections[0].trim()) {
    return 'The reactant side is empty.'
  }
  if (!sections[2].trim()) {
    return 'The product side is empty.'
  }
  return null
}

function App() {
  const structServiceProvider = useMemo(
    () => new StandaloneStructServiceProvider(),
    [],
  )
  const [ketcher, setKetcher] = useState<Ketcher | null>(null)
  const [reactionSmiles, setReactionSmiles] = useState('')
  const [message, setMessage] = useState<Message>({
    kind: 'idle',
    text: 'Draw reactants, add a reaction arrow, then draw the product.',
  })

  const generateSmiles = async () => {
    if (!ketcher) return

    try {
      const smiles = (await ketcher.getSmiles()).trim()
      const error = reactionError(smiles)
      setReactionSmiles(smiles)
      setMessage(
        error
          ? { kind: 'error', text: error }
          : { kind: 'success', text: 'Reaction SMILES generated locally.' },
      )
    } catch (error) {
      setMessage({
        kind: 'error',
        text:
          error instanceof Error
            ? error.message
            : 'Ketcher could not export this drawing.',
      })
    }
  }

  const loadExample = async () => {
    if (!ketcher) return
    try {
      await ketcher.setMolecule(EXAMPLE_REACTION)
      setReactionSmiles('')
      setMessage({
        kind: 'idle',
        text: 'Example loaded. Edit it or generate its reaction SMILES.',
      })
    } catch {
      setMessage({ kind: 'error', text: 'The example could not be loaded.' })
    }
  }

  const clearCanvas = async () => {
    if (!ketcher) return
    await ketcher.setMolecule('')
    setReactionSmiles('')
    setMessage({
      kind: 'idle',
      text: 'Canvas cleared. Draw a new reaction when ready.',
    })
  }

  const copySmiles = async () => {
    if (!reactionSmiles) return
    try {
      await navigator.clipboard.writeText(reactionSmiles)
      setMessage({ kind: 'success', text: 'Reaction SMILES copied.' })
    } catch {
      setMessage({
        kind: 'error',
        text: 'Clipboard access was blocked. Select and copy the text manually.',
      })
    }
  }

  return (
    <main className="app-shell">
      <header className="topbar">
        <div>
          <div className="eyebrow">LOCAL CHEMISTRY POC</div>
          <h1>Reaction Sketchpad</h1>
          <p>Draw a transformation and export it as reaction SMILES.</p>
        </div>
        <div className="local-badge" aria-label="Runs locally">
          <span className="status-dot" />
          Browser only
        </div>
      </header>

      <section className="workspace" aria-label="Reaction drawing workspace">
        <div className="editor-panel">
          <div className="panel-heading">
            <div>
              <span className="step-number">1</span>
              <h2>Draw the reaction</h2>
            </div>
            <div className="toolbar-actions">
              <button
                className="button secondary"
                type="button"
                onClick={loadExample}
                disabled={!ketcher}
              >
                Load example
              </button>
              <button
                className="button secondary"
                type="button"
                onClick={clearCanvas}
                disabled={!ketcher}
              >
                Clear
              </button>
            </div>
          </div>

          <div className="editor-frame">
            {!ketcher && <div className="editor-loading">Loading editor…</div>}
            <Editor
              staticResourcesUrl="/"
              structServiceProvider={structServiceProvider}
              onInit={setKetcher}
              errorHandler={(error) =>
                setMessage({ kind: 'error', text: String(error) })
              }
              disableMacromoleculesEditor
            />
          </div>
        </div>

        <aside className="output-panel">
          <div className="panel-heading compact">
            <div>
              <span className="step-number">2</span>
              <h2>Generate SMILES</h2>
            </div>
          </div>

          <p className="guidance">
            Use Ketcher&apos;s reaction arrow between the starting materials and
            product. Multiple molecules on one side are separated with dots.
          </p>

          <button
            className="button primary"
            type="button"
            onClick={generateSmiles}
            disabled={!ketcher}
          >
            Generate reaction SMILES
          </button>

          <label className="output-label" htmlFor="reaction-smiles">
            Reaction SMILES
          </label>
          <textarea
            id="reaction-smiles"
            value={reactionSmiles}
            readOnly
            placeholder="Your reaction SMILES will appear here…"
            spellCheck={false}
          />

          <button
            className="button copy"
            type="button"
            onClick={copySmiles}
            disabled={!reactionSmiles}
          >
            Copy SMILES
          </button>

          <div className={`message ${message.kind}`} role="status">
            <span className="message-icon" aria-hidden="true">
              {message.kind === 'success'
                ? '✓'
                : message.kind === 'error'
                  ? '!'
                  : 'i'}
            </span>
            <span>{message.text}</span>
          </div>

          <div className="format-note">
            <strong>Expected format</strong>
            <code>reactants&gt;&gt;products</code>
          </div>
        </aside>
      </section>

      <footer>
        Structures remain on this device. Server-side RDKit validation can be
        added in the next phase.
      </footer>
    </main>
  )
}

export default App
