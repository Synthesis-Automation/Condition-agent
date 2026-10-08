import { lazy, Suspense, useEffect, useId, useState } from 'react'
import type { Ketcher } from 'ketcher-core'
import { ReactionImage } from './ReactionImage'
import { DrawingEditorBoundary } from './DrawingEditorBoundary'

const KetcherCanvas = lazy(() => import('./KetcherCanvas'))

const EXAMPLE_REACTION =
  'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1'
const EXAMPLE_TARGET = 'Cc1ccnc(-c2ccccc2)c1'
const EXAMPLE_STARTING_MATERIALS = 'Brc1ccccc1.OB(O)c1ccccc1'

async function loadDrawing(instance: Ketcher, smiles: string): Promise<void> {
  // Ketcher 3.17 resolves setMolecule even on parser failure. Its event bus is
  // the failure signal; without this check the old canvas can be saved instead.
  let failed = false
  const onFailure = () => { failed = true }
  instance.eventBus.once('FAILURE', onFailure)
  try {
    await instance.setMolecule(smiles)
    if (failed) throw new Error('The drawing editor could not load these SMILES. Check the text or load a valid structure before saving.')
  } finally {
    instance.eventBus.removeListener('FAILURE', onFailure)
  }
}

interface ReactionEditorProps {
  value: string
  onChange: (value: string) => void
  onError: (message: string) => void
  allowMolecule?: boolean
  moleculeOnly?: boolean
  moleculePurpose?: 'target' | 'starting_materials' | 'fragment'
  queryFormat?: 'smiles' | 'smarts'
  disabled?: boolean
}

function formatError(smiles: string): string | null {
  const sections = smiles.split('>')
  if (sections.length !== 3) return 'Draw a reaction arrow between both sides.'
  if (!sections[0].trim()) return 'The reactant side is empty.'
  if (!sections[2].trim()) return 'The product side is empty.'
  return null
}

function inputFormatError(
  smiles: string,
  allowMolecule: boolean,
  moleculeOnly: boolean,
): string | null {
  if (moleculeOnly) {
    return smiles.includes('>')
      ? 'Enter molecules without a reaction arrow.'
      : null
  }
  if (allowMolecule && !smiles.includes('>')) return null
  return formatError(smiles)
}

interface DrawingDialogProps extends ReactionEditorProps {
  onClose: () => void
  fragment?: boolean
}

export function DrawingDialog({
  value,
  onChange,
  onError,
  allowMolecule = false,
  moleculeOnly = false,
  moleculePurpose = 'target',
  onClose,
  fragment = false,
}: DrawingDialogProps) {
  const isStartingMaterials = moleculeOnly && moleculePurpose === 'starting_materials'
  const moleculeLabel = fragment ? 'core fragment' : isStartingMaterials ? 'starting materials' : 'target molecule'
  const [ketcher, setKetcher] = useState<Ketcher | null>(null)
  const [draftSmiles, setDraftSmiles] = useState(value)
  const [draftEdited, setDraftEdited] = useState(false)
  const [editorError, setEditorError] = useState('')
  const [status, setStatus] = useState('Loading editor…')
  const [loading, setLoading] = useState(true)

  const reportError = (message: string) => {
    setEditorError(message)
    if (message) setStatus(message)
    onError(message)
  }

  useEffect(() => {
    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') onClose()
    }
    document.body.style.overflow = 'hidden'
    window.addEventListener('keydown', handleKeyDown)
    return () => {
      document.body.style.overflow = ''
      window.removeEventListener('keydown', handleKeyDown)
    }
  }, [onClose])

  const load = async (smiles: string) => {
    if (!ketcher) return
    if (!smiles.trim()) {
      reportError(`Enter ${moleculeOnly ? moleculeLabel : 'a reaction'} SMILES before loading it.`)
      return
    }
    setLoading(true)
    reportError('')
    setStatus('Loading drawing…')
    try {
      await loadDrawing(ketcher, smiles.trim())
      setDraftSmiles(smiles.trim())
      setDraftEdited(false)
      setStatus(`${fragment ? 'Fragment' : isStartingMaterials ? 'Starting materials' : moleculeOnly ? 'Target' : 'Reaction'} loaded into the drawing canvas.`)
      onError('')
    } catch (error) {
      reportError(error instanceof Error ? error.message : 'Ketcher could not load these SMILES.')
    } finally {
      setLoading(false)
    }
  }

  const clear = async () => {
    if (!ketcher) return
    setLoading(true)
    reportError('')
    try {
      await loadDrawing(ketcher, '')
      setDraftSmiles('')
      setDraftEdited(false)
      setStatus('Canvas cleared.')
    } catch (error) {
      reportError(error instanceof Error ? error.message : 'Could not clear the canvas.')
    } finally { setLoading(false) }
  }

  const finish = async () => {
    if (!ketcher || loading) return
    setLoading(true)
    reportError('')
    setStatus('Saving structure…')
    try {
      // Text edits are a separate draft until explicitly loaded. Never export
      // the stale canvas when the user finishes with pending SMILES changes.
      if (draftEdited) {
        const draft = draftSmiles.trim()
        const error = !draft ? 'Enter SMILES or load a drawing before saving.'
          : inputFormatError(draft, allowMolecule, moleculeOnly)
        if (error) { reportError(error); return }
        await loadDrawing(ketcher, draft)
        setDraftEdited(false)
      }
      const smiles = (await ketcher.getSmiles()).trim()
      const error = inputFormatError(smiles, allowMolecule, moleculeOnly)
        || (fragment && (!smiles || smiles.includes('.')) ? 'Draw one connected core fragment.' : null)
      if (error) {
        reportError(error)
        return
      }
      onChange(smiles)
      onError('')
      onClose()
    } catch (error) {
      reportError(error instanceof Error ? error.message : 'Could not export drawing.')
    } finally { setLoading(false) }
  }

  return (
    <div className="modal-backdrop drawing-backdrop" role="presentation">
      <section
        className="modal-card drawing-modal"
        role="dialog"
        aria-modal="true"
        aria-labelledby="drawing-title"
      >
        <div className="modal-heading drawing-heading">
          <div>
            <span className="eyebrow">{fragment ? 'FRAGMENT DRAWING' : isStartingMaterials ? 'STARTING MATERIALS' : moleculeOnly ? 'TARGET DRAWING' : 'REACTION DRAWING'}</span>
            <h2 id="drawing-title">
              {fragment ? 'Draw the core fragment' : isStartingMaterials ? 'Draw the starting materials' : moleculeOnly ? 'Draw the target molecule' : allowMolecule ? 'Draw a molecule or reaction' : 'Draw the transformation'}
            </h2>
            <p>
              {fragment
                ? 'Draw one connected core to search in reaction products. Keep the ring system and important substituents; do not add a reaction arrow.'
                : isStartingMaterials
                ? 'Draw every starting material as a separate molecular component. Do not add a reaction arrow or product.'
                : moleculeOnly
                ? 'Draw the product structure for single-step precursor generation.'
                : allowMolecule
                ? 'Draw one molecule, or place reactants and products around a reaction arrow.'
                : 'Place reactants and products on opposite sides of a reaction arrow.'}
            </p>
          </div>
          <div className="button-row">
            <button
              className="button quiet"
              type="button"
              onClick={() => void load(fragment ? 'c1ccc2c(c1)COc1ccccc1-2' : isStartingMaterials ? EXAMPLE_STARTING_MATERIALS : moleculeOnly ? EXAMPLE_TARGET : EXAMPLE_REACTION)}
              disabled={!ketcher || loading}
            >
              Load example
            </button>
            <button className="button quiet" type="button" onClick={clear} disabled={!ketcher || loading}>
              Clear
            </button>
            <button className="icon-button" type="button" onClick={onClose} aria-label="Close drawing editor">×</button>
          </div>
        </div>

        <div className="editor-frame drawing-editor-frame">
          <DrawingEditorBoundary onClose={onClose} onFailure={() => {
            setKetcher(null)
            setLoading(false)
            setStatus('Editor unavailable. Close this window to keep working with SMILES.')
          }}>
          {!ketcher && <div className="editor-loading">Loading editor…</div>}
          <Suspense fallback={null}><KetcherCanvas
            onInit={(instance) => {
              setKetcher(instance)
              if (
                value.trim()
                && inputFormatError(value.trim(), allowMolecule, moleculeOnly) === null
              ) {
                void loadDrawing(instance, value.trim())
                  .then(() => setStatus(fragment ? 'Existing fragment loaded.' : 'Existing reaction loaded.'))
                  .catch(() => reportError('Ketcher could not load the existing structure. Edit the SMILES or draw a replacement.'))
                  .finally(() => setLoading(false))
              } else {
                setStatus('Editor ready.')
                setLoading(false)
              }
            }}
            onError={reportError}
          /></Suspense>
          </DrawingEditorBoundary>
        </div>

        <div className="drawing-smiles-row">
          <label htmlFor="drawing-reaction-smiles">
            <span>{fragment ? 'Fragment SMILES' : isStartingMaterials ? 'Starting-material SMILES' : moleculeOnly ? 'Target molecule SMILES' : 'Reaction SMILES'}</span>
            <textarea
              id="drawing-reaction-smiles"
              aria-label={fragment ? 'Fragment SMILES' : isStartingMaterials ? 'Starting-material SMILES' : moleculeOnly ? 'Target molecule SMILES' : 'Reaction SMILES'}
              value={draftSmiles}
              disabled={loading}
              onChange={(event) => { setDraftSmiles(event.target.value); setDraftEdited(true); reportError('') }}
              placeholder={fragment ? 'connected core fragment' : isStartingMaterials ? 'starting.materials' : moleculeOnly ? 'target product' : 'reactants>>products'}
              spellCheck={false}
            />
          </label>
          <button
            className="button quiet"
            type="button"
            onClick={() => void load(draftSmiles)}
            disabled={!ketcher || loading || !draftSmiles.trim()}
          >
            Load SMILES
          </button>
        </div>

        {draftEdited && <p className="drawing-draft-notice">Edited SMILES have not been loaded into the canvas. “Use edited SMILES” will load and save this text.</p>}
        {editorError && <div className="alert error drawing-error" role="alert">{editorError}</div>}

        <div className="modal-actions drawing-actions">
          <span>{status}</span>
          <div className="button-row">
            <button className="button quiet" type="button" onClick={onClose}>Cancel</button>
            <button className="button primary" type="button" onClick={finish} disabled={!ketcher || loading}>
              {loading && ketcher ? 'Working…' : draftEdited ? 'Use edited SMILES' : 'Use drawing'}
            </button>
          </div>
        </div>
      </section>
    </div>
  )
}

export function ReactionEditor({
  value,
  onChange,
  onError,
  allowMolecule = false,
  moleculeOnly = false,
  moleculePurpose = 'target',
  queryFormat = 'smiles',
  disabled = false,
}: ReactionEditorProps) {
  const editorId = useId()
  const isStartingMaterials = moleculeOnly && moleculePurpose === 'starting_materials'
  const isFragment = moleculeOnly && moleculePurpose === 'fragment'
  const textOnly = queryFormat === 'smarts'
  const [open, setOpen] = useState(false)
  const normalizedValue = value.trim()
  const inputError = normalizedValue
    ? inputFormatError(normalizedValue, allowMolecule, moleculeOnly)
    : null
  const canPreview = Boolean(normalizedValue && inputError === null && !textOnly)
  const detectedKind = normalizedValue && !normalizedValue.includes('>')
    ? 'molecule'
    : 'reaction'

  return (
    <section className="editor-card reaction-paper" aria-labelledby={`${editorId}-title`}>
      <div className="section-heading reaction-paper-heading">
        <div>
          <span className="step-number">2</span>
          <div>
            <h2 id={`${editorId}-title`}>
              {isFragment ? 'Define the core fragment' : isStartingMaterials ? 'Define the starting materials' : moleculeOnly ? 'Define the target' : allowMolecule ? 'Define the structure' : 'Define the reaction'}
            </h2>
            <p>
              {isFragment
                ? textOnly ? 'Enter a connected SMARTS query; query features are preserved as text.' : 'Enter or draw one connected core to find its synthesis precedents.'
                : isStartingMaterials
                ? <>Enter or draw dot-separated starting-material SMILES without a reaction arrow.</>
                : moleculeOnly
                ? 'Enter or draw one product molecule for precursor generation.'
                : allowMolecule
                ? 'Enter or draw a molecule or complete reaction SMILES.'
                : <>Enter or draw complete <code>reactants&gt;&gt;products</code> SMILES.</>}
            </p>
          </div>
        </div>
        <div className="button-row">
          {value && (
            <button className="button quiet" type="button" disabled={disabled} onClick={() => onChange('')}>
              Clear
            </button>
          )}
          <button className="button secondary draw-button" type="button" disabled={disabled || textOnly} onClick={() => setOpen(true)}>
            {value ? 'Edit drawing' : 'Draw'}
          </button>
        </div>
      </div>

      <div className="reaction-main-input">
        <label htmlFor={`${editorId}-smiles`}>
          <span>{isFragment ? 'Core fragment' : isStartingMaterials ? 'Starting-material SMILES' : moleculeOnly ? 'Target molecule SMILES' : allowMolecule ? 'Molecule or reaction SMILES' : 'Reaction SMILES'}</span>
          <input
            id={`${editorId}-smiles`}
            type="text"
            value={value}
            disabled={disabled}
            maxLength={isFragment ? 2000 : undefined}
            onChange={(event) => {
              onChange(event.target.value)
              onError('')
            }}
            placeholder={isFragment ? textOnly ? 'Core SMARTS' : 'Core SMILES' : isStartingMaterials ? 'reactant1.reactant2' : moleculeOnly ? 'target product' : allowMolecule ? 'CCO or reactants>>products' : 'reactants>>products'}
            spellCheck={false}
          />
        </label>
      </div>

      {textOnly ? (
        <div className="reaction-paper-empty">
          <strong>SMARTS query</strong>
          <small>Edit SMARTS queries as text. Drawing is available in SMILES mode.</small>
        </div>
      ) : canPreview ? (
        <div className="reaction-paper-preview">
          <ReactionImage
            smiles={normalizedValue}
            label={`Current ${isFragment ? 'fragment' : detectedKind} drawing`}
            kind={detectedKind}
          />
        </div>
      ) : normalizedValue ? (
        <button className="reaction-paper-empty incomplete" type="button" disabled={disabled} onClick={() => setOpen(true)}>
          <span className="empty-reaction-mark">!</span>
          <strong>{isFragment ? 'Fragment input is not valid' : isStartingMaterials ? 'Starting-material input is not valid' : moleculeOnly ? 'Target input is not valid' : 'Reaction SMILES is not complete'}</strong>
          <small>{inputError ?? `Check the ${moleculeOnly ? 'target' : 'reaction'} text, or finish it in the drawing editor.`}</small>
        </button>
      ) : (
        <button className="reaction-paper-empty" type="button" disabled={disabled} onClick={() => setOpen(true)}>
          <span className="empty-reaction-mark">{moleculeOnly ? '⌬' : '→'}</span>
          <strong>{isFragment ? 'No fragment drawn yet' : isStartingMaterials ? 'No starting materials drawn yet' : moleculeOnly ? 'No target drawn yet' : 'No reaction drawn yet'}</strong>
          <small>Click to open the {isFragment ? 'fragment' : isStartingMaterials ? 'starting-material' : moleculeOnly ? 'target' : 'reaction'} drawing editor</small>
        </button>
      )}

      {open && (
        <DrawingDialog
          value={value}
          onChange={onChange}
          onError={onError}
          allowMolecule={allowMolecule}
          moleculeOnly={moleculeOnly}
          moleculePurpose={moleculePurpose}
          fragment={isFragment}
          onClose={() => setOpen(false)}
        />
      )}
    </section>
  )
}
