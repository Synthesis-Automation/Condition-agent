import { useEffect, useMemo, useState } from 'react'
import type { Ketcher } from 'ketcher-core'
import { Editor } from 'ketcher-react'
import { StandaloneStructServiceProvider } from 'ketcher-standalone'
import { ReactionImage } from './ReactionImage'

const EXAMPLE_REACTION =
  'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1'
const EXAMPLE_TARGET = 'Cc1ccnc(-c2ccccc2)c1'
const EXAMPLE_STARTING_MATERIALS = 'Brc1ccccc1.OB(O)c1ccccc1'

interface ReactionEditorProps {
  value: string
  onChange: (value: string) => void
  onError: (message: string) => void
  allowMolecule?: boolean
  moleculeOnly?: boolean
  moleculePurpose?: 'target' | 'starting_materials'
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
}

function DrawingDialog({
  value,
  onChange,
  onError,
  allowMolecule = false,
  moleculeOnly = false,
  moleculePurpose = 'target',
  onClose,
}: DrawingDialogProps) {
  const isStartingMaterials = moleculeOnly && moleculePurpose === 'starting_materials'
  const moleculeLabel = isStartingMaterials ? 'starting materials' : 'target molecule'
  const provider = useMemo(() => new StandaloneStructServiceProvider(), [])
  const [ketcher, setKetcher] = useState<Ketcher | null>(null)
  const [draftSmiles, setDraftSmiles] = useState(value)
  const [status, setStatus] = useState('Loading editor…')

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
      onError(`Enter ${moleculeOnly ? moleculeLabel : 'a reaction'} SMILES before loading it.`)
      return
    }
    try {
      await ketcher.setMolecule(smiles.trim())
      setDraftSmiles(smiles.trim())
      setStatus(`${isStartingMaterials ? 'Starting materials' : moleculeOnly ? 'Target' : 'Reaction'} loaded into the drawing canvas.`)
      onError('')
    } catch {
      onError('Ketcher could not load this reaction SMILES.')
    }
  }

  const clear = async () => {
    if (!ketcher) return
    await ketcher.setMolecule('')
    setDraftSmiles('')
    setStatus('Canvas cleared.')
    onError('')
  }

  const finish = async () => {
    if (!ketcher) return
    try {
      const smiles = (await ketcher.getSmiles()).trim()
      const error = inputFormatError(smiles, allowMolecule, moleculeOnly)
      if (error) {
        onError(error)
        setStatus(error)
        return
      }
      onChange(smiles)
      onError('')
      onClose()
    } catch (error) {
      onError(error instanceof Error ? error.message : 'Could not export drawing.')
    }
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
            <span className="eyebrow">{isStartingMaterials ? 'STARTING MATERIALS' : moleculeOnly ? 'TARGET DRAWING' : 'REACTION DRAWING'}</span>
            <h2 id="drawing-title">
              {isStartingMaterials ? 'Draw the starting materials' : moleculeOnly ? 'Draw the target molecule' : allowMolecule ? 'Draw a molecule or reaction' : 'Draw the transformation'}
            </h2>
            <p>
              {isStartingMaterials
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
              onClick={() => void load(isStartingMaterials ? EXAMPLE_STARTING_MATERIALS : moleculeOnly ? EXAMPLE_TARGET : EXAMPLE_REACTION)}
              disabled={!ketcher}
            >
              Load example
            </button>
            <button className="button quiet" type="button" onClick={clear} disabled={!ketcher}>
              Clear
            </button>
            <button className="icon-button" type="button" onClick={onClose} aria-label="Close drawing editor">×</button>
          </div>
        </div>

        <div className="editor-frame drawing-editor-frame">
          {!ketcher && <div className="editor-loading">Loading editor…</div>}
          <Editor
            staticResourcesUrl="/"
            structServiceProvider={provider}
            onInit={(instance) => {
              setKetcher(instance)
              if (
                value.trim()
                && inputFormatError(value.trim(), allowMolecule, moleculeOnly) === null
              ) {
                void instance
                  .setMolecule(value.trim())
                  .then(() => setStatus('Existing reaction loaded.'))
                  .catch(() => onError('Ketcher could not load the existing reaction.'))
              } else {
                setStatus('Editor ready.')
              }
            }}
            errorHandler={(error) => onError(String(error))}
            disableMacromoleculesEditor
          />
        </div>

        <div className="drawing-smiles-row">
          <label htmlFor="drawing-reaction-smiles">
            <span>{isStartingMaterials ? 'Starting-material SMILES' : moleculeOnly ? 'Target molecule SMILES' : 'Reaction SMILES'}</span>
            <textarea
              id="drawing-reaction-smiles"
              value={draftSmiles}
              onChange={(event) => setDraftSmiles(event.target.value)}
              placeholder={isStartingMaterials ? 'starting.materials' : moleculeOnly ? 'target product' : 'reactants>>products'}
              spellCheck={false}
            />
          </label>
          <button
            className="button quiet"
            type="button"
            onClick={() => void load(draftSmiles)}
            disabled={!ketcher || !draftSmiles.trim()}
          >
            Load SMILES
          </button>
        </div>

        <div className="modal-actions drawing-actions">
          <span>{status}</span>
          <div className="button-row">
            <button className="button quiet" type="button" onClick={onClose}>Cancel</button>
            <button className="button primary" type="button" onClick={finish} disabled={!ketcher}>
              Use drawing
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
}: ReactionEditorProps) {
  const isStartingMaterials = moleculeOnly && moleculePurpose === 'starting_materials'
  const [open, setOpen] = useState(false)
  const normalizedValue = value.trim()
  const inputError = normalizedValue
    ? inputFormatError(normalizedValue, allowMolecule, moleculeOnly)
    : null
  const canPreview = Boolean(normalizedValue && inputError === null)
  const detectedKind = normalizedValue && !normalizedValue.includes('>')
    ? 'molecule'
    : 'reaction'

  return (
    <section className="editor-card reaction-paper" aria-labelledby="editor-title">
      <div className="section-heading reaction-paper-heading">
        <div>
          <span className="step-number">2</span>
          <div>
            <h2 id="editor-title">
              {isStartingMaterials ? 'Define the starting materials' : moleculeOnly ? 'Define the target' : allowMolecule ? 'Define the structure' : 'Define the reaction'}
            </h2>
            <p>
              {isStartingMaterials
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
            <button className="button quiet" type="button" onClick={() => onChange('')}>
              Clear
            </button>
          )}
          <button className="button secondary draw-button" type="button" onClick={() => setOpen(true)}>
            {value ? 'Edit drawing' : 'Draw'}
          </button>
        </div>
      </div>

      <div className="reaction-main-input">
        <label htmlFor="main-reaction-smiles">
          <span>{isStartingMaterials ? 'Starting-material SMILES' : moleculeOnly ? 'Target molecule SMILES' : allowMolecule ? 'Molecule or reaction SMILES' : 'Reaction SMILES'}</span>
          <input
            id="main-reaction-smiles"
            type="text"
            value={value}
            onChange={(event) => {
              onChange(event.target.value)
              onError('')
            }}
            placeholder={isStartingMaterials ? 'reactant1.reactant2' : moleculeOnly ? 'target product' : allowMolecule ? 'CCO or reactants>>products' : 'reactants>>products'}
            spellCheck={false}
          />
        </label>
      </div>

      {canPreview ? (
        <div className="reaction-paper-preview">
          <ReactionImage
            smiles={normalizedValue}
            label={`Current ${detectedKind} drawing`}
            kind={detectedKind}
          />
        </div>
      ) : normalizedValue ? (
        <button className="reaction-paper-empty incomplete" type="button" onClick={() => setOpen(true)}>
          <span className="empty-reaction-mark">!</span>
          <strong>{isStartingMaterials ? 'Starting-material input is not valid' : moleculeOnly ? 'Target input is not valid' : 'Reaction SMILES is not complete'}</strong>
          <small>{inputError ?? `Check the ${moleculeOnly ? 'target' : 'reaction'} text, or finish it in the drawing editor.`}</small>
        </button>
      ) : (
        <button className="reaction-paper-empty" type="button" onClick={() => setOpen(true)}>
          <span className="empty-reaction-mark">→</span>
          <strong>{isStartingMaterials ? 'No starting materials drawn yet' : moleculeOnly ? 'No target drawn yet' : 'No reaction drawn yet'}</strong>
          <small>Click to open the {isStartingMaterials ? 'starting-material' : moleculeOnly ? 'target' : 'reaction'} drawing editor</small>
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
          onClose={() => setOpen(false)}
        />
      )}
    </section>
  )
}
