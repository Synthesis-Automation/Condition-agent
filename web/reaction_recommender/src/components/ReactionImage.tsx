import { useEffect, useState } from 'react'
import { api } from '../api/client'

interface ReactionImageProps {
  smiles: string
  label: string
  compact?: boolean
  kind?: 'molecule' | 'reaction'
  focusable?: boolean
}

export function ReactionImage({
  smiles,
  label,
  compact = false,
  kind = 'reaction',
  focusable = true,
}: ReactionImageProps) {
  const [source, setSource] = useState<string | null>(null)
  const [failed, setFailed] = useState(false)

  useEffect(() => {
    let active = true
    let objectUrl: string | null = null
    setFailed(false)
    setSource(null)
    if (!smiles) return () => undefined
    const render = kind === 'molecule' ? api.renderMolecule : api.renderReaction
    render(smiles, kind === 'molecule' ? 260 : compact ? 660 : 980, kind === 'molecule' ? 180 : compact ? 260 : 240)
      .then((blob) => {
        if (!active) return
        objectUrl = URL.createObjectURL(blob)
        setSource(objectUrl)
      })
      .catch(() => {
        if (active) setFailed(true)
      })
    return () => {
      active = false
      if (objectUrl) URL.revokeObjectURL(objectUrl)
    }
  }, [smiles, compact, kind])

  if (!smiles) return null
  return (
    <div className={`reaction-image scaled-structure ${compact ? 'compact' : ''}`} tabIndex={focusable ? 0 : undefined} role="region" aria-label={label}>
      {source ? <img src={source} alt={label} /> : <span>{failed ? 'Preview unavailable' : 'Rendering…'}</span>}
    </div>
  )
}
