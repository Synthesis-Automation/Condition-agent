import { useEffect, useMemo, useState } from 'react'
import type { FeatureAnalysisResult } from '../api/types'
import { ReactionImage } from './ReactionImage'
import { displayName } from './Results'

function valueText(value: unknown): string {
  if (Array.isArray(value)) return value.map(displayName).join(', ') || 'None'
  if (typeof value === 'number') {
    return Number.isInteger(value) ? String(value) : value.toFixed(3)
  }
  if (typeof value === 'boolean') return value ? 'Yes' : 'No'
  return displayName(value ?? 'Unavailable')
}

function Overview({ values }: { values: Record<string, unknown> }) {
  return (
    <dl className="feature-overview-grid">
      {Object.entries(values)
        .filter(([, value]) => value != null && value !== '')
        .map(([name, value]) => (
          <div key={name}>
            <dt>{displayName(name)}</dt>
            <dd>{valueText(value)}</dd>
          </div>
        ))}
    </dl>
  )
}

function InlineSvg({ source, label }: { source: string; label: string }) {
  const [url, setUrl] = useState('')
  useEffect(() => {
    const next = URL.createObjectURL(
      new Blob([source], { type: 'image/svg+xml' }),
    )
    setUrl(next)
    return () => URL.revokeObjectURL(next)
  }, [source])
  return (
    <div className="reaction-image feature-core-image">
      {url ? <img src={url} alt={label} /> : <span>Rendering…</span>}
    </div>
  )
}

function InlineNotes({ label, values }: { label: string; values: string[] }) {
  const notes = [...new Set(values.filter(Boolean))]
  if (!notes.length) return null
  return (
    <div className="feature-inline-notes">
      <strong>{label}</strong>
      <ul>{notes.map((value) => <li key={value}>{displayName(value)}</li>)}</ul>
    </div>
  )
}

export function FeatureResults({ result }: { result: FeatureAnalysisResult }) {
  const mapping = result.mapping ?? null
  const isReaction = result.input_kind === 'reaction'
  const mappingRows = useMemo(
    () => mapping
      ? Object.fromEntries(
          Object.entries(mapping).filter(([name, value]) => (
            ['status', 'mapper_confidence', 'structure_preserved', 'reactant_mapping_coverage', 'product_mapping_coverage'].includes(name)
            && value != null
          )),
        )
      : null,
    [mapping],
  )
  const overviewValues = useMemo(
    () => isReaction
      ? Object.fromEntries(
          Object.entries(result.overview).filter(([name]) => (
            !['signature_id', 'compatible_named_families'].includes(name)
          )),
        )
      : result.overview,
    [isReaction, result.overview],
  )
  const chemistryOverviewValues = useMemo(
    () => ({
      ...overviewValues,
      ...(mappingRows
        ? {
            atom_mapping: mappingRows.status,
            mapping_confidence: mappingRows.mapper_confidence,
            structure_preserved: mappingRows.structure_preserved,
            reactant_mapping_coverage: mappingRows.reactant_mapping_coverage,
            product_mapping_coverage: mappingRows.product_mapping_coverage,
          }
        : {}),
    }),
    [mappingRows, overviewValues],
  )

  return (
    <section className={`results-card feature-results ${isReaction ? 'reaction-features' : 'molecule-features'}`}>
      <div className="results-summary">
        <div>
          <span className="eyebrow">{isReaction ? 'REACTION ANALYSIS' : 'STRUCTURE ANALYSIS'}</span>
          <h2>{isReaction ? 'Reaction features' : 'Molecular features'}</h2>
        </div>
        <div className="metric-strip">
          <div><strong>{result.valid ? 'Valid' : 'Invalid'}</strong><span>graph</span></div>
          <div><strong>{isReaction ? result.reaction_core?.event_count ?? 0 : result.motifs.length}</strong><span>{isReaction ? 'core events' : 'motifs'}</span></div>
          <div><strong>{result.reactive_sites.length}</strong><span>reactive sites</span></div>
        </div>
      </div>

      {!isReaction && (
        <div className="feature-structure-panel">
          <h3>Structure</h3>
          <ReactionImage
            smiles={result.input_smiles}
            label="Analyzed molecule"
            kind="molecule"
          />
        </div>
      )}

      {!result.valid && <div className="alert error">{result.error ?? 'The structure is invalid.'}</div>}

      <div className="feature-layout">
        <section className="feature-panel">
          <h3>Chemistry overview</h3>
          <Overview values={chemistryOverviewValues} />
          <InlineNotes label="Analysis notes" values={result.warnings} />
        </section>

        {result.reaction_core && (
          <section className="feature-panel feature-panel-wide feature-core-panel">
            <div className="feature-panel-heading">
              <div><h3>Reaction core</h3><p>{displayName(result.reaction_core.evidence_status)} evidence · {(result.reaction_core.confidence * 100).toFixed(0)}% confidence</p></div>
              <span className={`quality-badge ${result.reaction_core.quality.status}`}>{displayName(result.reaction_core.quality.status)}</span>
            </div>
            {result.core_graphic_svg && <InlineSvg source={result.core_graphic_svg} label="Minimized reaction core" />}
            {result.core_graphic_error && <div className="alert error">{result.core_graphic_error}</div>}
            <div className="core-events">
              {result.reaction_core.events.map((event, index) => (
                <div key={event.event_id}>
                  <strong>Event {index + 1}</strong>
                  <ul>{event.edit_tokens.map((token) => <li key={token}>{displayName(token.replaceAll(':', ' · ').replace('>', ' → '))}</li>)}</ul>
                </div>
              ))}
            </div>
            <InlineNotes label="Core review" values={[
              ...result.reaction_core.quality.review_reasons,
              ...result.reaction_core.quality.blocking_reasons,
            ]} />
          </section>
        )}

        {result.partners.length > 0 && (
          <section className="feature-panel">
            <h3>Reaction partners and functional motifs <span className="feature-count">{result.partners.length}</span></h3>
            <div className="feature-card-list">
              {result.partners.map((partner, index) => {
                const motifs = result.motifs.filter((motif) => (
                  (!motif.side || motif.side === 'reactant')
                  && motif.component_index === partner.component_index
                ))
                return (
                  <article key={`${partner.component_index}-${partner.role}-${index}`}>
                    <strong>{partner.label || `Reactant ${partner.component_index + 1}`}</strong>
                    <span>
                      Reactant {partner.component_index + 1}
                      {partner.role ? ` · ${displayName(partner.role)} (${(partner.role_confidence * 100).toFixed(0)}%)` : ''}
                    </span>
                    <div className="partner-motifs">
                      <small>Functional motifs</small>
                      {motifs.length > 0
                        ? <ul>{motifs.map((motif, motifIndex) => <li key={`${motif.motif_id}-${motifIndex}`}>{motif.label || displayName(motif.motif_id)}</li>)}</ul>
                        : <em>None recognized</em>}
                    </div>
                  </article>
                )
              })}
            </div>
          </section>
        )}

        {(!isReaction || result.partners.length === 0) && result.motifs.length > 0 && (
          <section className="feature-panel">
            <h3>Functional motifs <span className="feature-count">{result.motifs.length}</span></h3>
            <div className="feature-card-list">
              {result.motifs.map((motif, index) => (
                <article key={`${motif.side ?? ''}-${motif.component_index}-${motif.motif_id}-${index}`}>
                  <strong>{motif.label || displayName(motif.motif_id)}</strong>
                  <span>{motif.side ? `${displayName(motif.side)} · ` : ''}component {motif.component_index + 1}</span>
                  <small>{motif.tags.map(displayName).join(', ')}</small>
                </article>
              ))}
            </div>
          </section>
        )}

        {result.reactive_sites.length > 0 && (
          <section className="feature-panel feature-panel-wide">
            <h3>Reactive-site hypotheses <span className="feature-count">{result.reactive_sites.length}</span></h3>
            <div className="feature-card-list site-list">
              {result.reactive_sites.map((site, index) => (
                <article key={`${site.side ?? ''}-${site.component_index}-${site.hypothesis_id}-${index}`}>
                  <strong>{site.label}</strong>
                  <span>{displayName(site.site_type)} · {displayName(site.availability)}</span>
                  <small>{site.side ? `${displayName(site.side)} · ` : ''}component {site.component_index + 1} · atoms {site.atom_indices.join(', ')}{site.context_kind ? ` · ${displayName(site.context_kind)}` : ''}</small>
                </article>
              ))}
            </div>
          </section>
        )}
      </div>

      <details className="trace-panel feature-json">
        <summary>Canonical nested analysis JSON</summary>
        <pre>{JSON.stringify(result.analysis, null, 2)}</pre>
      </details>
    </section>
  )
}
