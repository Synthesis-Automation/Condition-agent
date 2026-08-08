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

function FeatureWarnings({ values }: { values: string[] }) {
  if (!values.length) return null
  return (
    <div className="message-group caution feature-warning">
      <h4>Review notes</h4>
      <ul>{values.map((value) => <li key={value}>{displayName(value)}</li>)}</ul>
    </div>
  )
}

export function FeatureResults({ result }: { result: FeatureAnalysisResult }) {
  const mapping = result.mapping ?? null
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

  return (
    <section className="results-card feature-results">
      <div className="results-summary">
        <div>
          <span className="eyebrow">FEATURE ANALYSIS</span>
          <h2>{result.input_kind === 'reaction' ? 'Reaction features' : 'Molecular features'}</h2>
        </div>
        <div className="metric-strip">
          <div><strong>{displayName(result.input_kind)}</strong><span>detected input</span></div>
          <div><strong>{result.valid ? 'Valid' : 'Invalid'}</strong><span>graph status</span></div>
          <div><strong>{result.schema_version}</strong><span>schema</span></div>
        </div>
      </div>

      <div className="feature-structure-panel">
        <h3>Full structure</h3>
        <ReactionImage
          smiles={result.input_smiles}
          label={`Analyzed ${result.input_kind}`}
          kind={result.input_kind}
        />
      </div>

      {!result.valid && <div className="alert error">{result.error ?? 'The structure is invalid.'}</div>}
      <FeatureWarnings values={result.warnings} />

      <div className="feature-layout">
        <section className="feature-panel">
          <h3>Chemistry overview</h3>
          <Overview values={result.overview} />
        </section>

        {mappingRows && Object.keys(mappingRows).length > 0 && (
          <section className="feature-panel">
            <h3>Atom-mapping evidence</h3>
            <Overview values={mappingRows} />
          </section>
        )}

        {result.reaction_core && (
          <section className="feature-panel feature-panel-wide">
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
            <FeatureWarnings values={[
              ...result.reaction_core.quality.review_reasons,
              ...result.reaction_core.quality.blocking_reasons,
            ]} />
          </section>
        )}

        {result.partners.length > 0 && (
          <section className="feature-panel">
            <h3>Reaction partners</h3>
            <div className="feature-card-list">
              {result.partners.map((partner, index) => (
                <article key={`${partner.component_index}-${partner.role}-${index}`}>
                  <strong>{partner.label || `Component ${partner.component_index + 1}`}</strong>
                  <span>{displayName(partner.role ?? 'unassigned')} · {(partner.role_confidence * 100).toFixed(0)}%</span>
                  {partner.anchor_contexts.length > 0 && <small>{partner.anchor_contexts.map(displayName).join(', ')}</small>}
                </article>
              ))}
            </div>
          </section>
        )}

        {result.motifs.length > 0 && (
          <section className="feature-panel">
            <h3>Functional motifs</h3>
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
            <h3>Reactive-site hypotheses</h3>
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
