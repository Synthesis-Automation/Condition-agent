import { useMemo } from 'react'
import type { Ketcher } from 'ketcher-core'
import { Editor } from 'ketcher-react'
import { StandaloneStructServiceProvider } from 'ketcher-standalone'

/** Load the large drawing engine only when a chemist opens the editor. */
export default function KetcherCanvas({ onInit, onError }: {
  onInit: (instance: Ketcher) => void
  onError: (message: string) => void
}) {
  const provider = useMemo(() => new StandaloneStructServiceProvider(), [])
  return <Editor staticResourcesUrl="/" structServiceProvider={provider}
    onInit={onInit} errorHandler={error => onError(String(error))}
    disableMacromoleculesEditor />
}
