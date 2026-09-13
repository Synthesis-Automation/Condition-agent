import { Component, type ReactNode } from 'react'

interface Props {
  children: ReactNode
  onFailure: () => void
  onClose: () => void
}

/** Keep an editor import or render failure from unmounting the workbench. */
export class DrawingEditorBoundary extends Component<Props, { failed: boolean }> {
  state = { failed: false }

  static getDerivedStateFromError() {
    return { failed: true }
  }

  componentDidCatch() {
    this.props.onFailure()
  }

  render() {
    if (this.state.failed) {
      return <div className="editor-load-error" role="alert">
        <h3>The drawing editor could not load.</h3>
        <p>Close this editor and refresh the page before trying Draw again.
          You can also continue by entering SMILES directly.</p>
        <p>Your existing input is still available when you close this window.</p>
        <button className="button secondary" type="button" onClick={this.props.onClose}>
          Close editor
        </button>
      </div>
    }
    return this.props.children
  }
}
