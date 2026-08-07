import { createRoot } from 'react-dom/client'
import './ketcherRuntime'
import 'ketcher-react/dist/index.css'
import './styles.css'
import App from './App'

createRoot(document.getElementById('root')!).render(
  <App />,
)
