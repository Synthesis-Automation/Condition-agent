import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig(({ mode }) => ({
  plugins: [react()],
  define: {
    'process.env': {
      NODE_ENV: mode === 'production' ? 'production' : 'development',
    },
  },
  server: {
    host: '127.0.0.1',
    port: 5173,
  },
}))
