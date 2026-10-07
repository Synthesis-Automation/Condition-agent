import { defineConfig } from '@playwright/test'

process.env.WEBUI_TEST_PLANNER = '1'

export default defineConfig({
  testDir: './e2e', testMatch: 'interactive-planning.spec.ts',
  outputDir: '../../results/planner-playwright', workers: 1, timeout: 60_000,
  expect: { timeout: 15_000 },
  use: { baseURL: 'http://127.0.0.1:5184', channel: process.env.BROWSER_CHANNEL,
    viewport: { width: 1440, height: 1000 }, screenshot: 'only-on-failure' },
  webServer: { command: 'python -m tests.interactive_planner_server --port 5184',
    cwd: '../..', url: 'http://127.0.0.1:5184/api/v1/health', timeout: 60_000 },
})
