import { defineConfig } from '@playwright/test'

export default defineConfig({
  testDir: './e2e',
  outputDir: '../../results/webui-playwright',
  workers: 1,
  timeout: 60_000,
  expect: { timeout: 20_000 },
  use: {
    baseURL: process.env.WEBUI_TEST_URL ?? 'http://127.0.0.1:8000',
    channel: process.env.BROWSER_CHANNEL,
    viewport: { width: 1440, height: 1000 },
    screenshot: 'only-on-failure',
  },
})
