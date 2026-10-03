import { expect, test } from '@playwright/test'

test.skip(process.env.COMPOSITE_STRATEGY_SMOKE !== '1', 'Requires the local composite catalogue and operator library')

test('composite strategy displays both reactions, target-site evidence and physical cost', async ({ page }) => {
  test.setTimeout(180_000)
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await page.goto('/')
  await page.getByRole('radio', { name: 'Composite two-step strategies', exact: true }).check()
  await page.getByLabel('Top strategies', { exact: true }).fill('1')
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('CCOC(=O)c1cc(N)cc2c1OC(C)(C)C2')
  const responsePromise = page.waitForResponse(response => response.url().endsWith('/retrosynthesis/coupled-strategies'), { timeout: 150_000 })
  await page.getByRole('button', { name: 'Test two-step strategies', exact: true }).click()
  const response = await responsePromise
  expect(response.ok()).toBeTruthy()
  const result = (await response.json()).data
  expect(result.catalog_id).toMatch(/^COMPOSITECATALOG1:/)
  expect(result.actions.length).toBeGreaterThan(0)
  expect(result.actions[0].dependency.admitted).toBe(true)
  expect(result.actions[0].physical_steps).toHaveLength(2)
  expect(result.actions[0].physical_step_cost).toBe(2)
  await expect(page.getByRole('heading', { name: '1 composite strategy', exact: true })).toBeVisible()
  await expect(page.getByText('1 logical action · 2 physical steps', { exact: true })).toBeVisible()
  await expect(page.getByText('Target-site dependency', { exact: true })).toBeVisible()
  await expect(page.getByText('Conditions / one-pot execution', { exact: true })).toBeVisible()
  await expect(page.locator('.coupled-strategy-detail .route-step-card')).toHaveCount(2)
  await page.getByText('Strategy evidence', { exact: true }).click()
  await expect(page.getByText('Step 1 precedents', { exact: true })).toBeVisible()
  await expect(page.getByText('Step 2 precedents', { exact: true })).toBeVisible()
  expect(errors).toEqual([])
  await page.screenshot({ path: '../../results/webui-playwright/composite-strategies.png', fullPage: true })
})
