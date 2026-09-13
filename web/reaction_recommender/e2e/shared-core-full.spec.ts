import { expect, test } from '@playwright/test'

test.skip(process.env.SHARED_CORE_FULL !== '1', 'Requires completed Full and Compact shared-core artifacts')

test('larger libraries report their actual size and retrieve cyanation precedents', async ({ page }) => {
  test.setTimeout(120_000)
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  const capsResponse = await page.request.get('/api/v1/capabilities')
  expect(capsResponse.ok()).toBeTruthy()
  const caps = (await capsResponse.json()).data
  expect(caps.library_modes.full.row_count).toBeGreaterThan(500_000)
  expect(caps.library_modes.compact.row_count).toBeGreaterThan(90_000)
  expect(caps.library_modes.full.custom_index).toBe(false)
  await page.goto('/')
  const library = page.getByRole('combobox', { name: 'Precedent library' })
  await expect(library.locator('option:checked')).toHaveText(/Full · .* records/)
  await page.getByLabel('Reaction SMILES', { exact: true }).fill('Brc1c2c(cccc2)[nH]n1.N#C>>N#CC1=NNC2=CC=CC=C12')
  await page.getByText('Advanced options', { exact: true }).click()
  const mapping = page.getByRole('checkbox', { name: /Use RXNMapper/ })
  if (await mapping.isEnabled()) await mapping.uncheck()
  const responsePromise = page.waitForResponse(response => response.url().endsWith('/recommendations'), { timeout: 90_000 })
  await page.getByRole('button', { name: 'Recommend conditions', exact: true }).click()
  const response = await responsePromise
  expect(response.ok()).toBeTruthy()
  const result = (await response.json()).data
  expect(result.recommendation_mode).toBe('experimental_shared_core')
  expect(result.recommendations).toHaveLength(5)
  expect(result.recommendations.every((item: { match_namespace: string }) => item.match_namespace === 'shared_reaction_core.v2')).toBe(true)
  await expect(page.getByText(/L[12]: (shared transformation|broader core)/).first()).toBeVisible()
  await page.screenshot({ path: '../../results/shared_core_v2/workbench_full.png', fullPage: true })
  await page.getByLabel('Reaction SMILES', { exact: true }).fill('CC=O>>CCO')
  const reductionPromise = page.waitForResponse(response => response.url().endsWith('/recommendations'), { timeout: 90_000 })
  await page.getByRole('button', { name: 'Recommend conditions', exact: true }).click()
  const reductionResponse = await reductionPromise
  expect(reductionResponse.ok()).toBeTruthy()
  const reduction = (await reductionResponse.json()).data
  expect(reduction.recommendation_mode).toBe('experimental_shared_core')
  expect(reduction.recommendations.length).toBeGreaterThan(0)
  expect(reduction.recommendations.length).toBeLessThanOrEqual(5)
  await page.screenshot({ path: '../../results/shared_core_v2/workbench_reduction.png', fullPage: true })
  await library.selectOption('compact')
  await expect(library.locator('option:checked')).toHaveText(/Compact · .* records/)
  expect(errors).toEqual([])
})
