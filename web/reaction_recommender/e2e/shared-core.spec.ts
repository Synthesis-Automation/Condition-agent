import { expect, test } from '@playwright/test'

test.skip(process.env.SHARED_CORE_FIXTURE !== '1', 'Requires the shared-core development artifact')

test('workbench shows graph levels and product-side evidence without changing the query', async ({ page }) => {
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await page.goto('/')
  const query = 'C#N.IC1=NNC2=C1C=CC=C2>>N#CC1=NNC2=CC=CC=C12'
  await page.getByLabel('Reaction SMILES', { exact: true }).fill(query)
  await page.getByText('Advanced options', { exact: true }).click()
  const mapping = page.getByRole('checkbox', { name: /Use RXNMapper/ })
  if (await mapping.isEnabled()) await mapping.uncheck()
  await page.getByLabel('Independent evidence target').fill('4')
  const responsePromise = page.waitForResponse(response => response.url().endsWith('/recommendations'))
  await page.getByRole('button', { name: 'Recommend conditions', exact: true }).click()
  const response = await responsePromise
  expect(response.ok()).toBeTruthy()
  const result = (await response.json()).data
  expect(result.recommendation_mode).toBe('experimental_shared_core')
  expect(result.query_reaction_smiles).toBe(query)
  expect(result.recommendations).toHaveLength(3)
  expect(result.recommendations[0].match_namespace).toBe('shared_reaction_core.v1')
  await expect(page.getByText('L1: shared transformation', { exact: true }).first()).toBeVisible()
  await expect(page.getByText('Also found through product-side reaction matching.', { exact: true })).toBeVisible()
  await expect(page.getByText('Level 3 · L1: shared transformation', { exact: true })).toHaveCount(0)
  await page.screenshot({ path: '../../results/shared_core_validation/workbench.png', fullPage: true })
  await page.getByRole('combobox', { name: /Search scope/ }).selectOption('same_handle')
  const strictPromise = page.waitForResponse(response => response.url().endsWith('/recommendations'))
  await page.getByRole('button', { name: 'Recommend conditions', exact: true }).click()
  const strict = (await (await strictPromise).json()).data
  expect(strict.recommendations).toHaveLength(0)
  expect(errors).toEqual([])
})
