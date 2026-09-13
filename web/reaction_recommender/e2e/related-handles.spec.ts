import { expect, test } from '@playwright/test'

// Uses the one-bromide SQLite fixture documented in the validation report.
test.skip(process.env.RELATED_HANDLE_FIXTURE !== '1', 'Requires the controlled Br-only index')

test('simple priorities and scope control actual related-handle retrieval', async ({ page }) => {
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await page.goto('/')
  const scope = page.getByRole('combobox', { name: /Search scope/ })
  const priority = page.getByRole('combobox', { name: /Prioritize/ })
  await expect(scope).toHaveValue('automatic')
  await expect(priority.locator('option')).toHaveText([
    'Balanced', 'Strongest supporting evidence', 'Closest chemistry',
  ])
  await page.getByRole('button', { name: 'Advanced weights', exact: true }).click()
  await expect(page.getByRole('dialog', { name: 'Advanced ranking weights' })).toBeVisible()
  await expect(page.getByText('Evidence for your functional groups', { exact: true })).toBeVisible()
  await page.getByRole('button', { name: 'Cancel', exact: true }).click()
  await page.getByLabel('Reaction SMILES', { exact: true }).fill('Ic1ccccc1.CN>>CNc1ccccc1')
  await page.getByText('Advanced options', { exact: true }).click()
  const mapping = page.getByRole('checkbox', { name: /Use RXNMapper/ })
  if (await mapping.isEnabled()) await mapping.uncheck()
  await priority.selectOption('closest_chemistry')
  const responsePromise = page.waitForResponse(response => response.url().endsWith('/recommendations'))
  await page.getByRole('button', { name: 'Recommend conditions', exact: true }).click()
  const response = await responsePromise
  expect(response.request().postDataJSON().search_scope).toBe('automatic')
  expect(response.request().postDataJSON().ranking_preferences.profile_id).toBe('closest_chemistry')
  const result = (await response.json()).data
  expect(result.recommendations[0].match_level).toBe(3)
  await expect(page.getByText('Level 3 · Related handle', { exact: true })).toBeVisible()
  await expect(page.getByText('Query Ar-I; precedent Ar-Br', { exact: true })).toBeVisible()
  await scope.selectOption('same_handle')
  const strictPromise = page.waitForResponse(response => response.url().endsWith('/recommendations'))
  await page.getByRole('button', { name: 'Recommend conditions', exact: true }).click()
  const strict = await strictPromise
  expect(strict.request().postDataJSON().search_scope).toBe('same_handle')
  expect((await strict.json()).data.recommendations).toHaveLength(0)
  await expect(page.getByText('Level 3 · Related handle', { exact: true })).toBeHidden()
  expect(errors).toEqual([])
})
