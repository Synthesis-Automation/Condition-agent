import { expect, test } from '@playwright/test'

test.skip(process.env.RETRO_STRATEGY_SMOKE !== '1', 'Requires the local Compact operator library and Workbench server')

test('single-step results group precursor choices and preserve selection during condition loading', async ({ page }) => {
  test.setTimeout(180_000)
  const errors: string[] = []
  const conditionQueries: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  page.on('request', request => {
    if (request.url().endsWith('/retrosynthesis/conditions')) {
      conditionQueries.push(request.postDataJSON().reaction_smiles)
    }
  })
  await page.goto('/')
  await page.getByRole('radio', { name: 'Single-step retrosynthesis', exact: true }).check()
  await page.getByRole('combobox', { name: 'Operator library' }).selectOption('compact')
  await page.getByLabel('Top strategies', { exact: true }).fill('5')
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('CCN')
  await page.getByText('Advanced options', { exact: true }).click()
  await page.getByRole('checkbox', { name: /Independently replay/ }).uncheck()
  const responsePromise = page.waitForResponse(response => response.url().endsWith('/retrosynthesis'))
  await page.getByRole('button', { name: 'Plan one step', exact: true }).click()
  const response = await responsePromise
  expect(response.ok()).toBeTruthy()
  const result = (await response.json()).data
  expect(result.schema_version).toBe('2.0')
  expect(result.strategy_count).toBe(5)
  expect(new Set(result.strategies.map((strategy: { strategy_id: string }) => strategy.strategy_id)).size).toBe(5)
  expect(result.candidates).toBeUndefined()
  await expect(page.getByRole('heading', { name: '5 validated strategies', exact: true })).toBeVisible()
  const index = result.strategies.findIndex((strategy: { alternate_realizations: unknown[] }) => strategy.alternate_realizations.length > 0)
  expect(index).toBeGreaterThanOrEqual(0)
  await page.locator('.results-layout .table-scroll tbody tr').nth(index).click()
  const choice = page.getByRole('combobox', { name: 'Precursor choice' })
  await choice.selectOption('1')
  const alternate = result.strategies[index].alternate_realizations[0]
  await expect(page.locator('.retrosynthesis-detail .detail-title-row p')).toHaveText(alternate.precursor_smiles)
  await expect(page.getByText(/Done — 5 strategies, conditions loaded/)).toBeVisible({ timeout: 120_000 })
  await expect(choice).toHaveValue('1')
  await expect(page.locator('.retrosynthesis-detail .detail-title-row p')).toHaveText(alternate.precursor_smiles)
  expect(conditionQueries).toContain(alternate.condition_query_reaction_smiles || alternate.proposed_reaction_smiles)
  expect(conditionQueries.length).toBe(result.returned_realization_count)
  await page.getByText(/^Search coverage/).click()
  await expect(page.getByText('Strategies found / requested', { exact: true })).toBeVisible()
  expect(errors).toEqual([])
  await page.screenshot({ path: '../../results/single_step_strategy_stage/workbench.png', fullPage: true })
})
