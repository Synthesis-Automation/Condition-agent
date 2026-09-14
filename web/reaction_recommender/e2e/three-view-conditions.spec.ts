import { expect, test } from '@playwright/test'

test.skip(process.env.THREE_VIEW_CONDITIONS !== '1', 'Requires rebuilt Full and Compact projection indexes')

test('condition recommendations combine three graph views through one qualification path', async ({ page }) => {
  test.setTimeout(180_000)
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await page.goto('/')
  await page.getByRole('radio', { name: 'Condition recommendation', exact: true }).check()
  await page.getByText('Advanced options', { exact: true }).click()
  await page.getByRole('combobox', { name: /^Search scope/ }).selectOption('broad')
  const mapping = page.getByRole('checkbox', { name: /Use RXNMapper/ })
  if (await mapping.isEnabled()) await mapping.uncheck()
  await page.getByLabel('Reaction SMILES', { exact: true }).fill('Brc1c2c(cccc2)[nH]n1.N#C>>N#CC1=NNC2=CC=CC=C12')
  for (const mode of ['full', 'compact']) {
    await page.getByRole('combobox', { name: 'Precedent library', exact: true }).selectOption(mode)
    const pending = page.waitForResponse(response => response.url().endsWith('/recommendations'), { timeout: 90_000 })
    await page.getByRole('button', { name: 'Recommend conditions', exact: true }).click()
    const response = await pending
    expect(response.ok()).toBeTruthy()
    const result = (await response.json()).data
    expect(result.valid).toBe(true)
    expect(result.retrieval_definition_version).toContain('shared_core_retrieval.v3@3.0')
    expect(result.candidate_count).toBeLessThanOrEqual(768)
    const observations = result.shared_core_trace.filter((entry: { position?: number }) => entry.position !== undefined)
    expect(new Set(observations.map((entry: { position: number }) => entry.position)).size).toBe(result.candidate_count)
    const channels = new Set(observations.flatMap((entry: { channels: string[] }) => entry.channels))
    for (const channel of ['direct', 'reactant_side', 'product_side']) expect(channels.has(channel)).toBe(true)
    const matchedIndex = result.recommendations.findIndex((item: { candidate_channels: string[] }) => item.candidate_channels.includes('reactant_side'))
    expect(matchedIndex).toBeGreaterThanOrEqual(0)
    await page.locator('.results-layout .table-scroll tbody tr').nth(matchedIndex).click()
    await expect(page.getByText('Also found through reactant-side reaction matching.').first()).toBeVisible()
    await page.screenshot({ path: `../../results/three_view_conditions_20260914/workbench_${mode}.png`, fullPage: true })
  }
  expect(errors).toEqual([])
})
