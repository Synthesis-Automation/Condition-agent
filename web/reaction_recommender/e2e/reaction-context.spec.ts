import { expect, test } from '@playwright/test'
import { writeFileSync } from 'node:fs'

test.skip(process.env.REACTION_CONTEXT_LIVE !== '1', 'Requires local operator and condition libraries')

for (const mode of ['compact', 'full']) {
  test(`condition context runs graph planning separately (${mode})`, async ({ page }) => {
    test.setTimeout(300_000)
    const errors: string[] = []
    page.on('pageerror', error => errors.push(error.message))
    await page.goto('/')
    await page.getByRole('radio', { name: 'Condition recommendation', exact: true }).check()
    await page.getByRole('combobox', { name: 'Precedent library', exact: true }).selectOption(mode)
    await page.getByText('Advanced options', { exact: true }).click()
    const mapping = page.getByRole('checkbox', { name: /Use RXNMapper/ })
    if (await mapping.isEnabled()) await mapping.uncheck()
    await page.getByLabel('Reaction SMILES', { exact: true }).fill('Ic1ccccc1.N#C[Cu]>>N#Cc1ccccc1')
    const recommendation = page.waitForResponse(r => r.url().endsWith('/recommendations'), { timeout: 120_000 })
    await page.getByRole('button', { name: 'Recommend conditions', exact: true }).click()
    expect((await recommendation).ok()).toBeTruthy()
    const originalTable = await page.locator('.results-layout table').innerText()
    await expect(page.getByRole('button', { name: 'Explore reaction context', exact: true })).toBeHidden()
    await page.getByText('Reaction context: possible products and alternative precursors', { exact: true }).click()
    const context = page.waitForResponse(r => r.url().endsWith('/recommendations/context'), { timeout: 240_000 })
    await page.getByRole('button', { name: 'Explore reaction context', exact: true }).click()
    const response = await context
    expect(response.ok()).toBeTruthy()
    const result = (await response.json()).data
    writeFileSync(`../../results/reaction_context_20260915/${mode}.json`, JSON.stringify(result, null, 2))
    expect(result.advisory_only).toBe(true)
    expect(result.reactant_analysis.status).toBe('complete')
    expect(result.product_analysis.status).toBe('complete')
    expect(result.product_analysis.alternatives.length).toBeGreaterThan(0)
    expect(result.product_analysis.alternatives.length).toBeLessThanOrEqual(4)
    await expect(page.getByRole('heading', { name: 'Possible products from your reactants' })).toBeVisible()
    await expect(page.getByRole('heading', { name: 'Precursor choices for your product' })).toBeVisible()
    expect(await page.locator('.results-layout table').innerText()).toBe(originalTable)
    const choices = page.locator('[aria-label="Product-side analysis"] > details')
    await choices.first().locator('summary').first().click()
    await expect(choices.first().getByAltText('Proposed precursor alternative')).toBeVisible()
    await page.screenshot({ path: `../../results/reaction_context_20260915/${mode}.png`, fullPage: true })
    expect(errors).toEqual([])
  })
}
