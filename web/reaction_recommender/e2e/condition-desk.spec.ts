import { expect, test, type Page } from '@playwright/test'
import { readFile } from 'node:fs/promises'

const examples: { label: string; reaction_smiles: string }[] = JSON.parse(
  await readFile(new URL('../src/conditions/examples.json', import.meta.url), 'utf-8'),
)

test.skip(process.env.CONDITION_DESK_TEST !== '1', 'Requires the conditions-only frontend')

const reaction = 'Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1'

for (const [index, example] of examples.entries()) {
  test(`example picker loads ${example.label} and avoids an immediate repeat`, async ({ page }) => {
    await page.addInitScript(value => { Math.random = () => value }, (index + 0.1) / examples.length)
    await page.route('**/api/v1/capabilities', route => route.fulfill({ json: { data: { recommendation: true } } }))
    await page.goto('/')
    await expect(page).toHaveTitle('ZBS chemistry recommender')
    await expect(page.getByRole('heading', { name: 'ZBS chemistry recommender', exact: true })).toBeVisible()
    await expect(page.getByText('REACTION CHEMISTRY · CONDITION DESK', { exact: true })).toHaveCount(0)
    await page.getByRole('button', { name: 'Try an example', exact: true }).click()
    await expect(page.getByLabel('Reaction SMILES', { exact: true })).toHaveValue(example.reaction_smiles)
    await expect(page.getByText(example.label, { exact: true })).toBeVisible()
    await page.getByRole('button', { name: 'Try an example', exact: true }).click()
    await expect(page.getByLabel('Reaction SMILES', { exact: true })).not.toHaveValue(example.reaction_smiles)
  })
}

function fixture() {
  const evidence = (source: 'generic' | 'weak_label', rank: number) => ({
    source, recommendation_mode: source === 'generic' ? 'experimental_shared_core' : 'weak_label_screening',
    warnings: [], recommendation: {
      rank, support: 4, reference_support: 2, match_label: 'L0: detailed local reaction',
      match_namespace: 'shared_reaction_core.v2', evidence_relation: 'analogue_evidence',
      source_reaction_types: ['Suzuki-Miyaura'], source_row_numbers: [24, 31],
      explanation: [], match_details: ['Supplied molecular inputs differ from the precedent'],
    },
  })
  const options = Array.from({ length: 7 }, (_, i) => ({
    option_id: `option-${i + 1}`, rank: i + 1,
    evidence_kind: i === 6 ? 'weak_label' : 'structure_review',
    evidence_label: i === 6 ? 'Screening suggestion' : 'Broader structural match',
    resolved_recipe: {
      recipe_id: `RCR2:fixture-${i + 1}`,
      catalysts: [{ canonical_name: 'Pd(PPh3)4', identity_status: 'resolved' }],
      bases: [{ canonical_name: 'Potassium carbonate', identity_status: 'resolved' }],
      solvents: [{ canonical_name: i === 6 ? 'Ethanol' : '1,4-Dioxane', identity_status: 'resolved' }],
      temperature_c: 80, time_h: 12,
    },
    synthesis_protocol: { missing_required_fields: ['reaction_scale', 'addition_order'] },
    cautions: [], evidence: [evidence(i === 6 ? 'weak_label' : 'generic', i === 6 ? 1 : i + 1)],
  }))
  options[4].evidence.push(evidence('weak_label', 2))
  return {
    valid: true, query_reaction_smiles: reaction, shortlist_size: 3, schema_version: '1.0', warnings: [],
    recommendations: options,
    sources: [
      { source: 'generic', status: 'ok', message: '', result: { reaction_label: { text: 'HetAr–Cl + HetAr–B(OH)₂ → HetAr–HetAr' } } },
      { source: 'weak_label', status: 'ok', message: '', result: {} },
    ],
    automation_exports: Object.fromEntries(options.map(option => [option.option_id, { option_id: option.option_id, execution_ready: false }])),
  }
}

async function prepare(page: Page, data = fixture()) {
  await page.route('**/api/v1/capabilities', route => route.fulfill({ json: { data: { recommendation: true, weak_label_recommendation: true } } }))
  await page.route('**/api/v1/reactions/prepare', route => route.fulfill({ json: { data: { valid: true, completion_proposal: { requirements: [] } } } }))
  await page.route('**/api/v1/conditions/recommend', route => route.fulfill({ json: { data } }))
  await page.goto('/')
  await page.getByLabel('Reaction SMILES', { exact: true }).fill(reaction)
}

test('tabs preserve the viewed row and merged recipes export once', async ({ page }) => {
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await prepare(page)
  await page.getByLabel('Max recipes per tab').selectOption('20')
  await page.getByLabel('Reaction matches').selectOption('broad')
  const request = page.waitForRequest(request => request.url().endsWith('/conditions/recommend'))
  await page.getByRole('button', { name: 'Find conditions' }).click()
  expect((await request).postDataJSON()).toMatchObject({ top_k: 20, search_scope: 'broad' })
  const literature = page.getByRole('tab', { name: 'Literature-based 6', exact: true })
  const screening = page.getByRole('tab', { name: 'Screening suggestions 2', exact: true })
  const details = page.getByRole('region', { name: 'Selected suggestion details', exact: true })
  await expect(literature).toHaveAttribute('aria-selected', 'true')
  await expect(page.getByRole('tabpanel')).toHaveCount(1)
  await expect(page.getByRole('table').getByRole('row')).toHaveCount(7)
  await expect(details.getByRole('heading', { name: 'Precedent 1', exact: true })).toBeVisible()
  await page.getByRole('table').getByRole('row').nth(2).click()
  await expect(details.getByRole('heading', { name: 'Precedent 2', exact: true })).toBeVisible()
  await expect(page.getByRole('button', { name: 'View precedent 2', exact: true })).toHaveAttribute('aria-current', 'true')
  await expect(page.getByRole('region', { name: 'Reaction chemistry' })).toContainText('HetAr–Cl')
  await screening.click()
  await expect(page.getByRole('table').getByRole('row')).toHaveCount(3)
  await expect(details).toContainText('Ethanol')
  await expect(page.getByRole('tabpanel')).toContainText('Source reaction structures are not verified')
  await page.getByLabel('Select screen 2', { exact: true }).check()
  await page.getByRole('button', { name: 'View screen 2', exact: true }).click()
  await expect(details).toContainText('1,4-Dioxane')
  await expect(details).not.toContainText('Ethanol')
  await literature.click()
  await expect(details.getByRole('heading', { name: 'Precedent 2', exact: true })).toBeVisible()
  await expect(page.getByLabel('Select precedent 5', { exact: true })).toBeChecked()
  await expect(page.getByText('1 condition selected', { exact: true })).toBeVisible()
  const downloaded = page.waitForEvent('download')
  await page.getByRole('button', { name: 'Export screening set', exact: true }).click()
  const download = await downloaded
  const payload = JSON.parse(await readFile((await download.path())!, 'utf-8'))
  expect(payload.handoffs).toHaveLength(1)
  expect(payload.handoffs[0].option_id).toBe('option-5')
  expect(payload.execution_ready).toBe(false)
  const recipeDownloaded = page.waitForEvent('download')
  await details.getByRole('button', { name: 'Export recipe', exact: true }).click()
  const recipe = JSON.parse(await readFile((await (await recipeDownloaded).path())!, 'utf-8'))
  expect(recipe.option_id).toBe('option-2')
  await page.screenshot({ path: '../../results/zbs_tabs/desktop.png', fullPage: true })
  await screening.click()
  await expect(details.getByRole('heading', { name: 'Screen 2', exact: true })).toBeVisible()
  expect(errors).toEqual([])
})

test('screening remains usable when the structural library is unavailable', async ({ page }) => {
  const data = fixture()
  data.sources[0].status = 'unavailable'
  data.sources[0].message = 'The reaction library is unavailable.'
  data.recommendations = data.recommendations.filter(option => option.evidence_kind === 'weak_label')
  await prepare(page, data)
  await page.getByRole('button', { name: 'Find conditions' }).click()
  await expect(page.getByText('The reaction library is unavailable.', { exact: false })).toBeVisible()
  await expect(page.getByRole('tab', { name: 'Screening suggestions 1', exact: true })).toHaveAttribute('aria-selected', 'true')
  await expect(page.getByRole('table').getByRole('row')).toHaveCount(2)
  await page.getByLabel('Select screen 1', { exact: true }).check()
  await expect(page.getByRole('button', { name: 'Export screening set', exact: true })).toBeEnabled()
  await page.getByRole('tab', { name: 'Literature-based 0', exact: true }).click()
  await expect(page.getByRole('tabpanel')).toContainText('No qualifying literature precedents')
  await expect(page.getByRole('region', { name: 'Selected suggestion details' })).toHaveCount(0)
})

test('mobile layout fits and reports an unavailable screening source', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 })
  const data = fixture()
  data.sources[1].status = 'unavailable'
  data.sources[1].message = 'The screening library is unavailable.'
  data.recommendations = data.recommendations.filter(option => option.evidence_kind !== 'weak_label')
  data.recommendations.forEach(option => { option.evidence = option.evidence.filter(item => item.source === 'generic') })
  await prepare(page, data)
  await page.getByRole('button', { name: 'Find conditions' }).click()
  await expect(page.getByText('The screening library is unavailable.', { exact: false })).toBeVisible()
  await page.getByRole('button', { name: 'View precedent 2', exact: true }).click()
  await expect(page.getByRole('heading', { name: 'Precedent 2', exact: true })).toBeInViewport()
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true)
  await page.screenshot({ path: '../../results/zbs_tabs/mobile.png', fullPage: true })
  await page.getByRole('tab', { name: 'Screening suggestions 0', exact: true }).click()
  await expect(page.getByRole('tabpanel')).toContainText('No supported screening suggestions')
  await expect(page.getByRole('table')).toHaveCount(0)
  await expect(page.getByRole('region', { name: 'Selected suggestion details' })).toHaveCount(0)
})

test('tabs and row details work with the keyboard and reset for a new search', async ({ page }) => {
  await prepare(page)
  await page.getByRole('button', { name: 'Find conditions' }).click()
  const literature = page.getByRole('tab', { name: 'Literature-based 6', exact: true })
  const screening = page.getByRole('tab', { name: 'Screening suggestions 2', exact: true })
  await literature.focus()
  await page.keyboard.press('ArrowRight')
  await expect(screening).toBeFocused()
  await expect(screening).toHaveAttribute('aria-selected', 'true')
  await page.getByRole('button', { name: 'View screen 2', exact: true }).focus()
  await page.keyboard.press('Enter')
  await expect(page.getByRole('heading', { name: 'Screen 2', exact: true })).toBeVisible()
  await page.getByLabel('Select screen 2', { exact: true }).check()
  await screening.focus()
  await page.keyboard.press('Home')
  await expect(literature).toBeFocused()
  await page.keyboard.press('End')
  await expect(screening).toBeFocused()
  await page.getByRole('button', { name: 'Find conditions' }).click()
  await expect(literature).toHaveAttribute('aria-selected', 'true')
  await expect(page.getByRole('heading', { name: 'Precedent 1', exact: true })).toBeVisible()
  await expect(page.getByLabel('Select precedent 5', { exact: true })).not.toBeChecked()
  await expect(page.getByRole('button', { name: 'Export screening set', exact: true })).toHaveCount(0)
})

test('real libraries show structural chemistry and screening recipes', async ({ page }) => {
  test.skip(process.env.CONDITION_DESK_FULL !== '1', 'Requires the local Full and weak-label libraries')
  test.setTimeout(120_000)
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await page.goto('/')
  await page.getByLabel('Reaction SMILES', { exact: true }).fill(reaction)
  const response = page.waitForResponse(response => response.url().endsWith('/conditions/recommend'))
  await page.getByRole('button', { name: 'Find conditions' }).click()
  const result = (await (await response).json()).data
  expect(result.sources.map((source: { status: string }) => source.status)).toEqual(['ok', 'ok'])
  expect(result.sources[1].result.recommendation_mode).toBe('weak_label_screening')
  await expect(page.getByRole('region', { name: 'Reaction chemistry' })).toContainText('HetAr–Cl')
  const literatureCount = result.recommendations.filter((option: { evidence: { source: string }[] }) => option.evidence.some(item => item.source === 'generic')).length
  await expect(page.getByRole('table').getByRole('row')).toHaveCount(literatureCount + 1)
  const details = page.getByRole('region', { name: 'Selected suggestion details' })
  const tableBounds = await page.getByRole('table').boundingBox()
  const detailBounds = await details.boundingBox()
  expect(detailBounds!.x).toBeGreaterThan(tableBounds!.x + tableBounds!.width)
  await page.screenshot({ path: '../../results/zbs_tabs/real-literature.png', fullPage: true })
  await page.getByRole('tab', { name: 'Screening suggestions 10', exact: true }).click()
  await expect(page.getByRole('table').getByRole('row')).toHaveCount(11)
  await page.getByRole('button', { name: 'View screen 2', exact: true }).click()
  await expect(details.getByRole('heading', { name: 'Screen 2', exact: true })).toBeVisible()
  await page.screenshot({ path: '../../results/zbs_tabs/real-screening.png', fullPage: true })
  await page.setViewportSize({ width: 390, height: 844 })
  await page.getByRole('button', { name: 'View screen 3', exact: true }).click()
  await expect(page.getByRole('heading', { name: 'Screen 3', exact: true })).toBeInViewport()
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true)
  await page.screenshot({ path: '../../results/zbs_tabs/real-mobile.png', fullPage: true })
  expect(errors).toEqual([])
})
