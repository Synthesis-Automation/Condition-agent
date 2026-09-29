import { expect, test } from '@playwright/test'

const counted = (value: number, precision = 'exact') => ({ value, precision })
const base = {
  query: { expression: 'COC', query_format: 'smiles', topology: 'preserve_rings' },
  search_status: 'complete', stop_reason: null,
  counts: { products: counted(1), observations: counted(1), known_references: counted(1) },
  relationship_groups: { constructed: counted(1), unresolved: counted(1) },
  hits: [{ hit_id: 'hit-1', observation_id: 'obs-1', reference_id: 'ref-1',
    product_smiles: 'COC', relationships: ['constructed', 'unresolved'],
    citation_availability: 'identifier_present', procedure_availability: 'linked',
    procedure_match_scope: 'exact_observation', procedure_records_truncated: false,
    admission_tier: 'review', warnings: [], matches: [{ witnesses: ['formed bond'] }],
    record: { reaction_smiles: 'CBr.CO>>COC', reference_identity: { citation: 'Example patent' } },
    procedures: [{ link_scope: 'exact_observation', record: { text: {
      chunks: [{ text: 'Procedure text from source.' }], truncated: false,
    } } }],
  }], returned_count: 1, source_scope: 'full_source', source_coverage_complete: true,
  refinement_hints: [], limitations: ['Product presence is not evidence of construction.'],
  execution: { elapsed_seconds: 0.5 },
}

test.beforeEach(async ({ page }) => {
  await page.route('**/api/v1/capabilities', route => route.fulfill({ json: { data: { fragment_search: true } } }))
  await page.route('**/api/v1/ranking-profiles', route => route.fulfill({ json: { data: { profiles: [] } } }))
  await page.route('**/api/v1/forward-synthesis/condition-profiles', route => route.fulfill({ json: { data: {} } }))
  await page.route('**/api/v1/render/*', route => route.fulfill({ contentType: 'image/svg+xml',
    body: '<svg xmlns="http://www.w3.org/2000/svg" width="260" height="120"><text x="10" y="30">Structure preview</text></svg>' }))
})

test('fragment search renders evidence, procedures and exports full JSON', async ({ page }) => {
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await page.route('**/api/v1/fragments/search', async route => {
    expect(route.request().postDataJSON()).toEqual({ query: 'COC', query_format: 'smiles',
      topology: 'preserve_rings', limit: 5, timeout_seconds: 10 })
    await route.fulfill({ json: { data: base } })
  })
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment search', exact: true }).check()
  await expect(page.getByRole('button', { name: 'Recommend conditions', exact: true })).toBeHidden()
  await page.getByLabel('Core fragment', { exact: true }).fill('COC')
  await page.getByRole('button', { name: 'Search fragments', exact: true }).click()
  await expect(page.getByRole('heading', { name: '1. constructed · unresolved' })).toBeVisible()
  await expect(page.getByAltText('Precedent 1')).toBeVisible()
  await page.getByText('Procedure 1 · exact observation', { exact: true }).click()
  await expect(page.getByText(/Procedure text from source/)).toBeVisible()
  const download = page.waitForEvent('download')
  await page.getByRole('button', { name: 'Export JSON', exact: true }).click()
  expect((await download).suggestedFilename()).toBe('fragment_precedents.json')
  await page.getByLabel('Core fragment', { exact: true }).fill('CO')
  await expect(page.getByRole('heading', { name: 'Fragment precedents', exact: true })).toBeHidden()
  expect(errors).toEqual([])
})

test('broad and partial searches remain distinct from zero matches', async ({ page }) => {
  let calls = 0
  await page.route('**/api/v1/fragments/search', async route => {
    calls += 1
    await route.fulfill({ json: { data: { ...base, hits: [], returned_count: 0,
      search_status: calls === 1 ? 'too_broad' : 'partial',
      stop_reason: calls === 1 ? 'confirmed_product_limit' : 'deadline',
      counts: { products: counted(501, 'at_least'), observations: counted(0, 'at_least'), known_references: counted(0, 'at_least') },
      refinement_hints: ['Retain the complete ring system'],
    } } })
  })
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment search', exact: true }).check()
  await page.getByRole('button', { name: 'Cyclic ether example' }).click()
  await page.getByRole('button', { name: 'Search fragments', exact: true }).click()
  await expect(page.getByRole('heading', { name: 'Query too broad' })).toBeVisible()
  await expect(page.locator('.metric-strip').getByText('≥ 501', { exact: true })).toBeVisible()
  await page.getByRole('button', { name: 'Search fragments', exact: true }).click()
  await expect(page.getByRole('heading', { name: 'Partial search results' })).toBeVisible()
  await expect(page.getByText('No matching products in this index.')).toBeHidden()
})

test('server query errors are shown and inputs recover', async ({ page }) => {
  await page.route('**/api/v1/fragments/search', route => route.fulfill({ status: 422,
    json: { detail: { code: 'VALUEERROR', message: 'Invalid fragment' } } }))
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment search', exact: true }).check()
  await page.getByLabel('Core fragment', { exact: true }).fill('bad')
  await page.getByRole('button', { name: 'Search fragments', exact: true }).click()
  await expect(page.getByRole('alert')).toContainText('Invalid fragment')
  await expect(page.getByLabel('Core fragment', { exact: true })).toBeEnabled()
})

test('draw a fragment, reopen it, and search the exported structure', async ({ page }) => {
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  let submitted: Record<string, unknown> | null = null
  await page.route('**/api/v1/fragments/search', async route => {
    submitted = route.request().postDataJSON()
    await route.fulfill({ json: { data: base } })
  })
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment search', exact: true }).check()
  await page.getByRole('button', { name: 'Draw', exact: true }).click()
  const dialog = page.getByRole('dialog', { name: 'Draw the core fragment' })
  await expect(dialog.getByRole('button', { name: 'Load example', exact: true })).toBeEnabled()
  await dialog.getByRole('button', { name: 'Load example', exact: true }).click()
  await expect(dialog.getByText('Fragment loaded into the drawing canvas.', { exact: true })).toBeVisible()
  await dialog.getByRole('button', { name: 'Use drawing', exact: true }).click()
  await expect(dialog).toBeHidden()
  const input = page.getByLabel('Core fragment', { exact: true })
  const exported = await input.inputValue()
  expect(exported).toContain('O')
  expect(exported).not.toMatch(/[>.]/)
  await expect(page.getByAltText('Current fragment drawing', { exact: true })).toBeVisible()
  await page.getByRole('button', { name: 'Edit drawing', exact: true }).click()
  await expect(dialog.getByText('Existing fragment loaded.', { exact: true })).toBeVisible()
  await dialog.getByRole('button', { name: 'Clear', exact: true }).click()
  await dialog.getByRole('button', { name: 'Use drawing', exact: true }).click()
  await expect(dialog.getByText('Draw one connected core fragment.', { exact: true })).toBeVisible()
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click()
  await expect(input).toHaveValue(exported)
  await page.getByRole('button', { name: 'Search fragments', exact: true }).click()
  await expect(page.getByRole('heading', { name: 'Fragment precedents', exact: true })).toBeVisible()
  expect(submitted).toMatchObject({ query: exported, query_format: 'smiles' })
  expect(errors).toEqual([])
})

test('SMARTS remains text input and is never silently converted by drawing', async ({ page }) => {
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment search', exact: true }).check()
  await page.getByRole('combobox', { name: /Query format/ }).selectOption('smarts')
  const input = page.getByLabel('Core fragment', { exact: true })
  await input.fill('C[O,N]C')
  await expect(page.getByRole('button', { name: 'Edit drawing', exact: true })).toBeDisabled()
  await expect(page.getByText('Edit SMARTS queries as text. Drawing is available in SMILES mode.')).toBeVisible()
  await expect(input).toHaveValue('C[O,N]C')
  await expect(page.getByAltText('Current fragment drawing', { exact: true })).toBeHidden()
})

test('fragment mode shares the workbench layout and keeps its input on mode changes', async ({ page }) => {
  await page.goto('/')
  const modes = page.locator('.mode-switch')
  const original = await modes.boundingBox()
  await page.getByRole('radio', { name: 'Fragment search', exact: true }).check()
  const fragment = await modes.boundingBox()
  expect(fragment?.width).toBeCloseTo(original!.width, 0)
  await expect(page.locator('.analysis-options').getByLabel('Query format', { exact: true })).toBeVisible()
  await expect(page.locator('.editor-action-layout .reaction-paper').getByRole('heading', { name: 'Define the core fragment' })).toBeVisible()
  await expect(page.getByRole('heading', { name: 'Inspect fragment precedents' })).toBeVisible()
  await page.getByLabel('Core fragment', { exact: true }).fill('COC')
  await page.getByRole('radio', { name: 'Single-step retrosynthesis', exact: true }).check()
  await expect(page.getByRole('heading', { name: 'Define the target', exact: true })).toBeVisible()
  await page.getByRole('radio', { name: 'Fragment search', exact: true }).check()
  await expect(page.getByLabel('Core fragment', { exact: true })).toHaveValue('COC')
  await page.getByRole('button', { name: 'Clear', exact: true }).click()
  await expect(page.getByLabel('Core fragment', { exact: true })).toHaveValue('')
  await expect(page.getByRole('button', { name: 'Search fragments', exact: true })).toBeDisabled()
})
