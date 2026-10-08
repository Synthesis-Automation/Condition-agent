import { expect, test } from '@playwright/test'

const result = {
  target_smiles: 'CC(=O)N1CCCC1', search_status: 'complete', returned_count: 1,
  exact_target: { status: 'not_found' }, attempts: [{ level: 'core', decision: 'add_context' }],
  limitations: ['Discovery is not a validated route.'],
  execution: { elapsed_seconds: 1, library_loads: 1, candidate_reuses: 1 },
  hits: [{ observation_id: 'obs-1', reference_id: 'ref-1', matched_molecule_smiles: 'N1CCCC1',
    record: { reaction_smiles: 'NCCCC>>N1CCCC1', conditions: 'Fixture source conditions' }, procedures: [],
    discovery: { query: '[#7]1-[#6]-[#6]-[#6]-[#6]-1', core_relationship: 'constructed',
      explanation: 'The source constructs the selected ring core.', alignment_ambiguous: false,
      nitrogen_hydrogen_differences: [{ target_atom_id: 3, target_hydrogens: 0, source_hydrogen_counts: [1] }] },
  }],
}

test.beforeEach(async ({ page }) => {
  await page.route('**/api/v1/capabilities', route => route.fulfill({ json: { data: { fragment_search: true } } }))
  await page.route('**/api/v1/render/*', route => route.fulfill({ contentType: 'image/svg+xml', body: '<svg xmlns="http://www.w3.org/2000/svg"/>' }))
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
})

test('target-only default discovers and explains sources, then permits explicit refinement', async ({ page }) => {
  await expect(page.getByLabel('Workflow', { exact: true })).toHaveValue('discovery')
  await expect(page.getByLabel('Core fragment', { exact: true })).toBeHidden()
  await expect(page.getByLabel('Query format', { exact: true })).toBeHidden()
  let calls = 0
  await page.route('**/api/v1/fragments/discover', async route => {
    calls++
    expect(route.request().postDataJSON()).toEqual({ target_smiles: result.target_smiles })
    await route.fulfill({ json: { data: result } })
  })
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill(result.target_smiles)
  expect(calls).toBe(0)
  await page.getByRole('button', { name: 'Find synthesis precedents', exact: true }).click()
  await expect(page.getByText(/Nitrogen substitution differs/)).toBeVisible()
  await page.getByText('Inspect source and conditions', { exact: true }).click()
  await expect(page.getByText(/Fixture source conditions/)).toBeVisible()
  const download = page.waitForEvent('download')
  await page.getByRole('button', { name: 'Export search evidence' }).click()
  expect((await download).suggestedFilename()).toBe('synthesis_precedents.json')
  await page.getByText('Refine search or investigate transfer', { exact: true }).click()
  await page.getByRole('button', { name: 'Use this core in detailed research' }).click()
  await expect(page.getByLabel('Workflow', { exact: true })).toHaveValue('manual')
  await expect(page.getByLabel('Query format', { exact: true })).toHaveValue('smarts')
  await expect(page.getByLabel('Query topology', { exact: true })).toHaveValue('subgraph')
  await expect(page.getByLabel('Core fragment', { exact: true })).toHaveValue(result.hits[0].discovery.query)
})

test('failed discovery is visible and a target edit clears stale results', async ({ page }) => {
  await page.route('**/api/v1/fragments/discover', route => route.fulfill({ json: { data: result } }))
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill(result.target_smiles)
  await page.getByRole('button', { name: 'Find synthesis precedents', exact: true }).click()
  await expect(page.getByRole('heading', { name: 'Synthesis precedents' })).toBeVisible()
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('bad')
  await expect(page.getByRole('heading', { name: 'Synthesis precedents' })).toBeHidden()
  await page.route('**/api/v1/fragments/discover', route => route.fulfill({ status: 422, json: { detail: { message: 'Invalid target structure' } } }))
  await page.getByRole('button', { name: 'Find synthesis precedents', exact: true }).click()
  await expect(page.getByRole('alert')).toHaveText('Invalid target structure')
})
