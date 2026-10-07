import { expect, test } from '@playwright/test'

const arm = (precursor: string, focus = false) => ({
  status: 'completed', diagnostics: { validation_attempt_count: 1 },
  candidates: [{ precursor_smiles: precursor, proposed_reaction_smiles: `${precursor}>>COC`,
    forward_validation_status: 'verified_signature', abstraction_level: 'L2',
    precedent_reaction_ids: ['operator-source'], bond_focus_check: focus ? { status: 'verified' } : null }],
})

const result = {
  schema_version: 'fragment_guided_retrosynthesis.v1', target_smiles: 'COC', library_mode: 'full',
  policy: { query_limit: 3, max_focus_bonds: 3, require_complete_search: true },
  suggestions: { candidates: [{ candidate_id: 'fragment-1', kind: 'whole_target', query: 'COC',
    reasons: ['whole target'], cautions: [], target_highlight_svg: '<svg xmlns="http://www.w3.org/2000/svg" width="260" height="120"><text x="10" y="30">Target region</text></svg>' }] },
  searches: [{ candidate_id: 'fragment-1', artifact_ref: 'query:fragment-1', result: {
    search_status: 'complete', returned_count: 2, stop_reason: null,
  } }],
  guidance: { focus_bonds: [{ target_atom_ids: [0, 1], supports: [{ reaction_id: 'fragment-source', reference_id: 'ref-1' }] }],
    eligible_bond_count: 1, eligible_source_count: 1, source_records: [], exclusions: [] },
  transfers: { baseline: arm('CBr.CO'), guided: [{ target_atom_ids: [0, 1],
    direct_source_transfer: { status: 'no_admitted_source_operators', candidates: [], diagnostics: null },
    witness_directed_library: arm('CI.CO', true) }], compiled_source_template_count: 0,
    source_admissions: [{ reaction_id: 'fragment-source', status: 'rejected', reason: 'materialized_core_not_verified' }],
    comparison: { baseline_unique_precursor_count: 1, guided_unique_precursor_count: 1,
      additional_guided_precursor_sets: ['CI.CO'], baseline_work: { template_applications: 40, validation_attempts: 10 },
      guided_work: { template_applications: 80, validation_attempts: 12 } } },
  execution: { elapsed_seconds: 2.5 }, limitations: ['Experimental single-step proposals.'],
}

test.beforeEach(async ({ page }) => {
  await page.route('**/api/v1/capabilities', route => route.fulfill({ json: { data: {
    fragment_search: true, fragment_guided_retrosynthesis: true,
    retrosynthesis_library_modes: { compact: { library_available: true }, full: { library_available: true } },
  } } }))
  await page.route('**/api/v1/ranking-profiles', route => route.fulfill({ json: { data: { profiles: [] } } }))
  await page.route('**/api/v1/forward-synthesis/condition-profiles', route => route.fulfill({ json: { data: {} } }))
  await page.route('**/api/v1/render/*', route => route.fulfill({ contentType: 'image/svg+xml',
    body: '<svg xmlns="http://www.w3.org/2000/svg" width="260" height="120"><text x="10" y="30">Structure preview</text></svg>' }))
})

test('compares baseline and guided evidence, exports JSON, and resets on edits', async ({ page }) => {
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  await page.route('**/api/v1/retrosynthesis/fragment-guided', async route => {
    expect(route.request().postDataJSON()).toEqual({
      target_smiles: 'Fc(cn1)cc2c1c(c3ccccc3OC)n[nH]2', library_mode: 'full',
      query_limit: 3, max_focus_bonds: 3, top_k: 3,
    })
    await route.fulfill({ json: { data: result } })
  })
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  await page.getByLabel('Workflow', { exact: true }).selectOption('automatic')
  await expect(page.getByRole('button', { name: 'Recommend conditions', exact: true })).toBeHidden()
  await page.getByRole('button', { name: 'Fragment retro example' }).click()
  await page.getByRole('button', { name: 'Evaluate fragment-guided retro' }).click()
  await expect(page.getByRole('heading', { name: 'Fragment-guided comparison' })).toBeVisible()
  await expect(page.getByRole('heading', { name: 'Unrestricted baseline' })).toBeVisible()
  await expect(page.getByRole('heading', { name: 'Witness-directed library' })).toBeVisible()
  await expect(page.getByText('Additional to baseline', { exact: true })).toHaveCount(2)
  await expect(page.getByText('No source operators passed whole-reaction admission.')).toBeVisible()
  await expect(page.getByText(/Guided arms use more total search work/)).toBeVisible()
  await expect(page.getByAltText('Selected fragment 1 on target')).toBeVisible()
  const download = page.waitForEvent('download')
  await page.getByRole('button', { name: 'Export JSON', exact: true }).click()
  expect((await download).suggestedFilename()).toBe('fragment_guided_retrosynthesis.json')
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('CCNC')
  await expect(page.getByRole('heading', { name: 'Fragment-guided comparison' })).toBeHidden()
  expect(errors).toEqual([])
})

test('retains incomplete search status and shows absent guidance separately from empty matches', async ({ page }) => {
  await page.route('**/api/v1/retrosynthesis/fragment-guided', route => route.fulfill({ json: { data: {
    ...result,
    searches: [{ ...result.searches[0], result: { search_status: 'partial', stop_reason: 'deadline', returned_count: 0 } }],
    guidance: { ...result.guidance, focus_bonds: [], exclusions: [{ reason: 'incomplete_or_unvalidated_search' }] },
    transfers: { ...result.transfers, guided: [] },
  } } }))
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  await page.getByLabel('Workflow', { exact: true }).selectOption('automatic')
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('COC')
  await page.getByRole('button', { name: 'Evaluate fragment-guided retro' }).click()
  await expect(page.getByText('partial', { exact: true })).toBeVisible()
  await expect(page.getByText(/No eligible construction guidance/)).toBeVisible()
  await expect(page.getByRole('heading', { name: 'Unrestricted baseline' })).toBeVisible()
})

test('reports unavailable libraries and server errors while allowing recovery', async ({ page }) => {
  await page.route('**/api/v1/retrosynthesis/fragment-guided', route => route.fulfill({ status: 422,
    json: { detail: { code: 'VALUEERROR', message: 'Target must be one connected molecule' } } }))
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  await page.getByLabel('Workflow', { exact: true }).selectOption('automatic')
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('C.C')
  await page.getByRole('button', { name: 'Evaluate fragment-guided retro' }).click()
  await expect(page.getByRole('alert')).toContainText('one connected molecule')
  await expect(page.getByLabel('Target molecule SMILES', { exact: true })).toBeEnabled()
  await page.route('**/api/v1/capabilities', route => route.fulfill({ json: { data: {
    fragment_search: true, retrosynthesis_library_modes: { full: { library_available: false } },
  } } }))
  await page.reload()
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  await page.getByLabel('Workflow', { exact: true }).selectOption('automatic')
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('COC')
  await expect(page.getByRole('button', { name: 'Evaluate fragment-guided retro' })).toBeDisabled()
  await expect(page.getByText(/Requires a prepared fragment index/)).toBeVisible()
})

test('changing modes clears results and ignores an in-flight response', async ({ page }) => {
  let release!: () => void
  let completed!: () => void
  const gate = new Promise<void>(resolve => { release = resolve })
  const delivered = new Promise<void>(resolve => { completed = resolve })
  await page.route('**/api/v1/retrosynthesis/fragment-guided', async route => {
    await gate
    try { await route.fulfill({ json: { data: result } }) } finally { completed() }
  })
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  await page.getByLabel('Workflow', { exact: true }).selectOption('automatic')
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('COC')
  const requested = page.waitForRequest('**/api/v1/retrosynthesis/fragment-guided')
  await page.getByRole('button', { name: 'Evaluate fragment-guided retro' }).click()
  await requested
  await page.getByRole('radio', { name: 'Analyze reactions', exact: true }).check()
  release()
  await delivered
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  await page.getByLabel('Workflow', { exact: true }).selectOption('automatic')
  await expect(page.getByRole('heading', { name: 'Fragment-guided comparison' })).toBeHidden()
  await expect(page.getByRole('button', { name: 'Evaluate fragment-guided retro' })).toBeEnabled()
  await expect(page.getByRole('button', { name: 'Evaluating…' })).toBeHidden()
})
