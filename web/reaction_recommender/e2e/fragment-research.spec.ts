import { expect, test } from '@playwright/test'
import { readFile } from 'node:fs/promises'

const count = (value: number) => ({ value, precision: 'exact' })
const search = {
  query: { expression: 'COC', query_format: 'smiles', topology: 'preserve_rings' },
  search_status: 'complete', stop_reason: null, returned_count: 1,
  counts: { products: count(1), observations: count(1), known_references: count(1) },
  relationship_groups: { constructed: count(1) },
  hits: [{ hit_id: 'hit-1', observation_id: 'obs-1', reaction_id: 'source-1', reference_id: 'ref-1',
    product_smiles: 'COC', relationships: ['constructed'], warnings: [], admission_tier: 'review',
    citation_availability: 'identifier_present', procedure_availability: 'linked',
    procedure_match_scope: 'exact_observation', procedure_records_truncated: false,
    matches: [{ witnesses: ['formed bond'] }],
    record: { reaction_id: 'source-1', reaction_smiles: 'CBr.CO>>COC', reference_identity: { citation: 'Patent example' } },
    procedures: [{ link_scope: 'exact_observation', record: { text: { chunks: [{ text: 'Source procedure.' }], truncated: false } } }],
  }], source_scope: 'full_source', source_coverage_complete: true, refinement_hints: [],
  limitations: ['Presence does not establish construction.'], execution: { elapsed_seconds: 0.1 },
}
const transfer = {
  schema_version: 'fragment_guided_retrosynthesis.v1', target_smiles: 'CCOCC', library_mode: 'full',
  policy: {}, suggestions: { candidates: [] }, searches: [], query_search: search,
  guidance: { focus_bonds: [], eligible_bond_count: 0, eligible_source_count: 0, source_records: [], exclusions: [] },
  transfers: { baseline: { status: 'not_requested', candidates: [], diagnostics: null }, guided: [],
    compiled_source_template_count: 0, source_admissions: [],
    comparison: { baseline_requested: false, baseline_unique_precursor_count: 0, guided_unique_precursor_count: 0,
      additional_guided_precursor_sets: [], baseline_work: { template_applications: 0, validation_attempts: 0 },
      guided_work: { template_applications: 0, validation_attempts: 0 } } },
  source_comparisons: [{ observation_id: 'obs-1', reaction_id: 'source-1', comparison: {
    status: 'completed', core_atom_count: 3, right_coverage: 0.6, alignment_ambiguous: false,
  } }], execution: { elapsed_seconds: 0.2 }, limitations: ['Single-step hypotheses.'],
}

test.beforeEach(async ({ page }) => {
  await page.route('**/api/v1/capabilities', route => route.fulfill({ json: { data: {
    fragment_search: true, fragment_guided_retrosynthesis: true,
    retrosynthesis_library_modes: { compact: { library_available: true }, full: { library_available: true } },
  } } }))
  await page.route('**/api/v1/ranking-profiles', route => route.fulfill({ json: { data: { profiles: [] } } }))
  await page.route('**/api/v1/forward-synthesis/condition-profiles', route => route.fulfill({ json: { data: {} } }))
  await page.route('**/api/v1/render/*', route => route.fulfill({ contentType: 'image/svg+xml',
    body: '<svg xmlns="http://www.w3.org/2000/svg" width="260" height="120"><text x="10" y="30">Preview</text></svg>' }))
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  await expect(page.getByLabel('Workflow', { exact: true })).toHaveValue('manual')
  await page.getByLabel('Target molecule SMILES', { exact: true }).fill('CCOCC')
})

test('revises queries independently, selects sources, assesses and exports linked history', async ({ page }) => {
  const errors: string[] = []
  page.on('pageerror', error => errors.push(error.message))
  let searches = 0
  await page.route('**/api/v1/fragments/search', async route => {
    const request = route.request().postDataJSON()
    expect(request.target_smiles).toBe('CCOCC')
    expect(request.query).toBe(searches === 0 ? 'CCOCC' : 'COC')
    searches += 1
    await route.fulfill({ json: { data: searches === 1 ? { ...search, hits: [], returned_count: 0,
      counts: { products: count(0), observations: count(0), known_references: count(0) } } : search } })
  })
  await page.route('**/api/v1/fragments/suggest', async route => {
    expect(route.request().postDataJSON()).toEqual({ target_smiles: 'CCOCC' })
    await route.fulfill({ json: { data: { target_smiles: 'CCOCC', candidates: [{
      candidate_id: 'f1', kind: 'functional_region', query: 'COC', query_format: 'smiles',
      topology: 'preserve_rings', target_atom_ids: [1, 2, 3], reasons: ['Ether region'], cautions: [],
    }] } } })
  })
  await page.route('**/api/v1/retrosynthesis/fragment-transfer', async route => {
    expect(route.request().postDataJSON()).toEqual({ target_smiles: 'CCOCC', query: 'COC',
      query_format: 'smiles', topology: 'preserve_rings', limit: 10, timeout_seconds: 30,
      selected_observation_ids: ['obs-1'], library_mode: 'full', max_focus_bonds: 3, top_k: 3,
      include_baseline: false })
    await route.fulfill({ json: { data: transfer } })
  })
  await page.getByRole('button', { name: 'Use full target as query' }).click()
  await page.getByRole('button', { name: 'Search chosen fragment' }).click()
  await expect(page.getByText('No matching compounds on the selected reaction side.')).toBeVisible()
  await page.getByRole('button', { name: 'Suggest simpler queries' }).click()
  await page.getByRole('button', { name: 'Use region 1' }).click()
  await expect(page.getByLabel('Target molecule SMILES', { exact: true })).toHaveValue('CCOCC')
  await expect(page.getByLabel('Core fragment', { exact: true })).toHaveValue('COC')
  expect(searches).toBe(1)
  await page.getByLabel('Why this region or query revision?').fill('Omit peripheral alkyl context')
  await page.getByRole('button', { name: 'Search chosen fragment' }).click()
  await page.getByRole('checkbox', { name: /source-1/ }).check()
  await page.getByText(/Procedure 1/).click()
  await expect(page.locator('.fragment-hit').getByText('Source procedure.', { exact: false })).toBeVisible()
  await page.getByRole('button', { name: 'Assess selected precedents on target' }).click()
  await expect(page.getByRole('heading', { name: 'Selected precedent transfer assessment' })).toBeVisible()
  await expect(page.getByRole('heading', { name: 'Unrestricted baseline' })).toBeHidden()
  await expect(page.getByText(/Baseline comparison was not requested/)).toBeVisible()
  await expect(page.getByRole('heading', { name: 'Research history (3/20)' })).toBeVisible()
  const downloaded = page.waitForEvent('download')
  await page.getByRole('button', { name: 'Export research JSON' }).click()
  const file = await downloaded
  expect(file.suggestedFilename()).toBe('fragment_research_session.json')
  const saved = JSON.parse(await readFile((await file.path())!, 'utf8'))
  expect(saved.history.map((item: { parent_id: string | null }) => item.parent_id)).toEqual([null, 'attempt-1', 'attempt-2'])
  expect(saved.history[1].note).toBe('Omit peripheral alkyl context')
  expect(saved.history[2].request.selected_observation_ids).toEqual(['obs-1'])
  await page.getByLabel('Core fragment', { exact: true }).fill('CO')
  await expect(page.getByRole('heading', { name: 'Selected precedent transfer assessment' })).toBeHidden()
  await expect(page.getByRole('heading', { name: 'Research history (3/20)' })).toBeVisible()
  await page.getByRole('combobox', { name: 'Construction bonds', exact: true }).selectOption('5')
  await page.getByText(/attempt-3 .* transfer .* assessed/).click()
  await page.getByRole('button', { name: 'Restore attempt-3' }).click()
  await expect(page.getByRole('combobox', { name: 'Construction bonds', exact: true })).toHaveValue('3')
  await expect(page.getByRole('checkbox', { name: /source-1/ })).toBeChecked()
  await expect(page.getByLabel('Core fragment', { exact: true })).toHaveValue('COC')
  expect(errors).toEqual([])
})

test('records query errors, restores an attempt, and blocks partial evidence from transfer', async ({ page }) => {
  let calls = 0
  await page.route('**/api/v1/fragments/search', async route => {
    calls += 1
    if (calls === 1) await route.fulfill({ status: 422, json: { detail: { message: 'Query does not match target' } } })
    else await route.fulfill({ json: { data: { ...search, search_status: 'partial', stop_reason: 'deadline' } } })
  })
  await page.getByLabel('Core fragment', { exact: true }).fill('CN')
  await page.getByRole('button', { name: 'Search chosen fragment' }).click()
  await expect(page.getByRole('alert')).toContainText('Query does not match target')
  await page.getByLabel('Core fragment', { exact: true }).fill('COC')
  await page.getByRole('button', { name: 'Search chosen fragment' }).click()
  await expect(page.getByRole('heading', { name: 'Partial search results' })).toBeVisible()
  await expect(page.getByRole('checkbox', { name: /source-1/ })).toBeDisabled()
  await expect(page.getByRole('button', { name: 'Assess selected precedents on target' })).toBeDisabled()
  await page.getByText(/attempt-1 .* search .* error/).click()
  await page.getByRole('button', { name: 'Restore attempt-1' }).click()
  await expect(page.getByLabel('Core fragment', { exact: true })).toHaveValue('CN')
  await expect(page.getByRole('heading', { name: 'Research history (2/20)' })).toBeVisible()
})

test('SMARTS stays explicit and optional baseline is sent', async ({ page }) => {
  await page.getByLabel('Query format', { exact: true }).selectOption('smarts')
  await page.getByLabel('Query topology', { exact: true }).selectOption('subgraph')
  await page.getByLabel('Core fragment', { exact: true }).fill('C[O,N]C')
  await expect(page.getByRole('region', { name: 'Define the core fragment', exact: true }).getByRole('button', { name: 'Edit drawing', exact: true })).toBeDisabled()
  await expect(page.getByRole('button', { name: 'Suggest simpler queries' })).toBeDisabled()
  await page.route('**/api/v1/fragments/search', async route => {
    expect(route.request().postDataJSON()).toMatchObject({ target_smiles: 'CCOCC', query: 'C[O,N]C', query_format: 'smarts', topology: 'subgraph' })
    await route.fulfill({ json: { data: search } })
  })
  await page.getByRole('button', { name: 'Search chosen fragment' }).click()
  await page.getByRole('checkbox', { name: /source-1/ }).check()
  await page.getByText('Advanced options', { exact: true }).click()
  await page.getByRole('checkbox', { name: /Include bounded unrestricted baseline/ }).check()
  await page.route('**/api/v1/retrosynthesis/fragment-transfer', async route => {
    expect(route.request().postDataJSON().include_baseline).toBe(true)
    await route.fulfill({ json: { data: transfer } })
  })
  await page.getByRole('button', { name: 'Assess selected precedents on target' }).click()
  await expect(page.getByRole('heading', { name: 'Selected precedent transfer assessment' })).toBeVisible()
})

test('previews an explicit relaxation and preserves a rejected precedent investigation', async ({ page }) => {
  await page.route('**/api/v1/fragments/query-alternatives', async route => {
    expect(route.request().postDataJSON()).toEqual({ target_smiles: 'CCOCC', query: 'COC', query_format: 'smiles', topology: 'preserve_rings' })
    await route.fulfill({ json: { data: { target_smiles: 'CCOCC', parent_query: { expression: 'COC', query_id: 'q1' }, query_atoms: [], limitations: [],
      variants: [{ variant_id: 'v1', parent_query_id: 'q1', query: 'COC', query_format: 'smiles', topology: 'subgraph', relaxations: ['ring_boundary'], reason: 'Allow extra ring fusion.', alignment_ambiguous: true, target_alignments_truncated: false }] } } })
  })
  await page.route('**/api/v1/fragments/search', route => route.fulfill({ json: { data: search } }))
  await page.route('**/api/v1/fragments/investigate', async route => {
    expect(route.request().postDataJSON().observation_id).toBe('obs-1')
    expect(route.request().postDataJSON().topology).toBe('subgraph')
    await route.fulfill({ json: { data: { target_smiles: 'CCOCC', source: search.hits[0],
      comparison: { status: 'completed', alignment_ambiguous: true }, search_scope: { search_status: 'complete', stop_reason: null },
      transfer: { status: 'source_compilation_rejected', source_admissions: [{ reaction_id: 'source-1', status: 'rejected', reason: 'materialized_core_not_verified' }] }, limitations: [] } } })
  })
  await page.getByLabel('Core fragment', { exact: true }).fill('COC')
  await page.getByRole('button', { name: 'Preview query alternatives' }).click()
  await expect(page.getByRole('region', { name: 'Query alternatives' })).toContainText('ring boundary')
  await page.getByRole('button', { name: 'Use alternative 1' }).click()
  await expect(page.getByLabel('Query topology', { exact: true })).toHaveValue('subgraph')
  await page.getByRole('button', { name: 'Search chosen fragment' }).click()
  await page.getByRole('button', { name: 'Investigate this precedent' }).click()
  await expect(page.getByRole('region', { name: 'Precedent investigation' })).toContainText('operator compilation rejected it')
  await expect(page.getByRole('region', { name: 'Precedent investigation' })).toContainText('materialized_core_not_verified')
  await expect(page.getByRole('heading', { name: 'Research history (2/20)' })).toBeVisible()
  const downloadPromise = page.waitForEvent('download')
  await page.getByRole('button', { name: 'Export research JSON' }).click()
  const download = await downloadPromise
  const saved = JSON.parse(await readFile((await download.path())!, 'utf8'))
  expect(saved.history[0].query_revision.variant_id).toBe('v1')
  expect(saved.history[1].parent_id).toBe(saved.history[0].id)
  expect(saved.history[1].result.source.observation_id).toBe('obs-1')
  await page.getByLabel('Core fragment', { exact: true }).fill('CC')
  await expect(page.getByRole('region', { name: 'Precedent investigation' })).toBeHidden()
  await page.getByText(/attempt-2 .* investigation .* assessed/).click()
  await page.getByRole('button', { name: 'Restore attempt-2' }).click()
  await expect(page.getByRole('region', { name: 'Precedent investigation' })).toBeVisible()
})

test('leaving the mode aborts a pending revision and preserves completed history', async ({ page }) => {
  let calls = 0
  let release!: () => void
  let delivered!: () => void
  const gate = new Promise<void>(resolve => { release = resolve })
  const completed = new Promise<void>(resolve => { delivered = resolve })
  await page.route('**/api/v1/fragments/search', async route => {
    calls += 1
    if (calls === 2) await gate
    try { await route.fulfill({ json: { data: search } }) } finally { if (calls === 2) delivered() }
  })
  await page.getByLabel('Core fragment', { exact: true }).fill('COC')
  await page.getByRole('button', { name: 'Search chosen fragment' }).click()
  await expect(page.getByRole('heading', { name: 'Research history (1/20)' })).toBeVisible()
  await page.getByLabel('Core fragment', { exact: true }).fill('CO')
  const requested = page.waitForRequest('**/api/v1/fragments/search')
  await page.getByRole('button', { name: 'Search chosen fragment' }).click()
  await requested
  await page.getByRole('radio', { name: 'Analyze reactions', exact: true }).check()
  release(); await completed
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  await expect(page.getByLabel('Core fragment', { exact: true })).toHaveValue('CO')
  await expect(page.getByRole('heading', { name: 'Research history (1/20)' })).toBeVisible()
  await expect(page.getByRole('heading', { name: 'Fragment precedents', exact: true })).toBeHidden()
  await expect(page.getByRole('button', { name: 'Search chosen fragment' })).toBeEnabled()
})
