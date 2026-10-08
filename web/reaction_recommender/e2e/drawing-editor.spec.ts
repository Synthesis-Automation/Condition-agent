import { expect, test } from '@playwright/test'

const reaction = 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1'
const profile = process.env.WEBUI_TEST_PROFILE ?? 'recommendation_only'

test('fragment drawing saves edited SMILES and displays the changed structure', async ({ page }) => {
  test.skip(profile !== 'research_workbench', 'Fragment research is a workbench feature')
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  const input = page.getByLabel('Core fragment', { exact: true })
  const region = page.getByRole('region', { name: 'Define the core fragment', exact: true })
  await input.fill('CCO')
  await region.getByRole('button', { name: 'Edit drawing', exact: true }).click()
  const dialog = page.getByRole('dialog', { name: 'Draw the core fragment' })
  await expect(dialog.getByText('Existing fragment loaded.', { exact: true })).toBeVisible()
  await dialog.getByLabel('Fragment SMILES', { exact: true }).fill('c1ccccc1')
  await dialog.getByRole('button', { name: 'Use edited SMILES', exact: true }).click()
  await expect(dialog).toBeHidden()
  await expect(input).toHaveValue('c1ccccc1')
  await expect(region.getByRole('img', { name: 'Current fragment drawing', exact: true })).toBeVisible()
  await region.getByRole('button', { name: 'Edit drawing', exact: true }).click()
  await expect(dialog.getByText('Existing fragment loaded.', { exact: true })).toBeVisible()
  await dialog.getByLabel('Fragment SMILES', { exact: true }).fill('CCN')
  await dialog.getByRole('button', { name: 'Load SMILES', exact: true }).click()
  await expect(dialog.getByText('Fragment loaded into the drawing canvas.', { exact: true })).toBeVisible()
  await dialog.getByRole('button', { name: 'Use drawing', exact: true }).click()
  await expect(dialog).toBeHidden()
  await expect(input).toHaveValue('CCN')
  await expect(region.getByRole('img', { name: 'Current fragment drawing', exact: true })).toBeVisible()
})

test('fragment drawing shows validation errors in the dialog and preserves input on cancel', async ({ page }) => {
  test.skip(profile !== 'research_workbench', 'Fragment research is a workbench feature')
  await page.goto('/')
  await page.getByRole('radio', { name: 'Fragment-guided retro', exact: true }).check()
  const input = page.getByLabel('Core fragment', { exact: true })
  await input.fill('CCO')
  await page.getByRole('region', { name: 'Define the core fragment', exact: true }).getByRole('button', { name: 'Edit drawing', exact: true }).click()
  const dialog = page.getByRole('dialog', { name: 'Draw the core fragment' })
  await expect(dialog.getByText('Existing fragment loaded.', { exact: true })).toBeVisible()
  await dialog.getByLabel('Fragment SMILES', { exact: true }).fill('CC.O')
  await dialog.getByRole('button', { name: 'Use edited SMILES', exact: true }).click()
  await expect(dialog.getByRole('alert')).toHaveText('Draw one connected core fragment.')
  await expect(input).toHaveValue('CCO')
  await dialog.getByLabel('Fragment SMILES', { exact: true }).fill('invalid structure')
  await dialog.getByRole('button', { name: 'Use edited SMILES', exact: true }).click()
  await expect(dialog.getByRole('alert')).toBeVisible()
  await expect(dialog.getByRole('button', { name: 'Use edited SMILES', exact: true })).toBeEnabled()
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click()
  await expect(input).toHaveValue('CCO')
})

test.beforeEach(async ({ request }) => {
  const response = await request.get('/api/v1/health')
  expect((await response.json()).data.deployment_profile).toBe(profile)
})

for (const path of ['/', '/workbench.html']) {
  test(`${path}: load, export, and reopen a drawing`, async ({ page }) => {
    const errors: string[] = []
    page.on('pageerror', error => errors.push(error.message))
    await page.goto(path)
    await expect(page).toHaveURL(/\/$/)
    await expect(page).toHaveTitle(profile === 'research_workbench'
      ? 'Reaction Research Workbench' : /Condition Desk/)
    await page.getByRole('button', { name: 'Draw', exact: true }).click()
    const dialog = page.getByRole('dialog', { name: 'Draw the transformation' })
    await expect(dialog.getByRole('button', { name: 'Load example', exact: true })).toBeEnabled()
    await dialog.getByRole('button', { name: 'Load example', exact: true }).click()
    await expect(dialog.getByText('Reaction loaded into the drawing canvas.', { exact: true })).toBeVisible()
    await dialog.getByRole('button', { name: 'Use drawing', exact: true }).click()
    await expect(dialog).toBeHidden()
    const input = page.getByLabel('Reaction SMILES', { exact: true })
    await expect(input).toHaveValue(/.+>>.+/)
    const exported = await input.inputValue()
    await page.getByRole('button', { name: 'Edit drawing', exact: true }).click()
    await expect(dialog.getByText('Existing reaction loaded.', { exact: true })).toBeVisible()
    await dialog.getByRole('button', { name: 'Use drawing', exact: true }).click()
    await expect(input).toHaveValue(exported)
    expect(errors).toEqual([])
  })

  for (const failure of ['download', 'render']) {
    test(`${path}: ${failure} failure preserves the app and input`, async ({ page }) => {
      await page.route('**/assets/KetcherCanvas-*.js', route => failure === 'download'
        ? route.abort()
        : route.fulfill({
          contentType: 'application/javascript',
          body: 'export default function Editor() { throw new Error("Editor initialization failed"); }',
        }))
      await page.goto(path)
      const input = page.getByLabel('Reaction SMILES', { exact: true })
      await input.fill(reaction)
      await page.getByRole('button', { name: 'Edit drawing', exact: true }).click()
      await expect(page.getByRole('alert').filter({ hasText: 'The drawing editor could not load.' })).toBeVisible()
      await expect(page.getByRole('button', { name: 'Use drawing', exact: true })).toBeDisabled()
      await page.getByRole('button', { name: 'Close editor', exact: true }).click()
      await expect(page.getByRole('dialog')).toBeHidden()
      await expect(input).toHaveValue(reaction)
      await input.fill('CCO>>CC=O')
      await expect(input).toHaveValue('CCO>>CC=O')
    })
  }
}
