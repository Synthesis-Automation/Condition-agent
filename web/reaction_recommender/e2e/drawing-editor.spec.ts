import { expect, test } from '@playwright/test'

const reaction = 'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1'
const profile = process.env.WEBUI_TEST_PROFILE ?? 'recommendation_only'

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
