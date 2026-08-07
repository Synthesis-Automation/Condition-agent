import Raphael from 'raphael'

type LegacyRequire = (moduleName: string) => unknown

const browserGlobal = globalThis as typeof globalThis & {
  Raphael?: typeof Raphael
  require?: LegacyRequire
}

browserGlobal.Raphael ??= Raphael
browserGlobal.require ??= (moduleName: string) => {
  if (moduleName === 'raphael') return Raphael
  if (moduleName === 'acorn') return undefined
  throw new Error(`Unsupported browser-side require: ${moduleName}`)
}

