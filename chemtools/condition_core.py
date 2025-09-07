from typing import Dict, Any, List
from .contracts import Reagent
def parse_core(reagents: List[Reagent], text: str = "") -> Dict[str, Any]:
    metal=None; ligand=None
    def norm(n:str)->str:
        n=(n or '').lower()
        if 'phen' in n: return 'Phenanthroline'
        if 'dmeda' in n: return 'DMEDA'
        if 'tmeda' in n: return 'TMEDA'
        if 'xphos' in n: return 'XPhos'
        return (n or 'Ligand').title()
    for r in reagents:
        nm=(r.name or r.token or r.uid or '').lower()
        if r.role.upper() in ('CATALYST','METAL','METAL_SOURCE'):
            if 'cu' in nm: metal='Cu'
            if 'pd' in nm: metal='Pd'
        if r.role.upper()=='LIGAND':
            ligand=norm(nm)
    if not metal and ('copper' in text.lower() or 'cui' in text.lower()): metal='Cu'
    if not ligand and text: ligand=norm(text)
    core = f"{metal}/{ligand}" if metal and ligand else (f"{metal}/(no_ligand)" if metal else "Unknown/Unknown")
    return {"core": core,
            "metal_source_uid": next((r.uid for r in reagents if r.role.upper() in ('CATALYST','METAL','METAL_SOURCE')), None),
            "ligand_uid": next((r.uid for r in reagents if r.role.upper()=='LIGAND'), None),
            "precatalyst": False}
