"""
Notes tools for ChemCoworker.

Manages a multi-type chemistry knowledge base stored as markdown files:
  knowledge_base/notes/reactions/   — per-reaction-type notes (conditions, warnings, HTE data)
  knowledge_base/notes/mechanisms/  — mechanistic principle notes
  knowledge_base/notes/substrates/  — substrate class notes
  knowledge_base/notes/protocols/   — practical technique / procedure notes

An auto-generated _index.json (built by knowledge_base/notes/build_index.py) enables fast
faceted retrieval by type, bond, metal, and tags without reading every file.

Folder resolution (in priority order):
  1. CHEM_NOTES_PATH environment variable
  2. {project_root}/knowledge_base/notes/

Tools:
  - read_notes          : exact lookup by ID or reaction type across all note types
  - search_notes        : keyword + faceted search (type, bond_formed, metal, tags)
  - list_notes          : browse available notes by type
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional
import os

from ._helpers import _error, _success
from ._base import ToolPlugin

logger = logging.getLogger(__name__)

# Note type subfolders
_NOTE_SUBDIRS = ["reactions", "mechanisms", "substrates", "protocols"]


# ---------------------------------------------------------------------------
# Folder / path resolution
# ---------------------------------------------------------------------------

def get_notes_dir() -> Path:
    """Resolve the root notes folder."""
    env_path = os.getenv("CHEM_NOTES_PATH")
    if env_path:
        return Path(env_path)
    return Path(__file__).parent.parent.parent / "knowledge_base" / "notes"


def get_index_path() -> Path:
    return get_notes_dir() / "_index.json"


def _load_index() -> Dict[str, Any]:
    """Load _index.json if it exists. Returns {} on failure."""
    p = get_index_path()
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}


def get_notes_path(note_id: str, note_type: str = "reactions") -> Path:
    """Return the write-target path for a note file."""
    key = _normalize_reaction_type(note_id)
    return get_notes_dir() / note_type / f"{key}.md"


def _find_note_file(note_id: str, note_type: Optional[str] = None) -> Optional[Path]:
    """
    Find an existing note file by ID (or alias), searching subfolders.

    If note_type is given, only that subfolder is checked.
    Returns None if not found.
    """
    key = _normalize_reaction_type(note_id)
    notes_dir = get_notes_dir()

    # Fast path via index
    idx = _load_index()
    if idx.get("entries"):
        for entry in idx["entries"]:
            if entry.get("id") == key or key in entry.get("aliases", []):
                if note_type and entry.get("type") not in (note_type, note_type.rstrip("s")):
                    continue
                file_rel = entry.get("file", "")
                if file_rel:
                    p = notes_dir / file_rel
                    if p.exists():
                        return p

    # Fallback: scan subfolders
    check_subdirs = [note_type] if note_type else _NOTE_SUBDIRS
    for subdir in check_subdirs:
        p = notes_dir / subdir / f"{key}.md"
        if p.exists():
            return p

    return None


# ---------------------------------------------------------------------------
# Reaction type normalization (preserved for backward compatibility)
# ---------------------------------------------------------------------------

_ALIASES: Dict[str, str] = {
    # Suzuki
    "suzuki": "suzuki_miyaura",
    "suzuki coupling": "suzuki_miyaura",
    "suzuki-miyaura": "suzuki_miyaura",
    "suzuki miyaura": "suzuki_miyaura",
    # Buchwald-Hartwig
    "buchwald": "buchwald_hartwig_amination",
    "buchwald-hartwig": "buchwald_hartwig_amination",
    "buchwald hartwig": "buchwald_hartwig_amination",
    "bha": "buchwald_hartwig_amination",
    "buchwald_hartwig": "buchwald_hartwig_amination",
    "c-n coupling": "buchwald_hartwig_amination",
    # C-O coupling
    "c-o coupling": "c_o_coupling",
    "aryl etherification": "c_o_coupling",
    # Negishi
    "negishi": "negishi_coupling",
    # Heck
    "heck": "heck_reaction",
    # Kumada
    "kumada": "kumada_coupling",
    # Sonogashira
    "sonogashira": "sonogashira_coupling",
    # Chan-Lam
    "chan-lam": "chan_lam",
    "chan lam": "chan_lam",
    # SNAr
    "snar": "snar",
    "nucleophilic aromatic substitution": "snar",
    # Mitsunobu
    "mitsunobu": "mitsunobu",
    # Reductive amination
    "reductive amination": "reductive_amination",
    # Amide coupling
    "amide coupling": "amide_coupling",
    "peptide coupling": "amide_coupling",
    # Borylation
    "alkyl_borylation": "radical_borylation",
    "alkyl borylation": "radical_borylation",
    # Reduction
    "pinacolborane_reduction": "carboxylic_acid_to_aldehyde_reduction",
    "hbpin reduction": "carboxylic_acid_to_aldehyde_reduction",
}


def _normalize_reaction_type(raw: str) -> str:
    """Normalize a name to canonical snake_case."""
    cleaned = raw.strip().lower()
    if cleaned in _ALIASES:
        return _ALIASES[cleaned]
    return re.sub(r"[\s\-]+", "_", cleaned)


# ---------------------------------------------------------------------------
# Notes file I/O
# ---------------------------------------------------------------------------

def read_notes_text(note_id: str, note_type: Optional[str] = None) -> Optional[str]:
    """Read a note file by ID. Returns None if not found."""
    p = _find_note_file(note_id, note_type)
    if not p:
        return None
    text = p.read_text(encoding="utf-8").strip()
    return text if text else None


def append_notes(
    reaction_type: str,
    content: str,
    note_type: str = "reactions",
) -> Path:
    """Append extracted notes content to the notes file. Creates it if needed."""
    path = get_notes_path(reaction_type, note_type)
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.exists():
        existing = path.read_text(encoding="utf-8")
        path.write_text(
            existing.rstrip() + "\n\n---\n\n" + content.strip() + "\n",
            encoding="utf-8",
        )
    else:
        header = f"# {reaction_type.replace('_', ' ').title()} — {note_type.rstrip('s').title()} Notes\n\n"
        path.write_text(header + content.strip() + "\n", encoding="utf-8")

    return path


def _list_available(note_type: Optional[str] = None) -> List[str]:
    """Return available note IDs, optionally filtered by type."""
    idx = _load_index()
    if idx.get("entries"):
        entries = idx["entries"]
        if note_type:
            # match both "reaction" and "reactions"
            nt = note_type.rstrip("s")
            entries = [e for e in entries if e.get("type", "").rstrip("s") == nt]
        return [e["id"] for e in entries]

    # Fallback: scan filesystem
    notes_dir = get_notes_dir()
    ids = []
    subdirs = [note_type] if note_type else _NOTE_SUBDIRS
    for subdir in subdirs:
        d = notes_dir / subdir
        if d.exists():
            ids.extend(p.stem for p in sorted(d.glob("*.md"))
                       if not p.name.startswith(("README", "_", ".")))
    return ids


# ---------------------------------------------------------------------------
# Tool: read_notes
# ---------------------------------------------------------------------------

def _read_notes(
    id: str,
    note_type: Optional[str] = None,
    max_chars: int = 8000,
) -> Dict[str, Any]:
    """Load chemistry notes for any note type by ID.

    Searches reactions, mechanisms, substrates, and protocols.
    Use note_type to restrict to a specific subfolder.

    Args:
        id: Note ID or common alias (e.g. "suzuki_miyaura", "oxidative_addition",
            "aryl_halides", "suzuki", "buchwald"). Also accepts reaction type names.
        note_type: Optional filter — "reactions", "mechanisms", "substrates",
                   "protocols". If omitted, all subfolders are searched.
        max_chars: Maximum characters to return (default 8000).

    Returns:
        dict with id, note_type, has_notes, content, char_count, notes_file.
    """
    canonical = _normalize_reaction_type(id)
    note_path = _find_note_file(id, note_type)
    notes_text = None
    if note_path:
        text = note_path.read_text(encoding="utf-8").strip()
        notes_text = text if text else None

    if notes_text is None:
        available = _list_available(note_type)
        return _success({
            "id": canonical,
            "has_notes": False,
            "content": "",
            "char_count": 0,
            "notes_file": str(get_notes_path(canonical, note_type or "reactions")),
            "available": available,
            "note": (
                f"No notes found for '{canonical}'. "
                + (f"Available: {', '.join(available)}" if available else
                   "Notes folder is empty — run intake to extract notes from literature.")
            ),
        })

    total_chars = len(notes_text)
    display = notes_text
    if total_chars > max_chars:
        display = notes_text[:max_chars].rsplit("\n\n", 1)[0]
        display += f"\n\n[... truncated at {max_chars} chars; {total_chars} total]"

    resolved_type = note_path.parent.name if note_path else (note_type or "reactions")
    return _success({
        "id": canonical,
        "note_type": resolved_type,
        "has_notes": True,
        "content": display,
        "char_count": total_chars,
        "notes_file": str(note_path),
    })


read_notes_tool = ToolPlugin(
    name="read_notes",
    category="notes",
    description=(
        "Load chemistry notes by ID across all note types (reactions, mechanisms, "
        "substrates, protocols). Accepts reaction type names, mechanism names, "
        "substrate class names, or aliases like 'suzuki', 'oxidative_addition'. "
        "Use note_type='reactions' to restrict to reaction notes. "
        "Run in parallel with recommend_conditions when reaction type is confirmed."
    ),
    prerequisites=[],
    fn=_read_notes,
)

# Backward-compatible alias for existing agent prompts
read_reaction_notes_tool = ToolPlugin(
    name="read_reaction_notes",
    category="notes",
    description=(
        "Load accumulated chemistry notes for a reaction type (backward-compatible alias "
        "for read_notes with note_type='reactions'). Accepts 'suzuki_miyaura', 'buchwald_hartwig', etc."
    ),
    prerequisites=[],
    fn=lambda reaction_type, max_chars=6000: _read_notes(
        reaction_type, note_type="reactions", max_chars=max_chars
    ),
)


# ---------------------------------------------------------------------------
# Tool: search_notes
# ---------------------------------------------------------------------------

def _search_notes(
    query: str,
    note_types: Optional[List[str]] = None,
    tags: Optional[List[str]] = None,
    bond_formed: Optional[str] = None,
    metal: Optional[str] = None,
    top_k: int = 5,
) -> Dict[str, Any]:
    """Search all notes by keyword with optional faceted pre-filtering.

    Uses _index.json for fast faceted pre-filtering, then keyword-scores only
    the pre-filtered files. Falls back to scanning all files if no index.

    Args:
        query: Free-text search. Include catalysts, bond types, symptoms, etc.
               e.g. "copper sp3 coupling alkyl halide"
               e.g. "beta hydride elimination side reaction"
        note_types: Optional list to restrict note type, e.g. ["reactions"].
                    If omitted, all types are searched.
        tags: Optional list of tag strings to require — entry must match at least one.
              e.g. ["pd", "cross_coupling"]
        bond_formed: Optional bond type filter, e.g. "C-N" or "C-C".
        metal: Optional metal filter, e.g. "palladium" or "copper".
        top_k: Number of top results to return (default 5, max 10).

    Returns:
        dict with found, query, results (rank, score, id, note_type, source_file, tags, excerpt).
    """
    from .literature import _tokenize, _score_chunk

    notes_dir = get_notes_dir()
    if not notes_dir.exists():
        return _success({"found": 0, "query": query, "results": [], "notes_searched": []})

    # --- Faceted pre-filter via index ---
    idx = _load_index()
    candidate_files: List[Path] = []
    candidate_meta: Dict[str, Dict[str, Any]] = {}  # file_path_str → entry

    if idx.get("entries"):
        for entry in idx["entries"]:
            entry_id = entry.get("id", "")
            entry_type = entry.get("type", "")

            # note_types filter
            if note_types:
                nt_norm = [t.rstrip("s") for t in note_types]
                if entry_type.rstrip("s") not in nt_norm:
                    continue

            # bond_formed filter
            if bond_formed:
                bf_lower = bond_formed.lower()
                entry_bonds = [b.lower() for b in entry.get("bond_formed", [])]
                if not any(bf_lower in b for b in entry_bonds):
                    continue

            # metal filter
            if metal:
                m_lower = metal.lower()
                entry_metals = [m.lower() for m in entry.get("metal", [])]
                if not any(m_lower in m for m in entry_metals):
                    continue

            # tags filter (any match)
            if tags:
                entry_tags = [t.lower() for t in entry.get("tags", [])]
                req_lower = [t.lower() for t in tags]
                if not any(r in entry_tags for r in req_lower):
                    continue

            file_rel = entry.get("file", "")
            if file_rel:
                p = notes_dir / file_rel
                if p.exists():
                    candidate_files.append(p)
                    candidate_meta[str(p)] = entry
    else:
        # No index — scan all subfolders
        subdirs_to_search = note_types if note_types else _NOTE_SUBDIRS
        for subdir in subdirs_to_search:
            d = notes_dir / subdir
            if d.exists():
                for f in sorted(d.glob("*.md")):
                    if not f.name.startswith(("README", "_", ".")):
                        candidate_files.append(f)

    if not candidate_files:
        return _success({
            "found": 0, "query": query, "results": [],
            "notes_searched": [],
            "note": "No notes matched the faceted filters.",
        })

    query_tokens = _tokenize(query)
    if not query_tokens:
        return _error("Query is empty or contains only stopwords.")

    scored: List[tuple] = []

    for nf in candidate_files:
        try:
            text = nf.read_text(encoding="utf-8")
        except Exception:
            continue

        # Split into per-source sections
        sections = re.split(r"(?=## Source:)", text)
        # Also handle files without ## Source: headers (e.g. mechanism notes)
        if len(sections) <= 1:
            sections = [text]

        for section in sections:
            if len(section.strip()) < 50:
                continue

            tags_match = re.search(
                r"^tags\s*:\s*([^\n]+)", section, re.IGNORECASE | re.MULTILINE
            )
            tags_str = tags_match.group(1).strip() if tags_match else ""

            from .literature import _tokenize as tok
            tag_tokens = tok(tags_str) if tags_str else []
            tag_bonus = sum(2.0 for qt in query_tokens if qt in tag_tokens)

            body_score = _score_chunk(section, query_tokens)
            total_score = body_score + tag_bonus

            if total_score > 0:
                entry = candidate_meta.get(str(nf), {})
                note_type_resolved = entry.get("type") or nf.parent.name
                scored.append((total_score, section, nf, tags_str, note_type_resolved))

    scored.sort(key=lambda x: x[0], reverse=True)
    top_k = min(max(top_k, 1), 10)
    top = scored[:top_k]

    if not top:
        return _success({
            "found": 0,
            "query": query,
            "results": [],
            "notes_searched": [f.name for f in candidate_files],
            "note": f"No sections matched '{query}' across {len(candidate_files)} note file(s).",
        })

    results = []
    for rank, (score, section, nf, tags_str, note_type_res) in enumerate(top, 1):
        entry = candidate_meta.get(str(nf), {})
        excerpt = section.strip()
        if len(excerpt) > 800:
            excerpt = excerpt[:800].rsplit("\n", 1)[0] + "\n[...]"
        results.append({
            "rank": rank,
            "score": round(score, 3),
            "id": entry.get("id") or nf.stem,
            "note_type": note_type_res,
            "source_file": f"{nf.parent.name}/{nf.name}",
            "tags": tags_str,
            "excerpt": excerpt,
        })

    return _success({
        "found": len(results),
        "query": query,
        "results": results,
        "notes_searched": [f"{f.parent.name}/{f.name}" for f in candidate_files],
    })


search_notes_tool = ToolPlugin(
    name="search_notes",
    category="notes",
    description=(
        "Keyword + faceted search across ALL notes (reactions, mechanisms, substrates, protocols). "
        "Optional filters: note_types=['reactions'], bond_formed='C-N', metal='palladium', tags=['pd']. "
        "Use when reaction type is uncertain, for troubleshooting, or cross-cutting topics "
        "(e.g. 'copper catalyst sp3 coupling', 'beta-hydride elimination', 'oxidative addition palladium'). "
        "Tag-exact matches rank above body-text matches."
    ),
    prerequisites=[],
    fn=_search_notes,
)


# ---------------------------------------------------------------------------
# Tool: list_notes
# ---------------------------------------------------------------------------

def _list_notes(note_type: Optional[str] = None) -> Dict[str, Any]:
    """List all available notes, optionally filtered by type.

    Use this to discover what knowledge exists before deciding which
    read_notes or search_notes call to make.

    Args:
        note_type: Optional filter — "reactions", "mechanisms", "substrates",
                   "protocols". If omitted, all types are returned.

    Returns:
        dict with total_count and entries (id, type, title, tags, bond_formed, metal).
    """
    idx = _load_index()

    if idx.get("entries"):
        entries = idx["entries"]
        if note_type:
            nt = note_type.rstrip("s")
            entries = [e for e in entries if e.get("type", "").rstrip("s") == nt]
        summary = []
        for e in entries:
            item: Dict[str, Any] = {
                "id": e.get("id"),
                "type": e.get("type"),
                "title": e.get("title"),
            }
            if e.get("tags"):
                item["tags"] = e["tags"]
            if e.get("bond_formed"):
                item["bond_formed"] = e["bond_formed"]
            if e.get("metal"):
                item["metal"] = e["metal"]
            summary.append(item)
        return _success({"total_count": len(summary), "entries": summary})

    # Fallback if no index
    notes_dir = get_notes_dir()
    available: List[Dict[str, Any]] = []
    subdirs = [note_type] if note_type else _NOTE_SUBDIRS
    for subdir in subdirs:
        d = notes_dir / subdir
        if d.exists():
            for p in sorted(d.glob("*.md")):
                if not p.name.startswith(("README", "_", ".")):
                    available.append({"id": p.stem, "type": subdir.rstrip("s"), "file": f"{subdir}/{p.name}"})

    if not available:
        return _success({
            "total_count": 0,
            "entries": [],
            "note": "No notes found. Run 'python knowledge_base/notes/build_index.py' to build the index.",
        })

    return _success({"total_count": len(available), "entries": available})


list_notes_tool = ToolPlugin(
    name="list_notes",
    category="notes",
    description=(
        "List all available notes by type. Use to discover what knowledge exists "
        "before calling read_notes or search_notes. Optional filter: note_type='reactions'|'mechanisms'|'substrates'|'protocols'."
    ),
    prerequisites=[],
    fn=_list_notes,
)


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

NOTES_TOOLS = [
    read_notes_tool,
    read_reaction_notes_tool,
    search_notes_tool,
    list_notes_tool,
]
