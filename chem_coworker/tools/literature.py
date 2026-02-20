"""
Literature search and web fetch tools for ChemCoworker.

Searches a local folder of chemistry documents (papers, procedures, lab notes)
for relevant passages. No external API or embeddings required — pure keyword
scoring with chemistry-aware weighting.

Supported formats: .txt, .md, .pdf (requires pypdf or pdfminer.six)

Document folder (in priority order):
  1. CHEM_DOCS_PATH environment variable
  2. {project_root}/knowledge_base/sources/

Adding documents:
  - Drop .txt, .md, or .pdf files into the knowledge_base/sources/ folder
  - Optional: add a sidecar {filename}.meta.json with {"title": ..., "doi": ..., "year": ...}
  - Documents are cached in memory after first load; restart to pick up new files

Tools:
  - search_literature : keyword search over local chemistry documents
  - fetch_webpage     : fetch and extract text from a URL (optionally save to knowledge_base/sources/)
"""
from __future__ import annotations

import html as _html_module
import json
import logging
import math
import os
import re
import urllib.error
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ._helpers import _error, _success
from ._base import ToolPlugin

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Document folder resolution
# ---------------------------------------------------------------------------

def _get_docs_dir() -> Path:
    """Resolve the literature documents folder."""
    env_path = os.getenv("CHEM_DOCS_PATH")
    if env_path:
        return Path(env_path)
    # Default: {project_root}/knowledge_base/sources/
    return Path(__file__).parent.parent.parent / "knowledge_base" / "sources"


# ---------------------------------------------------------------------------
# Document loading
# ---------------------------------------------------------------------------

# Module-level cache: {filepath_str: [(chunk_text, metadata), ...]}
_DOC_CACHE: Dict[str, List[Tuple[str, Dict[str, Any]]]] = {}


def _load_text_file(path: Path) -> str:
    """Load a plain text or markdown file."""
    return path.read_text(encoding="utf-8", errors="replace")


def _extract_pdf_text(source: str, path: Optional[Path] = None, raw_bytes: Optional[bytes] = None) -> str:
    """Extract text from a PDF given either a file path or raw bytes."""
    import io

    # Try pypdf first (lightweight)
    try:
        import pypdf
        if raw_bytes is not None:
            reader = pypdf.PdfReader(io.BytesIO(raw_bytes))
        else:
            reader = pypdf.PdfReader(str(path))
        pages = [page.extract_text() or "" for page in reader.pages]
        return "\n\n".join(pages)
    except ImportError:
        pass
    except Exception as exc:
        logger.warning(f"[literature] pypdf failed on {source}: {exc}")

    # Try pdfminer.six
    try:
        from pdfminer.high_level import extract_text, extract_text_to_fp
        if raw_bytes is not None:
            from pdfminer.high_level import extract_text_to_fp
            buf = io.StringIO()
            extract_text_to_fp(io.BytesIO(raw_bytes), buf, output_type="text")
            return buf.getvalue()
        else:
            return extract_text(str(path))
    except ImportError:
        pass
    except Exception as exc:
        logger.warning(f"[literature] pdfminer failed on {source}: {exc}")

    logger.warning(f"[literature] Cannot read PDF {source}: install pypdf or pdfminer.six")
    return ""


def _load_pdf_file(path: Path) -> str:
    """Load a PDF file. Returns cached .txt sidecar text if available."""
    # Check for a pre-extracted text sidecar (avoids re-parsing on every run)
    txt_cache = path.with_suffix(path.suffix + ".txt")
    if txt_cache.exists():
        return txt_cache.read_text(encoding="utf-8", errors="replace")
    return _extract_pdf_text(source=path.name, path=path)


def _load_meta(path: Path) -> Dict[str, Any]:
    """Load optional sidecar metadata file: {path}.meta.json"""
    meta_path = path.with_suffix(path.suffix + ".meta.json")
    if meta_path.exists():
        try:
            return json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {}


def _chunk_text(text: str, source: str, meta: Dict[str, Any]) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Split text into paragraph-level chunks.

    Splits on double-newlines. Merges short paragraphs (<40 words) with the
    next one to avoid tiny isolated sentences. Caps chunks at ~400 words.
    """
    # Normalize line endings and collapse excessive blank lines
    text = re.sub(r"\r\n", "\n", text)
    text = re.sub(r"\n{3,}", "\n\n", text)

    raw_paras = [p.strip() for p in text.split("\n\n") if p.strip()]

    chunks: List[Tuple[str, Dict[str, Any]]] = []
    buffer = ""

    for para in raw_paras:
        word_count = len(para.split())
        if not buffer:
            buffer = para
        elif len(buffer.split()) + word_count <= 400:
            buffer = buffer + "\n\n" + para
        else:
            # Flush current buffer
            if len(buffer.split()) >= 20:  # skip very short chunks
                chunk_meta = {
                    "source": source,
                    "title": meta.get("title", source),
                    "doi": meta.get("doi", ""),
                    "year": meta.get("year", ""),
                }
                chunks.append((buffer, chunk_meta))
            buffer = para

    if buffer and len(buffer.split()) >= 20:
        chunk_meta = {
            "source": source,
            "title": meta.get("title", source),
            "doi": meta.get("doi", ""),
            "year": meta.get("year", ""),
        }
        chunks.append((buffer, chunk_meta))

    return chunks


def _load_document(path: Path) -> List[Tuple[str, Dict[str, Any]]]:
    """Load and chunk a single document. Results are cached by path."""
    cache_key = str(path)
    if cache_key in _DOC_CACHE:
        return _DOC_CACHE[cache_key]

    suffix = path.suffix.lower()
    if suffix == ".pdf":
        text = _load_pdf_file(path)
        # Persist extracted text so future runs skip PDF parsing entirely
        if text.strip():
            txt_cache = path.with_suffix(path.suffix + ".txt")
            if not txt_cache.exists():
                try:
                    txt_cache.write_text(text, encoding="utf-8")
                    logger.info(f"[literature] Saved text cache: {txt_cache.name}")
                except Exception as exc:
                    logger.debug(f"[literature] Could not write text cache: {exc}")
    elif suffix in (".txt", ".md"):
        text = _load_text_file(path)
    else:
        return []

    if not text.strip():
        return []

    meta = _load_meta(path)
    chunks = _chunk_text(text, source=path.name, meta=meta)
    _DOC_CACHE[cache_key] = chunks
    return chunks


def _load_all_documents(docs_dir: Path) -> List[Tuple[str, Dict[str, Any]]]:
    """Load all supported documents from the docs folder."""
    if not docs_dir.exists():
        return []

    all_chunks: List[Tuple[str, Dict[str, Any]]] = []
    supported = {".txt", ".md", ".pdf"}

    for path in sorted(docs_dir.iterdir()):
        if path.suffix.lower() in supported and not path.name.startswith("."):
            try:
                chunks = _load_document(path)
                all_chunks.extend(chunks)
            except Exception as exc:
                logger.warning(f"[literature] Failed to load {path.name}: {exc}")

    return all_chunks


# ---------------------------------------------------------------------------
# Keyword scoring
# ---------------------------------------------------------------------------

# Common English stopwords to ignore when building the query token set.
# Intentionally keeps chemistry-relevant words (yield, conditions, reaction, etc.)
_STOPWORDS = {
    "a", "an", "the", "and", "or", "but", "in", "on", "at", "to", "for",
    "of", "with", "by", "from", "as", "is", "was", "are", "were", "be",
    "been", "have", "has", "had", "do", "does", "did", "will", "would",
    "could", "should", "may", "might", "this", "that", "these", "those",
    "it", "its", "we", "our", "they", "their", "which", "who", "what",
    "when", "where", "how", "if", "not", "no", "can", "i", "you", "he",
    "she", "us", "them", "also", "both", "all", "any", "into", "through",
}


def _tokenize(text: str) -> List[str]:
    """Lowercase, strip punctuation, split into tokens, remove stopwords."""
    tokens = re.findall(r"[a-zA-Z0-9]+(?:['-][a-zA-Z0-9]+)*", text.lower())
    return [t for t in tokens if t not in _STOPWORDS and len(t) > 1]


def _score_chunk(chunk: str, query_tokens: List[str]) -> float:
    """
    Score a chunk against query tokens.

    Score = (weighted token hits) / sqrt(chunk_word_count)
    Exact phrase match in the chunk gets a 3x bonus.
    """
    if not query_tokens:
        return 0.0

    chunk_lower = chunk.lower()
    chunk_tokens = _tokenize(chunk)
    chunk_len = max(len(chunk_tokens), 1)

    # Count how many query tokens appear in the chunk
    hits = sum(1 for qt in query_tokens if qt in chunk_lower)

    if hits == 0:
        return 0.0

    # Bonus for exact multi-word phrase match
    query_phrase = " ".join(query_tokens[:4])  # use first 4 tokens as phrase
    phrase_bonus = 3.0 if query_phrase in chunk_lower and len(query_tokens) >= 2 else 1.0

    # TF-like normalization: more hits is better, but penalize very long chunks
    score = (hits * phrase_bonus) / math.sqrt(chunk_len)
    return score


# ---------------------------------------------------------------------------
# Main tool function
# ---------------------------------------------------------------------------

def _search_literature(query: str, top_k: int = 5) -> Dict[str, Any]:
    """Search local chemistry documents for passages relevant to the query.

    Performs keyword-based search over all documents in the literature folder.
    Returns the most relevant passages with source citations.

    Use this when:
    - You need experimental procedures or specific conditions from a paper
    - The user references a specific reaction class and you want precedent
    - HTE database conditions exist but you want literature context/rationale
    - The user has uploaded domain-specific documents or lab procedures

    Args:
        query: Free-text search query (e.g. "Suzuki coupling hindered aryl chloride").
               Include reaction name, substrate class, or specific reagents/conditions.
        top_k: Number of top passages to return (default 5, max 10).

    Returns:
        dict with:
          found: number of matching passages
          query: the search query used
          results: list of {rank, score, excerpt, source, title, doi, year}
          docs_folder: path to the documents folder that was searched
          total_docs: total number of documents indexed
    """
    docs_dir = _get_docs_dir()

    try:
        all_chunks = _load_all_documents(docs_dir)
    except Exception as exc:
        return _error(f"Failed to load documents from {docs_dir}: {exc}")

    if not all_chunks:
        return _success({
            "found": 0,
            "query": query,
            "results": [],
            "docs_folder": str(docs_dir),
            "total_docs": 0,
            "note": (
                f"No documents found in {docs_dir}. "
                "Add .txt, .md, or .pdf files to that folder and they will be searchable."
            ),
        })

    # Count unique source documents
    unique_sources = len({meta["source"] for _, meta in all_chunks})

    query_tokens = _tokenize(query)
    if not query_tokens:
        return _error("Query is empty or contains only stopwords.")

    # Score all chunks
    scored: List[Tuple[float, str, Dict[str, Any]]] = []
    for chunk_text, meta in all_chunks:
        score = _score_chunk(chunk_text, query_tokens)
        if score > 0:
            scored.append((score, chunk_text, meta))

    # Sort by score descending
    scored.sort(key=lambda x: x[0], reverse=True)

    top_k = min(top_k, 10)
    top = scored[:top_k]

    if not top:
        return _success({
            "found": 0,
            "query": query,
            "results": [],
            "docs_folder": str(docs_dir),
            "total_docs": unique_sources,
            "note": f"No passages matched '{query}' across {unique_sources} document(s).",
        })

    results = []
    for rank, (score, chunk_text, meta) in enumerate(top, 1):
        # Truncate excerpt to ~600 chars for readability in the LLM context
        excerpt = chunk_text.strip()
        if len(excerpt) > 600:
            excerpt = excerpt[:600].rsplit(" ", 1)[0] + " [...]"

        entry: Dict[str, Any] = {
            "rank": rank,
            "score": round(score, 3),
            "excerpt": excerpt,
            "source": meta["source"],
            "title": meta["title"],
        }
        if meta.get("doi"):
            entry["doi"] = meta["doi"]
        if meta.get("year"):
            entry["year"] = meta["year"]

        results.append(entry)

    return _success({
        "found": len(results),
        "query": query,
        "results": results,
        "docs_folder": str(docs_dir),
        "total_docs": unique_sources,
    })


# ---------------------------------------------------------------------------
# ToolPlugin registration
# ---------------------------------------------------------------------------

search_literature_tool = ToolPlugin(
    name="search_literature",
    category="literature",
    description=(
        "Search local chemistry papers, procedures, and lab notes for relevant "
        "experimental conditions, mechanisms, or substrate-scope precedents. "
        "Add documents to the literature/ folder to make them searchable."
    ),
    prerequisites=[],  # no dependencies — can run in parallel with G0
    fn=_search_literature,
)


# ---------------------------------------------------------------------------
# Public API for agent._describe_resources()
# ---------------------------------------------------------------------------

def describe_literature_folder() -> str:
    """Return a brief description of the literature folder for the resource context."""
    docs_dir = _get_docs_dir()
    if not docs_dir.exists():
        return f"  • Literature folder ({docs_dir}): empty — add .txt/.md/.pdf files to enable search"

    supported = {".txt", ".md", ".pdf"}
    doc_files = [
        p for p in docs_dir.iterdir()
        if p.suffix.lower() in supported and not p.name.startswith(".")
    ]
    if not doc_files:
        return f"  • Literature folder ({docs_dir}): empty — add .txt/.md/.pdf files to enable search"

    names = [p.name for p in sorted(doc_files)[:5]]
    more = f" and {len(doc_files) - 5} more" if len(doc_files) > 5 else ""
    return (
        f"  • Literature folder: {len(doc_files)} document(s) — "
        f"{', '.join(names)}{more} "
        f"(searchable via search_literature tool)"
    )


# ---------------------------------------------------------------------------
# Web fetch tool
# ---------------------------------------------------------------------------

class _TextExtractor(HTMLParser):
    """Minimal HTML → plain text extractor using stdlib only."""

    # Tags whose content we skip entirely
    _SKIP_TAGS = {"script", "style", "noscript", "head", "meta", "link", "nav",
                  "footer", "aside", "iframe"}
    # Tags that should introduce a newline break
    _BLOCK_TAGS = {"p", "br", "div", "li", "h1", "h2", "h3", "h4", "h5", "h6",
                   "tr", "td", "th", "section", "article", "blockquote", "pre"}

    def __init__(self) -> None:
        super().__init__()
        self._parts: List[str] = []
        self._skip_depth = 0

    def handle_starttag(self, tag: str, attrs: Any) -> None:
        if tag in self._SKIP_TAGS:
            self._skip_depth += 1
        elif tag in self._BLOCK_TAGS and self._skip_depth == 0:
            self._parts.append("\n")

    def handle_endtag(self, tag: str) -> None:
        if tag in self._SKIP_TAGS:
            self._skip_depth = max(0, self._skip_depth - 1)
        elif tag in self._BLOCK_TAGS and self._skip_depth == 0:
            self._parts.append("\n")

    def handle_data(self, data: str) -> None:
        if self._skip_depth == 0:
            self._parts.append(data)

    def get_text(self) -> str:
        text = "".join(self._parts)
        # Unescape HTML entities (&amp; → &, &lt; → <, etc.)
        text = _html_module.unescape(text)
        # Collapse runs of whitespace/blank lines
        text = re.sub(r"[ \t]+", " ", text)
        text = re.sub(r"\n{3,}", "\n\n", text)
        return text.strip()


def _fetch_webpage(
    url: str,
    save_as: str = "",
    max_chars: int = 12000,
) -> Dict[str, Any]:
    """Fetch and extract text from a URL (webpage or direct text/HTML link).

    Useful for reading chemistry paper abstracts, procedure pages, supplier
    product pages, or any online resource the user provides. Returns the
    extracted text directly for use in the current query.

    Optionally saves the fetched text to the literature/ folder so it becomes
    permanently searchable via search_literature.

    Args:
        url: The URL to fetch (http or https).
        save_as: Optional filename (e.g. "suzuki_review.txt") to save the
                 fetched content into the literature/ folder for future
                 searches. Leave empty to use the fetched text only for
                 this query without saving.
        max_chars: Maximum characters of extracted text to return (default 12000).
                   Web pages can be large; this keeps the LLM context manageable.

    Returns:
        dict with:
          url: the fetched URL
          title: page <title> if found
          text: extracted plain text (up to max_chars)
          char_count: total extracted characters (before truncation)
          saved_as: filename if saved to literature folder, else ""
    """
    if not url.startswith(("http://", "https://")):
        return _error(f"Invalid URL: must start with http:// or https://. Got: {url!r}")

    # Fetch
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "Mozilla/5.0 (compatible; ChemCoworker/1.0; chemistry research assistant)"
            )
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            content_type = resp.headers.get("Content-Type", "").lower()
            raw_bytes = resp.read(1_000_000)  # cap at ~1 MB
    except urllib.error.HTTPError as exc:
        return _error(f"HTTP {exc.code} fetching {url}: {exc.reason}")
    except urllib.error.URLError as exc:
        return _error(f"Network error fetching {url}: {exc.reason}")
    except Exception as exc:
        return _error(f"Failed to fetch {url}: {exc}")

    # Extract plain text — branch by content type
    is_pdf = "application/pdf" in content_type or url.lower().endswith(".pdf")
    is_html = (not is_pdf) and ("html" in content_type or raw_bytes[:5] == b"<!DOC" or raw_bytes[:5] == b"<html")
    page_title = ""

    if is_pdf:
        # Extract text from PDF bytes directly
        extracted = _extract_pdf_text(source=url, raw_bytes=raw_bytes)
        if not extracted.strip():
            return _error(
                f"Fetched a PDF from {url} but could not extract text. "
                "Install pypdf or pdfminer.six: pip install pypdf"
            )
    elif is_html:
        # Decode bytes → str for HTML parsing
        charset = "utf-8"
        if "charset=" in content_type:
            charset = content_type.split("charset=")[-1].split(";")[0].strip()
        try:
            raw_text = raw_bytes.decode(charset, errors="replace")
        except (LookupError, ValueError):
            raw_text = raw_bytes.decode("utf-8", errors="replace")

        title_match = re.search(r"<title[^>]*>(.*?)</title>", raw_text, re.IGNORECASE | re.DOTALL)
        if title_match:
            page_title = _html_module.unescape(title_match.group(1).strip())

        extractor = _TextExtractor()
        try:
            extractor.feed(raw_text)
            extracted = extractor.get_text()
        except Exception:
            extracted = re.sub(r"<[^>]+>", " ", raw_text)
            extracted = _html_module.unescape(extracted)
            extracted = re.sub(r"\s+", " ", extracted).strip()
    else:
        # Plain text / markdown / other
        charset = "utf-8"
        if "charset=" in content_type:
            charset = content_type.split("charset=")[-1].split(";")[0].strip()
        try:
            extracted = raw_bytes.decode(charset, errors="replace").strip()
        except (LookupError, ValueError):
            extracted = raw_bytes.decode("utf-8", errors="replace").strip()

    total_chars = len(extracted)

    # Truncate
    if len(extracted) > max_chars:
        # Try to truncate at a paragraph boundary
        cut = extracted[:max_chars].rsplit("\n\n", 1)[0]
        extracted_display = cut + f"\n\n[... truncated at {max_chars} chars; {total_chars} total]"
    else:
        extracted_display = extracted

    # Optionally save to literature folder
    saved_as = ""
    if save_as:
        save_as = save_as.strip()
        # Sanitize filename
        save_as = re.sub(r"[^\w\-. ]", "_", save_as)
        if not save_as.endswith((".txt", ".md")):
            save_as += ".txt"

        docs_dir = _get_docs_dir()
        docs_dir.mkdir(parents=True, exist_ok=True)
        save_path = docs_dir / save_as

        try:
            if is_pdf:
                # Save raw PDF bytes + pre-extracted text cache
                pdf_path = docs_dir / re.sub(r"\.txt$", ".pdf", save_as)
                pdf_path.write_bytes(raw_bytes)
                txt_cache = pdf_path.with_suffix(pdf_path.suffix + ".txt")
                txt_cache.write_text(extracted, encoding="utf-8")
                save_path = pdf_path  # report the .pdf as the saved file
                save_as = pdf_path.name
            else:
                save_path.write_text(extracted, encoding="utf-8")

            # Write meta file
            meta = {"title": page_title or save_as, "url": url}
            meta_path = save_path.with_suffix(save_path.suffix + ".meta.json")
            meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
            # Invalidate cache so next search picks up the new file
            _DOC_CACHE.pop(str(save_path), None)
            saved_as = save_path.name
            logger.info(f"[literature] Saved {url} → {save_path}")
        except Exception as exc:
            logger.warning(f"[literature] Could not save {save_as}: {exc}")

    return _success({
        "url": url,
        "title": page_title,
        "text": extracted_display,
        "char_count": total_chars,
        "saved_as": saved_as,
    })


fetch_webpage_tool = ToolPlugin(
    name="fetch_webpage",
    category="literature",
    description=(
        "Fetch and extract text from a URL (paper page, supplier site, protocol page). "
        "Use when the user provides a URL or when a specific web resource is needed. "
        "Optionally save to the literature/ folder with save_as='filename.txt'."
    ),
    prerequisites=[],  # can run in parallel with any group
    fn=_fetch_webpage,
)


# ---------------------------------------------------------------------------
# read_literature_source — read a specific saved source file
# ---------------------------------------------------------------------------

def _read_literature_source(filename: str, max_chars: int = 8000) -> Dict[str, Any]:
    """Read the full text of a specific source file from the literature/ folder.

    Use this when:
    - Notes cite a source file and you need the full experimental procedure,
      exact quantities, or detail that was intentionally omitted from the notes
    - The user asks "what is the exact procedure / amounts / workup?"
    - You want to verify or expand on a condensed note

    Lookup order:
      1. Exact filename match (e.g. "demo.aspx_prep_v102p0086.txt")
      2. Partial name match (filename is contained in a stored name or vice versa)
      3. URL match — checks .meta.json sidecars for a matching url field

    Args:
        filename: Filename in the literature/ folder, or a URL that was used
                  to fetch the file originally.
        max_chars: Maximum characters to return (default 8000). Use a higher
                   value (e.g. 20000) for full procedures.

    Returns:
        dict with: filename, text, char_count, truncated, available_files (if not found)
    """
    docs_dir = _get_docs_dir()
    if not docs_dir.exists():
        return _error(f"Literature folder not found: {docs_dir}")

    supported = {".txt", ".md", ".pdf"}
    candidates = [
        p for p in sorted(docs_dir.iterdir())
        if p.suffix.lower() in supported and not p.name.startswith(".")
    ]

    target_path: Optional[Path] = None
    query_lower = filename.strip().lower()

    # 1. Exact match
    for p in candidates:
        if p.name.lower() == query_lower:
            target_path = p
            break

    # 2. Partial name match (query is substring of filename or vice versa)
    if target_path is None:
        for p in candidates:
            if query_lower in p.name.lower() or p.name.lower() in query_lower:
                target_path = p
                break

    # 3. URL match via .meta.json sidecars
    if target_path is None and ("http://" in filename or "https://" in filename):
        for p in candidates:
            meta_path = p.with_suffix(p.suffix + ".meta.json")
            if meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text(encoding="utf-8"))
                    if meta.get("url", "").rstrip("/") == filename.rstrip("/"):
                        target_path = p
                        break
                except Exception:
                    pass

    if target_path is None:
        available = [p.name for p in candidates[:20]]
        return _error(
            f"No file matching {filename!r} in literature/. "
            f"Available ({len(candidates)} files): {', '.join(available)}"
            + (f" ... and {len(candidates) - 20} more" if len(candidates) > 20 else "")
        )

    # Load content
    suffix = target_path.suffix.lower()
    try:
        if suffix == ".pdf":
            text = _load_pdf_file(target_path)
        else:
            text = target_path.read_text(encoding="utf-8", errors="replace")
    except Exception as exc:
        return _error(f"Could not read {target_path.name}: {exc}")

    if not text.strip():
        return _error(f"File {target_path.name} is empty or could not be parsed.")

    total_chars = len(text)
    truncated = total_chars > max_chars
    display_text = text[:max_chars] + f"\n\n[... truncated at {max_chars} chars; {total_chars} total ...]" if truncated else text

    # Load sidecar metadata if available
    meta_path = target_path.with_suffix(target_path.suffix + ".meta.json")
    meta: Dict[str, Any] = {}
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:
            pass

    return _success({
        "filename": target_path.name,
        "title": meta.get("title", target_path.name),
        "url": meta.get("url", ""),
        "text": display_text,
        "char_count": total_chars,
        "truncated": truncated,
    })


read_literature_source_tool = ToolPlugin(
    name="read_literature_source",
    category="literature",
    description=(
        "Read the full saved text of a specific source document from the literature/ folder. "
        "Use when notes cite a source file and you need the complete experimental procedure, "
        "exact quantities, or detail intentionally omitted from the extracted notes. "
        "Accepts filename, partial name, or original URL."
    ),
    prerequisites=[],
    fn=_read_literature_source,
)


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

LITERATURE_TOOLS = [
    search_literature_tool,
    fetch_webpage_tool,
    read_literature_source_tool,
]
