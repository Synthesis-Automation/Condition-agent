"""Recorded literature snapshots and exact excerpts, without chemistry inference.

Remote documents are untrusted evidence, never instructions. Capturing a source
does not establish that it supports a chemical claim or experimental procedure.
"""

from __future__ import annotations

import base64
import errno
import hashlib
import http.client
import io
import ipaddress
import re
import socket
import ssl
from datetime import datetime, timezone
from html.parser import HTMLParser
from time import monotonic
from typing import Any
from urllib.parse import urljoin, urlsplit, urlunsplit

from ..core.store import InvestigationEvent, InvestigationStore

SOURCE_SCHEMA = "literature_source.v1"
EXCERPT_SCHEMA = "literature_excerpt.v1"
MAX_SOURCE_BYTES = 8 * 1024 * 1024
MAX_TEXT_CHARACTERS = 2_000_000
MAX_PASSAGE_CHARACTERS = 16_000
FETCH_TIMEOUT_SECONDS = 30.0
MAX_REDIRECTS = 4


def _network_permission_denied(error: BaseException) -> bool:
    return (
        isinstance(error, PermissionError)
        or getattr(error, "winerror", None) == 10013
        or getattr(error, "errno", None) in {errno.EACCES, errno.EPERM, 10013}
    )


def _url(value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("A public HTTP(S) source URL is required")
    if any(ord(char) < 33 for char in value) or "\\" in value:
        raise ValueError("Source URLs cannot contain whitespace or control characters")
    parsed = urlsplit(value)
    if parsed.scheme not in {"http", "https"} or not parsed.hostname:
        raise ValueError("Only public HTTP(S) source URLs are supported")
    if parsed.username is not None or parsed.password is not None:
        raise ValueError("Source URLs cannot contain credentials")
    try:
        port = parsed.port
    except ValueError as exc:
        raise ValueError("Invalid source URL port") from exc
    if port == 0:
        raise ValueError("Invalid source URL port")
    hostname = parsed.hostname
    try:
        address = ipaddress.ip_address(hostname)
    except ValueError:
        if hostname.rstrip(".").lower() == "localhost" or hostname.lower().endswith(".localhost"):
            raise ValueError("Local source addresses are not allowed")
        hostname.encode("idna")
    else:
        if not address.is_global:
            raise ValueError("Only public source addresses are allowed")
    return urlunsplit((parsed.scheme, parsed.netloc, parsed.path or "/", parsed.query, ""))


def _label(value: str | None, field: str) -> str | None:
    if value is not None and (not isinstance(value, str) or len(value) > 2000):
        raise ValueError(f"{field} must be a string of at most 2000 characters")
    return value.strip() if value else None


def _public_addresses(host: str, port: int) -> list[tuple[Any, ...]]:
    addresses = socket.getaddrinfo(host, port, type=socket.SOCK_STREAM)
    if not addresses or any(not ipaddress.ip_address(item[4][0]).is_global for item in addresses):
        raise ValueError("Source DNS must resolve exclusively to public addresses")
    return addresses


def _request(url: str, *, timeout: float) -> dict[str, Any]:
    """Fetch one hop through pinned public DNS addresses, without system proxies."""
    parsed = urlsplit(_url(url))
    host = (parsed.hostname or "").encode("idna").decode("ascii")
    port = parsed.port or (443 if parsed.scheme == "https" else 80)
    deadline = monotonic() + timeout
    addresses = _public_addresses(host, port)

    def connect(address: Any, timeout: float, source_address: Any = None) -> socket.socket:
        last_error: OSError | None = None
        for family, socktype, protocol, _, sockaddr in addresses:
            remaining = deadline - monotonic()
            if remaining <= 0:
                raise TimeoutError("Literature retrieval time budget exhausted")
            connection = socket.socket(family, socktype, protocol)
            try:
                connection.settimeout(min(timeout, remaining))
                connection.connect(sockaddr)
                return connection
            except OSError as exc:
                connection.close()
                if _network_permission_denied(exc):
                    # A blocked socket is not a DNS-address selection problem.
                    # Preserve the denial instead of retrying equivalent sockets.
                    raise
                last_error = exc
        raise last_error or OSError("No public source address reachable")

    connection: http.client.HTTPConnection
    if parsed.scheme == "https":
        connection = http.client.HTTPSConnection(
            host, port, timeout=timeout, context=ssl.create_default_context(),
        )
    else:
        connection = http.client.HTTPConnection(host, port, timeout=timeout)
    # HTTPConnection uses this hook; HTTPS still verifies TLS against the original hostname.
    connection._create_connection = connect
    try:
        connection.request("GET", urlunsplit(("", "", parsed.path, parsed.query, "")), headers={
            "User-Agent": "ScientificWorkspace/1.0 (literature evidence capture)",
            "Accept": "text/html,application/pdf,text/plain,application/xhtml+xml",
            "Accept-Encoding": "identity",
        })
        transport_socket = connection.sock
        response = connection.getresponse()
        headers = {key.lower(): value for key, value in response.getheaders()}
        if response.status in {301, 302, 303, 307, 308}:
            return {"status": response.status, "headers": headers, "body": b""}
        length = headers.get("content-length")
        if length and int(length) > MAX_SOURCE_BYTES:
            raise ValueError("Source exceeds the 8 MiB retrieval limit")
        chunks: list[bytes] = []
        received = 0
        while True:
            remaining = deadline - monotonic()
            if remaining <= 0:
                raise TimeoutError("Literature retrieval time budget exhausted")
            if transport_socket is not None:
                transport_socket.settimeout(remaining)
            chunk = response.read1(min(65536, MAX_SOURCE_BYTES + 1 - received))
            if not chunk:
                break
            received += len(chunk)
            if received > MAX_SOURCE_BYTES:
                raise ValueError("Source exceeds the 8 MiB retrieval limit")
            chunks.append(chunk)
        return {"status": response.status, "headers": headers, "body": b"".join(chunks)}
    finally:
        connection.close()


class _DocumentText(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.parts: list[str] = []
        self.title_parts: list[str] = []
        self.hidden: list[str] = []
        self.in_title = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in {"script", "style", "noscript", "template"}:
            self.hidden.append(tag)
        if self.hidden:
            return
        if tag == "title":
            self.in_title = True
        if tag in {"p", "div", "section", "article", "br", "li", "tr", "h1", "h2", "h3", "h4"}:
            self.parts.append("\n")
        if tag in {"td", "th"}:
            self.parts.append("\t")

    def handle_endtag(self, tag: str) -> None:
        if self.hidden:
            if tag == self.hidden[-1]:
                self.hidden.pop()
            return
        if tag == "title":
            self.in_title = False
        if tag in {"p", "div", "section", "article", "li", "tr", "h1", "h2", "h3", "h4"}:
            self.parts.append("\n")

    def handle_data(self, data: str) -> None:
        if not self.hidden:
            if self.in_title:
                self.title_parts.append(data)
            else:
                self.parts.append(re.sub(r"\s+", " ", data))


def _extract(body: bytes, content_type: str) -> dict[str, Any]:
    media_type = content_type.split(";", 1)[0].strip().lower()
    if media_type == "application/pdf" or body.startswith(b"%PDF-"):
        return _extract_pdf(body)
    if media_type not in {"text/html", "application/xhtml+xml", "text/plain"}:
        return {"status": "unsupported_media_type", "text": "", "pages": []}
    match = re.search(r"charset\s*=\s*[\"']?([^;\s\"']+)", content_type, re.I)
    encoding = match.group(1) if match else "utf-8"
    try:
        decoded = body.decode(encoding, errors="replace")
    except LookupError:
        return {"status": "unsupported_charset", "text": "", "pages": []}
    title = None
    if media_type in {"text/html", "application/xhtml+xml"}:
        parser = _DocumentText()
        parser.feed(decoded)
        text = "\n".join(line.strip() for line in "".join(parser.parts).splitlines() if line.strip())
        title = " ".join("".join(parser.title_parts).split()) or None
    else:
        text = decoded.replace("\r\n", "\n").replace("\r", "\n")
    truncated = len(text) > MAX_TEXT_CHARACTERS
    return {
        "status": "partial_text_limit" if truncated else ("completed" if text.strip() else "empty"),
        "text": text[:MAX_TEXT_CHARACTERS], "pages": [], "document_title": title,
        "encoding": encoding, "decoding_replacement_characters": decoded.count("\ufffd"),
        "method": "stdlib_html_parser" if media_type != "text/plain" else "text_decode",
        "limitations": ["Text extraction omits molecular drawings and does not validate source claims."],
    }


def _extract_pdf(body: bytes) -> dict[str, Any]:
    try:
        import pypdf
    except ImportError:
        return {"status": "pdf_parser_unavailable", "text": "", "pages": [],
                "limitations": ["Install requirements-web.txt for PDF text extraction; OCR is not provided."]}
    parser = {"name": "pypdf", "version": getattr(pypdf, "__version__", "unreported")}
    try:
        reader = pypdf.PdfReader(io.BytesIO(body), strict=False)
        if reader.is_encrypted:
            return {"status": "encrypted_pdf", "text": "", "pages": [], "parser": parser}
        pages, parts, cursor = [], [], 0
        truncated = len(reader.pages) > 200
        for index, page in enumerate(reader.pages):
            if index >= 200 or cursor >= MAX_TEXT_CHARACTERS:
                truncated = True
                break
            text = page.extract_text() or ""
            if len(text) > MAX_TEXT_CHARACTERS - cursor:
                text = text[:MAX_TEXT_CHARACTERS - cursor]
                truncated = True
            parts.append(text)
            pages.append({"page": index + 1, "start": cursor, "end": cursor + len(text)})
            cursor += len(text) + 1
        text = "\n".join(parts)
        return {
            "status": "partial_text_limit" if truncated else (
                "completed" if text.strip() else "no_extractable_text_ocr_required"
            ), "text": text, "pages": pages, "method": "pypdf_text_extraction", "parser": parser,
            "limitations": ["No OCR; molecular drawings, table structure, and reading order may be lost."],
        }
    except Exception as exc:
        return {"status": "pdf_extraction_failed", "text": "", "pages": [],
                "parser": parser,
                "error": f"{type(exc).__name__}: {str(exc)[:500]}"}


def _source_record(
    url: str, title: str | None, acquisition: str, reference_id: str | None = None,
) -> dict[str, Any]:
    record = {
        "schema_version": SOURCE_SCHEMA, "source_url": url, "title": title,
        "captured_at": datetime.now(timezone.utc).isoformat(), "acquisition": acquisition,
        "origin": "external_literature", "review_status": "unreviewed",
        "claim_support": "not_assessed", "experimental_verification": "not_performed",
    }
    if reference_id is not None:
        if not isinstance(reference_id, str) or not re.fullmatch(r"REF1:[0-9a-f]{64}", reference_id):
            raise ValueError("reference_id must be an explicit REF1 publication identity")
        record["reported_reference_id"] = reference_id
        record["bibliography_attribution"] = "agent_supplied_not_independently_verified"
    return record


def fetch_source(
    store: InvestigationStore, url: str, *, title: str | None = None, retry_network: bool = False,
    reference_id: str | None = None,
) -> InvestigationEvent:
    """Record a bounded public HTTP(S) snapshot, extracted text, and explicit failures."""
    current = _url(url)
    record = _source_record(current, _label(title, "title"), "http_fetch", reference_id)
    record.update({"retrieval_status": "failed", "redirects": [], "snapshot": None, "final_url": current,
                   "extraction": {"status": "not_attempted", "text": "", "pages": []}})
    if type(retry_network) is not bool:
        raise ValueError("retry_network must be a boolean")
    if not retry_network:
        for event in reversed(store.events()):
            if event.kind != "literature_source":
                continue
            previous = store.read_artifact(event.artifact_ref)
            if previous.get("acquisition") != "http_fetch":
                continue
            if previous.get("retrieval_status") == "completed":
                break
            if previous.get("error", {}).get("category") == "network_permission_denied":
                record.update({"error": dict(previous["error"]), "network_request_skipped": True,
                               "blocked_by_ref": event.artifact_ref,
                               "recovery": {**previous.get("recovery", {}), "source_url": current},
                               "text_sha256": hashlib.sha256(b"").hexdigest(), "text_characters": 0})
                return store.append("literature_source", record, evidence_refs=(event.artifact_ref,))
    deadline = monotonic() + FETCH_TIMEOUT_SECONDS
    try:
        for hop in range(MAX_REDIRECTS + 1):
            remaining = deadline - monotonic()
            if remaining <= 0:
                raise TimeoutError("Literature retrieval time budget exhausted")
            response = _request(current, timeout=remaining)
            status, headers, body = response["status"], response["headers"], response["body"]
            record.update({"final_url": current, "http_status": status})
            if status in {301, 302, 303, 307, 308}:
                if hop == MAX_REDIRECTS:
                    raise ValueError("Literature source exceeded redirect limit")
                if not headers.get("location"):
                    raise ValueError("Literature redirect has no destination")
                destination = _url(urljoin(current, headers["location"]))
                record["redirects"].append({"from": current, "to": destination, "status": status})
                current = destination
                continue
            if not 200 <= status < 300:
                raise ValueError(f"Literature source returned HTTP {status}")
            if len(body) > MAX_SOURCE_BYTES:
                raise ValueError("Source exceeds the 8 MiB retrieval limit")
            record["retrieval_status"] = "completed"
            record["snapshot"] = {
                "encoding": "base64", "content": base64.b64encode(body).decode("ascii"),
                "sha256": hashlib.sha256(body).hexdigest(), "bytes": len(body),
                "content_type": headers.get("content-type", ""),
                "content_encoding": headers.get("content-encoding", "identity"),
                "etag": headers.get("etag"), "last_modified": headers.get("last-modified"),
            }
            if headers.get("content-encoding", "identity").lower() not in {"", "identity"}:
                record["extraction"] = {"status": "unsupported_content_encoding", "text": "", "pages": []}
            else:
                record["extraction"] = _extract(body, headers.get("content-type", ""))
            record["title"] = record["title"] or record["extraction"].get("document_title")
            break
    except (OSError, ValueError, http.client.HTTPException) as exc:
        record["error"] = {"type": type(exc).__name__, "message": str(exc)[:1000]}
        permission_denied = _network_permission_denied(exc)
        record["error"]["category"] = (
            "network_permission_denied" if permission_denied else "source_retrieval_failed"
        )
        record["recovery"] = {
            "action": "capture_browser_source", "source_url": current,
            "instructions": (
                "If an available browser or web tool can open this public source, copy its "
                "visible text exactly and call w.capture_source(text, url=source_url, "
                "title=title, locator=locator). Keep this failed artifact as provenance. "
                "The captured text remains an unverified agent-supplied excerpt. "
                "If the source cannot be opened, report the evidence gap."
            ),
            "retry_guidance": (
                "Do not repeat direct fetches while this network-permission restriction persists. "
                "Use the available browser/web tool; do not change sandbox permissions."
                if permission_denied else
                "Retry only when there is a concrete reason the source or connection has changed."
            ),
        }
    text = record["extraction"]["text"]
    record["text_sha256"] = hashlib.sha256(text.encode("utf-8")).hexdigest()
    record["text_characters"] = len(text)
    return store.append("literature_source", record)


def _captured_record(
    store: InvestigationStore, text: str, *, url: str, title: str | None = None,
    locator: str | None = None, reference_id: str | None = None,
) -> dict[str, Any]:
    """Save browser-obtained text honestly as an unverified agent-supplied excerpt."""
    url = _url(url)
    if not isinstance(text, str) or not text.strip() or len(text) > MAX_TEXT_CHARACTERS:
        raise ValueError("Captured source text must be nonempty and at most 2000000 characters")
    record = _source_record(url, _label(title, "title"), "agent_supplied_excerpt", reference_id)
    record.update({
        "retrieval_status": "not_performed", "final_url": None,
        "reported_locator": _label(locator, "locator"), "snapshot": None,
        "text_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "text_characters": len(text),
        "extraction": {"status": "agent_supplied", "text": text, "pages": [],
                       "method": "agent_supplied_text", "limitations": [
                           "URL, completeness, transcription, and reported locator have not been HTTP-verified.",
                       ]},
    })
    if reference_id is not None:
        record["extraction"]["limitations"].append(
            "Publication identity is agent-attributed; independent bibliographic verification is absent."
        )
    return record


def capture_source(
    store: InvestigationStore, text: str, *, url: str, title: str | None = None,
    locator: str | None = None, reference_id: str | None = None,
) -> InvestigationEvent:
    """Save exact browser text with optional explicit, unverified publication attribution."""
    record = _captured_record(store, text, url=url, title=title, locator=locator, reference_id=reference_id)
    return store.append("literature_source", record)


def capture_source_file(
    store: InvestigationStore, file: str, *, url: str, title: str | None = None,
    locator: str | None = None, reference_id: str | None = None,
    text_path: list[str | int] | None = None,
) -> InvestigationEvent:
    """Import a bounded UTF-8 text or JSON browser export without console transcription.

    Files must be inside this investigation. A JSON string needs text_path=[];
    structured JSON needs literal keys/indices. Bytes and selected path are hashed
    as acquisition provenance, never promoted to an independently fetched source.
    """
    from pathlib import Path
    import json

    path = (store.root / file).resolve()
    if not path.is_relative_to(store.root) or not path.is_file():
        raise ValueError("Source export must be a file inside this investigation")
    if path.stat().st_size > MAX_SOURCE_BYTES:
        raise ValueError("Source export exceeds the 8 MiB limit")
    body = path.read_bytes()
    if len(body) > MAX_SOURCE_BYTES:
        raise ValueError("Source export exceeds the 8 MiB limit")
    text = body.decode("utf-8-sig")
    if text_path is None and path.suffix.lower() == ".json":
        decoded = json.loads(text)
        if not isinstance(decoded, str):
            raise ValueError("Structured JSON exports require an explicit text_path")
        text_path = []
    if text_path is not None:
        if (not isinstance(text_path, list) or len(text_path) > 12 or any(
                not (isinstance(part, str) and len(part) <= 80 or type(part) is int and part >= 0)
                for part in text_path)):
            raise ValueError("text_path must contain at most 12 literal keys or nonnegative indices")
        text = json.loads(text)
        for part in text_path:
            if isinstance(text, dict) and isinstance(part, str):
                text = text[part]
            elif isinstance(text, list) and type(part) is int:
                text = text[part]
            else:
                raise ValueError("text_path does not select a JSON text field")
    record = _captured_record(store, text, url=url, title=title, locator=locator, reference_id=reference_id)
    record["capture_file"] = {"path": path.relative_to(store.root).as_posix(),
                              "sha256": hashlib.sha256(body).hexdigest(), "bytes": len(body),
                              "text_path": text_path}
    return store.append("literature_source", record)


def _load_source(store: InvestigationStore, reference: str) -> dict[str, Any]:
    source = store.read_artifact(reference)
    if source.get("schema_version") != SOURCE_SCHEMA or not any(
        event.kind == "literature_source" and event.artifact_ref == reference for event in store.events()
    ):
        raise ValueError("Expected a recorded literature source reference")
    return source


def load_source_passage(store: InvestigationStore, reference: str) -> tuple[dict[str, Any], str, str]:
    """Verify source/excerpt lineage and return its literal captured text and root reference."""
    kinds = {event.artifact_ref: event.kind for event in store.events()}
    kind = kinds.get(reference)
    if kind not in {"literature_source", "literature_excerpt"}:
        raise ValueError("Expected a recorded literature source or exact excerpt")
    value = store.read_artifact(reference)
    root_ref = reference if kind == "literature_source" else value.get("source_ref")
    source = _load_source(store, root_ref)
    text = source.get("extraction", {}).get("text")
    if (not isinstance(text, str) or not text.strip() or not (
            source.get("acquisition") == "http_fetch" and source.get("retrieval_status") == "completed"
            or source.get("acquisition") == "agent_supplied_excerpt" and source.get("retrieval_status") == "not_performed")):
        raise ValueError("Captured source text is required; failed fetches are debugging records. Use capture_source.")
    if kind == "literature_excerpt":
        start, end = value.get("location", {}).get("start"), value.get("location", {}).get("end")
        if (value.get("schema_version") != EXCERPT_SCHEMA
                or type(start) is not int or type(end) is not int or not 0 <= start < end <= len(text)
                or value.get("text") != text[start:end]
                or value.get("verification") != "exact_match_to_captured_text"):
            raise ValueError("Literature excerpt must match its captured source exactly")
        text = text[start:end]
    return source, text, root_ref


def _location(source: dict[str, Any], start: int, end: int) -> dict[str, Any]:
    text = source["extraction"]["text"]
    return {
        "start": start, "end": end, "line_start": text.count("\n", 0, start) + 1,
        "line_end": text.count("\n", 0, max(start, end - 1)) + 1,
        "pages": [page["page"] for page in source["extraction"].get("pages", [])
                  if page["start"] < end and page["end"] > start],
        "coordinate_system": "zero_based_characters_end_exclusive_in_captured_text",
    }


def inspect_source(
    store: InvestigationStore, source_ref: str, *, query: str | None = None,
    offset: int = 0, limit: int = 4000,
) -> dict[str, Any]:
    """Read a bounded passage; optional case-insensitive search starts at offset."""
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= MAX_PASSAGE_CHARACTERS:
        raise ValueError("offset must be nonnegative and limit must be between 1 and 16000")
    if query is not None and (not isinstance(query, str) or not query.strip() or len(query) > 1000):
        raise ValueError("query must contain 1 to 1000 characters")
    source = _load_source(store, source_ref)
    text = source["extraction"]["text"]
    if offset > len(text):
        raise ValueError("offset is beyond the captured text")
    start = offset
    found = None
    match_start = None
    if query:
        match = re.compile(re.escape(query), re.IGNORECASE).search(text, offset)
        found = match is not None
        if match:
            match_start = match.start()
            start = max(offset, match.start() - min(240, limit // 4))
    end = min(start + limit, len(text)) if found is not False else start
    return {
        "source_ref": source_ref, "source_url": source["source_url"], "title": source["title"],
        "acquisition": source["acquisition"], "retrieval_status": source["retrieval_status"],
        "extraction_status": source["extraction"]["status"], "query_found": found,
        "query_match_start": match_start, "text": text[start:end],
        "location": _location(source, start, end), "total_characters": len(text),
        "next_offset": end if end < len(text) and found is not False else None,
        "reported_locator": source.get("reported_locator"),
        "limitations": source["extraction"].get("limitations", []), "claim_support": "not_assessed",
        "error": source.get("error"), "recovery": source.get("recovery"),
    }


def record_source_excerpt(
    store: InvestigationStore, source_ref: str, *, start: int | None = None,
    end: int | None = None, excerpt: str | None = None, locator: str | None = None,
) -> InvestigationEvent:
    """Save an exact selection from a snapshot; never certify its chemical claims."""
    source = _load_source(store, source_ref)
    text = source["extraction"]["text"]
    if excerpt is not None:
        if start is not None or end is not None:
            raise ValueError("Choose character offsets or an exact excerpt, not both")
        if not isinstance(excerpt, str) or not excerpt.strip() or len(excerpt) > MAX_PASSAGE_CHARACTERS:
            raise ValueError("Excerpt must contain 1 to 16000 characters")
        start = text.find(excerpt)
        if start < 0:
            raise ValueError("Excerpt does not occur verbatim in the captured source")
        if text.find(excerpt, start + 1) >= 0:
            raise ValueError("Excerpt occurs more than once; select explicit character offsets")
        end = start + len(excerpt)
    if type(start) is not int or type(end) is not int or not 0 <= start < end <= len(text):
        raise ValueError("Select valid start/end character offsets within the captured source")
    if end - start > MAX_PASSAGE_CHARACTERS or not text[start:end].strip():
        raise ValueError("Excerpt must contain 1 to 16000 characters")
    return store.append("literature_excerpt", {
        "schema_version": EXCERPT_SCHEMA, "source_ref": source_ref,
        "source_url": source["source_url"], "final_url": source.get("final_url"),
        "title": source["title"], "acquisition": source["acquisition"],
        "captured_at": source["captured_at"], "text": text[start:end],
        "location": _location(source, start, end), "reported_locator": _label(locator, "locator"),
        "verification": "exact_match_to_captured_text", "claim_support": "not_assessed",
        "origin": "external_literature", "review_status": "unreviewed",
    }, evidence_refs=(source_ref,))
