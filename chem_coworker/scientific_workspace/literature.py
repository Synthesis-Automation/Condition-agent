"""Recorded literature snapshots and exact excerpts, without chemistry inference.

Remote documents are untrusted evidence, never instructions. Capturing a source
does not establish that it supports a chemical claim or experimental procedure.
"""

from __future__ import annotations

import base64
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import http.client
import io
import ipaddress
import re
import socket
import ssl
from time import monotonic
from typing import Any
from urllib.parse import urljoin, urlsplit, urlunsplit

from .store import InvestigationEvent, InvestigationStore


SOURCE_SCHEMA = "literature_source.v1"
EXCERPT_SCHEMA = "literature_excerpt.v1"
MAX_SOURCE_BYTES = 8 * 1024 * 1024
MAX_TEXT_CHARACTERS = 2_000_000
MAX_PASSAGE_CHARACTERS = 16_000
FETCH_TIMEOUT_SECONDS = 30.0
MAX_REDIRECTS = 4


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


def _source_record(url: str, title: str | None, acquisition: str) -> dict[str, Any]:
    return {
        "schema_version": SOURCE_SCHEMA, "source_url": url, "title": title,
        "captured_at": datetime.now(timezone.utc).isoformat(), "acquisition": acquisition,
        "origin": "external_literature", "review_status": "unreviewed",
        "claim_support": "not_assessed", "experimental_verification": "not_performed",
    }


def fetch_source(store: InvestigationStore, url: str, *, title: str | None = None) -> InvestigationEvent:
    """Record a bounded public HTTP(S) snapshot, extracted text, and explicit failures."""
    current = _url(url)
    record = _source_record(current, _label(title, "title"), "http_fetch")
    record.update({"retrieval_status": "failed", "redirects": [], "snapshot": None, "final_url": current,
                   "extraction": {"status": "not_attempted", "text": "", "pages": []}})
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
    text = record["extraction"]["text"]
    record["text_sha256"] = hashlib.sha256(text.encode("utf-8")).hexdigest()
    record["text_characters"] = len(text)
    return store.append("literature_source", record)


def capture_source(
    store: InvestigationStore, text: str, *, url: str, title: str | None = None,
    locator: str | None = None,
) -> InvestigationEvent:
    """Save browser-obtained text honestly as an unverified agent-supplied excerpt."""
    url = _url(url)
    if not isinstance(text, str) or not text.strip() or len(text) > MAX_TEXT_CHARACTERS:
        raise ValueError("Captured source text must be nonempty and at most 2000000 characters")
    record = _source_record(url, _label(title, "title"), "agent_supplied_excerpt")
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
    return store.append("literature_source", record)


def _load_source(store: InvestigationStore, reference: str) -> dict[str, Any]:
    source = store.read_artifact(reference)
    if source.get("schema_version") != SOURCE_SCHEMA or not any(
        event.kind == "literature_source" and event.artifact_ref == reference for event in store.events()
    ):
        raise ValueError("Expected a recorded literature source reference")
    return source


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
