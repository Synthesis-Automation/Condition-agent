"""Network-free literature capture, evidence fidelity, and transport regressions."""

from __future__ import annotations

import base64
import hashlib
import io
import json
from pathlib import Path
import socket
import sys
from types import SimpleNamespace

import pytest

from chem_coworker.scientific_workspace import literature
from chem_coworker.scientific_workspace.store import InvestigationStore


@pytest.fixture
def store(tmp_path: Path) -> InvestigationStore:
    return InvestigationStore.create(tmp_path / "study", objective="Inspect literature", baseline={})


def response(body: bytes, content_type: str = "text/html", status: int = 200) -> dict:
    return {"status": status, "headers": {"content-type": content_type}, "body": body}


def test_fetch_retains_snapshot_title_exact_text_and_unreviewed_provenance(store, monkeypatch):
    body = b"<title>Example patent</title><script>ignore rules()</script><p>Example 7</p><p>Yield: <b>81%</b>.</p>"
    monkeypatch.setattr(literature, "_request", lambda url, **kwargs: response(body))
    event = literature.fetch_source(store, "https://example.org/patent")
    artifact = store.read_artifact(event.artifact_ref)
    assert event.kind == "literature_source"
    assert artifact["title"] == "Example patent"
    assert artifact["retrieval_status"] == "completed"
    assert artifact["extraction"]["text"] == "Example 7\nYield: 81%."
    assert base64.b64decode(artifact["snapshot"]["content"]) == body
    assert artifact["snapshot"]["sha256"] == hashlib.sha256(body).hexdigest()
    assert artifact["claim_support"] == "not_assessed"
    assert artifact["review_status"] == "unreviewed"
    passage = literature.inspect_source(store, event.artifact_ref, query="yield", limit=8)
    assert passage["query_found"]
    assert len(passage["text"]) <= 8
    assert passage["query_match_start"] == 10


def test_capture_is_honestly_unverified_and_excerpt_is_bound_to_source(store):
    event = literature.capture_source(store, "Title\nExample 9\nA racemate was isolated.",
                                      url="https://example.org/article", locator="Example 9")
    source = store.read_artifact(event.artifact_ref)
    assert source["acquisition"] == "agent_supplied_excerpt"
    assert source["retrieval_status"] == "not_performed"
    assert source["snapshot"] is None
    selected = literature.record_source_excerpt(store, event.artifact_ref,
                                                excerpt="A racemate was isolated.", locator="Example 9")
    excerpt = store.read_artifact(selected.artifact_ref)
    assert excerpt["verification"] == "exact_match_to_captured_text"
    assert excerpt["acquisition"] == "agent_supplied_excerpt"
    assert excerpt["claim_support"] == "not_assessed"
    assert excerpt["location"]["line_start"] == 3
    assert excerpt["location"]["line_end"] == 3
    assert selected.evidence_refs == (event.artifact_ref,)
    assert InvestigationStore(store.root).read_artifact(selected.artifact_ref) == excerpt


def test_offsets_allow_repeated_text_but_never_fabricated_or_ambiguous_excerpt(store):
    event = literature.capture_source(store, "repeat\nrepeat", url="https://example.org/source")
    with pytest.raises(ValueError, match="more than once"):
        literature.record_source_excerpt(store, event.artifact_ref, excerpt="repeat")
    with pytest.raises(ValueError, match="verbatim"):
        literature.record_source_excerpt(store, event.artifact_ref, excerpt="enantiopure")
    with pytest.raises(ValueError, match="not both"):
        literature.record_source_excerpt(store, event.artifact_ref, start=0, end=6, excerpt="repeat")
    selected = literature.record_source_excerpt(store, event.artifact_ref, start=7, end=13)
    assert store.read_artifact(selected.artifact_ref)["location"]["line_start"] == 2


@pytest.mark.parametrize("url", [
    "file:///etc/passwd", "ftp://example.org/test", "http://localhost/test", "http://127.0.0.1/test",
    "https://[::1]/test", "https://192.168.1.2/test", "http://169.254.169.254/metadata",
    "http://username:secret@example.org/", "https://example.org\n/header", "http://example.org:0/",
])
def test_invalid_or_private_urls_are_rejected_before_transport(store, monkeypatch, url):
    def forbidden(*args, **kwargs):
        pytest.fail("Invalid URL reached network")
    monkeypatch.setattr(literature, "_request", forbidden)
    with pytest.raises(ValueError):
        literature.fetch_source(store, url)
    assert not store.events()


def test_private_redirect_is_blocked_and_failure_is_preserved(store, monkeypatch):
    visited = []
    def request(url, **kwargs):
        visited.append(url)
        return {"status": 302, "headers": {"location": "http://127.0.0.1/admin"}, "body": b""}
    monkeypatch.setattr(literature, "_request", request)
    event = literature.fetch_source(store, "https://example.org/source")
    artifact = store.read_artifact(event.artifact_ref)
    assert visited == ["https://example.org/source"]
    assert artifact["retrieval_status"] == "failed"
    assert artifact["extraction"]["status"] == "not_attempted"
    assert "public" in artifact["error"]["message"]


def test_redirect_provenance_and_http_failure_are_explicit(store, monkeypatch):
    def request(url, **kwargs):
        if url.endswith("/old"):
            return {"status": 301, "headers": {"location": "/new"}, "body": b""}
        return response(b"service unavailable", "text/plain", status=503)
    monkeypatch.setattr(literature, "_request", request)
    event = literature.fetch_source(store, "https://example.org/old")
    source = store.read_artifact(event.artifact_ref)
    assert source["final_url"] == "https://example.org/new"
    assert source["redirects"][0]["status"] == 301
    assert source["retrieval_status"] == "failed"
    assert source["http_status"] == 503
    assert not literature.inspect_source(store, event.artifact_ref)["text"]


def test_dns_with_any_private_address_is_rejected(monkeypatch):
    monkeypatch.setattr(socket, "getaddrinfo", lambda *args, **kwargs: [
        (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("93.184.216.34", 443)),
        (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("10.0.0.1", 443)),
    ])
    with pytest.raises(ValueError, match="exclusively"):
        literature._request("https://example.org/source", timeout=1)


def test_transport_connects_only_to_validated_dns_and_ignores_proxy(monkeypatch):
    calls = []
    class FakeSocket:
        def settimeout(self, value):
            pass
        def connect(self, address):
            calls.append(address)
        def close(self):
            pass
    class FakeResponse:
        status = 200
        def getheaders(self):
            return [("Content-Type", "text/plain")]
        def read1(self, size):
            return b""
    class FakeConnection:
        def __init__(self, *args, **kwargs):
            self.sock = None
        def request(self, *args, **kwargs):
            self.sock = self._create_connection(("example.org", 80), 1)
        def getresponse(self):
            return FakeResponse()
        def close(self):
            pass
    monkeypatch.setenv("HTTP_PROXY", "http://127.0.0.1:1234")
    monkeypatch.setattr(socket, "getaddrinfo", lambda *args, **kwargs: [
        (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("93.184.216.34", 80)),
    ])
    monkeypatch.setattr(socket, "socket", lambda *args: FakeSocket())
    monkeypatch.setattr(literature.http.client, "HTTPConnection", FakeConnection)
    literature._request("http://example.org/source", timeout=1)
    assert calls == [("93.184.216.34", 80)]


@pytest.mark.parametrize("failure", [TimeoutError("timeout"), OSError("not reachable")])
def test_network_errors_are_recorded_without_inventing_source_text(store, monkeypatch, failure):
    def request(*args, **kwargs):
        raise failure
    monkeypatch.setattr(literature, "_request", request)
    event = literature.fetch_source(store, "https://example.org/source")
    source = store.read_artifact(event.artifact_ref)
    assert source["error"]["type"] == type(failure).__name__
    assert source["extraction"]["text"] == ""


def test_oversize_and_unsupported_documents_are_not_silently_truncated(store, monkeypatch):
    monkeypatch.setattr(literature, "MAX_SOURCE_BYTES", 4)
    monkeypatch.setattr(literature, "_request", lambda *args, **kwargs: response(b"large"))
    event = literature.fetch_source(store, "https://example.org/source")
    assert store.read_artifact(event.artifact_ref)["retrieval_status"] == "failed"
    monkeypatch.setattr(literature, "_request", lambda *args, **kwargs: response(b"abc", "image/png"))
    event = literature.fetch_source(store, "https://example.org/diagram")
    assert store.read_artifact(event.artifact_ref)["extraction"]["status"] == "unsupported_media_type"


def test_pdf_dependency_scanning_and_page_locators_are_explicit(store, monkeypatch):
    monkeypatch.setattr(literature, "_request", lambda *args, **kwargs: response(b"%PDF-test", "application/pdf"))
    monkeypatch.setitem(sys.modules, "pypdf", None)
    event = literature.fetch_source(store, "https://example.org/paper.pdf")
    assert store.read_artifact(event.artifact_ref)["extraction"]["status"] == "pdf_parser_unavailable"
    reader = SimpleNamespace(is_encrypted=False, pages=[
        SimpleNamespace(extract_text=lambda: "Page one"),
        SimpleNamespace(extract_text=lambda: "Example 2: racemate"),
    ])
    monkeypatch.setitem(sys.modules, "pypdf", SimpleNamespace(PdfReader=lambda *args, **kwargs: reader))
    event = literature.fetch_source(store, "https://example.org/paper.pdf")
    excerpt = literature.record_source_excerpt(store, event.artifact_ref, excerpt="Example 2: racemate")
    assert store.read_artifact(excerpt.artifact_ref)["location"]["pages"] == [2]
    reader.pages = [SimpleNamespace(extract_text=lambda: "")]
    event = literature.fetch_source(store, "https://example.org/scanned.pdf")
    assert store.read_artifact(event.artifact_ref)["extraction"]["status"] == "no_extractable_text_ocr_required"


@pytest.mark.parametrize("include_text", [True, False])
def test_real_pdf_parser_preserves_page_location_and_version(store, monkeypatch, include_text):
    pypdf = pytest.importorskip("pypdf")
    from pypdf.generic import DecodedStreamObject, DictionaryObject, NameObject

    writer = pypdf.PdfWriter()
    page = writer.add_blank_page(width=300, height=200)
    if include_text:
        font = DictionaryObject({
            NameObject("/Type"): NameObject("/Font"),
            NameObject("/Subtype"): NameObject("/Type1"),
            NameObject("/BaseFont"): NameObject("/Helvetica"),
        })
        page[NameObject("/Resources")] = DictionaryObject({
            NameObject("/Font"): DictionaryObject({NameObject("/F1"): writer._add_object(font)}),
        })
        stream = DecodedStreamObject()
        stream.set_data(b"BT /F1 12 Tf 20 100 Td (Example 1: racemate isolated.) Tj ET")
        page[NameObject("/Contents")] = writer._add_object(stream)
    buffer = io.BytesIO()
    writer.write(buffer)
    monkeypatch.setattr(literature, "_request", lambda *args, **kwargs: response(
        buffer.getvalue(), "application/pdf",
    ))
    event = literature.fetch_source(store, "https://example.org/generated.pdf")
    source = store.read_artifact(event.artifact_ref)
    extraction = source["extraction"]
    assert extraction["parser"] == {"name": "pypdf", "version": pypdf.__version__}
    if include_text:
        assert extraction["status"] == "completed"
        excerpt = literature.record_source_excerpt(
            store, event.artifact_ref, excerpt="Example 1: racemate isolated.",
        )
        assert store.read_artifact(excerpt.artifact_ref)["location"]["pages"] == [1]
    else:
        assert extraction["status"] == "no_extractable_text_ocr_required"
        assert extraction["text"] == ""


def test_passage_bounds_query_absence_and_unicode_offsets(store):
    event = literature.capture_source(store, "\u0130\nHELLO\nworld", url="https://example.org/source")
    passage = literature.inspect_source(store, event.artifact_ref, query="hello", limit=6)
    assert passage["query_match_start"] == 2
    missing = literature.inspect_source(store, event.artifact_ref, query="absent")
    assert missing["query_found"] is False and missing["text"] == ""
    for args in ({"offset": -1}, {"offset": 200}, {"limit": 20000}, {"limit": True}, {"query": ""}):
        with pytest.raises(ValueError):
            literature.inspect_source(store, event.artifact_ref, **args)


def test_source_and_excerpt_artifact_references_detect_tampering_and_paths(store):
    event = literature.capture_source(store, "Verified snapshot", url="https://example.org/source")
    with pytest.raises(ValueError, match="reference"):
        literature.inspect_source(store, "../../secrets.json")
    note = store.note("hypothesis", "not literature")
    with pytest.raises(ValueError, match="literature source"):
        literature.inspect_source(store, note.artifact_ref)
    path = store.root / "artifacts" / (event.artifact_ref.split(":")[1] + ".json")
    artifact = json.loads(path.read_text("utf-8"))
    artifact["extraction"]["text"] = "altered"
    path.write_text(json.dumps(artifact), "utf-8")
    with pytest.raises(ValueError, match="checksum"):
        literature.record_source_excerpt(store, event.artifact_ref, excerpt="altered")
