"""Saved-answer formatting, molecular SVGs and untrusted Markdown regressions."""

from __future__ import annotations

import base64
from copy import deepcopy
from html.parser import HTMLParser
import xml.etree.ElementTree as ET

from fastapi.testclient import TestClient

from app.web_api.main import create_app
from app.web_api.scientific_presentation import present_conversation, present_message


IDENTITY = "a" * 32
TARGET = "O=C(O)c1ccc2nc(-c3cc(Cl)cc(Cl)c3)oc2c1"
AMINOPHENOL = "COC(=O)c1ccc(N)c(O)c1"
ACID_CHLORIDE = "O=C(Cl)c1cc(Cl)cc(Cl)c1"
EXAMPLE = f"""Your target is **tafamidis free acid**.

Start from **methyl 4-amino-3-hydroxybenzoate**, `{AMINOPHENOL}`,
and **3,5-dichlorobenzoyl chloride**, `{ACID_CHLORIDE}`.

| Step | Transformation | Reported conditions | Reported yield |
| ---- | -------------- | ------------------- | -------------- |
| 1 | Acylation | NaHCO₃, dry THF; below 20 °C | 91.8% |
| 2 | Ring closure | p-TsOH·H₂O, toluene reflux | 92.4% |

[Patent examples](https://patents.google.com/patent/US12116352B2/en).
"""


class Tags(HTMLParser):
    def __init__(self, html: str) -> None:
        super().__init__()
        self.tags: list[tuple[str, dict[str, str | None]]] = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


def test_user_example_renders_table_emphasis_links_and_two_exact_structures() -> None:
    view = present_message(EXAMPLE, IDENTITY)
    tags = Tags(view["html"]).tags
    assert sum(tag == "table" for tag, _ in tags) == 1
    assert sum(tag == "th" for tag, _ in tags) == 4
    assert sum(tag == "td" for tag, _ in tags) == 8
    assert sum(tag == "strong" for tag, _ in tags) == 3
    link = next(attrs for tag, attrs in tags if tag == "a")
    assert link["href"] == "https://patents.google.com/patent/US12116352B2/en"
    assert link["rel"] == "noopener noreferrer"
    assert [item["smiles"] for item in view["structures"]] == [AMINOPHENOL, ACID_CHLORIDE]
    for item in view["structures"]:
        svg = ET.fromstring(base64.b64decode(item["image_url"].split(",", 1)[1]))
        assert svg.tag == "{http://www.w3.org/2000/svg}svg"
        assert len(svg.findall(".//{http://www.w3.org/2000/svg}path")) > 5


def test_plain_question_and_smiles_fence_are_drawn_without_changing_notation() -> None:
    question = present_message(f"how to make {TARGET}?", IDENTITY)
    assert [item["smiles"] for item in question["structures"]] == [TARGET]
    view = present_message("```smiles\n[NH4+].[Cl-]\n```\n\nAgain `[NH4+].[Cl-]`", IDENTITY)
    assert [item["smiles"] for item in view["structures"]] == ["[NH4+].[Cl-]"]


def test_invalid_smiles_prose_urls_and_reactions_do_not_create_fragment_drawings() -> None:
    text = "A reaction in THF at 20 °C. `C1CC` is invalid. `CCBr.N>>CCN` is a reaction. [URL](https://example.org/CCO)"
    view = present_message(text, IDENTITY)
    assert view["structures"] == []
    assert "C1CC" in view["html"]


def test_markdown_cannot_execute_html_load_images_or_link_to_executable_schemes() -> None:
    text = """<script>alert(1)</script><img src=x onerror=alert(2)>
<svg onload=alert(3)></svg>
[unsafe](javascript:alert%281%29) [data](data:text/html,evil)
[file](file:///etc/passwd) ![tracking](https://example.org/track.svg)
[relative](/api/v1/scientific/config)
"""
    view = present_message(text, IDENTITY)
    for tag, attrs in Tags(view["html"]).tags:
        assert tag not in {"script", "svg", "img", "iframe", "object"}
        assert not any(name.startswith("on") for name in attrs)
        if tag == "a":
            assert attrs["href"] == "#" or attrs["href"].startswith("https://")


def test_evidence_links_are_scoped_to_the_selected_conversation() -> None:
    reference = "sha256:" + "b" * 64
    view = present_message(f"[Call evidence]({reference})", IDENTITY)
    assert next(attrs["href"] for tag, attrs in Tags(view["html"]).tags if tag == "a") == (
        f"/api/v1/scientific/conversations/{IDENTITY}/artifacts/{reference}"
    )


def test_projection_preserves_saved_answers_and_scientific_evidence() -> None:
    conversation = {"id": IDENTITY, "turns": [{
        "question": f"how to make {TARGET}?",
        "answer": {"answer_markdown": EXAMPLE, "evidence_refs": ["sha256:" + "b" * 64]},
    }]}
    original = deepcopy(conversation)
    view = present_conversation(conversation)
    assert conversation == original
    assert view["turns"][0]["answer"] == original["turns"][0]["answer"]
    assert len(view["turns"][0]["question_presentation"]["structures"]) == 1


def test_saved_conversation_api_and_page_deliver_new_presentation_without_agent_calls() -> None:
    class SavedService:
        def get(self, identity):
            return {"id": identity, "turns": [{
                "question": f"how to make {TARGET}?", "status": "completed",
                "answer": {"answer_markdown": EXAMPLE},
            }]}

    client = TestClient(create_app(
        runtime=object(), scientific_service=SavedService(), recommendation_only=False,
    ), base_url="http://127.0.0.1")
    response = client.get(f"/api/v1/scientific/conversations/{IDENTITY}")
    assert response.status_code == 200
    turn = response.json()["turns"][0]
    assert "<table>" in turn["answer_presentation"]["html"]
    assert len(turn["answer_presentation"]["structures"]) == 2
    assert len(turn["question_presentation"]["structures"]) == 1
    page = client.get("/scientific").text
    assert 'id="question"' in page
    for name, media_type in (("chat.js", "text/javascript"), ("chat.css", "text/css")):
        url = f"/scientific/assets/{name}"
        assert url in page
        asset = client.get(url)
        assert asset.status_code == 200
        assert asset.headers["content-type"].startswith(media_type)
        assert asset.headers["cache-control"] == "no-store"
        assert client.get(url, headers={"origin": "https://attacker.example"}).status_code == 403
    assert "Download SVG" in client.get("/scientific/assets/chat.js").text
    assert client.get("/scientific/assets/scientific_chat.py").status_code == 404
