"""One launcher selects a matching frontend and API, including old page URLs."""

import sys

import pytest
from fastapi.testclient import TestClient

from app.web_api import __main__ as launcher
from app.web_api.main import app, create_app


@pytest.mark.parametrize("workbench", [False, True])
def test_launcher_selects_matching_profile(tmp_path, monkeypatch, workbench):
    (tmp_path / "index.html").write_text("conditions", encoding="utf-8")
    (tmp_path / "workbench.html").write_text("research", encoding="utf-8")
    monkeypatch.setattr(launcher, "DEFAULT_FRONTEND_DIST", tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        ["app.web_api", "--host", "127.0.0.1", "--port", "8123"]
        + (["--workbench"] if workbench else []),
    )
    captured = {}
    monkeypatch.setattr(
        launcher.uvicorn,
        "run",
        lambda application, **kwargs: captured.update(app=application, **kwargs),
    )
    launcher.main()
    assert captured["host"] == "127.0.0.1"
    assert captured["port"] == 8123
    assert captured["app"].state.deployment_profile == (
        "research_workbench" if workbench else "recommendation_only"
    )


@pytest.mark.parametrize("focused", [True, False])
def test_selected_page_matches_api_and_old_urls_redirect(tmp_path, focused):
    (tmp_path / "index.html").write_text("conditions", encoding="utf-8")
    (tmp_path / "workbench.html").write_text("research", encoding="utf-8")
    assets = tmp_path / "assets"
    assets.mkdir()
    (assets / "editor.js").write_text("export default {}", encoding="utf-8")
    client = TestClient(create_app(frontend_dist=tmp_path, recommendation_only=focused))
    response = client.get("/")
    assert response.text == ("conditions" if focused else "research")
    assert response.headers["cache-control"] == "no-store"
    for path in ("/index.html", "/workbench.html"):
        redirect = client.get(path, follow_redirects=False)
        assert redirect.status_code == 307
        assert redirect.headers["location"] == "/"
        assert client.get(path).text == response.text
    assert client.get("/WORKBENCH.HTML").status_code == 404
    assert client.get("/assets/editor.js").status_code == 200
    paths = client.get("/api/openapi.json").json()["paths"]
    assert ("/api/v1/features/analyze" in paths) is not focused
    assert "/api/v1/conditions/recommend" in paths


def test_asgi_default_and_missing_frontend_are_explicit(tmp_path):
    assert app.state.deployment_profile == "recommendation_only"
    client = TestClient(create_app(frontend_dist=tmp_path))
    response = client.get("/")
    assert response.status_code == 503
    assert "python -m app.web_api --build" in response.json()["detail"]
    assert (
        client.get("/api/v1/health").json()["data"]["deployment_profile"]
        == "recommendation_only"
    )
