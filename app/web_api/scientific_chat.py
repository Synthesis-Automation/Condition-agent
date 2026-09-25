"""Local browser adapter over the shared scientific conversation service."""

from __future__ import annotations

from pathlib import Path
import secrets
from typing import Any
from urllib.parse import urlsplit

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field

from chem_coworker.scientific_workspace.conversation import ConversationService

from .scientific_presentation import present_conversation


class ChatRequest(BaseModel):
    """Only user text and an optional saved conversation ID cross the browser boundary."""

    model_config = ConfigDict(extra="forbid")
    question: str = Field(min_length=1, max_length=20000)
    conversation_id: str | None = Field(default=None, pattern=r"^[0-9a-f]{32}$")


def create_scientific_router(service: ConversationService) -> APIRouter:
    """Expose an opt-in loopback-only UI; never accept executable paths from the browser."""
    router = APIRouter()
    token = secrets.token_urlsafe(32)

    def guard(request: Request, *, mutation: bool = False) -> None:
        hostname = request.url.hostname
        if hostname not in {"127.0.0.1", "localhost", "::1"}:
            raise HTTPException(403, "Scientific chat is only available on loopback")
        origin = request.headers.get("origin")
        if origin:
            parsed = urlsplit(origin)
            if (parsed.scheme, parsed.netloc) != (request.url.scheme, request.url.netloc):
                raise HTTPException(403, "Cross-origin scientific chat requests are not allowed")
        if mutation and not secrets.compare_digest(request.headers.get("x-scientific-token", ""), token):
            raise HTTPException(403, "Reload the scientific chat page to obtain its session token")

    def failure(exc: Exception) -> HTTPException:
        status = 404 if isinstance(exc, FileNotFoundError) else 422 if isinstance(exc, ValueError) else 409
        return HTTPException(status, str(exc))

    @router.get("/scientific", include_in_schema=False)
    def page(request: Request) -> FileResponse:
        guard(request)
        return FileResponse(Path(__file__).with_name("scientific_chat.html"), headers={"Cache-Control": "no-store"})

    @router.get("/scientific/assets/{name}", include_in_schema=False)
    def asset(name: str, request: Request) -> FileResponse:
        guard(request)
        files = {"chat.js": ("scientific_chat.js", "text/javascript"),
                 "chat.css": ("scientific_chat.css", "text/css")}
        if name not in files:
            raise HTTPException(404, "Unknown chat asset")
        filename, media_type = files[name]
        return FileResponse(Path(__file__).with_name(filename), media_type=media_type,
                            headers={"Cache-Control": "no-store"})

    @router.get("/api/v1/scientific/activity")
    def activity(request: Request) -> dict[str, Any]:
        """Expose actual worker-owned activity, even when another chat is selected."""
        guard(request)
        for row in service.list_conversations():
            try:
                conversation = service.get(row["id"])
            except FileNotFoundError:
                continue
            for turn in reversed(conversation["turns"]):
                if turn["status"] in {"queued", "preparing", "running"}:
                    return {"active": {"conversation_id": row["id"], "title": row["title"],
                                       **{key: turn.get(key) for key in (
                                           "id", "status", "created_at", "updated_at", "progress", "repair_attempts",
                                       )}}}
        return {"active": None}

    @router.get("/api/v1/scientific/config")
    def configuration(request: Request) -> dict[str, Any]:
        guard(request)
        return {**service.describe(), "token": token}

    @router.get("/api/v1/scientific/conversations")
    def conversations(request: Request) -> list[dict[str, Any]]:
        guard(request)
        return service.list_conversations()

    @router.post("/api/v1/scientific/turns", status_code=202)
    def submit(payload: ChatRequest, request: Request) -> dict[str, str]:
        guard(request, mutation=True)
        try:
            return service.submit(payload.question, payload.conversation_id)
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            raise failure(exc) from exc

    @router.get("/api/v1/scientific/conversations/{identity}")
    def conversation(identity: str, request: Request) -> dict[str, Any]:
        guard(request)
        try:
            return present_conversation(service.get(identity))
        except (ValueError, FileNotFoundError) as exc:
            raise failure(exc) from exc

    @router.post("/api/v1/scientific/conversations/{identity}/cancel")
    def cancel(identity: str, request: Request) -> dict[str, bool]:
        guard(request, mutation=True)
        try:
            return {"cancellation_requested": service.cancel(identity)}
        except ValueError as exc:
            raise failure(exc) from exc

    @router.get("/api/v1/scientific/conversations/{identity}/artifacts/{reference}")
    def artifact(identity: str, reference: str, request: Request) -> Any:
        guard(request)
        try:
            return service.artifact(identity, reference)
        except (ValueError, FileNotFoundError) as exc:
            raise failure(exc) from exc

    return router
