"""Shared provider boundary for schema-constrained assistance actions."""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Callable, Dict, Literal, Mapping, Protocol, Type

from pydantic import BaseModel, ConfigDict, Field, ValidationError


_PROMPT_PATH = (
    Path(__file__).resolve().parent.parent
    / "definitions"
    / "assistance_prompt.v1.txt"
)
_ALL_ACTIONS = Literal[
    "audit_target",
    "recommend_conditions",
    "inspect_condition_candidate",
    "compare_condition_candidates",
    "propose_clarification",
    "disconnect_target",
    "inspect_strategy",
    "compare_strategies",
    "inspect_strategy_conditions",
    "plan_routes",
    "inspect_route",
    "inspect_route_step",
    "search_step_precedents",
    "compare_routes",
    "apply_repair",
    "verify_strategy",
    "verify_route",
    "retry_route_search",
    "finish",
]


class ClaimPayload(BaseModel):
    """Strict provider payload for one evidence-linked advisory claim."""

    model_config = ConfigDict(extra="forbid")

    claim_type: str = Field(min_length=1, max_length=100)
    subject_id: str = Field(min_length=1, max_length=200)
    severity: Literal["info", "caution", "serious"]
    statement: str = Field(min_length=1, max_length=1200)
    evidence_ids: list[str] = Field(min_length=1)
    uncertainty: str = Field(max_length=800)
    suggested_user_action: str | None = Field(default=None, max_length=800)


class ActionArgumentsPayload(BaseModel):
    """Closed union-like argument object validated again per capability."""

    model_config = ConfigDict(extra="forbid")

    result_ref: str | None = None
    candidate_alias: str | None = None
    candidate_aliases: list[str] | None = None
    strategy_alias: str | None = None
    strategy_aliases: list[str] | None = None
    route_alias: str | None = None
    route_aliases: list[str] | None = None
    step_index: int | None = None
    question_id: str | None = None
    prompt: str | None = None
    constraint_owner: str | None = None
    constraint_kind: str | None = None
    reason: str | None = None
    terminal_status: Literal["completed", "abstained_insufficient_evidence"] | None = (
        None
    )
    stopping_reason: str | None = None
    search_depth_delta: int | None = None
    beam_width_delta: int | None = None
    max_expansions_delta: int | None = None
    proposal_id: str | None = Field(default=None, min_length=1, max_length=200)


class AssistanceActionPayload(BaseModel):
    """Exactly one requested action from a provider turn."""

    model_config = ConfigDict(extra="forbid")

    action_name: _ALL_ACTIONS
    arguments: ActionArgumentsPayload
    cited_evidence_ids: list[str]
    decision_summary: str = Field(min_length=1, max_length=500)
    claims: list[ClaimPayload]


@dataclass(frozen=True)
class AssistanceProviderSettings:
    """Provider configuration independent of domain request contracts."""

    provider: Literal["openai", "aliyun"] = "openai"
    model: str = "gpt-5.6-terra"
    reasoning_effort: Literal["none", "low", "medium", "high", "xhigh", "max"] = (
        "medium"
    )
    max_output_tokens: int = 8_000

    def __post_init__(self) -> None:
        if not self.model.strip():
            raise ValueError("provider model must not be empty")
        if self.max_output_tokens < 256 or self.max_output_tokens > 20_000:
            raise ValueError("max_output_tokens must be between 256 and 20000")

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "AssistanceProviderSettings":
        """Build settings while rejecting unknown provider keys."""

        allowed = {"provider", "model", "reasoning_effort", "max_output_tokens"}
        unknown = set(value) - allowed
        if unknown:
            raise ValueError(f"unknown provider settings: {sorted(unknown)!r}")
        return cls(**dict(value))


@dataclass(frozen=True)
class StructuredTransportResult:
    """A parsed typed payload and measured provider metadata."""

    payload: BaseModel
    response_id: str | None = None
    input_tokens: int = 0
    output_tokens: int = 0
    attempts: int = 1
    elapsed_ms: int = 0


@dataclass(frozen=True)
class AssistanceTransportResult:
    """A parsed action and measured provider metadata."""

    payload: AssistanceActionPayload
    response_id: str | None = None
    input_tokens: int = 0
    output_tokens: int = 0
    attempts: int = 1
    elapsed_ms: int = 0

    def __post_init__(self) -> None:
        if min(
            self.input_tokens,
            self.output_tokens,
            self.attempts,
            self.elapsed_ms,
        ) < 0 or self.attempts < 1:
            raise ValueError("transport usage values must be non-negative")


class AssistanceTransport(Protocol):
    """Injectable transport boundary used by the bounded controller."""

    def complete(
        self,
        packet: Mapping[str, Any],
        settings: AssistanceProviderSettings,
        *,
        max_retries: int,
    ) -> AssistanceTransportResult:
        """Return exactly one schema-validated action."""


class IncompleteAssistanceOutputError(ValueError):
    """A provider response ended without a valid action payload."""

    def __init__(
        self,
        message: str,
        *,
        input_tokens: int = 0,
        output_tokens: int = 0,
    ) -> None:
        super().__init__(message)
        self.input_tokens = input_tokens
        self.output_tokens = output_tokens


class OpenAICompatibleStructuredTransport:
    """Shared OpenAI Responses and compatible-chat structured transport."""

    def __init__(self, client_factory: Callable[..., Any] | None = None) -> None:
        self._client_factory = client_factory

    def complete(
        self,
        packet: Mapping[str, Any],
        settings: AssistanceProviderSettings,
        *,
        payload_model: Type[BaseModel],
        schema_name: str,
        instructions: str,
        max_retries: int,
    ) -> StructuredTransportResult:
        if max_retries < 0 or max_retries > 1:
            raise ValueError("provider retries must be between zero and one")
        attempts = 1 + max_retries
        failures = []
        consumed_input = 0
        consumed_output = 0
        started = time.perf_counter()
        active = settings
        for attempt in range(1, attempts + 1):
            try:
                if active.provider == "openai":
                    result = self._openai_response(
                        packet,
                        active,
                        payload_model=payload_model,
                        schema_name=schema_name,
                        instructions=instructions,
                    )
                else:
                    result = self._compatible_chat(
                        packet,
                        active,
                        payload_model=payload_model,
                        instructions=instructions,
                    )
                return replace(
                    result,
                    attempts=attempt,
                    input_tokens=result.input_tokens + consumed_input,
                    output_tokens=result.output_tokens + consumed_output,
                    elapsed_ms=int((time.perf_counter() - started) * 1000),
                )
            except (IncompleteAssistanceOutputError, ValidationError) as exc:
                failures.append(self._concise_failure(exc))
                consumed_input += int(getattr(exc, "input_tokens", 0) or 0)
                consumed_output += int(getattr(exc, "output_tokens", 0) or 0)
                active = replace(
                    settings,
                    reasoning_effort="low",
                    max_output_tokens=min(20_000, max(8_000, settings.max_output_tokens * 2)),
                )
        raise IncompleteAssistanceOutputError(
            "provider returned no valid assistance action "
            f"after {attempts} attempts: {'; '.join(failures)}",
            input_tokens=consumed_input,
            output_tokens=consumed_output,
        )

    def _openai_response(
        self,
        packet: Mapping[str, Any],
        settings: AssistanceProviderSettings,
        *,
        payload_model: Type[BaseModel],
        schema_name: str,
        instructions: str,
    ) -> StructuredTransportResult:
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            raise RuntimeError("OPENAI_API_KEY is not set")
        client = self._build_client(
            api_key=api_key,
            base_url=os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1"),
        )
        request: Dict[str, Any] = {
            "model": settings.model,
            "instructions": instructions,
            "input": json.dumps(packet, ensure_ascii=False, sort_keys=True),
            "text": {
                "format": {
                    "type": "json_schema",
                    "name": schema_name,
                    "strict": True,
                    "schema": payload_model.model_json_schema(),
                }
            },
            "max_output_tokens": settings.max_output_tokens,
            "store": False,
        }
        if settings.model.startswith(("gpt-5", "o3", "o4")):
            request["reasoning"] = {"effort": settings.reasoning_effort}
        response = client.responses.create(**request)
        usage = getattr(response, "usage", None)
        input_tokens = int(getattr(usage, "input_tokens", 0) or 0)
        output_tokens = int(getattr(usage, "output_tokens", 0) or 0)
        status = str(getattr(response, "status", "unknown") or "unknown")
        content = str(getattr(response, "output_text", "") or "").strip()
        if status == "incomplete" or not content:
            raise IncompleteAssistanceOutputError(
                f"OpenAI response contained no complete final action (status={status})",
                input_tokens=input_tokens,
                output_tokens=output_tokens,
            )
        try:
            payload = payload_model.model_validate_json(content)
        except ValidationError as exc:
            raise IncompleteAssistanceOutputError(
                "OpenAI returned invalid assistance JSON: " + self._concise_failure(exc),
                input_tokens=input_tokens,
                output_tokens=output_tokens,
            ) from exc
        return StructuredTransportResult(
            payload=payload,
            response_id=getattr(response, "id", None),
            input_tokens=input_tokens,
            output_tokens=output_tokens,
        )

    def _compatible_chat(
        self,
        packet: Mapping[str, Any],
        settings: AssistanceProviderSettings,
        *,
        payload_model: Type[BaseModel],
        instructions: str,
    ) -> StructuredTransportResult:
        api_key = os.getenv("ALIYUN_API_KEY")
        if not api_key:
            raise RuntimeError("ALIYUN_API_KEY is not set")
        client = self._build_client(
            api_key=api_key,
            base_url=os.getenv(
                "ALIYUN_BASE_URL",
                "https://dashscope.aliyuncs.com/compatible-mode/v1",
            ),
        )
        response = client.chat.completions.create(
            model=settings.model,
            messages=(
                {
                    "role": "system",
                    "content": instructions,
                },
                {
                    "role": "user",
                    "content": json.dumps(
                        {
                            "response_instruction": (
                                "Return one JSON object matching output_schema."
                            ),
                            "output_schema": payload_model.model_json_schema(),
                            "assistance_packet": packet,
                        },
                        ensure_ascii=False,
                        sort_keys=True,
                    ),
                },
            ),
            response_format={"type": "json_object"},
            max_tokens=settings.max_output_tokens,
        )
        usage = response.usage
        input_tokens = int(getattr(usage, "prompt_tokens", 0) or 0)
        output_tokens = int(getattr(usage, "completion_tokens", 0) or 0)
        content = str(response.choices[0].message.content or "").strip()
        if not content:
            choice = response.choices[0]
            reasoning_content = str(
                getattr(choice.message, "reasoning_content", "") or ""
            )
            raise IncompleteAssistanceOutputError(
                "compatible provider contained no final action JSON "
                f"(finish_reason={getattr(choice, 'finish_reason', None) or 'unknown'}, "
                f"reasoning_chars={len(reasoning_content)})",
                input_tokens=input_tokens,
                output_tokens=output_tokens,
            )
        try:
            payload = payload_model.model_validate_json(content)
        except ValidationError as exc:
            raise IncompleteAssistanceOutputError(
                "compatible provider returned invalid assistance JSON: "
                + self._concise_failure(exc),
                input_tokens=input_tokens,
                output_tokens=output_tokens,
            ) from exc
        return StructuredTransportResult(
            payload=payload,
            response_id=getattr(response, "id", None),
            input_tokens=input_tokens,
            output_tokens=output_tokens,
        )

    def _build_client(self, **kwargs: Any) -> Any:
        if self._client_factory is not None:
            return self._client_factory(**kwargs)
        from openai import OpenAI

        return OpenAI(**kwargs)

    @staticmethod
    def _concise_failure(exc: Exception) -> str:
        if isinstance(exc, ValidationError):
            errors = exc.errors(include_url=False)
            if errors:
                return "invalid JSON/schema: " + str(errors[0].get("msg") or "")
        return str(exc).replace("\n", " ")[:300]


class OpenAICompatibleAssistanceTransport:
    """Assistance-action specialization of the shared structured transport."""

    def __init__(self, client_factory: Callable[..., Any] | None = None) -> None:
        self._structured = OpenAICompatibleStructuredTransport(client_factory)

    def complete(
        self,
        packet: Mapping[str, Any],
        settings: AssistanceProviderSettings,
        *,
        max_retries: int,
    ) -> AssistanceTransportResult:
        result = self._structured.complete(
            packet,
            settings,
            payload_model=AssistanceActionPayload,
            schema_name="chem_coworker_assistance_action",
            instructions=_PROMPT_PATH.read_text(encoding="utf-8"),
            max_retries=max_retries,
        )
        if not isinstance(result.payload, AssistanceActionPayload):
            raise TypeError("structured transport returned the wrong payload type")
        return AssistanceTransportResult(
            payload=result.payload,
            response_id=result.response_id,
            input_tokens=result.input_tokens,
            output_tokens=result.output_tokens,
            attempts=result.attempts,
            elapsed_ms=result.elapsed_ms,
        )
