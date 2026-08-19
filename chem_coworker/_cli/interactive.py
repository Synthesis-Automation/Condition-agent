"""Interactive terminal session for repeated condition recommendations."""

from __future__ import annotations

import json
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, ContextManager, Optional, Protocol

from condition_recommender import available_ranking_profiles

from chem_coworker.contracts import (
    CompletionChoice,
    ConditionRequest,
    ConditionResponse,
    ConditionReviewSettings,
    MultistepRetrosynthesisRequest,
    MultistepRetrosynthesisResponse,
    RetrosynthesisRequest,
    RetrosynthesisResponse,
)

from .config import save_config
from .models import infer_provider, provider_model_set, selectable_models


InputFunction = Callable[[str], str]
OutputFunction = Callable[[str], None]
CoworkerResponse = (
    ConditionResponse | RetrosynthesisResponse | MultistepRetrosynthesisResponse
)
ResponseRenderer = Callable[[CoworkerResponse, bool], None]
StatusFunction = Callable[[str], ContextManager[object]]
ClearFunction = Callable[[], None]


class _InteractiveCoworker(Protocol):
    def prepare_reaction(self, reaction_smiles: str) -> dict[str, object]:
        """Return the reaction-completion proposal."""

    def recommend(self, request: ConditionRequest) -> ConditionResponse:
        """Return one condition recommendation response."""


class _InteractiveRetrosynthesisCoworker(Protocol):
    def disconnect(self, request: RetrosynthesisRequest) -> RetrosynthesisResponse:
        """Return one single-step retrosynthesis response."""


class _InteractiveMultistepCoworker(Protocol):
    def plan(
        self, request: MultistepRetrosynthesisRequest
    ) -> MultistepRetrosynthesisResponse:
        """Return one bounded multistep retrosynthesis response."""


@dataclass
class InteractiveSettings:
    """Mutable presentation and recommendation settings for one CLI session."""

    top_k: int = 5
    minimum_pool_size: Optional[int] = None
    ranking_profile: str = "default"
    unrestricted_fallback: bool = False
    as_json: bool = False
    review_mode: str = "auto"
    provider: str = "openai"
    model: str = "gpt-5.6-terra"
    reasoning_effort: str = "medium"
    review_candidates: int = 10
    review_max_tokens: int = 8_000
    apply_review_order: bool = True
    mode: str = "conditions"
    max_realizations_per_strategy: int = 3
    max_templates_to_apply: int = 500
    max_candidates_to_validate: int = 100
    include_l0: bool = True
    use_retro_context: bool = True
    include_retro_conditions: bool = True
    retro_condition_top_k: int = 3
    multistep_max_depth: int = 3
    multistep_per_step_top_k: int = 5
    multistep_beam_width: int = 20
    multistep_max_expansions: int = 12
    multistep_max_templates_to_apply: int = 40
    multistep_max_candidates_to_validate: int = 10
    strategic_guidance: str = ""

    def review_settings(self) -> ConditionReviewSettings:
        """Return the immutable review settings for one request."""

        return ConditionReviewSettings(
            mode=self.review_mode,  # type: ignore[arg-type]
            provider=self.provider,
            model=self.model,
            reasoning_effort=self.reasoning_effort,
            max_candidates=self.review_candidates,
            max_output_tokens=self.review_max_tokens,
            apply_order=self.apply_review_order,
        )


class InteractiveSession:
    """Small command loop around one preloaded condition coworker."""

    def __init__(
        self,
        coworker: _InteractiveCoworker,
        *,
        retrosynthesis_coworker: Optional[_InteractiveRetrosynthesisCoworker] = None,
        multistep_coworker: Optional[_InteractiveMultistepCoworker] = None,
        settings: Optional[InteractiveSettings] = None,
        input_fn: InputFunction = input,
        output_fn: OutputFunction = print,
        response_renderer: Optional[ResponseRenderer] = None,
        status_fn: Optional[StatusFunction] = None,
        clear_fn: Optional[ClearFunction] = None,
    ) -> None:
        self.coworker = coworker
        self.retrosynthesis_coworker = retrosynthesis_coworker
        self.multistep_coworker = multistep_coworker
        self.settings = settings or InteractiveSettings()
        self.input = input_fn
        self.output = output_fn
        self.response_renderer = response_renderer
        self.status = status_fn or (lambda _: nullcontext())
        self.clear = clear_fn
        self.last_response: Optional[CoworkerResponse] = None
        self._profiles = {
            str(item["profile_id"]): item for item in available_ranking_profiles()
        }
        self._models = selectable_models()

    def run(
        self,
        initial_reaction: Optional[str] = None,
        *,
        show_welcome: bool = True,
    ) -> int:
        """Run until the user quits or input closes."""

        if show_welcome:
            self._welcome()
        if initial_reaction:
            self.submit(initial_reaction)
        while True:
            try:
                prompt = (
                    "target> "
                    if self.settings.mode in {"retro", "multistep"}
                    else "reaction> "
                )
                value = self.input(prompt).strip()
            except EOFError:
                self.output("\nSession closed.")
                return 0
            except KeyboardInterrupt:
                self.output("\nUse /quit to exit.")
                continue
            if not value:
                continue
            if value.startswith("/"):
                if not self._command(value):
                    return 0
                continue
            self.submit(value)

    def submit(self, value: str) -> Optional[CoworkerResponse]:
        """Dispatch one molecule/reaction input according to the active mode."""

        if self.settings.mode == "multistep":
            return self.plan_multistep(value)
        if self.settings.mode == "retro":
            return self.disconnect(value)
        return self.recommend(value)

    def recommend(self, reaction_smiles: str) -> Optional[ConditionResponse]:
        """Prepare, optionally complete, and recommend one reaction."""

        try:
            with self.status("Analyzing reaction graph..."):
                proposal = self.coworker.prepare_reaction(reaction_smiles)
            choices = self._completion_choices(proposal)
            if choices is None:
                self.output("Recommendation cancelled.")
                return None
            status = "Filtering chemistry and ranking precedents..."
            if self.settings.review_mode != "off":
                status = "Ranking precedents and reviewing uncertainty..."
            with self.status(status):
                response = self.coworker.recommend(
                    ConditionRequest(
                        reaction_smiles=reaction_smiles,
                        top_k=self.settings.top_k,
                        minimum_pool_size=self.settings.minimum_pool_size,
                        unrestricted_fallback=self.settings.unrestricted_fallback,
                        ranking_profile=self.settings.ranking_profile,
                        completion_choices=choices,
                        review=self.settings.review_settings(),
                    )
                )
        except (RuntimeError, ValueError) as exc:
            self.output(f"Error: {exc}")
            return None

        self.last_response = response
        self.output("")
        self._render_response(response)
        self.output("")
        return response

    def plan_multistep(
        self, target_smiles: str
    ) -> Optional[MultistepRetrosynthesisResponse]:
        """Run bounded route search and optional route-level LLM review."""

        if self.multistep_coworker is None:
            self.output("Multistep retrosynthesis is unavailable.")
            return None
        status = "Searching and validating bounded multistep routes..."
        if self.settings.review_mode != "off":
            status = "Searching routes, checking conditions, and reviewing..."
        try:
            with self.status(status):
                response = self.multistep_coworker.plan(
                    MultistepRetrosynthesisRequest(
                        target_smiles=target_smiles,
                        top_k=min(10, self.settings.top_k),
                        max_depth=self.settings.multistep_max_depth,  # type: ignore[arg-type]
                        per_step_top_k=self.settings.multistep_per_step_top_k,
                        beam_width=self.settings.multistep_beam_width,
                        max_expansions=self.settings.multistep_max_expansions,
                        max_templates_to_apply=(
                            self.settings.multistep_max_templates_to_apply
                        ),
                        max_candidates_to_validate=(
                            self.settings.multistep_max_candidates_to_validate
                        ),
                        use_context=self.settings.use_retro_context,
                        include_l0=self.settings.include_l0,
                        include_conditions=self.settings.include_retro_conditions,
                        condition_top_k=self.settings.retro_condition_top_k,
                        strategic_guidance=self.settings.strategic_guidance,
                        review=self.settings.review_settings(),
                    )
                )
        except (RuntimeError, ValueError) as exc:
            self.output(f"Error: {exc}")
            return None
        self.last_response = response
        self.output("")
        self._render_response(response)
        self.output("")
        return response

    def disconnect(self, target_smiles: str) -> Optional[RetrosynthesisResponse]:
        """Run one deterministic one-step search and optional LLM review."""

        if self.retrosynthesis_coworker is None:
            self.output("Retrosynthesis is unavailable; check the startup warnings.")
            return None
        status = "Generating and forward-validating disconnections..."
        if self.settings.review_mode != "off":
            status = "Validating strategies, checking conditions, and reviewing..."
        try:
            with self.status(status):
                response = self.retrosynthesis_coworker.disconnect(
                    RetrosynthesisRequest(
                        target_smiles=target_smiles,
                        top_k=self.settings.top_k,
                        max_realizations_per_strategy=(
                            self.settings.max_realizations_per_strategy
                        ),
                        max_templates_to_apply=self.settings.max_templates_to_apply,
                        max_candidates_to_validate=(
                            self.settings.max_candidates_to_validate
                        ),
                        use_context=self.settings.use_retro_context,
                        include_l0=self.settings.include_l0,
                        include_conditions=self.settings.include_retro_conditions,
                        condition_top_k=self.settings.retro_condition_top_k,
                        condition_minimum_pool_size=self.settings.minimum_pool_size,
                        unrestricted_condition_fallback=(
                            self.settings.unrestricted_fallback
                        ),
                        review=self.settings.review_settings(),
                    )
                )
        except (RuntimeError, ValueError) as exc:
            self.output(f"Error: {exc}")
            return None
        self.last_response = response
        self.output("")
        self._render_response(response)
        self.output("")
        return response

    def _completion_choices(
        self,
        proposal: dict[str, object],
    ) -> Optional[tuple[CompletionChoice, ...]]:
        requirements = proposal.get("requirements") or ()
        if not isinstance(requirements, (list, tuple)) or not requirements:
            return ()

        self.output("This reaction appears to require missing source confirmation.")
        choices = []
        for raw_requirement in requirements:
            if not isinstance(raw_requirement, dict):
                raise ValueError("invalid completion requirement payload")
            requirement_id = str(raw_requirement.get("requirement_id") or "")
            fragment = (
                raw_requirement.get("rooted_fragment_smiles")
                or raw_requirement.get("canonical_fragment_smiles")
                or raw_requirement.get("fragment_key")
                or requirement_id
            )
            options = raw_requirement.get("options") or ()
            if not isinstance(options, (list, tuple)):
                raise ValueError("invalid completion options payload")
            self.output(f"Missing source for {fragment}:")
            for index, option in enumerate(options, start=1):
                if not isinstance(option, dict):
                    continue
                label = option.get("display_name") or option.get("option_id")
                kind = str(option.get("option_kind") or "option").replace("_", " ")
                self.output(f"  {index}. {label} [{kind}]")
            self.output("  c. Enter a custom source identifier")
            self.output("  q. Cancel this recommendation")

            while True:
                try:
                    selected = self.input("source> ").strip()
                except (EOFError, KeyboardInterrupt):
                    self.output("")
                    return None
                if selected.casefold() == "q":
                    return None
                if selected.casefold() == "c":
                    custom = self.input("custom source> ").strip()
                    if custom:
                        choices.append(
                            CompletionChoice(
                                requirement_id=requirement_id,
                                custom_identifier=custom,
                            )
                        )
                        break
                    self.output("Custom source must not be empty.")
                    continue
                try:
                    option = options[int(selected) - 1]
                except (ValueError, IndexError):
                    self.output("Choose an option number, c, or q.")
                    continue
                if not isinstance(option, dict) or not option.get("option_id"):
                    self.output("That option is unavailable.")
                    continue
                choices.append(
                    CompletionChoice(
                        requirement_id=requirement_id,
                        option_id=str(option["option_id"]),
                    )
                )
                break
        return tuple(choices)

    def _command(self, value: str) -> bool:
        command, _, argument = value.partition(" ")
        name = command.casefold()
        argument = argument.strip()
        if name in {"/quit", "/exit", "/q"}:
            self.output("Goodbye.")
            return False
        if name in {"/help", "/h", "/?"}:
            self._help()
        elif name == "/settings":
            self._show_settings()
        elif name == "/mode":
            self._set_mode(argument)
        elif name == "/top-k":
            self._set_top_k(argument)
        elif name == "/minimum":
            self._set_minimum(argument)
        elif name == "/profile":
            self._set_profile(argument)
        elif name == "/model":
            self._set_model(argument)
        elif name == "/provider":
            self._set_provider(argument)
        elif name == "/review":
            self._set_review_mode(argument)
        elif name == "/reasoning":
            self._set_reasoning(argument)
        elif name == "/review-limit":
            self._set_review_limit(argument)
        elif name == "/max-tokens":
            self._set_max_tokens(argument)
        elif name == "/review-order":
            self._set_review_order(argument)
        elif name == "/retro-conditions":
            self._set_retro_conditions(argument)
        elif name == "/realizations":
            self._set_realizations(argument)
        elif name == "/validate-limit":
            self._set_validate_limit(argument)
        elif name == "/depth":
            self._set_depth(argument)
        elif name == "/per-step":
            self._set_per_step(argument)
        elif name == "/beam":
            self._set_beam(argument)
        elif name == "/expansions":
            self._set_expansions(argument)
        elif name == "/guidance":
            self._set_guidance(argument)
        elif name == "/json":
            self._set_json(argument)
        elif name == "/last":
            self._show_last()
        elif name == "/save":
            self._save(argument)
        elif name == "/clear":
            self._clear()
        else:
            self.output(f"Unknown command: {command}. Use /help.")
        return True

    def _welcome(self) -> None:
        self.output("Chem Coworker")
        self.output("Conditions plus graph-validated one- and multistep retrosynthesis.")
        self.output("Enter reaction SMILES, or select a retrosynthesis mode.")
        self.output("")

    def _help(self) -> None:
        self.output("Commands:")
        self.output("  /mode conditions|retro|multistep  Select the workflow")
        self.output("  /settings           Show current settings")
        self.output("  /top-k N            Set recommendation count (1-50)")
        self.output("  /minimum N|auto     Set minimum retrieval pool size")
        self.output("  /profile [ID]       List or select a ranking profile")
        self.output("  /model [NAME]       List or select the LLM review model")
        self.output("  /provider [NAME]    List or select openai/aliyun")
        self.output("  /review MODE        Set off, auto, or always")
        self.output("  /reasoning LEVEL    Set none/low/medium/high/xhigh/max")
        self.output("  /review-limit N     Review the first 1-10 recipes/strategies")
        self.output("  /max-tokens N       Set review output limit (256-20000)")
        self.output("  /review-order on|off Apply advisory presentation ordering")
        self.output("  /retro-conditions on|off Recommend conditions for routes")
        self.output("  /realizations N      Retain 1-10 precursor variants per strategy")
        self.output("  /validate-limit N    Bound forward-validation attempts")
        self.output("  /depth 2|3           Set maximum multistep route depth")
        self.output("  /per-step N          Set candidates retained per expansion")
        self.output("  /beam N              Set multistep frontier beam width")
        self.output("  /expansions N        Set multistep state-expansion budget")
        self.output("  /guidance TEXT|clear Set advisory route preferences")
        self.output("  /json on|off        Toggle full JSON output")
        self.output("  /last               Show the previous result again")
        self.output("  /save PATH          Save the previous result as JSON")
        self.output("  /clear              Clear the terminal")
        self.output("  /quit               Exit")

    def _show_settings(self) -> None:
        minimum = self.settings.minimum_pool_size or "auto"
        self.output(f"mode: {self.settings.mode}")
        self.output(f"top-k: {self.settings.top_k}")
        self.output(f"minimum pool: {minimum}")
        self.output(f"ranking profile: {self.settings.ranking_profile}")
        self.output(f"LLM review: {self.settings.review_mode}")
        self.output(f"model: {self.settings.model} ({self.settings.provider})")
        self.output(f"reasoning effort: {self.settings.reasoning_effort}")
        self.output(f"review candidate limit: {self.settings.review_candidates}")
        self.output(f"review max tokens: {self.settings.review_max_tokens}")
        self.output(
            "apply review order: "
            + ("on" if self.settings.apply_review_order else "off")
        )
        self.output(f"JSON output: {'on' if self.settings.as_json else 'off'}")
        self.output(
            "unrestricted fallback: "
            + ("on" if self.settings.unrestricted_fallback else "off")
        )
        self.output(
            "retrosynthesis conditions: "
            + ("on" if self.settings.include_retro_conditions else "off")
        )
        self.output(
            "retrosynthesis realizations per strategy: "
            f"{self.settings.max_realizations_per_strategy}"
        )
        self.output(
            "retrosynthesis validation limit: "
            f"{self.settings.max_candidates_to_validate}"
        )
        self.output(f"multistep max depth: {self.settings.multistep_max_depth}")
        self.output(
            f"multistep per-step top-k: {self.settings.multistep_per_step_top_k}"
        )
        self.output(f"multistep beam width: {self.settings.multistep_beam_width}")
        self.output(
            f"multistep expansion budget: {self.settings.multistep_max_expansions}"
        )
        self.output(
            "multistep template/validation budgets: "
            f"{self.settings.multistep_max_templates_to_apply}/"
            f"{self.settings.multistep_max_candidates_to_validate}"
        )
        self.output(
            "route guidance: " + (self.settings.strategic_guidance or "none")
        )

    def _set_mode(self, argument: str) -> None:
        aliases = {
            "condition": "conditions",
            "conditions": "conditions",
            "retro": "retro",
            "retrosynthesis": "retro",
            "multi": "multistep",
            "route": "multistep",
            "multistep": "multistep",
        }
        mode = aliases.get(argument.casefold())
        if mode is None:
            self.output("Usage: /mode conditions|retro|multistep.")
            return
        if mode == "retro" and self.retrosynthesis_coworker is None:
            self.output("Retrosynthesis is unavailable; check the startup warnings.")
            return
        if mode == "multistep" and self.multistep_coworker is None:
            self.output("Multistep retrosynthesis is unavailable.")
            return
        self.settings.mode = mode
        expected = (
            "target SMILES"
            if mode in {"retro", "multistep"}
            else "reaction SMILES"
        )
        self.output(f"mode set to {mode}; enter {expected}.")

    def _set_top_k(self, argument: str) -> None:
        try:
            value = int(argument)
        except ValueError:
            self.output("Usage: /top-k N, where N is 1-50.")
            return
        if value < 1 or value > 50:
            self.output("top-k must be between 1 and 50.")
            return
        self.settings.top_k = value
        self.output(f"top-k set to {value}.")

    def _set_minimum(self, argument: str) -> None:
        if argument.casefold() in {"auto", "none"}:
            self.settings.minimum_pool_size = None
            self.output("minimum pool set to auto.")
            return
        try:
            value = int(argument)
        except ValueError:
            self.output("Usage: /minimum N|auto.")
            return
        if value < 1:
            self.output("minimum pool must be positive.")
            return
        self.settings.minimum_pool_size = value
        self.output(f"minimum pool set to {value}.")

    def _set_profile(self, argument: str) -> None:
        if not argument:
            self.output("Ranking profiles:")
            for profile_id, profile in self._profiles.items():
                marker = " *" if profile_id == self.settings.ranking_profile else ""
                self.output(
                    f"  {profile_id}{marker}: "
                    f"{profile.get('description') or profile.get('label') or ''}"
                )
            return
        if argument not in self._profiles:
            self.output(f"Unknown ranking profile: {argument}. Use /profile to list.")
            return
        self.settings.ranking_profile = argument
        self.output(f"ranking profile set to {argument}.")

    def _set_json(self, argument: str) -> None:
        values = {"on": True, "off": False}
        selected = values.get(argument.casefold())
        if selected is None:
            self.output("Usage: /json on|off.")
            return
        self.settings.as_json = selected
        self.output(f"JSON output {'enabled' if selected else 'disabled'}.")

    def _set_model(self, argument: str) -> None:
        if not argument:
            self.output("Review models:")
            for item in self._models:
                marker = " *" if item["name"] == self.settings.model else ""
                self.output(f"  {item['name']} ({item['provider']}){marker}")
            return
        known = {item["name"] for item in self._models}
        if argument not in known:
            self.output(f"Unknown model: {argument}. Use /model to list.")
            return
        self.settings.model = argument
        self.settings.provider = infer_provider(argument)
        self._persist_review_settings()
        self.output(f"review model set to {argument} ({self.settings.provider}).")

    def _set_provider(self, argument: str) -> None:
        if not argument:
            self.output("Providers: openai, aliyun")
            return
        provider = argument.casefold()
        if provider not in {"openai", "aliyun"}:
            self.output("Usage: /provider openai|aliyun.")
            return
        self.settings.provider = provider
        provider_models = provider_model_set(provider)
        if self.settings.model not in provider_models:
            self.settings.model = next(
                item["name"] for item in self._models if item["provider"] == provider
            )
        self._persist_review_settings()
        self.output(f"review provider set to {provider}; model {self.settings.model}.")

    def _set_review_mode(self, argument: str) -> None:
        mode = argument.casefold()
        if mode not in {"off", "auto", "always"}:
            self.output("Usage: /review off|auto|always.")
            return
        self.settings.review_mode = mode
        self._persist_review_settings()
        self.output(f"LLM review mode set to {mode}.")

    def _set_reasoning(self, argument: str) -> None:
        effort = argument.casefold()
        if effort not in {"none", "low", "medium", "high", "xhigh", "max"}:
            self.output("Usage: /reasoning none|low|medium|high|xhigh|max.")
            return
        self.settings.reasoning_effort = effort
        self._persist_review_settings()
        self.output(f"reasoning effort set to {effort}.")

    def _set_review_limit(self, argument: str) -> None:
        try:
            value = int(argument)
        except ValueError:
            self.output("Usage: /review-limit N, where N is 1-10.")
            return
        if value < 1 or value > 10:
            self.output("review limit must be between 1 and 10.")
            return
        self.settings.review_candidates = value
        self._persist_review_settings()
        self.output(f"review candidate limit set to {value}.")

    def _set_max_tokens(self, argument: str) -> None:
        try:
            value = int(argument)
        except ValueError:
            self.output("Usage: /max-tokens N, where N is 256-20000.")
            return
        if value < 256 or value > 20_000:
            self.output("review max tokens must be between 256 and 20000.")
            return
        self.settings.review_max_tokens = value
        self._persist_review_settings()
        self.output(f"review max tokens set to {value}.")

    def _set_review_order(self, argument: str) -> None:
        values = {"on": True, "off": False}
        selected = values.get(argument.casefold())
        if selected is None:
            self.output("Usage: /review-order on|off.")
            return
        self.settings.apply_review_order = selected
        self._persist_review_settings()
        self.output(
            "review presentation ordering " + ("enabled." if selected else "disabled.")
        )

    def _set_retro_conditions(self, argument: str) -> None:
        values = {"on": True, "off": False}
        selected = values.get(argument.casefold())
        if selected is None:
            self.output("Usage: /retro-conditions on|off.")
            return
        self.settings.include_retro_conditions = selected
        self.output(
            "retrosynthesis condition lookup "
            + ("enabled." if selected else "disabled.")
        )

    def _set_realizations(self, argument: str) -> None:
        try:
            value = int(argument)
        except ValueError:
            self.output("Usage: /realizations N, where N is 1-10.")
            return
        if value < 1 or value > 10:
            self.output("realizations must be between 1 and 10.")
            return
        self.settings.max_realizations_per_strategy = value
        self.output(f"retrosynthesis realizations set to {value}.")

    def _set_validate_limit(self, argument: str) -> None:
        try:
            value = int(argument)
        except ValueError:
            self.output("Usage: /validate-limit N, where N is positive.")
            return
        if value < 1:
            self.output("validation limit must be positive.")
            return
        self.settings.max_candidates_to_validate = value
        self.output(f"retrosynthesis validation limit set to {value}.")

    def _set_depth(self, argument: str) -> None:
        try:
            value = int(argument)
        except ValueError:
            value = 0
        if value not in {2, 3}:
            self.output("Usage: /depth 2|3.")
            return
        self.settings.multistep_max_depth = value
        self.output(f"multistep maximum depth set to {value}.")

    def _set_per_step(self, argument: str) -> None:
        self._set_positive_multistep(argument, "per-step", "multistep_per_step_top_k")

    def _set_beam(self, argument: str) -> None:
        self._set_positive_multistep(argument, "beam", "multistep_beam_width")

    def _set_expansions(self, argument: str) -> None:
        self._set_positive_multistep(
            argument, "expansions", "multistep_max_expansions"
        )

    def _set_positive_multistep(
        self, argument: str, label: str, attribute: str
    ) -> None:
        try:
            value = int(argument)
        except ValueError:
            value = 0
        if value < 1:
            self.output(f"Usage: /{label} N, where N is positive.")
            return
        setattr(self.settings, attribute, value)
        self.output(f"multistep {label} set to {value}.")

    def _set_guidance(self, argument: str) -> None:
        if argument.casefold() in {"clear", "none", "off"}:
            self.settings.strategic_guidance = ""
            self.output("route guidance cleared.")
            return
        if not argument:
            self.output(
                "route guidance: " + (self.settings.strategic_guidance or "none")
            )
            return
        if len(argument) > 2_000:
            self.output("route guidance must be at most 2000 characters.")
            return
        self.settings.strategic_guidance = argument
        self.output("route guidance updated; it affects advisory LLM ranking only.")

    def _persist_review_settings(self) -> None:
        try:
            save_config(
                self.settings.model,
                self.settings.provider,
                review_mode=self.settings.review_mode,
                reasoning_effort=self.settings.reasoning_effort,
                review_candidates=self.settings.review_candidates,
                review_max_tokens=self.settings.review_max_tokens,
                apply_review_order=self.settings.apply_review_order,
            )
        except OSError as exc:
            self.output(f"Warning: could not save LLM settings: {exc}")

    def _show_last(self) -> None:
        if self.last_response is None:
            self.output("No recommendation has been run yet.")
            return
        self._render_response(self.last_response)

    def _render_response(self, response: CoworkerResponse) -> None:
        if self.response_renderer is not None:
            self.response_renderer(response, self.settings.as_json)
            return
        self.output(
            json.dumps(response.to_dict(), ensure_ascii=False, indent=2)
            if self.settings.as_json
            else response.answer
        )

    def _clear(self) -> None:
        if self.clear is not None:
            self.clear()
        else:
            self.output("\n" * 40)

    def _save(self, argument: str) -> None:
        if self.last_response is None:
            self.output("No recommendation has been run yet.")
            return
        if not argument:
            self.output("Usage: /save PATH.")
            return
        target = Path(argument).expanduser()
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(
                json.dumps(
                    self.last_response.to_dict(),
                    ensure_ascii=False,
                    indent=2,
                ),
                encoding="utf-8",
            )
        except OSError as exc:
            self.output(f"Could not save result: {exc}")
            return
        self.output(f"Saved {target.resolve()}.")


def run_interactive(
    coworker: _InteractiveCoworker,
    *,
    retrosynthesis_coworker: Optional[_InteractiveRetrosynthesisCoworker] = None,
    multistep_coworker: Optional[_InteractiveMultistepCoworker] = None,
    settings: Optional[InteractiveSettings] = None,
    initial_reaction: Optional[str] = None,
    enhanced: Optional[bool] = None,
    persistent_history: bool = True,
) -> int:
    """Run the enhanced TTY interface or its deterministic plain fallback."""

    if enhanced is not False:
        try:
            from .terminal_ui import can_use_terminal_ui, run_terminal_ui

            if enhanced is True or can_use_terminal_ui():
                return run_terminal_ui(
                    coworker,
                    retrosynthesis_coworker=retrosynthesis_coworker,
                    multistep_coworker=multistep_coworker,
                    settings=settings,
                    initial_reaction=initial_reaction,
                    persistent_history=persistent_history,
                )
        except ImportError:
            if enhanced is True:
                print("Enhanced terminal dependencies unavailable; using plain mode.")
    return InteractiveSession(
        coworker,
        retrosynthesis_coworker=retrosynthesis_coworker,
        multistep_coworker=multistep_coworker,
        settings=settings,
    ).run(initial_reaction)


__all__ = ["InteractiveSession", "InteractiveSettings", "run_interactive"]
