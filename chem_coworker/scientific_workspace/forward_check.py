"""Optional, bounded forward evidence for one previously assessed route step."""

from __future__ import annotations

from contextlib import contextmanager
import json
import os
from pathlib import Path
import subprocess
import sys
from time import monotonic
from typing import Any, Iterator, TYPE_CHECKING
from uuid import uuid4

if TYPE_CHECKING:
    from .operations import ScientificOperations


SCHEMA_VERSION = "route_step_forward_investigation.v1"
_LIMITATIONS = [
    "Forward predictions are optional model evidence, not observed products or experimental feasibility.",
    "Timeouts, errors and absent predictions establish neither success nor impossibility.",
    "The saved route assessment is unchanged; this check addresses only the selected step.",
]


def _selected_step(
    operations: ScientificOperations, source_ref: str, step_id: str,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    from .route_investigation import _record

    if not isinstance(step_id, str) or not step_id.strip():
        raise ValueError("step_id must identify one saved route step")
    if not any(event.artifact_ref == source_ref and event.kind == "call"
               for event in operations.store.events()):
        raise ValueError("source_ref must identify a recorded route assessment call")
    record = _record(operations, source_ref)
    proposals = [step for step in record["proposal"]["steps"] if step["external_step_id"] == step_id]
    assessments = [step["assessment"] for step in record["assessment"]["step_assessments"]
                   if step["external_step_id"] == step_id]
    if len(proposals) != 1 or len(assessments) != 1:
        raise ValueError("step_id must identify exactly one previously assessed route step")
    assessment = assessments[0]
    matches = [item for item in assessment["operator_matches"]
               if item["match_id"] == assessment["selected_operator_match_id"]
               and item["match_level"] == "exact_operator_signature"]
    if len(matches) != 1 or not all(assessment.get(key) for key in (
        "canonical_target_smiles", "canonical_precursor_smiles",
    )):
        raise ValueError("Forward check requires a saved step with an exact admitted operator match")
    return proposals[0], assessment, matches[0]


def _worker_command(directory: Path) -> list[str]:
    return [sys.executable, "-m", __name__, str(directory)]


def _read_stages(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    stages = []
    for line in path.read_text("utf-8", errors="replace").splitlines():
        try:
            value = json.loads(line)
        except ValueError:
            continue  # Termination can interrupt the final diagnostic line.
        if isinstance(value, dict):
            stages.append(value)
    return stages


def assess_step_forward(
    operations: ScientificOperations, source_ref: str, step_id: str, question: str,
    timeout_seconds: int = 30,
) -> dict[str, Any]:
    """Run canonical forward assessment in a deadline-limited child, retaining diagnostics."""
    from .agent_runtime import _hidden_process_options, _stop_process
    from .store import canonical_bytes

    if not isinstance(question, str) or not question.strip() or len(question) > 2000:
        raise ValueError("question must be nonempty text of at most 2000 characters")
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= 30:
        raise ValueError("timeout_seconds must be an integer between 1 and 30")
    proposal, assessment, match = _selected_step(operations, source_ref, step_id)
    paths = {name: operations._path(name) for name in ("retro_library", "forward_library")}
    identities = operations.store.manifest["baseline"]["artifacts"]
    inputs = {name: dict(identities[name]) for name in paths}
    if any(not item.get("sha256") for item in inputs.values()):
        raise ValueError("Forward check requires content-pinned library artifacts")
    directory = operations.store.root / "diagnostics" / "forward_checks" / uuid4().hex
    directory.mkdir(parents=True)
    request = {
        "inputs": inputs, "starting_materials": assessment["canonical_precursor_smiles"],
        "intended_product": assessment["canonical_target_smiles"],
        "operator_hint": match["operator_id"], "recipe": proposal.get("proposed_conditions"),
    }
    (directory / "request.json").write_bytes(canonical_bytes(request))
    result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION, "source_ref": source_ref, "step_id": step_id,
        "question": question.strip(), "execution_status": "error", "assessment": None,
        "provenance": {
            "inputs": inputs, "source_assessment_id": assessment["assessment_id"],
            "selected_operator_match_id": match["match_id"], "operator_id": match["operator_id"],
            "operator_signature": match["operator_signature"], "template_ids": match["template_ids"],
        },
        "experimental_feasibility": "not_established", "limitations": list(_LIMITATIONS), "error": None,
    }
    environment = os.environ.copy()
    repository = str(Path(__file__).resolve().parents[2])
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, (repository, environment.get("PYTHONPATH"))))
    process = None
    cleanup_error = None
    started = monotonic()
    try:
        with (directory / "stdout.log").open("wb") as stdout, (directory / "stderr.log").open("wb") as stderr:
            process = subprocess.Popen(
                _worker_command(directory), cwd=repository, env=environment,
                stdout=stdout, stderr=stderr, start_new_session=(os.name != "nt"),
                **_hidden_process_options(),
            )
            process.wait(timeout=max(0, timeout_seconds - (monotonic() - started)))
        output = json.loads((directory / "result.json").read_text("utf-8"))
        if (not isinstance(output, dict) or output.get("execution_status") not in {"completed", "error"}
                or (output["execution_status"] == "completed" and
                    (process.returncode != 0 or not isinstance(output.get("assessment"), dict)))):
            raise ValueError("Forward worker did not return a valid assessment result")
        result.update({key: output.get(key) for key in ("execution_status", "assessment", "error")})
        result["provenance"].update(output.get("provenance", {}))
    except subprocess.TimeoutExpired:
        result.update(execution_status="timed_out", error={
            "type": "TimeoutExpired", "message": f"Forward check exceeded its {timeout_seconds}-second deadline.",
        })
    except KeyboardInterrupt:
        result.update(execution_status="cancelled", error={"type": "KeyboardInterrupt", "message": "Forward check interrupted."})
    except Exception as exc:
        result["error"] = {"type": type(exc).__name__, "message": str(exc)}
    finally:
        if process is not None and process.poll() is None:
            try:
                _stop_process(process)
            except (OSError, subprocess.SubprocessError) as exc:
                cleanup_error = {"type": type(exc).__name__, "message": str(exc)}
                # A process-exit race or failed platform tree cleanup must not
                # discard the scientific timeout/error and its stage records.
                try:
                    if process.poll() is None:
                        process.kill()
                    process.wait(timeout=5)
                except (OSError, subprocess.SubprocessError) as fallback:
                    cleanup_error["fallback_error"] = {"type": type(fallback).__name__, "message": str(fallback)}
    result["execution"] = {
        "timings": {"elapsed_seconds": round(monotonic() - started, 6), "timeout_seconds": timeout_seconds},
        "stages": _read_stages(directory / "stages.jsonl"),
        "diagnostics": {name: (directory / name).relative_to(operations.store.root).as_posix()
                        for name in ("request.json", "stages.jsonl", "stdout.log", "stderr.log", "result.json")
                        if (directory / name).exists()},
    }
    if cleanup_error is not None:
        result["execution"]["cleanup_error"] = cleanup_error
    return result


class _Stages:
    def __init__(self, directory: Path) -> None:
        self.path = directory / "stages.jsonl"
        self.started = monotonic()

    @contextmanager
    def run(self, stage: str) -> Iterator[None]:
        started = monotonic()
        self._write(stage, "running", started)
        try:
            yield
        except Exception:
            self._write(stage, "error", started)
            raise
        else:
            self._write(stage, "completed", started)

    def _write(self, stage: str, status: str, started: float) -> None:
        now = monotonic()
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({
                "stage": stage, "status": status, "elapsed_seconds": round(now - self.started, 6),
                "duration_seconds": round(now - started, 6),
            }) + "\n")
            handle.flush()


def _verify_inputs(inputs: dict[str, Any]) -> None:
    from .baseline import artifact_identity

    for name, expected in inputs.items():
        actual = artifact_identity(Path(expected["path"]))
        if any(actual.get(key) != expected.get(key) for key in ("status", "sha256", "size_bytes", "mtime_ns")):
            raise ValueError(f"Pinned {name} artifact changed; create a new scientific baseline")


def _run_worker(directory: Path) -> int:
    """Child entry point: load pinned inputs and call the public domain assessment once."""
    from .store import canonical_bytes

    stages = _Stages(directory)
    result: dict[str, Any] = {"execution_status": "error", "assessment": None, "error": None}
    try:
        with stages.run("verify_inputs"):
            request = json.loads((directory / "request.json").read_text("utf-8"))
            _verify_inputs(request["inputs"])
        with stages.run("load_retro_library"):
            from core_retrosynthesis.generic_library import load_generic_library

            generic = load_generic_library(request["inputs"]["retro_library"]["path"])
        with stages.run("load_forward_library"):
            from forward_synthesis import load_forward_library

            forward = load_forward_library(request["inputs"]["forward_library"]["path"])
        with stages.run("validate_library_source"):
            from forward_synthesis import validate_forward_library_source

            validate_forward_library_source(forward, generic)
            if not any(operator.operator_id == request["operator_hint"] for operator in forward.operators):
                raise ValueError("Selected saved operator is absent from the pinned forward library")
        result["provenance"] = {
            "input_hashes_verified": True, "library_source_verified": True,
            "forward_library_schema": forward.schema_version,
            "forward_library_definition_id": forward.definition_id,
            "source_library_definition_id": forward.source_library_definition_id,
        }
        with stages.run("assess_selected_step"):
            from core_retrosynthesis.external_proposal_assessment import load_external_proposal_admission_policy
            from forward_synthesis import assess_proposed_step

            limits = load_external_proposal_admission_policy().limits
            assessment = assess_proposed_step(
                request["starting_materials"], request["intended_product"], forward,
                recipe=request["recipe"], operator_hint=request["operator_hint"],
                top_k=limits.maximum_forward_products,
                max_operators_to_apply=limits.maximum_forward_operators,
            )
        with stages.run("verify_inputs_after_assessment"):
            _verify_inputs(request["inputs"])
        result.update(execution_status="completed", assessment=assessment.to_dict())
    except Exception as exc:
        result["error"] = {"type": type(exc).__name__, "message": str(exc)}
    (directory / "result.json").write_bytes(canonical_bytes(result))
    return 0 if result["execution_status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(_run_worker(Path(sys.argv[1])))
