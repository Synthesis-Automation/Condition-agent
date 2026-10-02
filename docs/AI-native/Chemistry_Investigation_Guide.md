# Chemistry investigation guide

Scientific investigations use the recorded Python/CLI workspace and the existing
agent tool loop. The agent chooses useful actions from the user's question and
available evidence. This document is an ownership map; the linked resources are
the instructions maintained by each layer.

| Responsibility | Instruction owner |
| --- | --- |
| Agent autonomy, fixed-baseline rules and uncertainty | [Core instructions](../../chem_coworker/scientific_workspace/agent_instructions/core.md) |
| Recorded calls, custom scripts and captured primary sources | [Workspace usage](../../chem_coworker/scientific_workspace/agent_instructions/workspace_usage.md) |
| Condition decisions and optional screening advice | [Conditions playbook](../../chem_coworker/scientific_workspace/task_playbooks/conditions.md) |
| Route decisions and optional specialized investigations | [Retrosynthesis playbook](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis.md) |
| Scientific attribution, step evidence and final self-review | [Answer authoring](../../chem_coworker/scientific_workspace/agent_instructions/answer_authoring.md) |
| Concise default answers and reaction-step display | [Presentation profile](../../chem_coworker/scientific_workspace/presentation/default.md) |

Main task playbooks name their optional supporting guides. The agent can skip,
reorder or replace suggestions; chemistry and evidence validators still apply.
Browser turns record the actual instructions independently of scientific identity.
See [task ownership and extension rules](Scientific_Workspace_Core.md#task-playbooks)
and the [workspace reference](readme.md) for implementation and tool details.

Evidence capture and deterministic attribution checks preserve traceability.
They do not replace chemical judgment, independent review or experimental validation.
