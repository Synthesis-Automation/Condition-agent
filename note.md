# My notes

python -m pip install --upgrade rdkit
python -m app.web_api

python -m app.web_api --workbench

Open http://127.0.0.1:8000/ (run only one web server).

Optional research tools: python -m app.web_api --workbench
Rebuild the browser files: python -m app.web_api --build

Shared-core v2 is the default on the actual Full/Compact libraries.
Restart the server; do not pass the small sample --index:

```powershell
python -m app.web_api --workbench
```

Validation and current limitations: docs/new/shared_reaction_core_v2_20260913.md
Default promotion and next steps: docs/new/shared_reaction_core_default_20260913.md
If previously set to '0', remove CONDITION_SHARED_CORE_EXPERIMENTAL from the shell environment.
Expanded review packet: results/shared_core_v2/chemist_review/review_packet.html
Disagreement review: results/shared_core_followup/disagreement_review/review_packet.html
Priority cases and rationale: results/shared_core_followup/disagreement_review/analyst_triage.md

cmd /d /c "cd /d C:\Git-softwares\Condition-agent && C:\Users\xubar\AppData\Local\Programs\Python\Python312\python.exe -m chem_coworker"

## 2-step tandem reaction

CC(C)(C)OC(=O)NCCN.O=c1oc2cc(Br)ccc2cc1-c1ccccc1>>NCCNc1ccc2cc(-c3ccccc3)c(=O)oc2c1

## double sites reaction

Brc1ccc(Br)cc1.CC(N)=O>>CC(=O)Nc1ccc(NC(C)=O)cc1
