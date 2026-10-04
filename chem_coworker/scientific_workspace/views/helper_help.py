"""Explicit public helper documentation for agents; no implementation-file discovery."""

from __future__ import annotations

HELPER_EXAMPLES = {
    "capture_source_file": "source = w.capture_source_file('browser.json', url=paper_url, text_path=['text'])",
    "capture_source": "source = w.capture_source(text, url=paper_url, locator='Example 1')",
    "capture_source_image": "image = w.capture_source_image('scheme.png', source_ref=source_ref, locator='Scheme 1, page 2')",
    "record_source_excerpt": "excerpt = w.record_source_excerpt(source_ref, start=120, end=640, locator='Example 1')",
    "inspect_source": "print(w.inspect_source(source_ref, query='Compound 8', limit=2000))",
    "batch_summary": "print(w.batch_summary([event.artifact_ref for event in events]))",
    "prepared_literature_reaction": "reaction = w.prepared_literature_reaction(preparation_ref)",
    "attach_literature_reaction": "draft = w.attach_literature_reaction(draft, 's1', preparation_ref)",
    "attach_recipe_check": "draft = w.attach_recipe_check(draft, 's1', recipe_check_ref)",
    "answer_template": "draft = w.answer_template('An evidence-linked proposal with explicit gaps.')",
    "answer_preflight": "print(w.answer_preflight(draft))",
    "finalize_answer": "receipt = w.finalize_answer(draft_path, draft, findings=findings)",
    "run_python": "event = w.run_python('audit.py', parameters, evidence_refs=(source_ref,))",
}
