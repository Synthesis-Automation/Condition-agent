from chem_coworker.classifier import TaskClassifier, TaskType
from chem_coworker.workflow import WORKFLOW_REGISTRY


def test_forward_synthesis_product_prediction_with_two_reactant_smiles() -> None:
    clf = TaskClassifier()
    q = (
        "Predict the major product of "
        "Brc1ccccc1 and Nc1ccccc1 under Pd catalysis"
    )
    result = clf.classify(q)
    assert result.task_type == TaskType.FORWARD_SYNTHESIS
    assert result.task_type_str == "forward_synthesis"
    assert WORKFLOW_REGISTRY.get_for_task(result.task_type_str).name == "forward_synthesis"


def test_forward_synthesis_explicit_forward_keyword_without_reaction_arrow() -> None:
    clf = TaskClassifier()
    q = "Forward synthesis prediction: reactant A anisidine, reactant B 4-bromopyridine"
    result = clf.classify(q)
    assert result.task_type == TaskType.FORWARD_SYNTHESIS


def test_reaction_smiles_stays_non_forward_even_if_predictive_language_present() -> None:
    clf = TaskClassifier()
    q = (
        "What will happen in this reaction "
        "Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1 ?"
    )
    result = clf.classify(q)
    assert result.has_reaction is True
    assert result.task_type != TaskType.FORWARD_SYNTHESIS


def test_retrosynthesis_target_only_still_routes_to_retrosynthesis() -> None:
    clf = TaskClassifier()
    q = "How do I synthesize c1ccccc1CO from simple starting materials?"
    result = clf.classify(q)
    assert result.task_type == TaskType.RETROSYNTHESIS


def test_reaction_smiles_with_synthesize_word_still_analyze_not_retro() -> None:
    clf = TaskClassifier()
    q = "How to synthesize via Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1"
    result = clf.classify(q)
    assert result.has_reaction is True
    assert result.task_type == TaskType.ANALYZE
