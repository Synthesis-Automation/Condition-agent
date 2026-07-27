from pathlib import Path

import pytest

from visualization import (
    RenderStyle,
    render_molecule_image,
    render_molecule_image_bytes,
    render_reaction_image_bytes,
)


def test_molecule_renderer_supports_in_memory_png_and_file_output(
    tmp_path: Path,
) -> None:
    drawing = render_molecule_image_bytes("CCO", size=(240, 160))

    assert drawing.startswith(b"\x89PNG\r\n\x1a\n")

    destination = tmp_path / "nested" / "ethanol.svg"
    result = render_molecule_image(
        "CCO",
        destination,
        size=(240, 160),
        image_format="svg",
    )

    assert result == destination
    assert "<svg" in destination.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "reaction_smiles",
    (
        "CCO>>CC=O",
        "CCO>O>CC=O",
    ),
)
def test_reaction_renderer_supports_both_reaction_smiles_forms(
    reaction_smiles: str,
) -> None:
    drawing = render_reaction_image_bytes(
        reaction_smiles,
        size=(480, 180),
    )

    assert drawing.startswith(b"\x89PNG\r\n\x1a\n")


def test_reaction_renderer_supports_vector_output() -> None:
    drawing = render_reaction_image_bytes(
        "CCO>>CC=O",
        size=(480, 180),
        image_format="svg",
    )

    assert b"<svg" in drawing[:512]


@pytest.mark.parametrize(
    "reaction_smiles",
    (
        "",
        "CCO",
        ">>CCO",
        "CCO>>",
        "CCO>>CCO>>CCO",
    ),
)
def test_reaction_renderer_rejects_invalid_reaction_smiles(
    reaction_smiles: str,
) -> None:
    with pytest.raises(ValueError):
        render_reaction_image_bytes(reaction_smiles)


@pytest.mark.parametrize("size", ((0, 100), (100, -1), (100.5, 100)))
def test_render_style_rejects_invalid_sizes(size) -> None:
    with pytest.raises(ValueError):
        RenderStyle(size=size).validated_size()
