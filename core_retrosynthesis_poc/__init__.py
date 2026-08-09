"""Reaction-core-derived one-step retrosynthesis proof of concept."""

from .compiler import CompilationResult, compile_core_templates
from .comparison import run_comparison, split_by_reference
from .ensemble import EnsembleCandidate, disconnect_ensemble
from .generic_compiler import GenericCompilationResult, compile_generic_templates
from .generic_library import (
    build_generic_library,
    load_generic_library,
    save_generic_library,
)
from .generic_models import (
    GenericCoreTemplate,
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
)
from .generic_search import disconnect_generic_target
from .html_report import render_comparison_html, write_comparison_html
from .library import build_library, load_library, save_library
from .models import (
    CenterReactivityContext,
    CoreDisconnectionCandidate,
    CoreLibraryBuildReport,
    CoreTemplate,
    CoreTemplateLibrary,
    CoreTemplatePrecedent,
    TemplateContext,
)
from .search import disconnect_target

__all__ = [
    "CenterReactivityContext",
    "CompilationResult",
    "CoreDisconnectionCandidate",
    "CoreLibraryBuildReport",
    "CoreTemplate",
    "CoreTemplateLibrary",
    "CoreTemplatePrecedent",
    "EnsembleCandidate",
    "GenericCompilationResult",
    "GenericCoreTemplate",
    "GenericDisconnectionCandidate",
    "GenericTemplateLibrary",
    "TemplateContext",
    "build_library",
    "build_generic_library",
    "compile_generic_templates",
    "compile_core_templates",
    "disconnect_target",
    "disconnect_ensemble",
    "disconnect_generic_target",
    "load_library",
    "load_generic_library",
    "run_comparison",
    "render_comparison_html",
    "save_library",
    "save_generic_library",
    "split_by_reference",
    "write_comparison_html",
]
