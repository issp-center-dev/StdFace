"""Output file generation subpackage.

This package contains shared modules for writing Expert-mode definition files.
Solver-specific writers have moved to ``stdface.solvers``.
"""
from __future__ import annotations

from .interaction_writer import build_interactions  # noqa: F401
