"""Backward-compatible shim for the Wannier90 export functions.

The implementation moved to :mod:`stdface.writer.wannier90_writer` (B3).
This module re-exports the public entry points for any code that still
imports them from the old location.
"""
from __future__ import annotations

from ...writer.wannier90_writer import (  # noqa: F401
    export_geometry,
    export_interaction,
)
