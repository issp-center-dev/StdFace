"""Built-in solver plugins for StdFace.

Importing this package auto-registers all built-in solver plugins.
"""
from __future__ import annotations

from . import hphi, mvmc, uhf, hwave  # noqa: F401 — auto-register plugins
