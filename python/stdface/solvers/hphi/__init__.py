"""HPhi solver plugin package.

Importing this package auto-registers the HPhi plugin.
"""
from __future__ import annotations

from ._plugin import HPhiPlugin  # noqa: F401 — public API

__all__ = ["HPhiPlugin"]
