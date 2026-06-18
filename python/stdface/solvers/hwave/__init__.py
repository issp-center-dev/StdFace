"""H-wave solver plugin package.

Importing this package auto-registers the H-wave plugins (UHFR / UHFK).
"""
from __future__ import annotations

from ._plugin import UHFRPlugin, UHFKPlugin  # noqa: F401 — public API

__all__ = ["UHFRPlugin", "UHFKPlugin"]
