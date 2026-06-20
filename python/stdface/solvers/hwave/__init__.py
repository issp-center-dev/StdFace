"""H-wave solver plugin package.

Importing this package auto-registers :class:`HWavePlugin` under
``HWAVE`` / ``UHFR`` / ``UHFK``.
"""
from __future__ import annotations

from ._plugin import HWavePlugin  # noqa: F401 — public API

__all__ = ["HWavePlugin"]
