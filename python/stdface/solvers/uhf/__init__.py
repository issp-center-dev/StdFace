"""UHF solver plugin package.

Importing this package auto-registers the UHF plugin.
"""
from __future__ import annotations

from ._plugin import UHFPlugin  # noqa: F401 — public API

__all__ = ["UHFPlugin"]
