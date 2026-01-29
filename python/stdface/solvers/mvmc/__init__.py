"""mVMC solver plugin package.

Importing this package auto-registers the mVMC plugin.
"""
from __future__ import annotations

from ._plugin import MVMCPlugin  # noqa: F401 — public API

__all__ = ["MVMCPlugin"]
