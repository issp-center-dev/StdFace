"""StdFace: Standard-mode input generator for HPhi / mVMC / UHF / H-wave.

This package provides tools for generating Expert-mode definition files
from a simplified Standard-mode input file.
"""
from __future__ import annotations

import logging

# Library convention: attach a no-op handler so importing applications that do
# not configure logging see no output (and no "No handlers" warning).
logging.getLogger("stdface").addHandler(logging.NullHandler())
