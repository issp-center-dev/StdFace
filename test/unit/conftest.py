"""Pytest configuration for StdFace unit tests.

Adds the python/ directory to sys.path so that the ``stdface`` package
can be imported (e.g., ``from stdface.core.stdface_vals import StdIntList``).
"""
from __future__ import annotations

import os
import sys

# Add python/ directory to the front of sys.path
_src_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "python")
_src_dir = os.path.abspath(_src_dir)
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)
