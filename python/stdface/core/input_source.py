"""Input-source abstraction for Standard-mode parameters.

This module defines :class:`InputSource`, an abstract interface that returns
a flat ``{keyword: value}`` dictionary, and concrete implementations that
read the parameters from different formats:

* :class:`StanFileSource` -- the legacy ``stan.in`` keyword-value format.
* :class:`TOMLSource`     -- a TOML file (with bare-key aliases for primes).
* :class:`JSONSource`     -- a JSON file.
* :class:`DictSource`     -- an in-memory dictionary.

Higher layers (``stdface_main`` today, ``generate()`` in D5) consume the
dictionary and apply it to a :class:`~stdface.core.stdface_vals.StdIntList`.

License
-------
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""
from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from .keyword_parser import trim_space_quote

logger = logging.getLogger(__name__)


class InputSource(ABC):
    """Abstract source of Standard-mode keyword/value pairs."""

    @abstractmethod
    def load(self) -> dict[str, Any]:
        """Return the input parameters as a ``{keyword: value}`` dict."""


class StanFileSource(InputSource):
    """Read parameters from a legacy ``stan.in`` keyword-value file.

    Each non-blank, non-comment line must contain ``keyword = value``.
    Keywords are lower-cased.  A keyword appearing twice raises
    :class:`ValueError` (duplicate), matching the behaviour of the original
    incremental parser.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def load(self) -> dict[str, str]:
        try:
            with open(self.path, "r") as fp:
                raw_lines = fp.readlines()
        except OSError as exc:
            logger.error("Cannot open input file: %s", self.path)
            raise FileNotFoundError(str(self.path)) from exc

        data: dict[str, str] = {}
        for raw_line in raw_lines:
            line = trim_space_quote(raw_line)
            if line.startswith("//") or line == "":
                continue

            parts = line.split("=", 1)
            if len(parts) < 2:
                msg = '"=" is NOT found.'
                logger.error(msg)
                raise ValueError(msg)

            keyword = parts[0].lower()
            value = parts[1]
            if keyword in data:
                msg = f"Keyword {keyword} is duplicated."
                logger.error(msg)
                raise ValueError(msg)
            data[keyword] = value
        return data


# ---------------------------------------------------------------------------
#  TOML bare-key aliases
# ---------------------------------------------------------------------------
#
# ``stan.in`` uses primed keywords (``t'``, ``J0''`` ...) that are not valid
# TOML bare keys.  Bare aliases replace each prime with the letter ``p`` so
# users can write ``tp`` / ``t0pp`` instead of quoting ``"t'"`` / ``"t0''"``.


def _build_bare_aliases() -> dict[str, str]:
    bases = ["t", "t0", "t1", "t2", "j", "j0", "j1", "j2", "v", "v0", "v1", "v2"]
    aliases: dict[str, str] = {}
    for base in bases:
        aliases[base + "p"] = base + "'"
        aliases[base + "pp"] = base + "''"
    return aliases


_TOML_BARE_ALIASES: dict[str, str] = _build_bare_aliases()
"""Maps bare TOML keys to their primed Standard-mode keyword."""


def _normalize_keys(data: dict[str, Any]) -> dict[str, Any]:
    """Translate bare TOML aliases (``tp`` -> ``t'``) into primed keywords.

    The alias lookup is case-insensitive; non-alias keys are passed through
    unchanged (the downstream parser lower-cases them anyway).
    """
    return {_TOML_BARE_ALIASES.get(k.lower(), k): v for k, v in data.items()}


class TOMLSource(InputSource):
    """Read parameters from a TOML file (with bare-key alias translation)."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def load(self) -> dict[str, Any]:
        try:
            import tomllib  # Python >= 3.11
        except ModuleNotFoundError:
            import tomli as tomllib  # Python 3.10 backport

        try:
            with open(self.path, "rb") as fp:
                raw = tomllib.load(fp)
        except OSError as exc:
            logger.error("Cannot open input file: %s", self.path)
            raise FileNotFoundError(str(self.path)) from exc
        return _normalize_keys(raw)


class JSONSource(InputSource):
    """Read parameters from a JSON file."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def load(self) -> dict[str, Any]:
        import json

        try:
            with open(self.path, "r") as fp:
                return json.load(fp)
        except OSError as exc:
            logger.error("Cannot open input file: %s", self.path)
            raise FileNotFoundError(str(self.path)) from exc


class DictSource(InputSource):
    """Wrap an in-memory dictionary as an input source."""

    def __init__(self, data: dict[str, Any]) -> None:
        self._data = data

    def load(self) -> dict[str, Any]:
        return dict(self._data)
