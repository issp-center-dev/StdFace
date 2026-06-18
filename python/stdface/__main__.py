"""Command-line entry point for the StdFace standard-mode input generator.

This module is the Python translation of ``dry.c``.  It can be invoked as::

    python -m stdface stan.in              # default solver: HPhi
    python -m stdface stan.in --solver mVMC

or, equivalently, via the installed console script (if packaged).

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

import argparse
import logging
import sys

from .core.version import print_version
from .core.stdface_main import stdface_main

logger = logging.getLogger(__name__)


def _setup_logging(verbose: bool = False) -> None:
    """Configure logging for command-line use.

    Library users do not call this; they configure logging themselves
    (the package only attaches a :class:`logging.NullHandler`).

    Parameters
    ----------
    verbose : bool
        Emit DEBUG-level messages when True, otherwise INFO and above.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def main(argv: list[str] | None = None) -> int:
    """Parse command-line arguments and run the standard-mode generator.

    Parameters
    ----------
    argv : list of str or None
        Command-line arguments.  ``None`` means ``sys.argv[1:]``.

    Returns
    -------
    int
        Exit status (0 = success, 1 = usage error).
    """
    parser = argparse.ArgumentParser(
        prog="stdface",
        description="StdFace: standard-mode input generator for HPhi / mVMC / UHF / H-wave",
    )
    parser.add_argument(
        "-v", "--version",
        action="store_true",
        help="print version and exit",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default=None,
        help="standard-mode input file (e.g. stan.in)",
    )
    parser.add_argument(
        "--solver",
        choices=["HPhi", "mVMC", "UHF", "HWAVE"],
        default="HPhi",
        help="target solver (default: HPhi)",
    )

    args = parser.parse_args(argv)

    if args.version:
        print_version()
        return 0

    if args.input_file is None:
        print_version()
        parser.print_usage()
        return 1

    _setup_logging()
    try:
        stdface_main(args.input_file, solver=args.solver)
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1
    except ValueError as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
