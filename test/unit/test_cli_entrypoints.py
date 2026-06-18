"""Tests for package entry points (CLI and ``python/__main__.py`` shim).

Aims to cover lines reported missing in ``cov*.log`` for ``python/__main__.py``
and ``stdface/__main__.py``.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PYTHONPATH = str(REPO_ROOT / "python")


def _run(
    args: list[str],
    *,
    cwd: str | None = None,
    **kwargs,
) -> subprocess.CompletedProcess:
    env = {**os.environ, "PYTHONPATH": PYTHONPATH}
    return subprocess.run(
        [sys.executable, *args],
        cwd=cwd if cwd is not None else str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        **kwargs,
    )


class TestStdfaceMainModule:
    """``python -m stdface`` and ``stdface.__main__.main``."""

    def test_version_flag_exits_zero(self):
        r = _run(["-m", "stdface", "-v"])
        assert r.returncode == 0
        assert r.stdout or r.stderr

    def test_missing_input_prints_usage_exits_one(self):
        r = _run(["-m", "stdface"])
        assert r.returncode == 1

    def test_main_invokes_stdface_on_valid_input(self, tmp_path):
        """Covers ``main()`` success path (``stdface_main`` + return 0)."""
        stan = tmp_path / "stan.in"
        stan.write_text(
            "model = hubbard\n"
            "lattice = chain\n"
            "L = 4\n"
            "nelec = 4\n"
            "method = lanczos\n"
        )
        r = _run(["-m", "stdface", str(stan)], cwd=str(tmp_path))
        assert r.returncode == 0
        assert (tmp_path / "namelist.def").is_file()

    def test_main_function_version(self):
        from stdface.__main__ import main

        assert main(["-v"]) == 0

    def test_main_function_no_input_file(self):
        from stdface.__main__ import main

        assert main([]) == 1

    def test_main_function_success(self, tmp_path):
        stan = tmp_path / "stan.in"
        stan.write_text(
            "model = hubbard\n"
            "lattice = chain\n"
            "L = 4\n"
            "nelec = 4\n"
            "method = lanczos\n"
        )
        from stdface.__main__ import main

        assert main([str(stan)]) == 0

    def test_main_function_invalid_input_returns_one(self, tmp_path, caplog):
        """ValueError from the run is caught and reported (exit code 1)."""
        stan = tmp_path / "stan.in"
        stan.write_text(
            "model = hubbard\n"
            "lattice = __no_such_lattice__\n"
            "L = 4\n"
            "nelec = 4\n"
            "method = lanczos\n"
        )
        from stdface.__main__ import main

        assert main([str(stan)]) == 1


class TestPythonPackageMainShim:
    """``python python/__main__.py`` delegates to ``stdface.__main__.main``."""

    def test_shim_version_invocation(self):
        main_py = REPO_ROOT / "python" / "__main__.py"
        r = _run([str(main_py), "-v"])
        assert r.returncode == 0
