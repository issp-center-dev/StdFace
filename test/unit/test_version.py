"""Unit tests for version module.

Tests for the Python translation of version.h.
"""
from __future__ import annotations

import stdface.core.version as version


class TestVersionConstants:
    """Tests for version constants."""

    def test_major_version(self):
        """VERSION_MAJOR should be 0."""
        assert version.VERSION_MAJOR == 0

    def test_minor_version(self):
        """VERSION_MINOR should be 5."""
        assert version.VERSION_MINOR == 5

    def test_patch_version(self):
        """VERSION_PATCH should be 0."""
        assert version.VERSION_PATCH == 0

    def test_prerelease_empty(self):
        """VERSION_PRERELEASE should be empty string."""
        assert version.VERSION_PRERELEASE == ""


class TestPrintVersion:
    """Tests for print_version function."""

    def test_print_version_output(self, capsys):
        """print_version should print correct version string."""
        version.print_version()
        captured = capsys.readouterr()
        assert captured.out.strip() == "StdFace version 0.5.0"

    def test_print_version_with_prerelease(self, capsys, monkeypatch):
        """print_version should append prerelease if set."""
        monkeypatch.setattr(version, "VERSION_PRERELEASE", "beta.1")
        version.print_version()
        captured = capsys.readouterr()
        assert captured.out.strip() == "StdFace version 0.5.0-beta.1"
