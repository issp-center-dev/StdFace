"""Unit tests for param_check module.

Tests for the parameter validation, printing, and program-exit utilities
extracted from ``stdface_model_util``.
"""
from __future__ import annotations


import numpy as np
import pytest

from stdface.core.param_check import (
    print_val_d,
    print_val_dd,
    print_val_c,
    print_val_i,
    not_used_d,
    not_used_j,
    not_used_i,
    required_val_i,
)
NaN_d = float("nan")


class TestPrintValD:
    """Tests for print_val_d."""

    def test_returns_default_when_nan(self):
        """Test that NaN is replaced by the default value."""
        result = print_val_d("test", NaN_d, 1.5)
        assert result == 1.5

    def test_returns_value_when_set(self):
        """Test that a set value is returned unchanged."""
        result = print_val_d("test", 2.5, 1.5)
        assert result == 2.5

    def test_prints_default_tag(self, caplog):
        """Test that the DEFAULT VALUE tag is logged for NaN."""
        import logging

        with caplog.at_level(logging.INFO, logger="stdface.core.param_check"):
            print_val_d("myvar", NaN_d, 3.0)
        assert "DEFAULT VALUE IS USED" in caplog.text
        assert "myvar" in caplog.text

    def test_no_default_tag_when_set(self, caplog):
        """Test that no DEFAULT tag is logged for a set value."""
        import logging

        with caplog.at_level(logging.INFO, logger="stdface.core.param_check"):
            print_val_d("myvar", 3.0, 0.0)
        assert "DEFAULT VALUE IS USED" not in caplog.text


class TestPrintValDD:
    """Tests for print_val_dd."""

    def test_primary_default(self):
        """Test that val0 is used when val is NaN and val0 is set."""
        result = print_val_dd("test", NaN_d, 2.0, 3.0)
        assert result == 2.0

    def test_secondary_default(self):
        """Test that val1 is used when both val and val0 are NaN."""
        result = print_val_dd("test", NaN_d, NaN_d, 4.0)
        assert result == 4.0

    def test_returns_value_when_set(self):
        """Test that a set value is returned unchanged."""
        result = print_val_dd("test", 1.0, 2.0, 3.0)
        assert result == 1.0


class TestPrintValC:
    """Tests for print_val_c."""

    def test_returns_default_when_nan(self):
        """Test that NaN is replaced by the default complex value."""
        result = print_val_c("test", complex(NaN_d, 0.0), complex(1.0, 2.0))
        assert result == complex(1.0, 2.0)

    def test_returns_value_when_set(self):
        """Test that a set complex value is returned unchanged."""
        result = print_val_c("test", complex(3.0, 4.0), complex(1.0, 2.0))
        assert result == complex(3.0, 4.0)


class TestPrintValI:
    """Tests for print_val_i."""

    def test_returns_default_when_sentinel(self):
        """Test that sentinel is replaced by the default value."""
        result = print_val_i("test", None, 10)
        assert result == 10

    def test_returns_value_when_set(self):
        """Test that a set value is returned unchanged."""
        result = print_val_i("test", 5, 10)
        assert result == 5

    def test_prints_default_tag(self, caplog):
        """Test that the DEFAULT VALUE tag is logged for sentinel."""
        import logging

        with caplog.at_level(logging.INFO, logger="stdface.core.param_check"):
            print_val_i("myvar", None, 7)
        assert "DEFAULT VALUE IS USED" in caplog.text


class TestNotUsedD:
    """Tests for not_used_d."""

    def test_no_exit_when_nan(self):
        """Test that NaN values do not trigger exit."""
        not_used_d("test", NaN_d)  # should not raise

    def test_exits_when_set(self):
        """Test that a set value raises ValueError."""
        with pytest.raises(ValueError):
            not_used_d("test", 1.0)

    def test_handles_complex_nan(self):
        """Test that complex NaN (real part NaN) does not trigger exit."""
        not_used_d("test", complex(NaN_d, 0.0))  # should not raise

    def test_exits_for_complex_set(self):
        """Test that a set complex value (real part set) raises ValueError."""
        with pytest.raises(ValueError):
            not_used_d("test", complex(1.0, 0.0))


class TestNotUsedJ:
    """Tests for not_used_j."""

    def test_no_exit_when_all_nan(self):
        """Test that all-NaN J values do not trigger exit."""
        J = np.full((3, 3), NaN_d)
        not_used_j("J0", NaN_d, J)  # should not raise

    def test_exits_when_scalar_set(self):
        """Test that a set scalar JAll raises ValueError."""
        J = np.full((3, 3), NaN_d)
        with pytest.raises(ValueError):
            not_used_j("J0", 1.0, J)

    def test_exits_when_matrix_element_set(self):
        """Test that a set matrix element raises ValueError."""
        J = np.full((3, 3), NaN_d)
        J[1, 1] = 0.5
        with pytest.raises(ValueError):
            not_used_j("J0", NaN_d, J)


class TestNotUsedI:
    """Tests for not_used_i."""

    def test_no_exit_when_sentinel(self):
        """Test that sentinel value does not trigger exit."""
        not_used_i("test", None)  # should not raise

    def test_exits_when_set(self):
        """Test that a set value raises ValueError."""
        with pytest.raises(ValueError):
            not_used_i("test", 5)


class TestRequiredValI:
    """Tests for required_val_i."""

    def test_exits_when_missing(self):
        """Test that missing (sentinel) value raises ValueError."""
        with pytest.raises(ValueError):
            required_val_i("test", None)

    def test_no_exit_when_present(self):
        """Test that a present value does not raise."""
        required_val_i("test", 5)  # should not raise

    def test_prints_value(self, caplog):
        """Test that the value is logged when present."""
        import logging

        with caplog.at_level(logging.INFO, logger="stdface.core.param_check"):
            required_val_i("myvar", 42)
        assert "myvar" in caplog.text
        assert "42" in caplog.text
