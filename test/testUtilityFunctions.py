# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

"""Unit tests for utility helper modules."""

import math
import unittest

import numpy as np

from pyfem.util.matrixUtils import skew
from pyfem.util.utilFunctions import (
    macauley,
    one_minus_cos_over_x2,
    sign,
    sin_over_x,
)


class TestMatrixUtils(unittest.TestCase):
    """Tests for matrix utility helpers."""

    def test_skew_returns_expected_matrix(self) -> None:
        """Test that skew builds the correct skew-symmetric matrix."""
        vec = np.array([1.0, 2.0, 3.0])
        expected = np.array(
            [
                [0.0, -3.0, 2.0],
                [3.0, 0.0, -1.0],
                [-2.0, 1.0, 0.0],
            ]
        )

        np.testing.assert_allclose(skew(vec), expected)

    def test_skew_is_antisymmetric(self) -> None:
        """Test that the skew matrix is anti-symmetric."""
        vec = np.array([0.5, -1.25, 4.0])
        mat = skew(vec)

        np.testing.assert_allclose(mat.transpose(), -mat)

    def test_skew_reproduces_cross_product(self) -> None:
        """Test that skew(a) @ b equals cross(a, b)."""
        a_vec = np.array([1.0, -2.0, 0.5])
        b_vec = np.array([0.25, 4.0, -3.0])

        np.testing.assert_allclose(skew(a_vec) @ b_vec, np.cross(a_vec, b_vec))


class TestUtilFunctions(unittest.TestCase):
    """Tests for scalar utility helpers."""

    def test_macauley(self) -> None:
        """Test Macaulay bracket evaluation."""
        self.assertEqual(macauley(2.5), 2.5)
        self.assertEqual(macauley(0.0), 0.0)
        self.assertEqual(macauley(-2.5), 0)

    def test_sign(self) -> None:
        """Test sign evaluation."""
        self.assertEqual(sign(4.0), 1.0)
        self.assertEqual(sign(0.0), 1.0)
        self.assertEqual(sign(-4.0), -1.0)

    def test_sin_over_x_matches_limit_near_zero(self) -> None:
        """Test the series-based limit handling near zero."""
        self.assertAlmostEqual(sin_over_x(0.0), 1.0)
        self.assertAlmostEqual(sin_over_x(1.0e-14), 1.0, places=12)

    def test_sin_over_x_matches_direct_evaluation(self) -> None:
        """Test sin_over_x against the direct analytical expression."""
        x_val = 0.75
        self.assertAlmostEqual(sin_over_x(x_val), math.sin(x_val) / x_val, places=14)

    def test_one_minus_cos_over_x2_matches_limit_near_zero(self) -> None:
        """Test the series-based limit handling near zero."""
        self.assertAlmostEqual(one_minus_cos_over_x2(0.0), 0.5)
        self.assertAlmostEqual(one_minus_cos_over_x2(1.0e-14), 0.5, places=12)

    def test_one_minus_cos_over_x2_matches_direct_evaluation(self) -> None:
        """Test one_minus_cos_over_x2 against the direct analytical expression."""
        x_val = 0.75
        expected = (1.0 - math.cos(x_val)) / (x_val * x_val)
        self.assertAlmostEqual(one_minus_cos_over_x2(x_val), expected, places=14)


if __name__ == "__main__":
    unittest.main()
