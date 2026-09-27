"""Tests for the solid-like shell element."""

import unittest

import numpy as np

from pyfem.elements.SLS import SLS
from pyfem.elements.SLSkinematic import SLSkinematic
from pyfem.elements.SLSutils import SLSparameters
from pyfem.elements.SolidLikeShell import SolidLikeShell


class TestSolidLikeShell(unittest.TestCase):
    def test_hmat_uses_linear_through_thickness_interpolation(self) -> None:
        kinematic = SLSkinematic(SLSparameters(8))
        sdat = type("ShapeData", (), {"h": np.array([0.25, 0.25, 0.25, 0.25])})()

        hmat = kinematic.getHmat(sdat, 0.4)

        np.testing.assert_allclose(hmat[0, 0], 0.5 * (1.0 - 0.4) * 0.25)
        np.testing.assert_allclose(hmat[0, 12], 0.5 * (1.0 + 0.4) * 0.25)
        np.testing.assert_allclose(hmat[:, :12].sum(axis=1), 0.5 * (1.0 - 0.4))
        np.testing.assert_allclose(hmat[:, 12:].sum(axis=1), 0.5 * (1.0 + 0.4))

    def test_new_and_legacy_element_names_are_available(self) -> None:
        self.assertTrue(issubclass(SLS, SolidLikeShell))
        self.assertEqual(SolidLikeShell.__name__, "SolidLikeShell")


if __name__ == "__main__":
    unittest.main()
