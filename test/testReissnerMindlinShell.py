# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

"""Unit tests for the Reissner-Mindlin shell element."""

import unittest

import numpy as np

from pyfem.elements.ReissnerMindlinShell import ReissnerMindlinShell
from pyfem.util.dataStructures import Properties, elementData, solverStatus


class TestReissnerMindlinShellLaminate(unittest.TestCase):
    """Tests laminate support in the Reissner-Mindlin shell element."""

    def setUp(self) -> None:
        """Create a flat Quad4 shell patch with a symmetric composite laminate."""
        self.props = Properties()
        self.props.rank = 3
        self.props.solverStat = solverStatus()
        self.props.materials = ["UD"]
        self.props.layers = ["ply0_bot", "ply90_bot", "ply90_top", "ply0_top"]
        self.props.shearCorrection = 5.0 / 6.0
        self.props.drillingScale = 1.0e-6

        self.props.UD = Properties(
            {
                "E1": 135.0e3,
                "E2": 10.0e3,
                "nu12": 0.3,
                "G12": 5.0e3,
                "G13": 4.0e3,
                "G23": 3.8e3,
                "rho": 1.6e-9,
            }
        )

        self.props.ply0_bot = Properties(
            {"material": "UD", "theta": 0.0, "thickness": 0.125}
        )
        self.props.ply90_bot = Properties(
            {"material": "UD", "theta": 90.0, "thickness": 0.125}
        )
        self.props.ply90_top = Properties(
            {"material": "UD", "theta": 90.0, "thickness": 0.125}
        )
        self.props.ply0_top = Properties(
            {"material": "UD", "theta": 0.0, "thickness": 0.125}
        )

        self.coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )

        self.shell = ReissnerMindlinShell([0, 1, 2, 3], self.props)

    def get_element_data(self) -> elementData:
        """Create zero-state element data for the shell patch."""
        data = elementData(np.zeros(24), np.zeros(24))
        data.coords = self.coords
        return data

    def test_multilayer_laminate_is_initialized(self) -> None:
        """Test that a laminate stack-up is stored correctly."""
        self.assertEqual(self.shell.material.layerCount(), 4)
        self.assertAlmostEqual(self.shell.material.thick, 0.5)

        zeta_points = list(self.shell.iterateLayers())
        self.assertEqual(len(zeta_points), 4)
        self.assertAlmostEqual(
            sum(weight for _, zeta_data in zeta_points for _, weight in zeta_data),
            self.shell.material.thick,
        )

    def test_layer_matrices_follow_layer_orientation(self) -> None:
        """Test that different ply orientations produce different layer matrices."""
        c0 = self.shell.getLayerMatrix(self.shell.material.layers[0])
        c90 = self.shell.getLayerMatrix(self.shell.material.layers[1])

        self.assertFalse(np.allclose(c0, c90))

    def test_tangent_stiffness_assembles_for_composite_shell(self) -> None:
        """Test tangent stiffness assembly for a multilayer shell patch."""
        data = self.get_element_data()

        self.shell.getTangentStiffness(data)

        np.testing.assert_allclose(data.fint, np.zeros_like(data.fint))
        np.testing.assert_allclose(data.stiff, data.stiff.transpose())
        self.assertGreater(np.linalg.norm(data.stiff), 0.0)

    def test_mass_matrix_uses_laminate_inertia(self) -> None:
        """Test mass assembly for a multilayer shell patch."""
        data = self.get_element_data()

        self.shell.getMassMatrix(data)

        np.testing.assert_allclose(data.mass, data.mass.transpose())
        self.assertGreater(np.trace(data.mass), 0.0)
        self.assertGreater(data.lumped, 0.0)


if __name__ == "__main__":
    unittest.main()
