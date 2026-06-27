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

    def rotate_coords(self, coords: np.ndarray, rotation: np.ndarray) -> np.ndarray:
        """Rotate nodal coordinates with an active 3D rotation."""
        return coords @ rotation.transpose()

    def rotate_state(self, state: np.ndarray, rotation: np.ndarray) -> np.ndarray:
        """Rotate translational and rotational DOFs with the same 3D rotation."""
        rotated = state.copy()

        for i_nod in range(4):
            base = 6 * i_nod
            rotated[base : base + 3] = state[base : base + 3] @ rotation.transpose()
            rotated[base + 3 : base + 6] = state[base + 3 : base + 6] @ rotation.transpose()

        return rotated

    def get_dof_rotation(self, rotation: np.ndarray) -> np.ndarray:
        """Return the block rotation acting on shell DOFs."""
        transform = np.zeros((24, 24))

        for i_nod in range(4):
            base = 6 * i_nod
            transform[base : base + 3, base : base + 3] = rotation
            transform[base + 3 : base + 6, base + 3 : base + 6] = rotation

        return transform

    def get_element_data(self) -> elementData:
        """Create zero-state element data for the shell patch."""
        data = elementData(np.zeros(24), np.zeros(24))
        data.coords = self.coords
        return data

    def test_multilayer_laminate_is_initialized(self) -> None:
        """Test that a laminate stack-up is stored correctly."""
        self.assertEqual(self.shell.material.layerCount(), 4)
        self.assertAlmostEqual(self.shell.material.thick, 0.5)
        self.assertTrue(self.shell.reducedShearIntegration)

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

    def test_tangent_stiffness_covaries_under_rigid_rotation(self) -> None:
        """Test that the local-frame shell formulation is objective under rotation."""
        data = self.get_element_data()
        data.state = np.array(
            [
                1.0e-3,
                -2.0e-3,
                5.0e-4,
                3.0e-2,
                -1.0e-2,
                2.0e-2,
                -1.5e-3,
                7.0e-4,
                -8.0e-4,
                -1.0e-2,
                1.5e-2,
                -2.5e-2,
                8.0e-4,
                1.2e-3,
                -6.0e-4,
                1.2e-2,
                1.0e-2,
                -1.0e-2,
                -9.0e-4,
                -1.1e-3,
                1.4e-3,
                -2.0e-2,
                5.0e-3,
                1.8e-2,
            ]
        )

        self.shell.getTangentStiffness(data)

        axis = np.array([1.0, -2.0, 0.5])
        axis /= np.linalg.norm(axis)
        angle = 0.37
        skew = np.array(
            [
                [0.0, -axis[2], axis[1]],
                [axis[2], 0.0, -axis[0]],
                [-axis[1], axis[0], 0.0],
            ]
        )
        rotation = (
            np.eye(3)
            + np.sin(angle) * skew
            + (1.0 - np.cos(angle)) * (skew @ skew)
        )
        transform = self.get_dof_rotation(rotation)

        rotated_data = self.get_element_data()
        rotated_data.coords = self.rotate_coords(self.coords, rotation)
        rotated_data.state = self.rotate_state(data.state, rotation)

        rotated_shell = ReissnerMindlinShell([0, 1, 2, 3], self.props)
        rotated_shell.getTangentStiffness(rotated_data)

        np.testing.assert_allclose(rotated_data.fint, transform @ data.fint, rtol=1.0e-8, atol=1.0e-8)
        np.testing.assert_allclose(
            rotated_data.stiff,
            transform @ data.stiff @ transform.transpose(),
            rtol=1.0e-8,
            atol=1.0e-8,
        )


if __name__ == "__main__":
    unittest.main()
