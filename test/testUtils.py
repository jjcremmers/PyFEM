"""
Unit tests for PyFEM utilities and shape functions.

These tests exercise the material models, laminate helpers, and finite element
shape functions. They are focused on numerical checks and basic API behaviour;
no production code is modified by these tests.

(c) Joris Remmers (2025)
"""

import unittest, os
import numpy as np
from typing import Any
import random

from pyfem.util.shapeFunctions import (
    shapeData, elemShapeData, 
    getShapeLine2, getShapeLine3,
    getShapeTria3, getShapeQuad4, getShapeTria6, getShapeQuad8, getShapeQuad9,
    getShapeTetra4, getShapePyramid5, getShapePrism6, getShapePrism18, getShapeHexa8,
    getElemType, tria_scheme, tetra_scheme, pyramid_scheme,
    getIntegrationPoints, calcWeightandDerivatives,
    getElemShapeData, getShapeData
)
from pyfem.util.itemList import itemList
from pyfem.util.BezierShapeFunctions import (
    getBezierLine4, calcWeight, getElemBezierData
)
from pyfem.util.BaseModule import BaseModule



class testShapeDataContainers(unittest.TestCase):
    """Tests for shapeData and elemShapeData container classes."""

    def test_shapeData_instantiation(self) -> None:
        """Test that shapeData can be instantiated."""
        sData = shapeData()
        self.assertIsNotNone(sData)

    def test_shapeData_attributes(self) -> None:
        """Test that shapeData attributes can be assigned."""
        sData = shapeData()
        sData.h = np.array([1.0, 0.0])
        sData.dhdxi = np.array([[1.0], [0.0]])
        self.assertIsNotNone(sData.h)
        self.assertIsNotNone(sData.dhdxi)

    def test_elemShapeData_initialization(self) -> None:
        """Test elemShapeData initializes with empty list."""
        elemData = elemShapeData()
        self.assertEqual(len(elemData), 0)
        self.assertEqual(len(elemData.sData), 0)

    def test_elemShapeData_iteration(self) -> None:
        """Test iteration over shape data."""
        elemData = elemShapeData()
        sData1 = shapeData()
        sData2 = shapeData()
        elemData.sData.append(sData1)
        elemData.sData.append(sData2)

        count = 0
        for sData in elemData:
            count += 1
        self.assertEqual(count, 2)

    def test_elemShapeData_length(self) -> None:
        """Test __len__ method of elemShapeData."""
        elemData = elemShapeData()
        self.assertEqual(len(elemData), 0)

        elemData.sData.append(shapeData())
        self.assertEqual(len(elemData), 1)


class testLine2ShapeFunction(unittest.TestCase):
    """Tests for Line2 (linear line element) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity: sum of shape functions = 1."""
        xi = 0.5
        sData = getShapeLine2(xi)
        self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10)

    def test_values_at_nodes(self) -> None:
        """Test shape functions equal 1 at own node, 0 at other."""
        sData_node0 = getShapeLine2(-1.0)
        self.assertAlmostEqual(sData_node0.h[0], 1.0, places=10)
        self.assertAlmostEqual(sData_node0.h[1], 0.0, places=10)

        sData_node1 = getShapeLine2(1.0)
        self.assertAlmostEqual(sData_node1.h[0], 0.0, places=10)
        self.assertAlmostEqual(sData_node1.h[1], 1.0, places=10)

    def test_derivatives(self) -> None:
        """Test shape function derivatives."""
        xi = 0.0
        sData = getShapeLine2(xi)
        self.assertAlmostEqual(sData.dhdxi[0, 0], -0.5, places=10)
        self.assertAlmostEqual(sData.dhdxi[1, 0], 0.5, places=10)

    def test_invalid_input(self) -> None:
        """Test error handling for invalid input."""
        with self.assertRaises(NotImplementedError):
            getShapeLine2(np.array([0.5, 0.5]))

    def test_centerpoint(self) -> None:
        """Test shape functions at element center."""
        xi = 0.0
        sData = getShapeLine2(xi)
        self.assertAlmostEqual(sData.h[0], 0.5, places=10)
        self.assertAlmostEqual(sData.h[1], 0.5, places=10)


class testLine3ShapeFunction(unittest.TestCase):
    """Tests for Line3 (quadratic line element) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity at multiple points."""
        for xi in [-1.0, -0.5, 0.0, 0.5, 1.0]:
            sData = getShapeLine3(xi)
            self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10, 
                                 msg=f"Failed at xi={xi}")

    def test_values_at_nodes(self) -> None:
        """Test shape functions equal 1 at own node, 0 at others."""
        # Corner node at xi=-1
        sData = getShapeLine3(-1.0)
        self.assertAlmostEqual(sData.h[0], 1.0, places=10)
        self.assertAlmostEqual(sData.h[1], 0.0, places=10)
        self.assertAlmostEqual(sData.h[2], 0.0, places=10)

        # Mid-side node at xi=0
        sData = getShapeLine3(0.0)
        self.assertAlmostEqual(sData.h[0], 0.0, places=10)
        self.assertAlmostEqual(sData.h[1], 1.0, places=10)
        self.assertAlmostEqual(sData.h[2], 0.0, places=10)

        # Corner node at xi=1
        sData = getShapeLine3(1.0)
        self.assertAlmostEqual(sData.h[0], 0.0, places=10)
        self.assertAlmostEqual(sData.h[1], 0.0, places=10)
        self.assertAlmostEqual(sData.h[2], 1.0, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapeLine3(0.5)
        self.assertEqual(sData.h.shape, (3,))
        self.assertEqual(sData.dhdxi.shape, (1, 3))

    def test_invalid_input(self) -> None:
        """Test error handling for invalid input."""
        with self.assertRaises(NotImplementedError):
            getShapeLine3(np.array([0.5]))


class testTria3ShapeFunction(unittest.TestCase):
    """Tests for Tria3 (linear triangle) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity at random points."""
        for _ in range(5):
            xi0 = random.random() * 0.8
            xi1 = random.random() * (0.8 - xi0)
            xi = np.array([xi0, xi1])
            sData = getShapeTria3(xi)
            self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10)

    def test_values_at_nodes(self) -> None:
        """Test Kronecker delta property at nodes."""
        nodes = [np.array([0., 0.]), np.array([1., 0.]), np.array([0., 1.])]
        for i, node in enumerate(nodes):
            sData = getShapeTria3(node)
            for j in range(3):
                if i == j:
                    self.assertAlmostEqual(sData.h[j], 1.0, places=10)
                else:
                    self.assertAlmostEqual(sData.h[j], 0.0, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapeTria3(np.array([0.5, 0.3]))
        self.assertEqual(sData.h.shape, (3,))
        self.assertEqual(sData.dhdxi.shape, (3, 2))

    def test_invalid_dimension(self) -> None:
        """Test error handling for wrong dimension."""
        with self.assertRaises(NotImplementedError):
            getShapeTria3(np.array([0.5]))


class testQuad4ShapeFunction(unittest.TestCase):
    """Tests for Quad4 (bilinear quadrilateral) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity at various points."""
        test_points = [
            np.array([-1.0, -1.0]), np.array([0.0, 0.0]), np.array([1.0, 1.0]),
            np.array([-0.5, 0.5]), np.array([0.7, -0.3])
        ]
        for xi in test_points:
            sData = getShapeQuad4(xi)
            self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10, 
                                 msg=f"Failed at {xi}")

    def test_corner_nodes(self) -> None:
        """Test shape functions at corner nodes."""
        corners = [
            (np.array([-1., -1.]), 0),
            (np.array([1., -1.]), 1),
            (np.array([1., 1.]), 2),
            (np.array([-1., 1.]), 3)
        ]
        for xi, node_idx in corners:
            sData = getShapeQuad4(xi)
            for j in range(4):
                if j == node_idx:
                    self.assertAlmostEqual(sData.h[j], 1.0, places=10)
                else:
                    self.assertAlmostEqual(sData.h[j], 0.0, places=10)

    def test_at_center(self) -> None:
        """Test shape functions at element center."""
        sData = getShapeQuad4(np.array([0.0, 0.0]))
        for h in sData.h:
            self.assertAlmostEqual(h, 0.25, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapeQuad4(np.array([0.5, -0.5]))
        self.assertEqual(sData.h.shape, (4,))
        self.assertEqual(sData.dhdxi.shape, (4, 2))


class testTria6ShapeFunction(unittest.TestCase):
    """Tests for Tria6 (quadratic triangle) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity."""
        test_points = [np.array([0.3, 0.3]), np.array([0.1, 0.1]), np.array([0.5, 0.3])]
        for xi in test_points:
            if np.sum(xi) <= 1.0:
                sData = getShapeTria6(xi)
                self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10)

    def test_corner_nodes(self) -> None:
        """Test at corner nodes."""
        corners = [np.array([0., 0.]), np.array([1., 0.]), np.array([0., 1.])]
        for i, xi in enumerate(corners):
            sData = getShapeTria6(xi)
            self.assertAlmostEqual(sData.h[i], 1.0, places=10)
            for j in range(6):
                if j != i:
                    self.assertAlmostEqual(sData.h[j], 0.0, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapeTria6(np.array([0.2, 0.3]))
        self.assertEqual(sData.h.shape, (6,))
        self.assertEqual(sData.dhdxi.shape, (6, 2))


class testQuad8ShapeFunction(unittest.TestCase):
    """Tests for Quad8 (serendipity quadrilateral) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity."""
        test_points = [
            np.array([-1., -1.]), np.array([0., 0.]), np.array([1., 1.]),
            np.array([-0.5, 0.3]), np.array([0.7, -0.8])
        ]
        for xi in test_points:
            sData = getShapeQuad8(xi)
            self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10, 
                                 msg=f"Failed at {xi}")

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapeQuad8(np.array([0.0, 0.0]))
        self.assertEqual(sData.h.shape, (8,))
        self.assertEqual(sData.dhdxi.shape, (8, 2))


class testQuad9ShapeFunction(unittest.TestCase):
    """Tests for Quad9 (biquadratic quadrilateral) shape functions."""
    
    def test_partition_of_unity(self) -> None:
        """Test partition of unity."""
        test_points = [np.array([0., 0.]), np.array([-1., -1.]), np.array([0.3, -0.7])]
        for xi in test_points:
            sData = getShapeQuad9(xi)
            self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapeQuad9(np.array([0.0, 0.0]))
        self.assertEqual(sData.h.shape, (9,))
        self.assertEqual(sData.dhdxi.shape, (9, 2))

    def test_composed_from_line3(self) -> None:
        """Test that Quad9 is tensor product of Line3."""
        xi = np.array([0.3, -0.5])
        sData = getShapeQuad9(xi)

        s0 = getShapeLine3(xi[0])
        s1 = getShapeLine3(xi[1])

        nodeMap = np.array([[0, 1, 2], [7, 8, 3], [6, 5, 4]])
        for i in range(3):
            for j in range(3):
                iNod = nodeMap[i, j]
                expected = s0.h[i] * s1.h[j]
                self.assertAlmostEqual(sData.h[iNod], expected, places=10)


class testTetra4ShapeFunction(unittest.TestCase):
    """Tests for Tetra4 (linear tetrahedron) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity."""
        test_points = [
            np.array([0., 0., 0.]),
            np.array([0.25, 0.25, 0.25]),
            np.array([0.5, 0.3, 0.1])
        ]
        for xi in test_points:
            if np.sum(xi) <= 1.0:
                sData = getShapeTetra4(xi)
                self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10)

    def test_corner_nodes(self) -> None:
        """Test at corner nodes."""
        corners = [
            np.array([0., 0., 0.]),
            np.array([1., 0., 0.]),
            np.array([0., 1., 0.]),
            np.array([0., 0., 1.])
        ]
        for i, xi in enumerate(corners):
            sData = getShapeTetra4(xi)
            self.assertAlmostEqual(sData.h[i], 1.0, places=10)
            for j in range(4):
                if j != i:
                    self.assertAlmostEqual(sData.h[j], 0.0, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapeTetra4(np.array([0.2, 0.3, 0.1]))
        self.assertEqual(sData.h.shape, (4,))
        self.assertEqual(sData.dhdxi.shape, (4, 3))


class testPyramid5ShapeFunction(unittest.TestCase):
    """Tests for Pyramid5 (pyramid element) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity."""
        test_points = [
            np.array([0., 0., -1.]),
            np.array([0., 0., 0.]),
            np.array([0.3, -0.2, -0.5])
        ]
        for xi in test_points:
            sData = getShapePyramid5(xi)
            self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapePyramid5(np.array([0., 0., 0.]))
        self.assertEqual(sData.h.shape, (5,))
        self.assertEqual(sData.dhdxi.shape, (5, 3))


class testPrism6ShapeFunction(unittest.TestCase):
    """Tests for Prism6 (linear prism) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity."""
        test_points = [
            np.array([0.3, 0.3, -1.]),
            np.array([0.5, 0.3, 0.]),
            np.array([0.2, 0.2, 1.])
        ]
        for xi in test_points:
            if np.sum(xi[:2]) <= 1.0:
                sData = getShapePrism6(xi)
                self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10)

    def test_tensor_product(self) -> None:
        """Test that Prism6 is tensor product of Line2 and Tria3."""
        xi = np.array([0.3, 0.2, 0.5])
        sData = getShapePrism6(xi)

        sDataLine2 = getShapeLine2(xi[2])
        sDataTria3 = getShapeTria3(xi[:2])

        for i in range(3):
            for j in range(2):
                expected = sDataLine2.h[j] * sDataTria3.h[i]
                self.assertAlmostEqual(sData.h[i * 2 + j], expected, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapePrism6(np.array([0.2, 0.3, 0.1]))
        self.assertEqual(sData.h.shape, (6,))
        self.assertEqual(sData.dhdxi.shape, (6, 3))


class testPrism18ShapeFunction(unittest.TestCase):
    """Tests for Prism18 (quadratic prism) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity."""
        xi = np.array([0.2, 0.2, 0.0])
        sData = getShapePrism18(xi)
        self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10)

    def test_tensor_product(self) -> None:
        """Test that Prism18 is tensor product of Line3 and Tria6."""
        xi = np.array([0.2, 0.3, 0.4])
        sData = getShapePrism18(xi)

        sDataLine3 = getShapeLine3(xi[2])
        sDataTria6 = getShapeTria6(xi[:2])

        for i in range(6):
            for j in range(3):
                expected = sDataLine3.h[j] * sDataTria6.h[i]
                self.assertAlmostEqual(sData.h[i * 3 + j], expected, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapePrism18(np.array([0.2, 0.3, 0.1]))
        self.assertEqual(sData.h.shape, (18,))
        self.assertEqual(sData.dhdxi.shape, (18, 3))


class testHexa8ShapeFunction(unittest.TestCase):
    """Tests for Hexa8 (linear hexahedron) shape functions."""

    def test_partition_of_unity(self) -> None:
        """Test partition of unity at various points."""
        test_points = [
            np.array([-1., -1., -1.]),
            np.array([0., 0., 0.]),
            np.array([1., 1., 1.]),
            np.array([-0.3, 0.5, -0.7])
        ]
        for xi in test_points:
            sData = getShapeHexa8(xi)
            self.assertAlmostEqual(np.sum(sData.h), 1.0, places=10, 
                                 msg=f"Failed at {xi}")

    def test_corner_nodes(self) -> None:
        """Test at corner nodes."""
        corners = [
            (np.array([-1., -1., -1.]), 0),
            (np.array([1., -1., -1.]), 1),
            (np.array([1., 1., -1.]), 2),
            (np.array([-1., 1., -1.]), 3),
            (np.array([-1., -1., 1.]), 4),
            (np.array([1., -1., 1.]), 5),
            (np.array([1., 1., 1.]), 6),
            (np.array([-1., 1., 1.]), 7)
        ]
        for xi, node_idx in corners:
            sData = getShapeHexa8(xi)
            self.assertAlmostEqual(sData.h[node_idx], 1.0, places=10)
            for j in range(8):
                if j != node_idx:
                    self.assertAlmostEqual(sData.h[j], 0.0, places=10)

    def test_at_center(self) -> None:
        """Test shape functions at element center."""
        sData = getShapeHexa8(np.array([0., 0., 0.]))
        for h in sData.h:
            self.assertAlmostEqual(h, 0.125, places=10)

    def test_shape_sizes(self) -> None:
        """Test output array shapes."""
        sData = getShapeHexa8(np.array([0.3, -0.2, 0.1]))
        self.assertEqual(sData.h.shape, (8,))
        self.assertEqual(sData.dhdxi.shape, (8, 3))

    def test_invalid_dimension(self) -> None:
        """Test error handling for wrong dimension."""
        with self.assertRaises(NotImplementedError):
            getShapeHexa8(np.array([0.5, 0.5]))


class testElemType(unittest.TestCase):
    """Tests for element type determination."""

    def test_getElemType_line2(self) -> None:
        """Test Line2 element detection."""
        elemCoords = np.array([[0.], [1.]])
        self.assertEqual(getElemType(elemCoords), "Line2")

    def test_getElemType_line3(self) -> None:
        """Test Line3 element detection."""
        elemCoords = np.array([[0.], [1.], [0.5]])
        self.assertEqual(getElemType(elemCoords), "Line3")

    def test_getElemType_tria3(self) -> None:
        """Test Tria3 element detection."""
        elemCoords = np.array([[0., 0.], [1., 0.], [0., 1.]])
        self.assertEqual(getElemType(elemCoords), "Tria3")

    def test_getElemType_quad4(self) -> None:
        """Test Quad4 element detection."""
        elemCoords = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
        self.assertEqual(getElemType(elemCoords), "Quad4")

    def test_getElemType_tria6(self) -> None:
        """Test Tria6 element detection."""
        elemCoords = np.array([
            [0., 0.], [1., 0.], [0., 1.],
            [0.5, 0.], [0.5, 0.5], [0., 0.5]
        ])
        self.assertEqual(getElemType(elemCoords), "Tria6")

    def test_getElemType_tetra4(self) -> None:
        """Test Tetra4 element detection."""
        elemCoords = np.array([
            [0., 0., 0.], [1., 0., 0.],
            [0., 1., 0.], [0., 0., 1.]
        ])
        self.assertEqual(getElemType(elemCoords), "Tetra4")

    def test_getElemType_hexa8(self) -> None:
        """Test Hexa8 element detection."""
        elemCoords = np.array([
            [0., 0., 0.], [1., 0., 0.],
            [1., 1., 0.], [0., 1., 0.],
            [0., 0., 1.], [1., 0., 1.],
            [1., 1., 1.], [0., 1., 1.]
        ])
        self.assertEqual(getElemType(elemCoords), "Hexa8")

    def test_getElemType_unknown(self) -> None:
        """Test error for unknown element type."""
        elemCoords = np.array([[0., 0.], [1., 0.]])
        with self.assertRaises(NotImplementedError):
            getElemType(elemCoords)

    def test_getElemType_invalid_rank(self) -> None:
        """Test error for invalid rank."""
        elemCoords = np.array([[0., 0., 0., 0.], [1., 0., 0., 0.]])
        with self.assertRaises(NotImplementedError):
            getElemType(elemCoords)


class testIntegrationSchemes(unittest.TestCase):
    """Tests for integration schemes."""

    def test_tria_scheme_order1(self) -> None:
        """Test triangular integration scheme order 1."""
        xi, weight = tria_scheme(1)
        self.assertEqual(len(xi), 1)
        self.assertEqual(len(weight), 1)
        self.assertAlmostEqual(weight[0], 0.5, places=10)
        self.assertAlmostEqual(xi[0][0], 1.0/3.0, places=10)
        self.assertAlmostEqual(xi[0][1], 1.0/3.0, places=10)

    def test_tria_scheme_order3(self) -> None:
        """Test triangular integration scheme order 3."""
        xi, weight = tria_scheme(3)
        self.assertEqual(len(xi), 3)
        self.assertEqual(len(weight), 3)
        self.assertAlmostEqual(np.sum(weight), 0.5, places=10)

    def test_tria_scheme_order7(self) -> None:
        """Test triangular integration scheme order 7."""
        xi, weight = tria_scheme(7)
        self.assertEqual(len(xi), 7)
        self.assertEqual(len(weight), 7)
        self.assertAlmostEqual(np.sum(weight), 1.0, places=10)

    def test_tria_scheme_invalid_order(self) -> None:
        """Test error for invalid order."""
        with self.assertRaises(NotImplementedError):
            tria_scheme(5)

    def test_tetra_scheme_order1(self) -> None:
        """Test tetrahedral integration scheme order 1."""
        xi, weight = tetra_scheme(1)
        self.assertEqual(len(xi), 1)
        self.assertEqual(len(weight), 1)

    def test_tetra_scheme_invalid_order(self) -> None:
        """Test error for invalid order."""
        with self.assertRaises(NotImplementedError):
            tetra_scheme(2)

    def test_pyramid_scheme_order1(self) -> None:
        """Test pyramid integration scheme order 1."""
        xi, weight = pyramid_scheme(1)
        self.assertEqual(len(xi), 1)
        self.assertEqual(len(weight), 1)
        self.assertAlmostEqual(xi[0][2], -0.5, places=10)

    def test_pyramid_scheme_invalid_order(self) -> None:
        """Test error for invalid order."""
        with self.assertRaises(NotImplementedError):
            pyramid_scheme(2)


class testIntegrationPoints(unittest.TestCase):
    """Tests for integration point generation."""

    def test_getIntegrationPoints_line2_standard(self) -> None:
        """Test standard integration for Line2."""
        xi, weight = getIntegrationPoints('Line2', 0, 'Gauss')
        self.assertEqual(len(xi), 2)
        self.assertEqual(len(weight), 2)
        self.assertAlmostEqual(np.sum(weight), 2.0, places=10)

    def test_getIntegrationPoints_line3_standard(self) -> None:
        """Test standard integration for Line3."""
        xi, weight = getIntegrationPoints('Line3', 0, 'Gauss')
        self.assertEqual(len(xi), 3)
        self.assertEqual(len(weight), 3)

    def test_getIntegrationPoints_quad4_standard(self) -> None:
        """Test standard integration for Quad4."""
        xi, weight = getIntegrationPoints('Quad4', 0, 'Gauss')
        self.assertEqual(len(xi), 4)
        self.assertEqual(len(weight), 4)

    def test_getIntegrationPoints_quad8_standard(self) -> None:
        """Test standard integration for Quad8."""
        xi, weight = getIntegrationPoints('Quad8', 0, 'Gauss')
        self.assertEqual(len(xi), 9)

    def test_getIntegrationPoints_tetra4_standard(self) -> None:
        """Test standard integration for Tetra4."""
        xi, weight = getIntegrationPoints('Tetra4', 0, 'Gauss')
        self.assertEqual(len(xi), 1)

    def test_getIntegrationPoints_hexa8_standard(self) -> None:
        """Test standard integration for Hexa8."""
        xi, weight = getIntegrationPoints('Hexa8', 0, 'Gauss')
        self.assertEqual(len(xi), 8)

    def test_getIntegrationPoints_order_adjustment(self) -> None:
        """Test order adjustment parameter."""
        xi1, w1 = getIntegrationPoints('Line2', 0, 'Gauss')
        xi2, w2 = getIntegrationPoints('Line2', 1, 'Gauss')
        self.assertGreater(len(xi2), len(xi1))

    def test_getIntegrationPoints_unknown_element(self) -> None:
        """Test error for unknown element type."""
        with self.assertRaises(NotImplementedError):
            getIntegrationPoints('Unknown', 0, 'Gauss')


class testCalcWeightandDerivatives(unittest.TestCase):
    """Tests for weight and derivatives calculation."""

    def test_calcWeightandDerivatives_line_element(self) -> None:
        """Test weight and derivatives calculation for line element."""
        elemCoords = np.array([[0.], [1.]])
        sData = getShapeLine2(0.5)
        weight = 1.0

        calcWeightandDerivatives(elemCoords, sData, weight)

        self.assertTrue(hasattr(sData, 'weight'))
        self.assertTrue(hasattr(sData, 'dhdx'))
        self.assertGreater(sData.weight, 0)

    def test_calcWeightandDerivatives_2d_element(self) -> None:
        """Test weight and derivatives calculation for 2D element."""
        elemCoords = np.array([[0., 0.], [1., 0.], [0., 1.]])
        sData = getShapeTria3(np.array([1./3., 1./3.]))
        weight = 0.5

        calcWeightandDerivatives(elemCoords, sData, weight)

        self.assertTrue(hasattr(sData, 'weight'))
        self.assertTrue(hasattr(sData, 'dhdx'))
        self.assertEqual(sData.dhdx.shape, (3, 2))

    def test_calcWeightandDerivatives_3d_element(self) -> None:
        """Test weight and derivatives calculation for 3D element."""
        elemCoords = np.array([
            [0., 0., 0.], [1., 0., 0.],
            [1., 1., 0.], [0., 1., 0.],
            [0., 0., 1.], [1., 0., 1.],
            [1., 1., 1.], [0., 1., 1.]
        ])
        sData = getShapeHexa8(np.array([0., 0., 0.]))
        weight = 1.0

        calcWeightandDerivatives(elemCoords, sData, weight)

        self.assertTrue(hasattr(sData, 'weight'))
        self.assertTrue(hasattr(sData, 'dhdx'))


class testElemShapeData(unittest.TestCase):
    """Tests for element shape data generation."""

    def test_getElemShapeData_line2(self) -> None:
        """Test getting shape data for Line2 element."""
        elemCoords = np.array([[0.], [2.]])
        elemData = getElemShapeData(elemCoords, order=0)

        self.assertGreater(len(elemData), 0)
        for sData in elemData:
            self.assertTrue(hasattr(sData, 'h'))
            self.assertTrue(hasattr(sData, 'dhdx'))
            self.assertTrue(hasattr(sData, 'x'))
            self.assertTrue(hasattr(sData, 'weight'))

    def test_getElemShapeData_quad4(self) -> None:
        """Test getting shape data for Quad4 element."""
        elemCoords = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
        elemData = getElemShapeData(elemCoords, order=0)

        self.assertEqual(len(elemData), 4)
        for sData in elemData:
            self.assertEqual(sData.h.shape[0], 4)
            self.assertEqual(sData.dhdx.shape, (4, 2))

    def test_getElemShapeData_tetra4(self) -> None:
        """Test getting shape data for Tetra4 element."""
        elemCoords = np.array([
            [0., 0., 0.], [1., 0., 0.],
            [0., 1., 0.], [0., 0., 1.]
        ])
        elemData = getElemShapeData(elemCoords, order=0)

        self.assertGreater(len(elemData), 0)
        for sData in elemData:
            self.assertEqual(sData.h.shape[0], 4)
            self.assertEqual(sData.dhdx.shape, (4, 3))

    def test_getElemShapeData_physical_coordinates(self) -> None:
        """Test that physical coordinates are computed."""
        elemCoords = np.array([[0., 0.], [2., 0.], [2., 2.], [0., 2.]])
        elemData = getElemShapeData(elemCoords, order=0)

        for sData in elemData:
            self.assertTrue(all(sData.x >= 0.))
            self.assertTrue(all(sData.x <= 2.))

    def test_getElemShapeData_order_adjustment(self) -> None:
        """Test order adjustment parameter."""
        elemCoords = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])

        elemData1 = getElemShapeData(elemCoords, order=0)
        elemData2 = getElemShapeData(elemCoords, order=1)

        self.assertGreater(len(elemData2), len(elemData1))

    def test_getElemShapeData_explicit_element_type(self) -> None:
        """Test explicit element type specification."""
        elemCoords = np.array([[0., 0.], [1., 0.], [0., 1.]])
        elemData = getElemShapeData(elemCoords, order=0, elemType='Tria3')

        self.assertGreater(len(elemData), 0)
        self.assertEqual(elemData.sData[0].h.shape[0], 3)


class testShapeData(unittest.TestCase):
    """Tests for shape data generation without physical coordinates."""

    def test_getShapeData_line2(self) -> None:
        """Test getting shape data for Line2 element."""
        shpData = getShapeData(order=0, elemType='Line2')

        self.assertGreater(len(shpData), 0)
        for sData in shpData:
            self.assertTrue(hasattr(sData, 'h'))
            self.assertTrue(hasattr(sData, 'dhdx'))
            self.assertTrue(hasattr(sData, 'weight'))

    def test_getShapeData_quad4(self) -> None:
        """Test getting shape data for Quad4 element."""
        shpData = getShapeData(order=0, elemType='Quad4')

        self.assertEqual(len(shpData), 4)
        for sData in shpData:
            self.assertEqual(sData.h.shape[0], 4)

    def test_getShapeData_tria3(self) -> None:
        """Test getting shape data for Tria3 element."""
        shpData = getShapeData(order=0, elemType='Tria3')

        self.assertEqual(len(shpData), 1)
        self.assertEqual(shpData.sData[0].h.shape[0], 3)


class testItemList(unittest.TestCase):
    """Tests for itemList container class."""

    def setUp(self) -> None:
        """Create a fresh itemList for each test."""
        self.items = itemList()

    def test_itemList_instantiation(self) -> None:
        """Test that itemList can be instantiated."""
        items = itemList()
        self.assertIsInstance(items, dict)
        self.assertEqual(len(items), 0)

    def test_itemList_add_single_item(self) -> None:
        """Test adding a single item to itemList."""
        self.items.add(1, "first_item")
        self.assertEqual(len(self.items), 1)
        self.assertIn(1, self.items)
        self.assertEqual(self.items[1], "first_item")

    def test_itemList_add_multiple_items(self) -> None:
        """Test adding multiple items to itemList."""
        self.items.add(10, "item_10")
        self.items.add(20, "item_20")
        self.items.add(30, "item_30")
        
        self.assertEqual(len(self.items), 3)
        self.assertEqual(self.items[10], "item_10")
        self.assertEqual(self.items[20], "item_20")
        self.assertEqual(self.items[30], "item_30")

    def test_itemList_add_various_types(self) -> None:
        """Test adding items of various types."""
        self.items.add(1, 42)
        self.items.add(2, 3.14)
        self.items.add(3, [1, 2, 3])
        self.items.add(4, {"key": "value"})
        self.items.add(5, np.array([1.0, 2.0]))
        
        self.assertEqual(self.items[1], 42)
        self.assertEqual(self.items[2], 3.14)
        self.assertEqual(self.items[3], [1, 2, 3])
        self.assertEqual(self.items[4]["key"], "value")
        np.testing.assert_array_equal(self.items[5], np.array([1.0, 2.0]))

    def test_itemList_add_duplicate_ID_raises_error(self) -> None:
        """Test that adding an item with duplicate ID raises RuntimeError."""
        self.items.add(1, "first")
        
        with self.assertRaises(RuntimeError) as context:
            self.items.add(1, "second")
        
        self.assertIn("already exists", str(context.exception))
        self.assertIn("1", str(context.exception))

    def test_itemList_get_single_item_by_int(self) -> None:
        """Test retrieving a single item by ID using get()."""
        self.items.add(5, "value_five")
        result = self.items.get(5)
        self.assertEqual(result, "value_five")

    def test_itemList_get_multiple_items_by_list(self) -> None:
        """Test retrieving multiple items by list of IDs using get()."""
        self.items.add(1, "one")
        self.items.add(2, "two")
        self.items.add(3, "three")
        
        result = self.items.get([1, 3])
        self.assertEqual(result, ["one", "three"])

    def test_itemList_get_all_items_by_list(self) -> None:
        """Test retrieving all items by list of all IDs."""
        self.items.add(10, "a")
        self.items.add(20, "b")
        self.items.add(30, "c")
        
        result = self.items.get([10, 20, 30])
        self.assertEqual(result, ["a", "b", "c"])

    def test_itemList_get_invalid_ID_raises_error(self) -> None:
        """Test that get() raises KeyError for non-existent ID."""
        self.items.add(1, "value")
        
        with self.assertRaises(KeyError):
            self.items.get(999)

    def test_itemList_get_invalid_argument_type(self) -> None:
        """Test that get() raises RuntimeError for invalid argument type."""
        self.items.add(1, "value")
        
        with self.assertRaises(RuntimeError) as context:
            self.items.get("invalid")
        
        self.assertIn("illegal argument", str(context.exception))

    def test_itemList_getIndices_default_all(self) -> None:
        """Test getIndices() with default argument returns all IDs."""
        self.items.add(100, "a")
        self.items.add(200, "b")
        self.items.add(300, "c")
        
        result = self.items.getIndices()
        self.assertEqual(set(result), {100, 200, 300})

    def test_itemList_getIndices_explicit_all(self) -> None:
        """Test getIndices(-1) explicitly returns all IDs."""
        self.items.add(10, "x")
        self.items.add(20, "y")
        
        result = self.items.getIndices(-1)
        self.assertEqual(set(result), {10, 20})

    def test_itemList_getIndices_single_ID(self) -> None:
        """Test getIndices() with single ID returns index."""
        self.items.add(5, "first")
        self.items.add(10, "second")
        self.items.add(15, "third")
        
        # The index of ID 10 should be its position in the key sequence
        index = self.items.getIndices(10)
        self.assertEqual(index, 1)

    def test_itemList_getIndices_list_of_IDs(self) -> None:
        """Test getIndices() with list of IDs returns list of indices."""
        self.items.add(1, "a")
        self.items.add(2, "b")
        self.items.add(3, "c")
        self.items.add(4, "d")
        
        indices = self.items.getIndices([1, 3])
        self.assertEqual(indices, [0, 2])

    def test_itemList_getIndices_invalid_argument_type(self) -> None:
        """Test that getIndices() raises RuntimeError for invalid type."""
        self.items.add(1, "value")
        
        with self.assertRaises(RuntimeError) as context:
            self.items.getIndices("invalid")
        
        self.assertIn("illegal argument", str(context.exception))

    def test_itemList_getIndices_invalid_ID_in_list(self) -> None:
        """Test that getIndices() raises ValueError for non-existent ID."""
        self.items.add(1, "value")
        
        with self.assertRaises(ValueError):
            self.items.getIndices([1, 999])

    def test_itemList_findID_first_item(self) -> None:
        """Test findID() returns correct ID for first index."""
        self.items.add(100, "first")
        self.items.add(200, "second")
        
        ID = self.items.findID(0)
        self.assertEqual(ID, 100)

    def test_itemList_findID_middle_item(self) -> None:
        """Test findID() returns correct ID for middle index."""
        self.items.add(10, "a")
        self.items.add(20, "b")
        self.items.add(30, "c")
        
        ID = self.items.findID(1)
        self.assertEqual(ID, 20)

    def test_itemList_findID_last_item(self) -> None:
        """Test findID() returns correct ID for last index."""
        self.items.add(5, "x")
        self.items.add(10, "y")
        self.items.add(15, "z")
        
        ID = self.items.findID(2)
        self.assertEqual(ID, 15)

    def test_itemList_findID_out_of_range(self) -> None:
        """Test that findID() raises IndexError for out of range index."""
        self.items.add(1, "value")
        
        with self.assertRaises(IndexError):
            self.items.findID(99)

    def test_itemList_findID_negative_index(self) -> None:
        """Test findID() with negative index (Python list behavior)."""
        self.items.add(1, "first")
        self.items.add(2, "second")
        self.items.add(3, "third")
        
        # Negative indices should work like Python lists
        ID = self.items.findID(-1)
        self.assertEqual(ID, 3)

    def test_itemList_bidirectional_mapping(self) -> None:
        """Test that findID and getIndices are inverse operations."""
        self.items.add(100, "a")
        self.items.add(200, "b")
        self.items.add(300, "c")
        
        # findID(i) should give ID, getIndices(ID) should give i
        for i in range(len(self.items)):
            ID = self.items.findID(i)
            index = self.items.getIndices(ID)
            self.assertEqual(i, index)

    def test_itemList_large_IDs(self) -> None:
        """Test itemList with large ID values."""
        large_id = 1000000
        self.items.add(large_id, "big_id_value")
        
        self.assertEqual(self.items.get(large_id), "big_id_value")
        self.assertEqual(self.items.findID(0), large_id)

    def test_itemList_negative_IDs(self) -> None:
        """Test itemList with negative ID values."""
        self.items.add(-1, "negative_one")
        self.items.add(-10, "negative_ten")
        self.items.add(0, "zero")
        
        self.assertEqual(self.items.get(-1), "negative_one")
        self.assertEqual(self.items.get(-10), "negative_ten")
        self.assertEqual(len(self.items), 3)

    def test_itemList_get_and_add_cycle(self) -> None:
        """Test repeated add and get operations."""
        for i in range(10):
            self.items.add(i, f"value_{i}")
        
        for i in range(10):
            self.assertEqual(self.items.get(i), f"value_{i}")

    def test_itemList_mixed_operations(self) -> None:
        """Test mixed operations with itemList."""
        # Add items
        self.items.add(5, "five")
        self.items.add(10, "ten")
        self.items.add(15, "fifteen")
        
        # Get items
        self.assertEqual(self.items.get(10), "ten")
        
        # Get indices
        indices = self.items.getIndices([5, 15])
        self.assertEqual(len(indices), 2)
        
        # Find IDs
        self.assertEqual(self.items.findID(0), 5)
        self.assertEqual(self.items.findID(2), 15)

    def test_itemList_with_numpy_arrays(self) -> None:
        """Test itemList storing numpy arrays as items."""
        arr1 = np.array([1.0, 2.0, 3.0])
        arr2 = np.array([[1, 2], [3, 4]])
        
        self.items.add(1, arr1)
        self.items.add(2, arr2)
        
        np.testing.assert_array_equal(self.items.get(1), arr1)
        np.testing.assert_array_equal(self.items.get(2), arr2)

    def test_itemList_with_custom_objects(self) -> None:
        """Test itemList storing custom objects."""
        class SimpleObject:
            def __init__(self, name: str):
                self.name = name
        
        obj1 = SimpleObject("object1")
        obj2 = SimpleObject("object2")
        
        self.items.add(1, obj1)
        self.items.add(2, obj2)
        
        self.assertEqual(self.items.get(1).name, "object1")
        self.assertEqual(self.items.get(2).name, "object2")


class testBezierShapeFunctions(unittest.TestCase):
    """Tests for Bezier shape function module."""

    def test_getBezierLine4_parametric_center(self) -> None:
        """Test getBezierLine4 at parametric center xi=0."""
        # Simple straight line from (0,0) to (1,0)
        C = np.array([[0.0, 0.0, 1.0, 1.0],
                      [0.0, 0.0, 0.0, 0.0]])
        
        sData = getBezierLine4(0.0, C)
        
        self.assertIsNotNone(sData)
        self.assertTrue(hasattr(sData, 'h'))
        self.assertTrue(hasattr(sData, 'dhdxi'))
        self.assertEqual(sData.h.shape, (2,))
        self.assertEqual(sData.dhdxi.shape, (2, 1))
        self.assertEqual(sData.xi, 0.0)

    def test_getBezierLine4_parametric_left(self) -> None:
        """Test getBezierLine4 at left boundary xi=-1."""
        C = np.array([[0.0, 1.0, 2.0, 3.0],
                      [0.0, 0.0, 0.0, 0.0]])
        
        sData = getBezierLine4(-1.0, C)
        
        # At xi=-1, first basis function B[0] = 1, others = 0
        # So h should equal first column of C
        expected_h = C[:, 0]
        np.testing.assert_array_almost_equal(sData.h, expected_h)

    def test_getBezierLine4_parametric_right(self) -> None:
        """Test getBezierLine4 at right boundary xi=1."""
        C = np.array([[0.0, 1.0, 2.0, 3.0],
                      [0.0, 0.0, 0.0, 0.0]])
        
        sData = getBezierLine4(1.0, C)
        
        # At xi=1, last basis function B[3] = 1, others = 0
        # So h should equal last column of C
        expected_h = C[:, 3]
        np.testing.assert_array_almost_equal(sData.h, expected_h)

    def test_getBezierLine4_2D_curve(self) -> None:
        """Test getBezierLine4 with 2D curved path."""
        # Quadratic Bezier curve approximated by cubic
        C = np.array([[0.0, 0.33, 0.67, 1.0],
                      [0.0, 0.5, 0.5, 0.0]])
        
        sData = getBezierLine4(0.5, C)
        
        # Check that point is reasonable
        self.assertTrue(0.0 <= sData.h[0] <= 1.0)
        self.assertTrue(0.0 <= sData.h[1] <= 0.5)

    def test_getBezierLine4_3D_control_points(self) -> None:
        """Test getBezierLine4 with 3D control points."""
        C = np.array([[0.0, 1.0, 2.0, 3.0],
                      [0.0, 1.0, 1.0, 0.0],
                      [0.0, 0.0, 0.5, 1.0]])
        
        sData = getBezierLine4(0.0, C)
        
        self.assertEqual(sData.h.shape, (3,))
        self.assertEqual(sData.dhdxi.shape, (3, 1))

    def test_getBezierLine4_invalid_control_points(self) -> None:
        """Test getBezierLine4 raises error with wrong number of control points."""
        # Should have 4 control points, provide 3
        C = np.array([[0.0, 1.0, 2.0],
                      [0.0, 0.0, 0.0]])
        
        with self.assertRaises(NotImplementedError):
            getBezierLine4(0.0, C)

    def test_getBezierLine4_derivatives_finite_difference(self) -> None:
        """Test derivatives match finite difference approximation."""
        C = np.array([[0.0, 0.25, 0.75, 1.0],
                      [0.0, 0.5, 0.5, 0.0]])
        
        xi = 0.3
        h = 1e-8
        
        # Get derivatives analytically
        sData = getBezierLine4(xi, C)
        dhdxi_analytical = sData.dhdxi.flatten()
        
        # Compute finite difference approximation
        sData_plus = getBezierLine4(xi + h, C)
        sData_minus = getBezierLine4(xi - h, C)
        dhdxi_fd = (sData_plus.h - sData_minus.h) / (2 * h)
        
        np.testing.assert_array_almost_equal(dhdxi_analytical, dhdxi_fd, decimal=5)

    def test_calcWeight_square_jacobian_identity(self) -> None:
        """Test calcWeight with identity Jacobian."""
        jac = np.eye(2)
        weight = calcWeight(jac)
        
        # Determinant of identity matrix is 1
        self.assertAlmostEqual(weight, 1.0)

    def test_calcWeight_square_jacobian_scaling(self) -> None:
        """Test calcWeight with scaled Jacobian."""
        jac = np.array([[2.0, 0.0],
                        [0.0, 3.0]])
        weight = calcWeight(jac)
        
        # Determinant should be 6
        self.assertAlmostEqual(weight, 6.0)

    def test_calcWeight_square_jacobian_3D(self) -> None:
        """Test calcWeight with 3D Jacobian."""
        jac = np.array([[1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0]])
        weight = calcWeight(jac)
        
        self.assertAlmostEqual(weight, 1.0)

    def test_calcWeight_rectangular_jacobian_1d(self) -> None:
        """Test calcWeight with 1D curve Jacobian (1, 2 shape)."""
        jac = np.array([[1.0, 1.0]])
        weight = calcWeight(jac)
        
        # Should be sqrt(1^2 + 1^2) = sqrt(2)
        expected = np.sqrt(2.0)
        self.assertAlmostEqual(weight, expected)

    def test_calcWeight_rectangular_jacobian_arc_length(self) -> None:
        """Test calcWeight gives arc length differential for 1D curve."""
        # Parametric line with slope 3
        jac = np.array([[3.0, 4.0]])
        weight = calcWeight(jac)
        
        # Arc length differential: sqrt(3^2 + 4^2) = 5
        self.assertAlmostEqual(weight, 5.0)

    def test_getElemBezierData_basic(self) -> None:
        """Test getElemBezierData returns elemShapeData."""
        C = np.array([[0.0, 0.25, 0.75, 1.0],
                      [0.0, 0.5, 0.5, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0]])
        
        elemData = getElemBezierData(elemCoords, C, order=2, elemType='Line4')
        
        self.assertIsInstance(elemData, elemShapeData)
        self.assertGreater(len(elemData.sData), 0)

    def test_getElemBezierData_integration_points(self) -> None:
        """Test that getElemBezierData generates correct number of integration points."""
        C = np.array([[0.0, 0.5, 1.0, 1.5],
                      [0.0, 0.0, 0.0, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0]])
        
        # Order 2 should give 2 integration points
        elemData = getElemBezierData(elemCoords, C, order=2, elemType='Line4')
        self.assertEqual(len(elemData.sData), 5)
        
        # Order 4 should give 4 integration points
        elemData = getElemBezierData(elemCoords, C, order=4, elemType='Line4')
        self.assertEqual(len(elemData.sData), 7)

    def test_getElemBezierData_shape_data_attributes(self) -> None:
        """Test that each integration point has required attributes."""
        C = np.array([[0.0, 0.33, 0.67, 1.0],
                      [0.0, 0.5, 0.5, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0]])
        
        elemData = getElemBezierData(elemCoords, C, order=2, elemType='Line4')
        
        for sData in elemData.sData:
            self.assertTrue(hasattr(sData, 'h'))
            self.assertTrue(hasattr(sData, 'dhdxi'))
            self.assertTrue(hasattr(sData, 'weight'))
            self.assertTrue(hasattr(sData, 'xi'))
            self.assertGreater(sData.weight, 0.0)

    def test_getElemBezierData_positive_weights(self) -> None:
        """Test that integration weights are positive."""
        C = np.array([[0.0, 0.5, 1.0, 1.5],
                      [0.0, 0.0, 0.0, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0]])
        
        elemData = getElemBezierData(elemCoords, C, order=3, elemType='Line4')
        
        for sData in elemData.sData:
            self.assertGreater(sData.weight, 0.0)

    def test_getElemBezierData_weight_sum_approximates_length(self) -> None:
        """Test that sum of weights approximates curve length."""
        # Straight line from 0 to 1
        C = np.array([[0.0, 0.33, 0.67, 1.0],
                      [0.0, 0.0, 0.0, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0]])
        
        elemData = getElemBezierData(elemCoords, C, order=4, elemType='Line4')
        
        total_weight = sum(sData.weight for sData in elemData.sData)
        
        # For a straight line from 0 to 1 in x-direction, length is 1
        # Integration should approximate this
        self.assertAlmostEqual(total_weight, 1.0, places=1)

    def test_getElemBezierData_invalid_element_type(self) -> None:
        """Test that invalid element type raises NotImplementedError."""
        C = np.array([[0.0, 0.5, 1.0, 1.5],
                      [0.0, 0.0, 0.0, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0]])
        
        with self.assertRaises(NotImplementedError):
            getElemBezierData(elemCoords, C, order=2, elemType='InvalidType')

    def test_getElemBezierData_curved_path(self) -> None:
        """Test getElemBezierData with curved path."""
        # Arc-like curve
        C = np.array([[0.0, 0.5, 0.5, 1.0],
                      [0.0, 0.5, 0.5, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0]])
        
        elemData = getElemBezierData(elemCoords, C, order=3, elemType='Line4')
        
        self.assertGreater(len(elemData.sData), 0)
        for sData in elemData.sData:
            # Point should be within bounding box of control points
            self.assertTrue(0.0 <= sData.h[0] <= 1.0)
            self.assertTrue(0.0 <= sData.h[1] <= 0.5)

    def test_getElemBezierData_gauss_method(self) -> None:
        """Test getElemBezierData with Gauss quadrature method."""
        C = np.array([[0.0, 0.33, 0.67, 1.0],
                      [0.0, 0.0, 0.0, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0]])
        
        elemData = getElemBezierData(elemCoords, C, order=2, method="Gauss", elemType='Line4')
        
        self.assertEqual(len(elemData.sData), 5)
        for sData in elemData.sData:
            self.assertGreater(sData.weight, 0.0)

    def test_getElemBezierData_3D_physical_space(self) -> None:
        """Test getElemBezierData with 3D physical space."""
        C = np.array([[0.0, 0.33, 0.67, 1.0],
                      [0.0, 0.5, 0.5, 0.0],
                      [0.0, 0.0, 0.1, 0.0]])
        elemCoords = np.array([[0.0, 1.0],
                               [0.0, 0.0],
                               [0.0, 0.0]])
        
        elemData = getElemBezierData(elemCoords, C, order=2, elemType='Line4')
        
        for sData in elemData.sData:
            self.assertEqual(sData.h.shape, (3,))

  



class testBaseModule(unittest.TestCase):
    """Unit tests for the `BaseModule` class in `pyfem.util.BaseModule`."""

    def test_instantiation_and_type(self) -> None:
        props = type('Props', (), {})()
        module = BaseModule(props)
        self.assertIsInstance(module, BaseModule)
        self.assertEqual(module.type, 'BaseModule')

    '''
    def test_solver_detection_by_name(self) -> None:
        # Create a derived class whose name contains 'Solver' to trigger detection
        class MySolver(BaseModule):
            pass

        MySolver.__name__ = 'NonlinearSolver'
        module = MySolver(type('Props', (), {})())
        self.assertTrue(module.isSolver)
    '''
    def test_currentModule_property_detection(self) -> None:
        props = type('Props', (), {})()
        props.currentModule = 'solver'
        module = BaseModule(props)
        self.assertTrue(module.isSolver)
    '''
    def test_nested_property_loading_single_level(self) -> None:
        props = type('Props', (), {})()
        props.currentModule = 'myconf'

        # Prepare an iterable of (name, value) pairs to simulate property list
        module_props = type('ModuleProps', (), {})()
        module_props_list = [('alpha', 1), ('beta', 'x')]
        module_props.__iter__ = lambda self: iter(module_props_list)
        props.myconf = module_props

        module = BaseModule(props)
        self.assertEqual(module.alpha, 1)
        self.assertEqual(module.beta, 'x')
    '''
    
    def test_writeHeader_and_writeFooter_timing(self) -> None:
        props = type('Props', (), {})()
        module = BaseModule(props)
        module.isSolver = True

        # writeHeader should set t0
        module.writeHeader(cycle=2)
        self.assertTrue(hasattr(module, 't0'))

        # Prepare global data object with startTime
        globdat = type('Global', (), {})()
        globdat.startTime = module.t0 - 0.5

        # writeFooter should not raise and uses globdat.startTime
        module.writeFooter(globdat)


class TestVtkUtils(unittest.TestCase):
    """Unit tests for pyfem.util.vtkUtils.

    Tests are skipped if the `vtk` package is not available in the
    execution environment.
    """

    def setUp(self) -> None:
        try:
            import importlib
            self.vtkUtils = importlib.import_module('pyfem.util.vtkUtils')
            import vtk  # ensure vtk is importable
            self.vtk = vtk
        except Exception as e:
            self.skipTest(f"VTK not available: {e}")

    def test_setCellNodes_sets_ids(self) -> None:
        cell = self.vtk.vtkLine()
        # Ensure the point id array has size 2
        cell.GetPointIds().SetNumberOfIds(2)
        self.vtkUtils.setCellNodes(cell, [10, 20])
        self.assertEqual(cell.GetPointIds().GetId(0), 10)
        self.assertEqual(cell.GetPointIds().GetId(1), 20)

    def test_insert2Dcontinuum_inserts_cells(self) -> None:
        grid = self.vtk.vtkUnstructuredGrid()
        n0 = grid.GetNumberOfCells()
        # triangle
        self.vtkUtils.insert2Dcontinuum(grid, [0, 1, 2])
        self.assertEqual(grid.GetNumberOfCells(), n0 + 1)
        # line
        self.vtkUtils.insert2Dcontinuum(grid, [3, 4])
        self.assertEqual(grid.GetNumberOfCells(), n0 + 2)
        # quad
        self.vtkUtils.insert2Dcontinuum(grid, [5, 6, 7, 8])
        self.assertEqual(grid.GetNumberOfCells(), n0 + 3)

    def test_insert3Dcontinuum_tetra_inserts_cell(self) -> None:
        grid = self.vtk.vtkUnstructuredGrid()
        n0 = grid.GetNumberOfCells()
        self.vtkUtils.insert3Dcontinuum(grid, [0, 1, 2, 3])
        self.assertEqual(grid.GetNumberOfCells(), n0 + 1)

    def test_insertBeam_and_shell(self) -> None:
        grid = self.vtk.vtkUnstructuredGrid()
        n0 = grid.GetNumberOfCells()
        # beam 2-node
        self.vtkUtils.insertBeam(grid, [0, 1])
        # beam 3-node
        self.vtkUtils.insertBeam(grid, [2, 3, 4])
        # shell triangle
        self.vtkUtils.insertShell(grid, [5, 6, 7])
        # shell quad
        self.vtkUtils.insertShell(grid, [8, 9, 10, 11])
        self.assertEqual(grid.GetNumberOfCells(), n0 + 4)



if __name__ == '__main__':
    unittest.main()
