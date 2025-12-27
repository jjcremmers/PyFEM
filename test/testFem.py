################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by                #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since        #
#  then augmented and maintained by Joris J.C. Remmers.                        #
#  All rights reserved.                                                        #
#                                                                              #
#  A github repository, with the most up to date version of the code,          #
#  can be found here:                                                          #
#     https://github.com/jjcremmers/PyFEM/                                     #
#     https://pyfem.readthedocs.io/                                            #	
#                                                                              #
#  The original code can be downloaded from the web-site:                      #
#     http://www.wiley.com/go/deborst                                          #
#                                                                              #
#  The code is open source and intended for educational and scientific         #
#  purposes only. If you use PyFEM in your research, the developers would      #
#  be grateful if you could cite the book.                                     #    
#                                                                              #
#  Disclaimer:                                                                 #
#  The authors reserve all rights but do not guarantee that the code is        #
#  free from errors. Furthermore, the authors shall not be liable in any       #
#  event caused by the use of the program.                                     #
################################################################################

"""
Unit tests for pyfem.fem.NodeSet module.

Tests the NodeSet class functionality including:
- Node creation and retrieval
- Node coordinate operations
- Node group management
- Input/output operations
"""

import unittest
import numpy as np
from pyfem.fem.NodeSet import NodeSet


class TestNodeSet(unittest.TestCase):
    """Test suite for NodeSet class."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.nodes = NodeSet()
        
    def test_initialization(self):
        """Test NodeSet initialization."""
        self.assertEqual(self.nodes.rank, -1)
        self.assertEqual(len(self.nodes.groups), 0)
        self.assertEqual(len(self.nodes), 0)

    def test_add_nodes_2d(self):
        """Test adding 2D nodes."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.add(2, [1.0, 1.0])
        
        self.assertEqual(len(self.nodes), 3)
        self.assertEqual(self.nodes.rank, 2)
        
    def test_add_nodes_3d(self):
        """Test adding 3D nodes."""
        self.nodes.add(0, [0.0, 0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0, 0.0])
        self.nodes.add(2, [1.0, 1.0, 0.0])
        
        self.assertEqual(len(self.nodes), 3)
        self.assertEqual(self.nodes.rank, 3)

    def test_getNodeCoords_single_node(self):
        """Test getting coordinates for a single node."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 2.0])
        self.nodes.add(2, [3.0, 4.0])
        
        coords = self.nodes.getNodeCoords(1)
        np.testing.assert_array_equal(coords, [1.0, 2.0])

    def test_getNodeCoords_multiple_nodes(self):
        """Test getting coordinates for multiple nodes."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 2.0])
        self.nodes.add(2, [3.0, 4.0])
        
        coords = self.nodes.getNodeCoords([0, 2])
        expected = np.array([[0.0, 0.0], [3.0, 4.0]])
        np.testing.assert_array_equal(coords, expected)

    def test_getNodeCoords_node_group(self):
        """Test getting coordinates for a node group."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.add(2, [1.0, 1.0])
        
        self.nodes.addToGroup('Left', 0)
        self.nodes.addToGroup('Left', 2)
        
        coords = self.nodes.getNodeCoords('Left')
        self.assertEqual(len(coords), 2)

    def test_getNodeCoords_invalid_group(self):
        """Test getting coordinates for non-existent group raises error."""
        self.nodes.add(0, [0.0, 0.0])
        
        with self.assertRaises(SystemExit):
            self.nodes.getNodeCoords('NonExistent')

    def test_getNodeIDs(self):
        """Test getting node IDs from a group."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.add(2, [1.0, 1.0])
        
        self.nodes.addToGroup('Bottom', 0)
        self.nodes.addToGroup('Bottom', 1)
        
        nodeIDs = self.nodes.getNodeIDs('Bottom')
        self.assertEqual(set(nodeIDs), {0, 1})

    def test_getNodeIDs_invalid_group(self):
        """Test getting node IDs from non-existent group raises error."""
        with self.assertRaises(SystemExit):
            self.nodes.getNodeIDs('NonExistent')

    def test_addToGroup_new_group(self):
        """Test adding nodes to a new group."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        
        self.nodes.addToGroup('Top', 0)
        self.nodes.addToGroup('Top', 1)
        
        self.assertIn('Top', self.nodes.groups)
        self.assertEqual(len(self.nodes.groups['Top']), 2)

    def test_addToGroup_existing_group(self):
        """Test adding nodes to an existing group."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.add(2, [2.0, 0.0])
        
        self.nodes.addToGroup('Right', 0)
        self.nodes.addToGroup('Right', 1)
        self.nodes.addToGroup('Right', 2)
        
        self.assertEqual(len(self.nodes.groups['Right']), 3)

    def test_addToGroup_with_string_id(self):
        """Test adding nodes with string ID to group."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.addToGroup('Group1', '0')
        
        self.assertIn(0, self.nodes.groups['Group1'])

    def test_repr(self):
        """Test string representation of NodeSet."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.addToGroup('Test', 0)
        
        repr_str = repr(self.nodes)
        self.assertIn('Number of nodes', repr_str)
        self.assertIn('Test', repr_str)

    def test_multiple_groups(self):
        """Test managing multiple node groups."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.add(2, [1.0, 1.0])
        self.nodes.add(3, [0.0, 1.0])
        
        self.nodes.addToGroup('Left', 0)
        self.nodes.addToGroup('Left', 3)
        self.nodes.addToGroup('Right', 1)
        self.nodes.addToGroup('Right', 2)
        self.nodes.addToGroup('Bottom', 0)
        self.nodes.addToGroup('Bottom', 1)
        self.nodes.addToGroup('Top', 2)
        self.nodes.addToGroup('Top', 3)
        
        self.assertEqual(len(self.nodes.groups), 4)
        self.assertEqual(len(self.nodes.groups['Left']), 2)
        self.assertEqual(len(self.nodes.groups['Right']), 2)
        self.assertEqual(len(self.nodes.groups['Bottom']), 2)
        self.assertEqual(len(self.nodes.groups['Top']), 2)

    def test_numpy_integer_types(self):
        """Test that numpy integer types work for node IDs."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        
        # Test with np.int64
        coords = self.nodes.getNodeCoords(np.int64(1))
        np.testing.assert_array_equal(coords, [1.0, 0.0])
        
        # Test with np.int32
        coords = self.nodes.getNodeCoords(np.int32(0))
        np.testing.assert_array_equal(coords, [0.0, 0.0])

    def test_coordinates_are_numpy_arrays(self):
        """Test that returned coordinates are numpy arrays."""
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 2.0])
        
        coords = self.nodes.getNodeCoords(1)
        self.assertIsInstance(coords, np.ndarray)
        
        coords_multiple = self.nodes.getNodeCoords([0, 1])
        self.assertIsInstance(coords_multiple, np.ndarray)


class TestNodeSetEdgeCases(unittest.TestCase):
    """Test edge cases and error handling for NodeSet."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.nodes = NodeSet()

    def test_empty_nodeset(self):
        """Test operations on empty NodeSet."""
        self.assertEqual(len(self.nodes), 0)
        self.assertEqual(len(self.nodes.groups), 0)

    def test_empty_group(self):
        """Test empty node group."""
        self.nodes.groups['Empty'] = []
        nodeIDs = self.nodes.getNodeIDs('Empty')
        self.assertEqual(len(nodeIDs), 0)

    def test_large_node_ids(self):
        """Test handling of large node IDs."""
        large_id = 1000000
        self.nodes.add(large_id, [1.0, 2.0])
        coords = self.nodes.getNodeCoords(large_id)
        np.testing.assert_array_equal(coords, [1.0, 2.0])


if __name__ == '__main__':
    unittest.main()
