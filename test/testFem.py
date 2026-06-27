# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

"""
Unit tests for pyfem.fem.NodeSet and pyfem.fem.ElementSet modules.

Tests the NodeSet class functionality including:
- Node creation and retrieval
- Node coordinate operations
- Node group management
- Input/output operations

Tests the ElementSet class functionality including:
- Element creation and retrieval
- Element group management
- Element family operations
- Iterator operations
"""

import unittest
import numpy as np
from pyfem.fem.NodeSet import NodeSet
from pyfem.fem.ElementSet import ElementSet
from pyfem.util.dataStructures import Properties


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
        self.nodes = NodeSet()

        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.add(2, [1.0, 1.0])
        
        self.assertEqual(len(self.nodes), 3)
        self.assertEqual(self.nodes.getRank(), 2)
        
    def test_add_nodes_3d(self):
        """Test adding 3D nodes."""
        self.nodes = NodeSet()

        self.nodes.add(0, [0.0, 0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0, 0.0])
        self.nodes.add(2, [1.0, 1.0, 0.0])
        
        self.assertEqual(len(self.nodes), 3)
        self.assertEqual(self.nodes.getRank(), 3)

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
        self.nodes.add(5, [1.0, 0.0])
        self.nodes.add(6, [1.0, 0.0])
        
        coords = self.nodes.getNodeCoords([1,5,6])
        np.testing.assert_array_equal(coords[0,:], [1.0, 0.0])

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


class TestElementSet(unittest.TestCase):
    """Test suite for ElementSet class."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create a simple 2D nodeset
        self.nodes = NodeSet()
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.add(2, [1.0, 1.0])
        self.nodes.add(3, [0.0, 1.0])
        
        # Create a properties object with a mock model
        self.props = Properties()
        
    def test_initialization(self):
        """Test ElementSet initialization."""
        elements = ElementSet(self.nodes, self.props)
        
        self.assertEqual(len(elements), 0)
        self.assertEqual(len(elements.groups), 0)
        self.assertIsNotNone(elements.solverStat)

    def test_get_dof_types_empty(self):
        """Test getDofTypes on empty element set."""
        elements = ElementSet(self.nodes, self.props)
        dofTypes = elements.getDofTypes()
        
        self.assertEqual(len(dofTypes), 0)

    def test_add_to_group_new(self):
        """Test adding element to a new group."""
        elements = ElementSet(self.nodes, self.props)
        
        elements.addToGroup('Group1', 0)
        
        self.assertIn('Group1', elements.groups)
        self.assertIn(0, elements.groups['Group1'])
        self.assertEqual(len(elements.groups['Group1']), 1)

    def test_add_to_group_existing(self):
        """Test adding elements to an existing group."""
        elements = ElementSet(self.nodes, self.props)
        
        elements.addToGroup('Material1', 0)
        elements.addToGroup('Material1', 1)
        elements.addToGroup('Material1', 2)
        
        self.assertEqual(len(elements.groups['Material1']), 3)
        self.assertIn(0, elements.groups['Material1'])
        self.assertIn(1, elements.groups['Material1'])
        self.assertIn(2, elements.groups['Material1'])

    def test_add_group(self):
        """Test creating/replacing a group with specified IDs."""
        elements = ElementSet(self.nodes, self.props)
        
        elements.addGroup('TestGroup', [0, 1, 2])
        
        self.assertIn('TestGroup', elements.groups)
        self.assertEqual(elements.groups['TestGroup'], [0, 1, 2])
        
        # Test replacing an existing group
        elements.addGroup('TestGroup', [3, 4])
        self.assertEqual(elements.groups['TestGroup'], [3, 4])

    def test_iter_group_names(self):
        """Test iterating over group names."""
        elements = ElementSet(self.nodes, self.props)
        
        elements.addToGroup('Group1', 0)
        elements.addToGroup('Group2', 1)
        elements.addToGroup('Group3', 2)
        
        groupNames = elements.iterGroupNames()
        
        self.assertEqual(len(groupNames), 3)
        self.assertIn('Group1', groupNames)
        self.assertIn('Group2', groupNames)
        self.assertIn('Group3', groupNames)

    def test_element_group_count_single(self):
        """Test counting elements in a single group."""
        elements = ElementSet(self.nodes, self.props)
        
        elements.addGroup('Material1', [0, 1, 2, 3])
        
        count = elements.elementGroupCount('Material1')
        self.assertEqual(count, 4)

    def test_element_group_count_multiple(self):
        """Test counting elements across multiple groups."""
        elements = ElementSet(self.nodes, self.props)
        
        elements.addGroup('Material1', [0, 1])
        elements.addGroup('Material2', [2, 3, 4])
        
        count = elements.elementGroupCount(['Material1', 'Material2'])
        self.assertEqual(count, 5)

    def test_repr(self):
        """Test string representation of ElementSet."""
        elements = ElementSet(self.nodes, self.props)
        
        elements.addToGroup('Material1', 0)
        elements.addToGroup('Material1', 1)
        elements.addToGroup('Material2', 2)
        
        repr_str = repr(elements)
        
        self.assertIn('Number of elements', repr_str)
        self.assertIn('Material1', repr_str)
        self.assertIn('Material2', repr_str)

    def test_multiple_groups(self):
        """Test managing multiple element groups."""
        elements = ElementSet(self.nodes, self.props)
        
        elements.addToGroup('Steel', 0)
        elements.addToGroup('Steel', 1)
        elements.addToGroup('Concrete', 2)
        elements.addToGroup('Concrete', 3)
        elements.addToGroup('Concrete', 4)
        elements.addToGroup('Wood', 5)
        
        self.assertEqual(len(elements.groups), 3)
        self.assertEqual(len(elements.groups['Steel']), 2)
        self.assertEqual(len(elements.groups['Concrete']), 3)
        self.assertEqual(len(elements.groups['Wood']), 1)


class TestElementSetEdgeCases(unittest.TestCase):
    """Test edge cases and error handling for ElementSet."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.nodes = NodeSet()
        self.nodes.add(0, [0.0, 0.0])
        self.nodes.add(1, [1.0, 0.0])
        self.nodes.add(2, [1.0, 1.0])
        
        self.props = Properties()

    def test_empty_elementset(self):
        """Test operations on empty ElementSet."""
        elements = ElementSet(self.nodes, self.props)
        
        self.assertEqual(len(elements), 0)
        self.assertEqual(len(elements.groups), 0)
        
        # Test iteration
        element_list = list(elements)
        self.assertEqual(len(element_list), 0)

    def test_empty_group_iteration(self):
        """Test iterating over an empty group."""
        elements = ElementSet(self.nodes, self.props)
        elements.groups['EmptyGroup'] = []
        
        count = elements.elementGroupCount('EmptyGroup')
        self.assertEqual(count, 0)

    def test_commit_history_empty(self):
        """Test commitHistory on empty element set."""
        elements = ElementSet(self.nodes, self.props)
        
        # Should not raise an error
        elements.commitHistory()


if __name__ == '__main__':
    unittest.main()
