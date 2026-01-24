import unittest
import numpy as np
from pyfem.fem.Assembly import MatrixBuilder

class TestMatrixBuilder(unittest.TestCase):
    def setUp(self):
        self.nDofs = 4
        self.mb = MatrixBuilder(self.nDofs)

    def test_initialization(self):
        self.assertEqual(self.mb.nDofs, self.nDofs)
        np.testing.assert_array_equal(self.mb.B, np.zeros(self.nDofs))
        self.assertEqual(self.mb.c, 0.0)
        self.assertEqual(len(self.mb.val), 0)
        self.assertEqual(len(self.mb.row), 0)
        self.assertEqual(len(self.mb.col), 0)

    def test_append_and_getMatrix(self):
        a = np.eye(self.nDofs)
        dofs = np.arange(self.nDofs)
        self.mb.append(a, dofs)
        mat = self.mb.getMatrix()
        np.testing.assert_array_equal(mat.toarray(), a)

    def test_clear(self):
        a = np.eye(self.nDofs)
        dofs = np.arange(self.nDofs)
        self.mb.append(a, dofs)
        self.mb.B[:] = 1.0
        self.mb.c = 5.0
        self.mb.clear()
        np.testing.assert_array_equal(self.mb.B, np.zeros(self.nDofs))
        self.assertEqual(self.mb.c, 0.0)
        self.assertEqual(len(self.mb.val), 0)
        self.assertEqual(len(self.mb.row), 0)
        self.assertEqual(len(self.mb.col), 0)

    def test_B_and_c_accumulation(self):
        self.mb.B[1] += 2.5
        self.mb.c += 3.5
        self.assertEqual(self.mb.B[1], 2.5)
        self.assertEqual(self.mb.c, 3.5)

if __name__ == '__main__':
    unittest.main()
