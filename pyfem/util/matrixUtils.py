# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

"""Matrix-related utility helpers."""

from numpy import zeros


def skew(a_vec):
    """Return the skew-symmetric matrix associated with a vector."""
    mat = zeros(shape=(3, 3))

    mat[0, 1] = -a_vec[2]
    mat[0, 2] = a_vec[1]
    mat[1, 0] = a_vec[2]
    mat[1, 2] = -a_vec[0]
    mat[2, 0] = -a_vec[1]
    mat[2, 1] = a_vec[0]

    return mat
